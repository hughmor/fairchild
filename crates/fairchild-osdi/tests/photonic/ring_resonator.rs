/// Ring resonator wavelength sweep — Phase 2 milestone.
///
/// Circuit topology:
///   Xlaser      laser_re laser_im                          cw_laser
///   Xcoupler    laser_re laser_im  ring_fb_re ring_fb_im
///               through_re through_im  ring_in_re ring_in_im  directional_coupler
///   Xring       ring_in_re ring_in_im  ring_fb_re ring_fb_im  waveguide (ring round-trip)
///   Xpd         through_re through_im  ph_a 0                 photodetector
///   Rload       ph_a 0  1k
///
/// The ring waveguide closes the feedback loop.  NR solves the algebraic loop
/// directly at each wavelength to find the steady-state SVEA envelope.
///
/// Physics parameters (L_ring = 100 µm, n_g = 4.2, kappa = 0.1, alpha = 2 dB/cm):
///   FSR   = λ² / (n_g · L_ring) ≈ 5.72 nm
///   λ_res = n_g · L_ring / m    where m = round(n_g · L_ring / λ_center)
///
/// Validation: resonance dip in simulated V(ph_a) must be within 0.1 nm of the
/// CMT analytical resonance wavelength.
use std::path::{Path, PathBuf};
use std::sync::Arc;

use fairchild_core::{
    build_devices, dc_op_nr_with_devices, CircuitTopology, Device, DeviceRegistry, SimContext,
};
use fairchild_osdi::{OsdiDevice, OsdiLibrary};
use fairchild_parser::parse_spice;

// Ring resonator physical parameters
const L_RING_UM: f64 = 100.0;
const N_G: f64 = 4.2;
const KAPPA_0: f64 = 0.1; // coupler cross-coupling power fraction
const ALPHA_DB_CM: f64 = 2.0; // waveguide propagation loss
const POWER_MW: f64 = 1.0; // laser power
const R_LOAD: f64 = 1e3; // load resistance

fn model_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(format!("../../../legacy/va-models/build/{name}.osdi"))
}

fn skip_if_missing(path: &Path) -> bool {
    if !path.exists() {
        eprintln!(
            "Skipping: {} not found. Compile legacy/va-models first with openvaf-r.",
            path.display()
        );
        true
    } else {
        false
    }
}

/// Build the ring resonator netlist string with given wavelength.
///
/// 3-wire topology: each optical port carries (re, im, lambda).
/// All lambda terminals connect to the shared `wl` net; the laser drives it.
fn ring_resonator_netlist(wavelength_nm: f64) -> String {
    format!(
        "* Ring resonator (wavelength {wavelength_nm:.3} nm)\n\
         Xlaser     laser_re laser_im wl                                        cw_laser \
               power_mW={POWER_MW} wavelength_nm={wavelength_nm:.6}\n\
         Xcoupler   laser_re laser_im wl  ring_fb_re ring_fb_im wl  \
               through_re through_im wl  ring_in_re ring_in_im wl  directional_coupler \
               kappa_0={KAPPA_0} wavelength_nm={wavelength_nm:.6}\n\
         Xring      ring_in_re ring_in_im wl  ring_fb_re ring_fb_im wl  waveguide \
               L_um={L_RING_UM} n_g={N_G} alpha_dB_cm={ALPHA_DB_CM} wavelength_nm={wavelength_nm:.6}\n\
         Xpd        through_re through_im wl  ph_a 0  photodetector\n\
         Rload      ph_a 0  1k\n\
         .optical   laser_re laser_im  ring_in_re ring_in_im  ring_fb_re ring_fb_im  \
               through_re through_im  wl\n\
         .op\n\
         .end\n"
    )
}

// ─── CMT analytical model ─────────────────────────────────────────────────────

/// Coupled-mode theory power transmission at the through port.
/// r = sqrt(1-kappa), a = round-trip amplitude, phi = round-trip phase.
fn cmt_transmission(wavelength_m: f64) -> f64 {
    let r = (1.0 - KAPPA_0).sqrt();
    let alpha_lin = ALPHA_DB_CM * 1e2 / 8.685_895; // dB/cm → Np/m
    let l_ring_m = L_RING_UM * 1e-6;
    let a = (-alpha_lin * l_ring_m / 2.0).exp(); // round-trip amplitude (one full loop = 2*half)

    let beta = 2.0 * std::f64::consts::PI * N_G / wavelength_m;
    let phi = beta * l_ring_m;

    (r * r - 2.0 * r * a * phi.cos() + a * a) / (1.0 - 2.0 * r * a * phi.cos() + r * r * a * a)
}

/// Find the resonance wavelength nearest to `lambda_center_m` using CMT.
fn cmt_resonance_nearest(lambda_center_m: f64) -> f64 {
    let l_ring_m = L_RING_UM * 1e-6;
    let m = (N_G * l_ring_m / lambda_center_m).round() as u64;
    N_G * l_ring_m / m as f64
}

// ─── Diagnostic tests ─────────────────────────────────────────────────────────

/// Low-level probe of the access() function to understand the OSDI id mapping.
/// Prints n_inst, n_model, all param names, and what offset access(MODEL|j) returns.
/// A debugging harness, not a test: it prints the descriptor's parameter table
/// and asserts nothing. `#[ignore]`d for the same reason as
/// `setup_runs_without_the_device_layer`.
#[test]
#[ignore = "diagnostic harness: prints, asserts nothing"]
fn access_ptr_diagnostic() {
    let path = model_path("cw_laser");
    if skip_if_missing(&path) {
        return;
    }

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen"));

    // Inspect the descriptor fields directly.
    let desc = lib.descriptors().next().expect("no descriptors");
    let n_total = desc.num_params as usize;
    let n_inst = desc.num_instance_params as usize;
    let n_opvar = desc.num_opvars as usize;
    let model_size = desc.model_size as usize;
    eprintln!("Descriptor: num_params={n_total}, num_instance_params={n_inst}, num_opvars={n_opvar}, model_size={model_size}");

    // Print all param/opvar names and their flags.
    if !desc.param_opvar.is_null() {
        let entries = unsafe { std::slice::from_raw_parts(desc.param_opvar, n_total + n_opvar) };
        for (idx, p) in entries.iter().enumerate() {
            let kind = if idx < n_inst {
                "INST"
            } else if idx < n_total {
                "MODEL"
            } else {
                "OPVAR"
            };
            let name = if p.name.is_null() {
                "<null>"
            } else {
                unsafe { std::ffi::CStr::from_ptr(*p.name).to_str().unwrap_or("?") }
            };
            eprintln!(
                "  param_opvar[{idx}] {kind}: name={name:20} flags={:#010x}",
                p.flags
            );
        }
    }

    // Build a device, setup, then probe access() for each model param id.
    let mut dev = OsdiDevice::from_library(Arc::clone(&lib), 0).expect("descriptor 0");
    let ctx = SimContext::default();
    dev.setup_model(&ctx);
    dev.setup_instance(&[Some(0), Some(1), Some(2)], &ctx);

    let model_raw = dev.model_ptr_raw();
    let n_model_words = model_size.div_ceil(8);
    eprintln!("\nModel buffer after setup (as f64):");
    for i in 0..n_model_words.min(12) {
        let v = unsafe { *((model_raw as *const f64).add(i)) };
        eprintln!("  [{}] offset {:3}: {}", i, i * 8, v);
    }

    eprintln!(
        "\naccess(null, model_ptr, PARA_KIND_MODEL|j, READ) for j=0..{}:",
        n_total + 2
    );
    use fairchild_osdi::ffi::{ACCESS_FLAG_READ, PARA_KIND_MODEL};
    if let Some(access_fn) = desc.access {
        for j in 0u32..(n_total as u32 + 2) {
            let id = PARA_KIND_MODEL | j;
            let ptr = unsafe { access_fn(std::ptr::null_mut(), model_raw, id, ACCESS_FLAG_READ) };
            if ptr.is_null() {
                eprintln!("  j={j}: -> null");
            } else {
                let off = unsafe { (ptr as *const u8).offset_from(model_raw as *const u8) };
                let val = unsafe { *(ptr as *const f64) };
                eprintln!("  j={j}: model_ptr+{off} = {val}");
            }
        }
    }

    // Check if set_real_param now works correctly.
    let ok = dev.set_real_param("power_mW", 4.0);
    eprintln!("\nset_real_param('power_mW', 4.0) → {ok}");
    let probe = dev.probe_model_param("power_mW");
    eprintln!("probe_model_param('power_mW') after set → {:?}", probe);

    let model_f64s: Vec<f64> = (0..n_model_words.min(12))
        .map(|i| unsafe { *((model_raw as *const f64).add(i)) })
        .collect();
    eprintln!("Model buffer after set_real_param: {:?}", model_f64s);

    // --- Node/Jacobian topology (cw_laser) ---
    eprintln!(
        "\nNode topology: num_nodes={}, num_terminals={}",
        desc.num_nodes, desc.num_terminals
    );
    let n_nodes_total = desc.num_nodes as usize;
    if !desc.nodes.is_null() {
        let nodes = unsafe { std::slice::from_raw_parts(desc.nodes, n_nodes_total) };
        for (i, node) in nodes.iter().enumerate() {
            let name = if node.name.is_null() {
                "<null>"
            } else {
                unsafe { std::ffi::CStr::from_ptr(node.name).to_str().unwrap_or("?") }
            };
            eprintln!(
                "  node[{i}] name={name:20} is_flow={} resist_res_off={}",
                node.is_flow, node.resist_residual_off
            );
        }
    }

    eprintln!(
        "\nJacobian: num_jacobian_entries={}, num_resistive={}",
        desc.num_jacobian_entries, desc.num_resistive_jacobian_entries
    );
    let n_jac = desc.num_jacobian_entries as usize;
    let n_resist = desc.num_resistive_jacobian_entries as usize;
    if !desc.jacobian_entries.is_null() && n_jac > 0 {
        let entries = unsafe { std::slice::from_raw_parts(desc.jacobian_entries, n_jac) };
        for (i, e) in entries.iter().enumerate() {
            let kind = if i < n_resist { "RESIST" } else { "REACT" };
            eprintln!(
                "  entry[{i}] {kind}: node_1={} node_2={} react_ptr_off={}",
                e.nodes.node_1, e.nodes.node_2, e.react_ptr_off
            );
        }
    }
}

/// Verify that set_real_param correctly updates a model parameter in the OSDI
/// instance so eval() picks up the new value. Tests power_mW: 1 mW → 4 mW.
/// Expected: V(laser_re) ≈ sqrt(power_W) since G_src >> R_load^-1.
#[test]
fn set_real_param_verification() {
    let path = model_path("cw_laser");
    if skip_if_missing(&path) {
        return;
    }

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen"));
    let mut registry = DeviceRegistry::new();
    lib.register_into(&mut registry);

    let netlist_str = "* laser set_real_param test\n\
                       Xlaser laser_re laser_im wl cw_laser\n\
                       Rre laser_re 0 1\n\
                       Rim laser_im 0 1\n\
                       .optical laser_re laser_im wl\n\
                       .op\n\
                       .end\n";
    let netlist = parse_spice(netlist_str).expect("parse netlist");
    let mut topo = CircuitTopology::build(&netlist);
    let ctx = SimContext::default();

    // Default power_mW = 1.0 → V(laser_re) ≈ sqrt(1e-3) ≈ 0.03162 V
    {
        let mut devs = build_devices(&netlist, &mut topo, &ctx, &registry).expect("build");
        let r = dc_op_nr_with_devices(&netlist, &topo, &mut devs, &ctx).expect("dc op default");
        let v = r.node_voltage("laser_re").unwrap();
        let expected = (1e-3_f64).sqrt();
        eprintln!("[default power_mW=1] V(laser_re)={v:.6e}  expected≈{expected:.6e}");
        assert!(
            (v - expected).abs() < 0.01 * expected,
            "default power: got {v:.6e}, expected {expected:.6e}"
        );
    }

    // set_real_param power_mW = 4.0 → V(laser_re) ≈ sqrt(4e-3) ≈ 0.06325 V
    {
        let mut devs = build_devices(&netlist, &mut topo, &ctx, &registry).expect("build");
        let ok = devs[0].set_real_param("power_mW", 4.0);
        eprintln!("[set_real_param] power_mW=4.0 returned: {ok}");
        assert!(ok, "set_real_param('power_mW', 4.0) must return true");

        let r = dc_op_nr_with_devices(&netlist, &topo, &mut devs, &ctx).expect("dc op 4mW");
        let v = r.node_voltage("laser_re").unwrap();
        let expected = (4e-3_f64).sqrt();
        eprintln!("[4mW] V(laser_re)={v:.6e}  expected≈{expected:.6e}");
        assert!(
            (v - expected).abs() < 0.01 * expected,
            "power_mW=4.0: got {v:.6e}, expected {expected:.6e}"
        );
    }
}

/// Single DC operating point at 1550 nm with all physical parameters
/// explicitly injected via set_real_param.  Prints every node voltage so
/// the circuit state can be inspected manually.  The only hard assertion is
/// that the laser amplitude is correct (|A_laser|² ≈ power_mW * 1e-3).
#[test]
fn ring_single_point_diagnostic() {
    let paths = [
        model_path("cw_laser"),
        model_path("waveguide"),
        model_path("directional_coupler"),
        model_path("photodetector"),
    ];
    for p in &paths {
        if skip_if_missing(p) {
            return;
        }
    }

    let libs: Vec<Arc<OsdiLibrary>> = paths
        .iter()
        .map(|p| Arc::new(unsafe { OsdiLibrary::open(p) }.expect("dlopen")))
        .collect();

    let mut registry = DeviceRegistry::new();
    for lib in &libs {
        lib.register_into(&mut registry);
    }

    let netlist_str = ring_resonator_netlist(1550.0);
    let netlist = parse_spice(&netlist_str).expect("parse netlist");
    let mut topo = CircuitTopology::build(&netlist);
    let ctx = SimContext::default();

    // build_devices element order: [laser=0, coupler=1, waveguide=2, PD=3]
    let mut devices = build_devices(&netlist, &mut topo, &ctx, &registry).expect("build_devices");

    let r_laser_pwr = devices[0].set_real_param("power_mW", POWER_MW);
    let r_laser_wl = devices[0].set_real_param("wavelength_nm", 1550.0);
    let r_kappa = devices[1].set_real_param("kappa_0", KAPPA_0);
    let r_coupler_wl = devices[1].set_real_param("wavelength_nm", 1550.0);
    let r_l_um = devices[2].set_real_param("L_um", L_RING_UM);
    let r_n_g = devices[2].set_real_param("n_g", N_G);
    let r_alpha = devices[2].set_real_param("alpha_dB_cm", ALPHA_DB_CM);
    let r_wg_wl = devices[2].set_real_param("wavelength_nm", 1550.0);

    eprintln!("set_real_param results:");
    eprintln!("  [laser]  power_mW={POWER_MW}       → {r_laser_pwr}");
    eprintln!("  [laser]  wavelength_nm=1550 → {r_laser_wl}");
    eprintln!("  [coupler] kappa_0={KAPPA_0}      → {r_kappa}");
    eprintln!("  [coupler] wavelength_nm=1550 → {r_coupler_wl}");
    eprintln!("  [wg]     L_um={L_RING_UM}       → {r_l_um}");
    eprintln!("  [wg]     n_g={N_G}           → {r_n_g}");
    eprintln!("  [wg]     alpha_dB_cm={ALPHA_DB_CM}  → {r_alpha}");
    eprintln!("  [wg]     wavelength_nm=1550 → {r_wg_wl}");

    let result = dc_op_nr_with_devices(&netlist, &topo, &mut devices, &ctx)
        .unwrap_or_else(|e| panic!("DC OP failed at 1550 nm: {e:?}"));

    eprintln!("\nAll node voltages at 1550.0 nm (iters={}):", result.iters);
    for (name, v) in result.all_voltages() {
        eprintln!("  V({name:20}) = {v:+.6e}");
    }

    let v_laser_re = result.node_voltage("laser_re").unwrap_or(0.0);
    let v_laser_im = result.node_voltage("laser_im").unwrap_or(0.0);
    let v_through_re = result.node_voltage("through_re").unwrap_or(0.0);
    let v_through_im = result.node_voltage("through_im").unwrap_or(0.0);
    let v_ring_fb_re = result.node_voltage("ring_fb_re").unwrap_or(0.0);
    let v_ring_fb_im = result.node_voltage("ring_fb_im").unwrap_or(0.0);
    let v_ph_a = result.node_voltage("ph_a").unwrap_or(0.0);

    let p_laser = v_laser_re.powi(2) + v_laser_im.powi(2);
    let p_through = v_through_re.powi(2) + v_through_im.powi(2);
    let p_ring_fb = v_ring_fb_re.powi(2) + v_ring_fb_im.powi(2);
    let expected_p_in = POWER_MW * 1e-3;
    let expected_v_ph_no_ring = expected_p_in * (1.0 - KAPPA_0) * R_LOAD;

    eprintln!("\nDerived quantities:");
    eprintln!("  P_laser   = {p_laser:.6e}  (expected ≈ {expected_p_in:.6e} W)");
    eprintln!("  P_through = {p_through:.6e}");
    eprintln!("  P_ring_fb = {p_ring_fb:.6e}");
    eprintln!("  V(ph_a)   = {v_ph_a:.6e}");
    eprintln!("  V(ph_a) expected (no ring, kappa=0.1) ≈ {expected_v_ph_no_ring:.6e}");

    assert!(
        (p_laser - expected_p_in).abs() < 0.05 * expected_p_in,
        "Laser power wrong: P_laser={p_laser:.4e} W, expected {expected_p_in:.4e} W"
    );
}

// ─── Sweep test ───────────────────────────────────────────────────────────────

#[test]
fn ring_resonator_wavelength_sweep() {
    let paths = [
        model_path("cw_laser"),
        model_path("waveguide"),
        model_path("directional_coupler"),
        model_path("photodetector"),
    ];
    for p in &paths {
        if skip_if_missing(p) {
            return;
        }
    }

    // Load all model libraries once.
    let libs: Vec<Arc<OsdiLibrary>> = paths
        .iter()
        .map(|p| Arc::new(unsafe { OsdiLibrary::open(p) }.expect("dlopen failed")))
        .collect();

    let mut registry = DeviceRegistry::new();
    for lib in &libs {
        lib.register_into(&mut registry);
    }

    // Wavelength sweep: 1544..=1558 nm in steps of 0.14 nm → 101 points.
    let n_points = 101usize;
    let wl_start = 1544.0_f64;
    let wl_end = 1558.0_f64;
    let wavelengths: Vec<f64> = (0..n_points)
        .map(|i| wl_start + (wl_end - wl_start) * i as f64 / (n_points - 1) as f64)
        .collect();

    let mut sweep_results: Vec<(f64, f64)> = Vec::with_capacity(n_points);

    // Parse one netlist per sweep point (different wavelength_nm in the X-element params).
    // A single netlist reference is enough since the topology is identical at every wavelength;
    // we build fresh devices per point and call set_real_param for correctness.
    let base_netlist_str = ring_resonator_netlist(1550.0);
    let base_netlist = parse_spice(&base_netlist_str).expect("parse netlist");
    let mut topo = CircuitTopology::build(&base_netlist);
    let ctx = SimContext::default();

    for &wl_nm in &wavelengths {
        let mut devices =
            build_devices(&base_netlist, &mut topo, &ctx, &registry).expect("build_devices");

        // Inject wavelength into every optical OSDI device.
        for dev in &mut devices {
            dev.set_real_param("wavelength_nm", wl_nm);
        }

        let result = dc_op_nr_with_devices(&base_netlist, &topo, &mut devices, &ctx)
            .unwrap_or_else(|e| panic!("NR failed at λ={wl_nm:.3} nm: {e:?}"));

        let v_ph = result.node_voltage("ph_a").unwrap();
        sweep_results.push((wl_nm, v_ph));
    }

    // Write sweep to CSV for Python plotting.
    let csv_path = std::env::temp_dir().join("ring_resonator_sweep.csv");
    {
        let mut f = std::fs::File::create(&csv_path).expect("create csv");
        use std::io::Write;
        writeln!(f, "wavelength_nm,V_ph_a_V,T_cmt").unwrap();
        for &(wl, v) in &sweep_results {
            let t_cmt = cmt_transmission(wl * 1e-9);
            writeln!(f, "{wl:.6},{v:.8e},{t_cmt:.8e}").unwrap();
        }
    }
    println!("Sweep CSV written to {}", csv_path.display());

    // Find the simulated transmission minimum and maximum.
    let v_off_resonance = sweep_results
        .iter()
        .map(|&(_, v)| v)
        .fold(f64::NEG_INFINITY, f64::max);
    let (sim_res_nm, v_min) = sweep_results
        .iter()
        .copied()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();

    // CMT reference: find resonance nearest to wherever the simulation found its dip.
    let cmt_res_nm = cmt_resonance_nearest(sim_res_nm * 1e-9) * 1e9;

    println!(
        "Ring resonator sweep summary:\n  Simulated resonance: {sim_res_nm:.3} nm  (V_min = {v_min:.4e} V)\n  CMT resonance:       {cmt_res_nm:.3} nm\n  delta_lambda = {:.4} nm  (tolerance: 0.1 nm)\n  Off-resonance V:     {v_off_resonance:.4e} V",
        (sim_res_nm - cmt_res_nm).abs()
    );

    // Sanity: minimum should be at least 2.5% below off-resonance peak.
    // (Nearest sweep point is ~0.02 nm from resonance, FWHM≈0.098 nm → apparent T_min≈0.930.)
    assert!(
        v_min < v_off_resonance * 0.975,
        "No resonance dip detected: V_min={v_min:.4e} V, off-res={v_off_resonance:.4e} V"
    );

    // Main validation: simulated resonance within 0.1 nm of CMT.
    assert!(
        (sim_res_nm - cmt_res_nm).abs() < 0.1,
        "Resonance mismatch: sim={sim_res_nm:.4} nm, CMT={cmt_res_nm:.4} nm, \
         Δ={:.4} nm > 0.1 nm tolerance",
        (sim_res_nm - cmt_res_nm).abs()
    );
}
