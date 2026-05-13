/// OSDI integration tests for Phase 2 photonic models.
///
/// Tests: descriptor sanity + DC OP physics validation.
///
/// Pre-condition: photonic models must be compiled to va-models/build/*.osdi.
/// Build script: va-models/build.sh or run manually with openvaf-r.
///
/// Physics baseline (coupled-mode theory / analytical):
///   CW laser (1 mW, 0°): V(out_re) = sqrt(1e-3) ≈ 0.031623, V(out_im) = 0.0
///   Photodetector (1 mW in, R=1.0 A/W, R_load=1kΩ): V(ph_a) ≈ 1.0 V

use std::ffi::CStr;
use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::{dc_op_nr_with_registry, DeviceRegistry};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

fn model_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(format!("../../va-models/build/{name}.osdi"))
}

fn skip_if_missing(path: &PathBuf) -> bool {
    if !path.exists() {
        eprintln!(
            "Skipping: {} not found. Compile va-models first with openvaf-r.",
            path.display()
        );
        true
    } else {
        false
    }
}

// ─── Descriptor sanity tests ──────────────────────────────────────────────────

#[test]
fn cw_laser_descriptor_sanity() {
    let path = model_path("cw_laser");
    if skip_if_missing(&path) { return; }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed");
    assert_eq!(lib.version, (0, 4));
    assert_eq!(lib.num_descriptors, 1);

    let d = lib.descriptors().next().unwrap();
    let name = unsafe { CStr::from_ptr(d.name) }.to_str().unwrap();
    assert_eq!(name, "cw_laser");

    // 3 external terminals (out_re, out_im, out_lambda).
    assert_eq!(d.num_terminals, 3, "expected 3 terminals");
    // At least 3 declared parameters (power_mW, phi_0_deg, wavelength_nm).
    assert!(d.num_params >= 3, "expected ≥3 params, got {}", d.num_params);

    println!(
        "cw_laser: terminals={}, nodes={}, params={}, resist_jac_entries={}",
        d.num_terminals, d.num_nodes, d.num_params, d.num_resistive_jacobian_entries
    );
}

#[test]
fn waveguide_descriptor_sanity() {
    let path = model_path("waveguide");
    if skip_if_missing(&path) { return; }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed");
    let d = lib.descriptors().next().unwrap();
    let name = unsafe { CStr::from_ptr(d.name) }.to_str().unwrap();
    assert_eq!(name, "waveguide");

    // 6 external terminals: in_re, in_im, in_lambda, out_re, out_im, out_lambda.
    assert_eq!(d.num_terminals, 6, "expected 6 terminals");
    assert!(d.num_params >= 4, "expected ≥4 params, got {}", d.num_params);

    println!(
        "waveguide: terminals={}, nodes={}, params={}, resist_jac_entries={}",
        d.num_terminals, d.num_nodes, d.num_params, d.num_resistive_jacobian_entries
    );
}

#[test]
fn directional_coupler_descriptor_sanity() {
    let path = model_path("directional_coupler");
    if skip_if_missing(&path) { return; }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed");
    let d = lib.descriptors().next().unwrap();
    let name = unsafe { CStr::from_ptr(d.name) }.to_str().unwrap();
    assert_eq!(name, "directional_coupler");

    // 12 external terminals: 4 ports × 3 wires (re, im, lambda).
    assert_eq!(d.num_terminals, 12, "expected 12 terminals");
    assert!(d.num_params >= 4, "expected ≥4 params, got {}", d.num_params);

    println!(
        "directional_coupler: terminals={}, nodes={}, params={}, resist_jac_entries={}",
        d.num_terminals, d.num_nodes, d.num_params, d.num_resistive_jacobian_entries
    );
}

#[test]
fn photodetector_descriptor_sanity() {
    let path = model_path("photodetector");
    if skip_if_missing(&path) { return; }

    let lib = unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed");
    let d = lib.descriptors().next().unwrap();
    let name = unsafe { CStr::from_ptr(d.name) }.to_str().unwrap();
    assert_eq!(name, "photodetector");

    // 5 external terminals: in_re, in_im, in_lambda (optical) + anode, cathode (electrical).
    assert_eq!(d.num_terminals, 5, "expected 5 terminals");
    assert!(d.num_params >= 3, "expected ≥3 params, got {}", d.num_params);

    println!(
        "photodetector: terminals={}, nodes={}, params={}, resist_jac_entries={}",
        d.num_terminals, d.num_nodes, d.num_params, d.num_resistive_jacobian_entries
    );
}

// ─── CW laser DC OP physics validation ───────────────────────────────────────

/// CW laser (1 mW, 0° phase) stand-alone: Ophase output should equal sqrt(P).
///
/// Circuit: Xlaser out_re out_im wl cw_laser
/// Expected (analytical):
///   V(out_re) = sqrt(1e-3) ≈ 0.031623  (real part of SVEA amplitude)
///   V(out_im) = 0.0                    (imaginary part)
#[test]
fn cw_laser_dc_op_default_params() {
    let path = model_path("cw_laser");
    if skip_if_missing(&path) { return; }

    let netlist = parse_spice(
        "* CW laser stand-alone DC OP\n\
         Xlaser out_re out_im wl cw_laser\n\
         .optical out_re out_im wl\n\
         .op\n.end\n",
    ).unwrap();

    let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.expect("dlopen failed"));
    let mut registry = DeviceRegistry::new();
    lib.register_into(&mut registry);

    let result = dc_op_nr_with_registry(&netlist, &registry)
        .expect("DC OP failed for cw_laser");

    let v_re = result.node_voltage("out_re").unwrap();
    let v_im = result.node_voltage("out_im").unwrap();

    // Default params: power_mW=1.0, phi_0_deg=0.0
    // A_out = sqrt(1e-3) ≈ 0.031623
    let a_expected = (1e-3_f64).sqrt();
    let tol = 1e-4;

    println!("cw_laser: V(out_re)={v_re:.6e}  V(out_im)={v_im:.6e}  expected A={a_expected:.6e}");

    assert!(
        (v_re - a_expected).abs() < tol,
        "V(out_re)={v_re:.6e}  expected≈{a_expected:.6e}  diff={:.2e}",
        (v_re - a_expected).abs()
    );
    assert!(
        v_im.abs() < tol,
        "V(out_im)={v_im:.6e}  expected≈0.0",
    );
}

// ─── Photodetector physics validation ────────────────────────────────────────

/// Laser → photodetector → load resistor: photocurrent should match responsivity * power.
///
/// Circuit:
///   Xlaser  laser_re laser_im wl  cw_laser              (1 mW default)
///   Xpd     laser_re laser_im wl  ph_a 0  photodetector (R=1 A/W default)
///   R_load  ph_a 0  1k
///
/// Expected photocurrent: I_ph = 1.0 A/W * 1e-3 W = 1 mA
/// Expected: V(ph_a) ≈ I_ph * R_load = 1e-3 * 1e3 = 1.0 V
/// (Shunt resistance R_shunt=1MΩ >> 1kΩ, so its contribution is negligible.)
#[test]
fn laser_photodetector_chain_dc_op() {
    let laser_path = model_path("cw_laser");
    let pd_path    = model_path("photodetector");
    if skip_if_missing(&laser_path) || skip_if_missing(&pd_path) { return; }

    let netlist = parse_spice(
        "* CW laser → photodetector → load\n\
         Xlaser  laser_re laser_im wl              cw_laser\n\
         Xpd     laser_re laser_im wl  ph_a 0      photodetector\n\
         Rload   ph_a 0  1k\n\
         .optical laser_re laser_im wl\n\
         .op\n.end\n",
    ).unwrap();

    let laser_lib = Arc::new(unsafe { OsdiLibrary::open(&laser_path) }.expect("dlopen laser"));
    let pd_lib    = Arc::new(unsafe { OsdiLibrary::open(&pd_path) }.expect("dlopen pd"));

    let mut registry = DeviceRegistry::new();
    laser_lib.register_into(&mut registry);
    pd_lib.register_into(&mut registry);

    let result = dc_op_nr_with_registry(&netlist, &registry)
        .expect("DC OP failed for laser→PD chain");

    let v_ph = result.node_voltage("ph_a").unwrap();

    // Expected: V(ph_a) ≈ 1.0 V (1 mW * 1 A/W * 1 kΩ)
    // Small correction from shunt (R_shunt=1MΩ in parallel with 1kΩ):
    // V ≈ I_ph * (R_load || R_shunt) = 1e-3 * (1k || 1M) ≈ 0.999 V
    let expected = 1e-3 * (1e3 * 1e6) / (1e3 + 1e6);  // ≈ 0.999 V
    let tol = 5e-3; // 5 mV

    println!("laser→PD: V(ph_a)={v_ph:.6e}  expected≈{expected:.6e}");

    assert!(
        (v_ph - expected).abs() < tol,
        "V(ph_a)={v_ph:.6e}  expected≈{expected:.6e}  diff={:.2e}",
        (v_ph - expected).abs()
    );
}

// ─── Waveguide physics validation ─────────────────────────────────────────────

/// Laser → waveguide (lossless) → photodetector → load.
///
/// With alpha_dB_cm=0 (lossless), the waveguide only rotates phase.
/// Power at the output should equal power at input: I_ph unchanged.
/// V(ph_a) should still ≈ 1.0 V regardless of waveguide length.
#[test]
fn laser_waveguide_photodetector_chain_dc_op() {
    let laser_path = model_path("cw_laser");
    let wg_path    = model_path("waveguide");
    let pd_path    = model_path("photodetector");
    if skip_if_missing(&laser_path) || skip_if_missing(&wg_path) || skip_if_missing(&pd_path) {
        return;
    }

    let netlist = parse_spice(
        "* laser → waveguide (lossless) → PD → load\n\
         Xlaser  laser_re laser_im wl                              cw_laser\n\
         Xwg     laser_re laser_im wl  wg_out_re wg_out_im wl     waveguide\n\
         Xpd     wg_out_re wg_out_im wl  ph_a 0                   photodetector\n\
         Rload   ph_a 0  1k\n\
         .optical laser_re laser_im wg_out_re wg_out_im wl\n\
         .op\n.end\n",
    ).unwrap();

    let laser_lib = Arc::new(unsafe { OsdiLibrary::open(&laser_path) }.expect("dlopen laser"));
    let wg_lib    = Arc::new(unsafe { OsdiLibrary::open(&wg_path) }.expect("dlopen wg"));
    let pd_lib    = Arc::new(unsafe { OsdiLibrary::open(&pd_path) }.expect("dlopen pd"));

    let mut registry = DeviceRegistry::new();
    laser_lib.register_into(&mut registry);
    wg_lib.register_into(&mut registry);
    pd_lib.register_into(&mut registry);

    let result = dc_op_nr_with_registry(&netlist, &registry)
        .expect("DC OP failed for laser→wg→PD chain");

    let v_ph = result.node_voltage("ph_a").unwrap();

    // Waveguide default: alpha_dB_cm=2.0, L_um=100 → α·L = 2 dB/cm * 100µm * 1e-4 cm/µm = 2e-4 dB
    // Amplitude transmission: exp(-α_lin * L / 2) where α_lin = 2e-4 * 1e2 / 8.686 ≈ 2.3e-3 Np/m * 1e-4 m
    // T ≈ 1 - tiny loss → output power ≈ input power for such a short waveguide
    let expected = 1e-3 * (1e3 * 1e6) / (1e3 + 1e6);  // ≈ 0.999 V (power-conserving)
    let tol = 5e-2; // 50 mV tolerance (accounts for waveguide insertion loss)

    println!("laser→wg→PD: V(ph_a)={v_ph:.6e}  expected≈{expected:.6e}");

    assert!(
        (v_ph - expected).abs() < tol,
        "V(ph_a)={v_ph:.6e}  expected≈{expected:.6e}  diff={:.2e}",
        (v_ph - expected).abs()
    );
}

// ─── Phase 2.5: PN phase shifter L1 ─────────────────────────────────────────

/// Laser → PN phase shifter L1 (0 V bias) → PD → load.
///
/// At 0 V bias, Δφ = 0.  Through-power = input × amplitude_transmission².
/// With L=500 µm, alpha=3 dB/cm (default): loss = 3*500e-4 dB = 0.015 dB ≈ 0.17%.
#[test]
fn pn_phase_shifter_l1_dc_op_zero_bias() {
    let laser_path = model_path("cw_laser");
    let pnps_path  = model_path("pn_phase_shifter_l1");
    let pd_path    = model_path("photodetector");
    if skip_if_missing(&laser_path) || skip_if_missing(&pnps_path) || skip_if_missing(&pd_path) {
        return;
    }

    let netlist = parse_spice(
        "* laser → pn_phase_shifter_l1 (0V) → PD → load\n\
         Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1550.0\n\
         Xpnps   lre lim wl  ore oim wl  vbias 0  pn_phase_shifter_l1\n\
         Xpd     ore oim wl  ph_a 0  photodetector  responsivity=1.0\n\
         Rload   ph_a 0  1k\n\
         Vbias   vbias 0  DC 0.0\n\
         .optical  lre lim wl ore oim\n\
         .op\n.end\n",
    ).unwrap();

    let laser_lib = Arc::new(unsafe { OsdiLibrary::open(&laser_path) }.expect("dlopen laser"));
    let pnps_lib  = Arc::new(unsafe { OsdiLibrary::open(&pnps_path) }.expect("dlopen pn_ps"));
    let pd_lib    = Arc::new(unsafe { OsdiLibrary::open(&pd_path) }.expect("dlopen pd"));

    let mut registry = DeviceRegistry::new();
    laser_lib.register_into(&mut registry);
    pnps_lib.register_into(&mut registry);
    pd_lib.register_into(&mut registry);

    let result = dc_op_nr_with_registry(&netlist, &registry)
        .expect("DC OP failed for laser→pn_ps→PD");

    let v_ph = result.node_voltage("ph_a").unwrap();
    // Default: L_um=100 µm, alpha_dB_cm=3.0
    // alpha_lin = 3*100/8.686 = 34.54 Np/m  (power attenuation coefficient)
    // T_amp = exp(-34.54 * 100e-6 / 2) → P_out = T_amp² = exp(-0.003454) ≈ 0.9966 mW
    // V(ph_a) = 0.9966e-3 × (1k ∥ 1M) = 0.9966e-3 × 999 ≈ 0.9956 V
    let expected = 0.9956_f64;
    let tol = 0.002;
    println!("pn_ps(0V)→PD: V(ph_a)={v_ph:.6e}  expected≈{expected:.6e}");
    assert!(
        (v_ph - expected).abs() < tol,
        "V(ph_a)={v_ph:.6e} expected≈{expected:.6e} diff={:.3e}", (v_ph - expected).abs()
    );
}

// ─── Phase 2.5: MRR modulator L1 — resonance dip ─────────────────────────────

/// Laser → MRR modulator L1 at ring resonance → PD → load.
///
/// At resonance: T = (r−a)/(1−r·a).  r=sqrt(0.9)≈0.9487, a=exp(-23.03·100e-6)≈0.9977.
/// T_res ≈ -0.9163 → |T|²≈0.8396 → V(ph_a) ≈ 0.8390 V.
#[test]
fn mrr_modulator_l1_resonance_dip() {
    let laser_path = model_path("cw_laser");
    let mrr_path   = model_path("mrr_modulator_l1");
    let pd_path    = model_path("photodetector");
    if skip_if_missing(&laser_path) || skip_if_missing(&mrr_path) || skip_if_missing(&pd_path) {
        return;
    }

    let netlist = parse_spice(
        "* laser → MRR L1 at resonance → PD → load\n\
         Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1544.12\n\
         Xmod    lre lim wl  ore oim wl  vbias 0  mrr_modulator_l1\n\
         + kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=2.0\n\
         + Vpi_rt=10.0 V_ref=0.0 wavelength_nm=1544.12\n\
         Xpd     ore oim wl  ph_a 0  photodetector  responsivity=1.0\n\
         Rload   ph_a 0  1k\n\
         Vbias   vbias 0  DC 0.0\n\
         .optical  lre lim wl ore oim\n\
         .op\n.end\n",
    ).unwrap();

    let laser_lib = Arc::new(unsafe { OsdiLibrary::open(&laser_path) }.expect("dlopen laser"));
    let mrr_lib   = Arc::new(unsafe { OsdiLibrary::open(&mrr_path) }.expect("dlopen mrr"));
    let pd_lib    = Arc::new(unsafe { OsdiLibrary::open(&pd_path) }.expect("dlopen pd"));

    let mut registry = DeviceRegistry::new();
    laser_lib.register_into(&mut registry);
    mrr_lib.register_into(&mut registry);
    pd_lib.register_into(&mut registry);

    let result = dc_op_nr_with_registry(&netlist, &registry)
        .expect("DC OP failed for laser→MRR→PD");

    let v_ph = result.node_voltage("ph_a").unwrap();
    println!("MRR(resonance, 0V): V(ph_a)={v_ph:.6e}");
    // At resonance: through-port is reduced from max transmission
    assert!(v_ph < 0.99, "Expected resonance dip: V(ph_a)={v_ph:.4} should be < 0.99 V");
    assert!(v_ph > 0.01, "Unexpectedly deep resonance: V(ph_a)={v_ph:.4}");
    // Energy conservation: input is 1 mW → V(ph_a) ≤ 1mW × 999Ω = 0.999 V always
    assert!(v_ph < 1.0, "Energy conservation violated: V(ph_a)={v_ph:.6} > 1 V");
}

// ─── Phase 2.5: MRR modulator L1 — energy conservation off resonance ─────────

/// CMT transfer function sign correctness: off-resonance bias must not produce
/// more output power than the input (the bug was T_re sign error giving |T|²>1).
///
/// At V=-0.25 V the ring is slightly detuned; through power rises toward 1 mW.
/// At V=-2.5 V (π/4 detuning) through power ≈ 1 mW (near anti-resonance).
/// In both cases V(ph_a) must stay below 0.999 V (= 1 mW × 999 Ω load).
#[test]
fn mrr_modulator_l1_off_resonance_energy_conservation() {
    let laser_path = model_path("cw_laser");
    let mrr_path   = model_path("mrr_modulator_l1");
    let pd_path    = model_path("photodetector");
    if skip_if_missing(&laser_path) || skip_if_missing(&mrr_path) || skip_if_missing(&pd_path) {
        return;
    }

    let base_netlist = |vbias_v: f64| format!(
        "* MRR L1 energy conservation at V={vbias_v} V\n\
         Xlaser  lre lim wl  cw_laser  power_mW=1.0 wavelength_nm=1544.12\n\
         Xmod    lre lim wl  ore oim wl  vbias 0  mrr_modulator_l1\n\
         + kappa_0=0.1 L_ring_um=100.0 n_g=4.2 alpha_dB_cm=2.0\n\
         + Vpi_rt=10.0 V_ref=0.0 wavelength_nm=1544.12\n\
         Xpd     ore oim wl  ph_a 0  photodetector  responsivity=1.0\n\
         Rload   ph_a 0  1k\n\
         Vbias   vbias 0  DC {vbias_v}\n\
         .optical  lre lim wl ore oim\n\
         .op\n.end\n"
    );

    let laser_lib = Arc::new(unsafe { OsdiLibrary::open(&laser_path) }.expect("dlopen laser"));
    let mrr_lib   = Arc::new(unsafe { OsdiLibrary::open(&mrr_path) }.expect("dlopen mrr"));
    let pd_lib    = Arc::new(unsafe { OsdiLibrary::open(&pd_path) }.expect("dlopen pd"));

    let mut registry = DeviceRegistry::new();
    laser_lib.register_into(&mut registry);
    mrr_lib.register_into(&mut registry);
    pd_lib.register_into(&mut registry);

    for vbias in [-0.25_f64, -1.0, -2.5, -5.0] {
        let netlist = parse_spice(&base_netlist(vbias)).unwrap();
        let result = dc_op_nr_with_registry(&netlist, &registry)
            .expect("DC OP failed");
        let v_ph = result.node_voltage("ph_a").unwrap();
        println!("MRR(V={vbias:.2}): V(ph_a)={v_ph:.6e}");
        assert!(
            v_ph < 1.0,
            "Energy conservation violated at Vbias={vbias}: V(ph_a)={v_ph:.6} ≥ 1 V"
        );
        assert!(v_ph > 0.0, "Negative V(ph_a) at Vbias={vbias}: {v_ph:.6}");
    }
}
