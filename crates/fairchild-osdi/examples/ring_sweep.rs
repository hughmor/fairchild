//! Ring resonator wavelength sweep — runnable example.
//!
//! Reads all physical parameters from environment variables so the Python
//! driver script (`scripts/run_ring_sweep.py`) can tune them without
//! recompiling.  Writes sweep data as CSV to the path in RING_CSV_OUT
//! (default: ring_resonator_sweep.csv in the current directory).
//!
//! Env vars (all optional, defaults shown):
//!   RING_KAPPA_0       = 0.1       coupler power cross-coupling fraction
//!   RING_L_UM          = 100.0     ring circumference (µm)
//!   RING_N_G           = 4.2       group index
//!   RING_ALPHA_DB_CM   = 2.0       propagation loss (dB/cm)
//!   RING_POWER_MW      = 1.0       laser power (mW)
//!   RING_R_LOAD        = 1000.0    photodetector load resistance (Ω)
//!   RING_WL_START_NM   = 1544.0    sweep start wavelength (nm)
//!   RING_WL_END_NM     = 1558.0    sweep end wavelength (nm)
//!   RING_N_POINTS      = 101       number of wavelength points
//!   RING_CSV_OUT       = ring_resonator_sweep.csv  output CSV path
//!   RING_MODEL_DIR     = ../../../legacy/va-models/build     directory for .osdi files

use std::path::PathBuf;
use std::sync::Arc;

use fairchild_core::{
    build_devices, dc_op_nr_with_devices, CircuitTopology, DeviceRegistry, SimContext,
};
use fairchild_osdi::OsdiLibrary;
use fairchild_parser::parse_spice;

fn env_f64(key: &str, default: f64) -> f64 {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}
fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn cmt_transmission(
    wavelength_m: f64,
    kappa_0: f64,
    l_ring_m: f64,
    n_g: f64,
    alpha_db_cm: f64,
) -> f64 {
    let r = (1.0 - kappa_0).sqrt();
    let alpha_lin = alpha_db_cm * 1e2 / 8.685_895;
    let a = (-alpha_lin * l_ring_m / 2.0).exp();
    let beta = 2.0 * std::f64::consts::PI * n_g / wavelength_m;
    let phi = beta * l_ring_m;
    (r * r - 2.0 * r * a * phi.cos() + a * a) / (1.0 - 2.0 * r * a * phi.cos() + r * r * a * a)
}

fn cmt_resonance_nearest(lambda_center_m: f64, n_g: f64, l_ring_m: f64) -> f64 {
    let m = (n_g * l_ring_m / lambda_center_m).round() as u64;
    n_g * l_ring_m / m as f64
}

fn model_path(dir: &str, name: &str) -> PathBuf {
    PathBuf::from(dir).join(format!("{name}.osdi"))
}

fn main() {
    // ── Parameters ────────────────────────────────────────────────────────────
    let kappa_0 = env_f64("RING_KAPPA_0", 0.1);
    let l_ring_um = env_f64("RING_L_UM", 100.0);
    let n_g = env_f64("RING_N_G", 4.2);
    let alpha_db_cm = env_f64("RING_ALPHA_DB_CM", 2.0);
    let power_mw = env_f64("RING_POWER_MW", 1.0);
    let r_load = env_f64("RING_R_LOAD", 1000.0);
    let wl_start = env_f64("RING_WL_START_NM", 1544.0);
    let wl_end = env_f64("RING_WL_END_NM", 1558.0);
    let n_points = env_usize("RING_N_POINTS", 101);
    let csv_out =
        std::env::var("RING_CSV_OUT").unwrap_or_else(|_| "ring_resonator_sweep.csv".to_string());
    let model_dir = std::env::var("RING_MODEL_DIR").unwrap_or_else(|_| {
        // default: relative to this example's manifest dir at compile time
        format!(
            "{}/../../../legacy/va-models/build",
            env!("CARGO_MANIFEST_DIR")
        )
    });

    let l_ring_m = l_ring_um * 1e-6;

    eprintln!("Ring resonator sweep");
    eprintln!("  kappa_0={kappa_0}, L_ring={l_ring_um} µm, n_g={n_g}, alpha={alpha_db_cm} dB/cm");
    eprintln!("  P_in={power_mw} mW, R_load={r_load} Ω");
    eprintln!("  λ: {wl_start}–{wl_end} nm ({n_points} points)");
    eprintln!("  models: {model_dir}");
    eprintln!("  output: {csv_out}");

    // ── Load OSDI libraries ───────────────────────────────────────────────────
    let model_names = [
        "cw_laser",
        "directional_coupler",
        "waveguide",
        "photodetector",
    ];
    let mut libs = Vec::new();
    for name in &model_names {
        let path = model_path(&model_dir, name);
        if !path.exists() {
            eprintln!("error: model not found: {}", path.display());
            eprintln!("  Compile with: cd legacy/va-models && bash build.sh");
            std::process::exit(1);
        }
        let lib = Arc::new(unsafe { OsdiLibrary::open(&path) }.unwrap_or_else(|e| {
            eprintln!("error: dlopen {}: {e}", path.display());
            std::process::exit(1);
        }));
        libs.push(lib);
    }

    let mut registry = DeviceRegistry::new();
    for lib in &libs {
        lib.register_into(&mut registry);
    }

    // ── Build base netlist (topology is wavelength-independent) ───────────────
    // 3-wire topology: each optical port carries (re, im, lambda).
    // All lambda terminals connect to the shared `wl` net; the laser drives it.
    let base_netlist_str = format!(
        "* Ring resonator sweep\n\
         Xlaser    laser_re laser_im wl                                        cw_laser \
               power_mW={power_mw} wavelength_nm=1550.0\n\
         Xcoupler  laser_re laser_im wl  ring_fb_re ring_fb_im wl  \
               through_re through_im wl  ring_in_re ring_in_im wl  directional_coupler \
               kappa_0={kappa_0} wavelength_nm=1550.0\n\
         Xring     ring_in_re ring_in_im wl  ring_fb_re ring_fb_im wl  waveguide \
               L_um={l_ring_um} n_g={n_g} alpha_dB_cm={alpha_db_cm} wavelength_nm=1550.0\n\
         Xpd       through_re through_im wl  ph_a 0  photodetector\n\
         Rload     ph_a 0  {r_load}\n\
         .optical  laser_re laser_im  ring_in_re ring_in_im  ring_fb_re ring_fb_im  \
               through_re through_im  wl\n\
         .op\n\
         .end\n"
    );
    let base_netlist = parse_spice(&base_netlist_str).unwrap_or_else(|e| {
        eprintln!("error: netlist parse failed: {e}");
        std::process::exit(1);
    });
    let mut topo = CircuitTopology::build(&base_netlist);
    let ctx = SimContext::default();

    // ── Wavelength sweep ──────────────────────────────────────────────────────
    let wavelengths: Vec<f64> = (0..n_points)
        .map(|i| wl_start + (wl_end - wl_start) * i as f64 / (n_points - 1) as f64)
        .collect();

    let mut sweep: Vec<(f64, f64, f64)> = Vec::with_capacity(n_points); // (wl_nm, v_sim, t_cmt)

    for (idx, &wl_nm) in wavelengths.iter().enumerate() {
        let mut devices =
            build_devices(&base_netlist, &mut topo, &ctx, &registry).unwrap_or_else(|e| {
                eprintln!("error: build_devices at {wl_nm:.3} nm: {e}");
                std::process::exit(1);
            });
        for dev in &mut devices {
            dev.set_real_param("wavelength_nm", wl_nm);
        }
        let result = dc_op_nr_with_devices(&base_netlist, &topo, &mut devices, &ctx)
            .unwrap_or_else(|e| {
                eprintln!("error: NR failed at {wl_nm:.3} nm: {e}");
                std::process::exit(1);
            });
        let v_ph = result.node_voltage("ph_a").unwrap_or(0.0);
        let t_cmt = cmt_transmission(wl_nm * 1e-9, kappa_0, l_ring_m, n_g, alpha_db_cm);
        sweep.push((wl_nm, v_ph, t_cmt));

        if idx % 10 == 0 || idx == n_points - 1 {
            eprint!(
                "\r  [{}/{}] {:.3} nm  V(ph_a)={:.4} V   ",
                idx + 1,
                n_points,
                wl_nm,
                v_ph
            );
        }
    }
    eprintln!();

    // ── Write CSV ─────────────────────────────────────────────────────────────
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&csv_out).unwrap_or_else(|e| {
            eprintln!("error: cannot write {csv_out}: {e}");
            std::process::exit(1);
        });
        writeln!(f, "wavelength_nm,V_ph_a_V,T_cmt").unwrap();
        for &(wl, v, t) in &sweep {
            writeln!(f, "{wl:.6},{v:.8e},{t:.8e}").unwrap();
        }
    }

    // ── Summary ───────────────────────────────────────────────────────────────
    let v_max = sweep
        .iter()
        .map(|&(_, v, _)| v)
        .fold(f64::NEG_INFINITY, f64::max);
    let (sim_res_nm, v_min, _) = sweep
        .iter()
        .copied()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    let cmt_res_nm = cmt_resonance_nearest(sim_res_nm * 1e-9, n_g, l_ring_m) * 1e9;
    let dip_pct = (1.0 - v_min / v_max) * 100.0;

    println!("CSV written to {csv_out}");
    println!("Simulated resonance: {sim_res_nm:.3} nm  (V_min={v_min:.4e} V, dip={dip_pct:.1}%)");
    println!("CMT resonance:       {cmt_res_nm:.3} nm");
    println!("Δλ = {:.4} nm", (sim_res_nm - cmt_res_nm).abs());
}
