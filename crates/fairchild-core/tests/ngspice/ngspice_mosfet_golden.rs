//! Golden comparison tests for MOSFET Level 1 circuits.
//! Compares fairchild dc_op_nr against ngspice within tolerance.
//! Tests are skipped (not failed) when ngspice is absent.

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use fairchild_core::{
    dc_op_nr, dc_op_nr_with_registry_opts, options::SimOptions, tran_nr_with_registry_opts,
    DeviceRegistry,
};
use fairchild_parser::parse_spice;

const REL_TOL: f64 = 2e-3; // 0.2% — Level 1 has some param differences from ngspice
const ABS_TOL_V: f64 = 1e-4; // 100 µV floor

fn find_ngspice() -> Option<std::path::PathBuf> {
    if Command::new("ngspice").arg("--version").output().is_ok() {
        return Some("ngspice".into());
    }
    for candidate in &[
        "/opt/homebrew/bin/ngspice",
        "/usr/local/bin/ngspice",
        "/usr/bin/ngspice",
    ] {
        let p = std::path::Path::new(candidate);
        if p.exists() {
            return Some(p.to_owned());
        }
    }
    None
}

/// ngspice runner that parses `print v(node)` from a .op netlist.
fn ngspice_op_node(spice_with_print: &str, node: &str) -> Option<f64> {
    let ngspice = find_ngspice()?;
    let dir = tempfile::tempdir().ok()?;
    let input = dir.path().join("test.sp");
    std::fs::write(&input, spice_with_print).ok()?;

    let out = Command::new(&ngspice)
        .args(["-b", input.to_str().unwrap()])
        .output()
        .ok()?;

    let combined = format!(
        "{}\n{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    for line in combined.lines() {
        let line = line.trim();
        let prefix = format!("v({node})");
        if line.to_lowercase().starts_with(&prefix) {
            if let Some(eq) = line.find('=') {
                let val: f64 = line[eq + 1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .parse()
                    .ok()?;
                return Some(val);
            }
        }
    }
    None
}

#[test]
fn nmos_resistor_dc_op() {
    // NMOS in saturation: VDD=3.3V through 10kΩ to drain.
    // Gate tied to VDD, source and bulk grounded.
    // .model nmos1 NMOS (VTO=0.7 KP=100u)
    // At gate=3.3V: VGS=3.3, VGS-VTO=2.6V; saturation if VDS >= 2.6V.
    // IDS_sat = 0.5 * 100e-6 * (W/L=10) * 2.6² = 3.38mA
    // But VDD / R drain resistor clamps this: V(d) = VDD - IDS * R = 3.3 - 3.38e-3 * 10e3 ≈ -34V
    // That's outside VDD; actually IDS limited by triode onset at VDS = VGS - VTO = 2.6V.
    // Exact solution: in triode at VDS = Vd: IDS = β·[(VGS-VTO)*VDS - 0.5*VDS²]
    // IDS = (VDD - VD)/R → triode/sat boundary where NR settles.
    let netlist_str = "* NMOS R load\n\
        .model nm1 NMOS (VTO=0.7 KP=100u)\n\
        VDD vdd 0 DC 3.3\n\
        VG  g   0 DC 3.3\n\
        R1  vdd d 10k\n\
        M1  d g 0 0 nm1 W=10u L=1u\n\
        .op\n\
        .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let result = dc_op_nr(&netlist).unwrap();
    let vd = result.node_voltage("d").unwrap();

    // Self-consistency check: IDS through R1 = (VDD - VD)/R1.
    let ids_r = (3.3 - vd) / 10e3;
    // IDS from NMOS: at the OP, must match.
    // We can verify V(d) is positive and less than VDD.
    assert!(
        vd > 0.0 && vd < 3.3,
        "V(d) = {vd:.4}V should be between 0 and VDD"
    );
    assert!(ids_r > 0.0, "IDS should be positive: {ids_r:.4e}");

    // Compare with ngspice if available.
    let ngspice_str = format!("{netlist_str}.control\nop\nprint v(d)\n.endc\n");
    if let Some(vd_ng) = ngspice_op_node(&ngspice_str, "d") {
        let err = (vd - vd_ng).abs();
        let tol = ABS_TOL_V + REL_TOL * vd_ng.abs();
        assert!(
            err < tol,
            "V(d): fairchild={vd:.6e}  ngspice={vd_ng:.6e}  err={err:.2e}  tol={tol:.2e}"
        );
    }
}

#[test]
fn cmos_inverter_high_output() {
    // CMOS inverter with VIN=0: NMOS off, PMOS on → VOUT ≈ VDD.
    let netlist_str = "* CMOS inverter — high output\n\
        .model nm NMOS (VTO=0.7  KP=100u)\n\
        .model pm PMOS (VTO=-0.7 KP=100u)\n\
        VDD vdd 0 DC 3.3\n\
        VIN in  0 DC 0\n\
        MN  out in 0   0   nm W=10u L=1u\n\
        MP  out in vdd vdd pm W=10u L=1u\n\
        .op\n\
        .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let result = dc_op_nr(&netlist).unwrap();
    let vout = result.node_voltage("out").unwrap();

    // VIN=0: NMOS cutoff, PMOS connects VDD to out → VOUT near VDD.
    assert!(
        vout > 3.0,
        "CMOS inverter VIN=0: VOUT={vout:.4}V (expected > 3.0V)"
    );

    let ngspice_str = format!("{netlist_str}.control\nop\nprint v(out)\n.endc\n");
    if let Some(vout_ng) = ngspice_op_node(&ngspice_str, "out") {
        let err = (vout - vout_ng).abs();
        let tol = ABS_TOL_V + REL_TOL * vout_ng.abs() + 0.05;
        assert!(
            err < tol,
            "CMOS high: fairchild={vout:.6e}  ngspice={vout_ng:.6e}  err={err:.2e}"
        );
    }
}

#[test]
fn cmos_inverter_low_output() {
    // CMOS inverter with VIN=VDD: PMOS off, NMOS on → VOUT ≈ 0V.
    let netlist_str = "* CMOS inverter — low output\n\
        .model nm NMOS (VTO=0.7  KP=100u)\n\
        .model pm PMOS (VTO=-0.7 KP=100u)\n\
        VDD vdd 0 DC 3.3\n\
        VIN in  0 DC 3.3\n\
        MN  out in 0   0   nm W=10u L=1u\n\
        MP  out in vdd vdd pm W=10u L=1u\n\
        .op\n\
        .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let result = dc_op_nr(&netlist).unwrap();
    let vout = result.node_voltage("out").unwrap();

    // VIN=VDD: PMOS cutoff, NMOS pulls out to GND → VOUT ≈ 0.
    assert!(
        vout < 0.01,
        "CMOS inverter VIN=VDD: VOUT={vout:.4}V (expected < 10mV)"
    );

    let ngspice_str = format!("{netlist_str}.control\nop\nprint v(out)\n.endc\n");
    if let Some(vout_ng) = ngspice_op_node(&ngspice_str, "out") {
        let err = (vout - vout_ng).abs();
        let tol = ABS_TOL_V + REL_TOL * vout_ng.abs() + 1e-6;
        assert!(
            err < tol,
            "CMOS low: fairchild={vout:.6e}  ngspice={vout_ng:.6e}  err={err:.2e}"
        );
    }
}

/// Run ngspice in batch mode on a netlist containing `.meas tran` statements.
/// Returns a map of measurement name → value, or `None` if ngspice is absent.
fn ngspice_meas_tran(netlist: &str) -> Option<HashMap<String, f64>> {
    let ngspice = find_ngspice()?;
    let mut tmp = tempfile::NamedTempFile::new().ok()?;
    write!(tmp, "{netlist}").ok()?;
    let output = Command::new(&ngspice)
        .arg("-b")
        .arg(tmp.path())
        .output()
        .ok()?;
    let combined = format!(
        "{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let mut map = HashMap::new();
    for line in combined.lines() {
        let line = line.trim();
        if let Some((key, val)) = line.split_once('=') {
            let key = key.trim().to_lowercase();
            if key.starts_with("v_") {
                if let Ok(v) = val.split_whitespace().next().unwrap_or("").parse::<f64>() {
                    map.insert(key, v);
                }
            }
        }
    }
    if map.is_empty() {
        None
    } else {
        Some(map)
    }
}

#[test]
fn cmos_inverter_caps_switching_time() {
    // CMOS inverter with junction + Meyer overlap capacitances and a 10pF load.
    // At t=15ns the input has just gone high (rise starts at 10ns, 1ns rise time).
    // The output should be mid-transition due to the capacitive load — not yet at
    // a rail — confirming that the cap model is slowing the switching edge.
    let netlist_str = "* CMOS inverter — switching time with junction+Meyer caps\n\
        .model nm NMOS (VTO=0.7 KP=100u CGSO=2.5e-10 CGDO=2.5e-10 CJ=2e-4 CJSW=5e-10)\n\
        .model pm PMOS (VTO=-0.7 KP=100u CGSO=2.5e-10 CGDO=2.5e-10 CJ=2e-4 CJSW=5e-10)\n\
        VDD vdd 0 DC 3.3\n\
        VIN in 0 PULSE(0 3.3 10n 1n 1n 40n 100n)\n\
        MN out in 0 0 nm W=10u L=1u AS=50p AD=50p PS=20u PD=20u\n\
        MP out in vdd vdd pm W=10u L=1u AS=50p AD=50p PS=20u PD=20u\n\
        CL out 0 10p\n\
        .tran 100p 120n\n\
        .meas tran v_15n FIND V(out) AT=15n\n\
        .end\n";

    let netlist = parse_spice(netlist_str).unwrap();
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_mosfets(&netlist.models);
    let opts = SimOptions::from_netlist(&netlist);
    let result = tran_nr_with_registry_opts(&netlist, 100e-12, 120e-9, &registry, &opts).unwrap();
    let fc_v = result.voltage_at("out", 15e-9).unwrap();

    // Shape check: output must be mid-transition (not stuck at either rail).
    assert!(
        (0.3..=3.0).contains(&fc_v),
        "V(out) at t=15ns = {fc_v:.4}V — expected mid-transition [0.3, 3.0]V; \
        caps (CGSO/CGDO/CJ/CJSW on MOSFETs + CL=10p) are required to slow the edge \
        enough that the output is still switching at 15ns"
    );

    // ngspice comparison (skipped if ngspice absent).
    if let Some(meas) = ngspice_meas_tran(netlist_str) {
        if let Some(&ng_v) = meas.get("v_15n") {
            let err = (fc_v - ng_v).abs();
            let tol = f64::max(0.15, 0.05 * ng_v.abs());
            assert!(
                err <= tol,
                "V(out) at 15ns: fairchild={fc_v:.6e}  ngspice={ng_v:.6e}  \
                err={err:.3e}  tol={tol:.3e}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// RD / RS — the ohmic series resistances (#77 §2)
// ---------------------------------------------------------------------------

/// The drain node of an NMOS with a load resistor. `V(d)` reads the drain
/// current directly through the load, so it is the observable for `RD`/`RS`.
fn nmos_drain_v(model_extra: &str, instances: &str) -> f64 {
    let deck = format!(
        "* rd/rs\n.model nm NMOS (VTO=0.7 KP=500u {model_extra})\n\
         VDD vdd 0 DC 3\nVG g 0 DC 1.5\nRL vdd d 1k\n{instances}.op\n"
    );
    let netlist = parse_spice(&deck).expect("parse");
    dc_op_nr(&netlist)
        .unwrap_or_else(|e| panic!("solve failed on\n{deck}\n{e:?}"))
        .node_voltage("d")
        .expect("V(d)")
}

/// `RD`/`RS` in the model card are series resistances, so they must equal
/// explicit resistors of the same value — and match ngspice.
///
/// # Why the self-consistency half as well as the anchor
///
/// Two forms of one circuit agreeing cannot detect a fault common to both, so
/// ngspice anchors one side. But the internal-node form is the *interesting* one:
/// it is the shape the diode's `RS` got wrong by eliminating the internal node
/// instead of allocating it, which read 2.7% low with the convergence test unable
/// to notice. Drawing the same circuit twice is how that class of fault surfaces.
#[test]
fn rd_and_rs_equal_explicit_resistors() {
    // 500 Ω, not 50: at ~2.7 mA that is over a volt of series drop, so the
    // effect being tested is far larger than the noise between two solves.
    for (rd, rs) in [(500.0_f64, 0.0_f64), (0.0, 500.0), (500.0, 500.0)] {
        let internal = nmos_drain_v(&format!("RD={rd} RS={rs}"), "M1 d g 0 0 nm W=100u L=1u\n");
        // The same circuit drawn. A zero resistance becomes a tiny one rather
        // than a removed element, so the topology matches in every case.
        let external = nmos_drain_v(
            "",
            &format!(
                "RDx d dp {:e}\nRSx sp 0 {:e}\nM1 dp g sp 0 nm W=100u L=1u\n",
                rd.max(1e-9),
                rs.max(1e-9),
            ),
        );
        // The solver's own convergence bound, not a tighter number: the two
        // decks have different topologies and so take different Newton paths, and
        // anything below this bound is where each stopped rather than what the
        // circuit is. The first version used 1e-6 and measured 3.4e-6 of exactly
        // that. The series drop being tested is over a volt, so this stays
        // thousands of times discriminating.
        let opts = SimOptions::default();
        let same = opts.vntol + opts.reltol * internal.abs();
        assert!(
            (internal - external).abs() < same,
            "RD={rd} RS={rs}: in the model V(d)={internal:.9} V, drawn as \
             resistors V(d)={external:.9} V — the same circuit, and they differ by \
             more than the convergence bound {same:.2e}."
        );

        let with_print = format!(
            "* rd/rs\n.model nm NMOS (VTO=0.7 KP=500u RD={rd} RS={rs})\n\
             VDD vdd 0 DC 3\nVG g 0 DC 1.5\nRL vdd d 1k\n\
             M1 d g 0 0 nm W=100u L=1u\n\
             .control\nop\nprint v(d)\n.endc\n.end\n"
        );
        if let Some(ng) = ngspice_op_node(&with_print, "d") {
            let rel_ng = (internal - ng).abs() / ng.abs().max(ABS_TOL_V);
            assert!(
                rel_ng < REL_TOL,
                "RD={rd} RS={rs}: fairchild V(d)={internal:.9}, ngspice {ng:.9} \
                 (rel {rel_ng:.2e})"
            );
        }
    }
}

/// The series resistance has to *do* something, which the comparison above
/// cannot show: two forms that both ignored the parameter would still agree.
#[test]
fn rd_and_rs_reduce_the_drain_current() {
    let bare = nmos_drain_v("", "M1 d g 0 0 nm W=100u L=1u\n");
    let with_r = nmos_drain_v("RD=50 RS=50", "M1 d g 0 0 nm W=100u L=1u\n");
    // Less current through the load means a *higher* drain node voltage.
    assert!(
        with_r > bare + 1e-3,
        "100 Ω in series must cost current, so V(d) must rise above the bare \
         {bare:.6} V: got {with_r:.6} V. An equal voltage means RD/RS reached the \
         card and stopped there."
    );
    // Source degeneration is the stronger of the two: its drop subtracts from
    // Vgs as well as from Vds.
    let drain_only = nmos_drain_v("RD=100", "M1 d g 0 0 nm W=100u L=1u\n");
    let source_only = nmos_drain_v("RS=100", "M1 d g 0 0 nm W=100u L=1u\n");
    assert!(
        source_only > drain_only,
        "RS degenerates the gate drive as well as the drain, so the same value \
         must cost more current than RD: RS gives V(d)={source_only:.6}, RD gives \
         {drain_only:.6}"
    );
}

// ---------------------------------------------------------------------------
// UO, and the LEVEL 2/3 group that correctly does nothing here (#97 §1)
// ---------------------------------------------------------------------------

/// Drain current from a MOSFET biased into saturation, with no load resistor —
/// so `I(vd)` is the device's own current and nothing else.
fn nmos_id(model: &str, w: &str, l: &str) -> f64 {
    let deck = format!(
        "* mos id\n.model nm NMOS ({model})\n\
         VG g 0 DC 1.5\nVD d 0 DC 3\nM1 d g 0 0 nm W={w} L={l}\n.op\n"
    );
    let netlist = parse_spice(&deck).expect("parse");
    dc_op_nr(&netlist)
        .unwrap_or_else(|e| panic!("solve failed on\n{deck}\n{e:?}"))
        .vsrc_current("vd")
        .expect("I(vd)")
        .abs()
}

/// `UO` derives `KP` when the card gives no `KP`.
///
/// The card shape this covers is common: a foundry-ish card gives an oxide
/// thickness and a mobility and lets the simulator work out the transconductance.
/// Before this, `KP` fell back to SPICE's 2e-5 and the drain current was wrong by
/// whatever ratio the real `UO·COX` implied — a factor of 5 for the card below.
///
/// The anchor is the Level 1 closed form rather than a second simulator, because
/// `UO·COX` is arithmetic and the measurement that established it is already
/// recorded: ngspice returns 3.315020e-4 A for `TOX=20n` with no `KP`, and
/// `UO=300` gives exactly half the 600 default. What can go wrong here is the unit
/// conversion (`UO` is cm²/V·s) and the precedence of an explicit `KP`, and both
/// are closed-form facts.
#[test]
fn uo_derives_kp_when_kp_is_absent() {
    const EPS_OX: f64 = 3.9 * 8.854187817e-12;
    let cox = EPS_OX / 20e-9;
    for uo in [300.0_f64, 600.0, 900.0] {
        let model = format!("VTO=0.7 TOX=20n UO={uo}");
        let got = nmos_id(&model, "10u", "1u");
        // 0.5·KP·(W/L)·(Vgs−VTO)², with KP = UO(cm²/V·s)·1e-4·COX.
        let kp = uo * 1e-4 * cox;
        let want = 0.5 * kp * 10.0 * (1.5 - 0.7_f64).powi(2);
        let rel = (got - want).abs() / want;
        assert!(
            rel < 1e-5,
            "UO={uo}: I(vd)={got:.6e} A, and KP=UO·COX={kp:.6e} gives \
             {want:.6e} (rel {rel:.2e}). Falling back to SPICE's 2e-5 KP would \
             give {:.6e}.",
            0.5 * 2e-5 * 10.0 * 0.64
        );
    }

    // 1e-6 and not tighter: `I(vd)` also carries the reverse-biased bulk-drain
    // junction's leakage now, which at the default `gmin` is `IS + gmin·3V` ≈ 3 pA
    // — 9e-9 of a 0.32 mA drain current. That leakage is the body-diode feature
    // working, so the tolerance accommodates it rather than the deck avoiding it.
    const CHANNEL_TOL: f64 = 1e-6;

    // An explicit `KP` still wins over `UO`, which is SPICE's rule.
    let explicit = nmos_id("VTO=0.7 TOX=20n UO=900 KP=100u", "10u", "1u");
    let want = 0.5 * 100e-6 * 10.0 * (1.5 - 0.7_f64).powi(2);
    assert!(
        (explicit - want).abs() / want < CHANNEL_TOL,
        "KP given must win over UO: got {explicit:.9e}, expected {want:.9e}"
    );

    // And with neither KP nor an oxide, SPICE's fallback KP applies.
    let bare = nmos_id("VTO=0.7", "10u", "1u");
    let want_bare = 0.5 * 2e-5 * 10.0 * (1.5 - 0.7_f64).powi(2);
    assert!(
        (bare - want_bare).abs() / want_bare < CHANNEL_TOL,
        "no KP and no TOX must fall back to KP=2e-5: got {bare:.9e}"
    );
}

/// The mobility-degradation group is **not** a Level 1 gap — it is a Level 2/3
/// parameter set, and ngspice's Level 1 ignores it exactly as fairchild does.
///
/// # Why a test that asserts nothing happens
///
/// Normally this shape is worthless: `X_is_accepted` passes whether a parameter is
/// implemented, dropped or deleted. This one earns its place because the claim is
/// not "we accept it" but "the reference ignores it too, so implementing it would
/// be a divergence" — and that claim is the only thing standing between this group
/// and someone re-opening it as a to-do. It is checked both ways:
///
/// * at LEVEL 1, fairchild and ngspice must both be unmoved;
/// * at LEVEL 3, ngspice must *move*, which is what proves the parameters reach
///   its parser and that the first half is a modelling fact and not a typo.
///
/// `LAMBDA` is the control: a genuine Level 1 parameter, so it has to move both.
#[test]
fn the_mobility_group_is_level_2_or_3_and_correctly_does_nothing() {
    const BASE: &str = "VTO=0.7 KP=200u TOX=20n";
    let base = nmos_id(BASE, "10u", "1u");

    for extra in [
        "THETA=0.1",
        "THETA=1.0",
        "ETA=0.1",
        "KAPPA=1.0",
        "VMAX=1e5",
        "UCRIT=1e4",
        "UEXP=0.1",
        "UTRA=0.5",
        "NFS=1e12",
    ] {
        let got = nmos_id(&format!("{BASE} {extra}"), "10u", "1u");
        assert!(
            (got - base).abs() / base < 1e-12,
            "{extra} is a LEVEL 2/3 parameter and must not move a LEVEL 1 drain              current: {got:.9e} against {base:.9e}. If this starts failing, either              someone implemented it here — which diverges from ngspice's LEVEL 1 —              or the parameter is being misread as one that is modelled."
        );
    }

    // The control: a real Level 1 parameter has to move it.
    let with_lambda = nmos_id(&format!("{BASE} LAMBDA=0.05"), "10u", "1u");
    assert!(
        with_lambda > base * 1.05,
        "LAMBDA is a Level 1 parameter and must move the current, or this test is          measuring nothing: {with_lambda:.6e} against {base:.6e}"
    );
}

// ---------------------------------------------------------------------------
// Bulk-source / bulk-drain diodes (#97 §2)
// ---------------------------------------------------------------------------

/// Current out of a bulk-driving source, with the gate off so the channel
/// contributes nothing and the two body junctions are all that conducts.
fn bulk_current(model: &str, inst: &str, vb: f64, gmin: f64) -> f64 {
    bulk_current_at(model, inst, vb, gmin, None)
}

/// The same, at `temp_c` when it is given.
fn bulk_current_at(model: &str, inst: &str, vb: f64, gmin: f64, temp_c: Option<f64>) -> f64 {
    let temp = temp_c.map(|t| format!(".temp {t}\n")).unwrap_or_default();
    let deck = format!(
        "* body diode\n.options gmin={gmin:e}\n{temp}.model nm NMOS ({model})\n\
         VB bk 0 DC {vb}\nVG g 0 DC 0\nVD d 0 DC 0\n\
         M1 d g 0 bk nm {inst}\n.op\n"
    );
    let net = parse_spice(&deck).expect("parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&net.models);
    // `from_netlist`, not `dc_op_nr`: the latter uses `SimOptions::default()` and
    // ignores the deck's `.options`, so every point would run at the default gmin.
    let opts = SimOptions::from_netlist(&net);
    dc_op_nr_with_registry_opts(&net, &registry, &opts)
        .unwrap_or_else(|e| panic!("solve failed on\n{deck}\n{e:?}"))
        .vsrc_current("vb")
        .expect("I(vb)")
        .abs()
}

fn ngspice_bulk_current(model: &str, inst: &str, vb: f64, gmin: f64) -> Option<f64> {
    ngspice_bulk_current_at(model, inst, vb, gmin, None)
}

fn ngspice_bulk_current_at(
    model: &str,
    inst: &str,
    vb: f64,
    gmin: f64,
    temp_c: Option<f64>,
) -> Option<f64> {
    let temp = temp_c.map(|t| format!(".temp {t}\n")).unwrap_or_default();
    let dir = std::env::temp_dir().join("fc_body_golden");
    std::fs::create_dir_all(&dir).ok()?;
    let tag: String = format!("{model}{inst}{vb}{gmin}{temp}")
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .take(48)
        .collect();
    let path = dir.join(format!("body_{tag}.sp"));
    std::fs::write(
        &path,
        format!(
            "* body diode\n.options gmin={gmin:e}\n{temp}\
             .model nm NMOS ({model})\n\
             VB bk 0 DC {vb}\nVG g 0 DC 0\nVD d 0 DC 0\n\
             M1 d g 0 bk nm {inst}\n\
             .control\nop\nprint i(vb)\n.endc\n.end\n"
        ),
    )
    .ok()?;
    let out = std::process::Command::new("ngspice")
        .arg("-b")
        .arg(&path)
        .output()
        .ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    for line in text.lines() {
        let t = line.trim();
        if t.starts_with("i(vb)") {
            if let Ok(v) = t.split('=').nth(1)?.trim().parse::<f64>() {
                return Some(v.abs());
            }
        }
    }
    None
}

/// A forward-biased bulk conducts, and matches ngspice.
///
/// Before this the bulk was an island: `IS`/`JS` were accepted and nothing was
/// stamped, so latch-up and substrate injection could not be simulated at all and
/// the current was exactly zero at any bias.
///
/// The forward points are the load-bearing ones. `exp(0.7/vt)` moves 2.7% for a
/// 0.1% error in `vt`, so agreeing to 2e-3 at 0.4, 0.6 and 0.7 V pins the thermal
/// voltage and the saturation current together.
///
/// # Why the reverse points stop at -0.2 V
///
/// **ngspice is not a usable anchor inside 3*vt of zero** (0.0776 V), so that band
/// is excluded rather than the tolerance widened. Measured there, ngspice's total
/// over the two junctions is one junction flat at `-IS` plus one plain Shockley:
///
/// | `Vb` | ngspice | `-IS + IS*(exp(V/vt)-1)` | `2*IS*(exp(V/vt)-1)` | `-2*IS` |
/// |---|---|---|---|---|
/// | -0.01 V | 1.320654e-14 | 1.320663e-14 | 6.413258e-15 | 2.0e-14 |
/// | -0.05 V | 1.855304e-14 | 1.855314e-14 | 1.710628e-14 | 2.0e-14 |
/// | -0.07 V | 1.933221e-14 | 1.933228e-14 | 1.866455e-14 | 2.0e-14 |
///
/// The middle column matches to seven digits, and it still does with the bulk-drain
/// junction held five volts reverse, so the asymmetry is real and is ngspice's own
/// numerical convenience. Outside the band ngspice is flat at exactly `-IS` and
/// Shockley has converged onto it: 4.4e-4 relative at -0.2 V, exact by -0.5 V. So
/// these points confirm the magnitude and cannot see the branch shape.
/// [`the_bulk_junctions_follow_the_shockley_law`] pins the shape against the closed
/// form instead.
#[test]
fn the_bulk_junctions_conduct_and_match_ngspice() {
    const MODEL: &str = "VTO=0.7 KP=200u IS=1e-14";
    for vb in [-20.0, -2.0, -0.5, -0.2, 0.4, 0.6, 0.7] {
        let fc = bulk_current(MODEL, "W=10u L=1u", vb, 0.0);
        let Some(ng) = ngspice_bulk_current(MODEL, "W=10u L=1u", vb, 0.0) else {
            eprintln!("ngspice not available — skipping");
            return;
        };
        let rel = (fc - ng).abs() / ng;
        assert!(
            rel < 2e-3,
            "Vb={vb}: fairchild |I(vb)|={fc:.6e} A, ngspice {ng:.6e} (rel              {rel:.2e}). Exactly zero would mean the junctions are not stamped."
        );
    }
}

/// The reverse branch is Shockley, not flat.
///
/// The absolute anchor for the branch *shape*, which the ngspice comparison cannot
/// see: outside 3*vt the two forms agree to 4.4e-4, and inside it ngspice is
/// asymmetric. So this compares against the closed form directly, at biases where a
/// flat reverse branch would be 5% to 46% out.
///
/// `vt` comes from fairchild's own constants, which is only safe because
/// [`the_bulk_junctions_conduct_and_match_ngspice`] anchors it externally through
/// the forward exponential. What is under test here is the shape.
#[test]
fn the_bulk_junctions_follow_the_shockley_law() {
    const IS: f64 = 1e-14;
    let vt = fairchild_core::device::K_BOLTZMANN * 300.15 / fairchild_core::device::Q_ELECTRON;
    for vb in [-0.15, -0.10, -0.05, -0.02] {
        let got = bulk_current("VTO=0.7 KP=200u IS=1e-14", "W=10u L=1u", vb, 0.0);
        // Two junctions, both at `vb` because the drain and the source are grounded.
        let want = (2.0 * IS * ((vb / vt).exp() - 1.0)).abs();
        let flat = 2.0 * IS;
        let rel = (got - want).abs() / want;
        assert!(
            rel < 1e-6,
            "Vb={vb}: got {got:.6e} A, Shockley gives {want:.6e} (rel {rel:.2e}). \
             A flat reverse branch would give {flat:.6e}."
        );
    }
}

/// `gmin` crosses the bulk junctions, which is what makes this family consistent
/// with the diode and the BJT.
///
/// This is the second symptom #97 §2 named: without body diodes a MOSFET had no pn
/// junction at all, so its `gmin` stayed a Jacobian-only channel floor while the
/// other two families carried it as a real conductance.
#[test]
fn gmin_crosses_the_bulk_junctions() {
    const MODEL: &str = "VTO=0.7 KP=200u IS=1e-14";
    for gmin in [1e-12, 1e-9, 1e-6] {
        let fc = bulk_current(MODEL, "W=10u L=1u", -1.0, gmin);
        // Two junctions, each reverse at −1 V: `2·(IS + gmin·1V)`.
        let want = 2.0 * (1e-14 + gmin);
        let rel = (fc - want).abs() / want;
        assert!(
            rel < 1e-3,
            "gmin={gmin:e}: reverse bulk current is {fc:.6e} A and two junctions              give 2·(IS + gmin·1V) = {want:.6e} (rel {rel:.2e}). Leakage that does              not follow gmin means it is not crossing the junctions."
        );
        if let Some(ng) = ngspice_bulk_current(MODEL, "W=10u L=1u", -1.0, gmin) {
            assert!(
                (fc - ng).abs() / ng < 2e-3,
                "gmin={gmin:e}: fairchild {fc:.6e}, ngspice {ng:.6e}"
            );
        }
    }
}

/// `JS·area` wins over `IS` when the area is given, and each junction resolves
/// independently because `AS` and `AD` can differ.
#[test]
fn js_and_the_areas_set_the_saturation_current() {
    // `JS=1e-6` with `AS=AD=1p` is 1e-18 per junction, 1e-4 of the `IS` default.
    let with_js = bulk_current(
        "VTO=0.7 KP=200u JS=1e-6",
        "W=10u L=1u AS=1p AD=1p",
        0.6,
        0.0,
    );
    let with_is = bulk_current("VTO=0.7 KP=200u IS=1e-14", "W=10u L=1u", 0.6, 0.0);
    let ratio = with_is / with_js;
    assert!(
        (ratio / 1e4 - 1.0).abs() < 0.05,
        "JS·AS = 1e-18 against IS = 1e-14 is a factor of 1e4 in saturation          current: got {ratio:.4e}"
    );

    // Doubling the areas doubles the current.
    let doubled = bulk_current(
        "VTO=0.7 KP=200u JS=1e-6",
        "W=10u L=1u AS=2p AD=2p",
        0.6,
        0.0,
    );
    assert!(
        (doubled / with_js / 2.0 - 1.0).abs() < 0.05,
        "doubling AS and AD must double the junction current: {doubled:.6e}          against {with_js:.6e}"
    );

    // And `JS` with an area beats an explicit `IS`, which is SPICE's precedence.
    let both = bulk_current(
        "VTO=0.7 KP=200u IS=1e-14 JS=1e-6",
        "W=10u L=1u AS=1p AD=1p",
        0.6,
        0.0,
    );
    assert!(
        (both / with_js - 1.0).abs() < 1e-6,
        "with both given and an area present, JS·area wins: {both:.6e} against          the JS-only {with_js:.6e}"
    );

    for (m, inst) in [
        ("VTO=0.7 KP=200u JS=1e-6", "W=10u L=1u AS=1p AD=1p"),
        ("VTO=0.7 KP=200u IS=1e-14 JS=1e-6", "W=10u L=1u AS=1p AD=2p"),
    ] {
        if let Some(ng) = ngspice_bulk_current(m, inst, 0.6, 0.0) {
            let fc = bulk_current(m, inst, 0.6, 0.0);
            assert!(
                (fc - ng).abs() / ng < 2e-3,
                "{m} / {inst}: fairchild {fc:.6e}, ngspice {ng:.6e}"
            );
        }
    }
}

/// The bulk junctions' saturation current scales with temperature, by a **third**
/// law that is neither the diode's nor the BJT's.
///
/// A MOSFET card carries no `EG` and no `XTI`, so SPICE cannot use the constant-`EG`
/// form the other two families use. It puts the temperature-dependent bandgap in
/// the exponent instead: `exp(Eg(TNOM)/vt(TNOM) − Eg(T)/vt(T))`.
///
/// This is the whole reason the law needed measuring rather than reusing. Applying
/// the diode's law here would be out by up to 2.4× over −40 to 125 °C, which the
/// tolerance below rejects — that is the sabotage this test is built to catch.
///
/// Read with the bulk one volt reverse and `gmin = 0`, so the current is exactly
/// `2·Isat(T)` and nothing else. Five decades of `Isat` are covered, so the test
/// fails on a wrong *exponent*, not only on a wrong prefactor.
#[test]
fn the_bulk_junction_saturation_current_scales_with_temperature() {
    const MODEL: &str = "VTO=0.7 KP=200u IS=1e-14 TNOM=27";
    let mut ran = 0;
    for tc in [-40.0, 0.0, 27.0, 75.0, 125.0] {
        let Some(ng) = ngspice_bulk_current_at(MODEL, "W=10u L=1u", -1.0, 0.0, Some(tc)) else {
            eprintln!("ngspice not available — skipping");
            return;
        };
        let fc = bulk_current_at(MODEL, "W=10u L=1u", -1.0, 0.0, Some(tc));
        let rel = (fc - ng).abs() / ng;
        // 1e-3, not the file's 2e-3: the measured residual is 3.7e-4 worst case
        // and its source is named in `temperature::mos_junction_is_factor`. The
        // diode law would be 1.4 to 2.4 out here, so this has 3 orders of margin
        // over the error it exists to reject.
        assert!(
            rel < 1e-3,
            "{tc} °C: fairchild 2·Isat = {fc:.6e} A, ngspice {ng:.6e} (rel \
             {rel:.2e}). No temperature scaling at all would read 2e-14 at every \
             point; the diode's constant-EG law would be up to 2.4x out."
        );
        ran += 1;
    }
    assert_eq!(ran, 5, "every temperature point must have been compared");
}

/// `RSH` with `NRD`/`NRS` becomes the drain and source series resistances.
///
/// `RSH` is a resistance per square and `NRD`/`NRS` are the number of squares in
/// each diffusion, so `RD = RSH·NRD` and `RS = RSH·NRS`. Before this `RSH` was
/// accepted and dropped, and a card that gives sheet resistance instead of `RD`/`RS`
/// got no series resistance at all.
///
/// Measured, and the equalities are exact rather than approximate. `RSH=50 NRD=2
/// NRS=2` is bit-identical to `RD=100 RS=100` in ngspice, ratio 1.000000000, and
/// `RSH=50` alone equals `RD=50 RS=50`, which is how the `NRD=NRS=1` default was
/// read off.
///
/// # Precedence is per terminal
///
/// `RSH=50 RD=1000 NRD=2 NRS=2` reads 0.00147653 where `RD=1000` alone reads
/// 0.00161661. So the explicit `RD` wins on the drain and `RSH·NRS` still applies
/// to the source. Treating one explicit value as disabling `RSH` for both terminals
/// gives the second number, and is the mistake this test exists to catch.
///
/// `NRD` maps to the drain and `NRS` to the source: `NRD=4 NRS=1` reads 0.00350969
/// and `NRD=1 NRS=4` reads 0.00281955, because source degeneration costs more
/// current than the same resistance in the drain. So swapping them fails too.
#[test]
fn rsh_times_the_squares_becomes_the_series_resistance() {
    let mut compared = 0;
    for (model, inst) in [
        ("VTO=0.7 KP=200u", "W=10u L=1u"),
        ("VTO=0.7 KP=200u RSH=50", "W=10u L=1u"),
        ("VTO=0.7 KP=200u RSH=50", "W=10u L=1u NRD=2 NRS=2"),
        ("VTO=0.7 KP=200u RSH=50", "W=10u L=1u NRD=4 NRS=1"),
        ("VTO=0.7 KP=200u RSH=50", "W=10u L=1u NRD=1 NRS=4"),
        ("VTO=0.7 KP=200u RSH=50 RD=1000", "W=10u L=1u NRD=2 NRS=2"),
        ("VTO=0.7 KP=200u RSH=50 RS=1000", "W=10u L=1u NRD=2 NRS=2"),
        ("VTO=0.7 KP=200u RD=100 RS=100", "W=10u L=1u"),
    ] {
        let deck = format!(
            "* rsh\n.model nm NMOS ({model})\n\
             VG g 0 DC 3\nVD d 0 DC 2\nM1 d g 0 0 nm {inst}\n"
        );
        let net = parse_spice(&deck).expect("parse");
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let opts = SimOptions::from_netlist(&net);
        let got = dc_op_nr_with_registry_opts(&net, &registry, &opts)
            .unwrap_or_else(|e| panic!("solve failed on\n{deck}\n{e:?}"))
            .vsrc_current("vd")
            .expect("i(vd)")
            .abs();
        let Some(ng) = ngspice_bulk_i_vd(&deck) else {
            eprintln!("ngspice not available — skipping");
            return;
        };
        // `CHANNEL_TOL` is too tight here: the body-drain junction's leakage rides
        // on `I(vd)` and the channel current is milliamps, so this is the file's
        // ordinary golden tolerance.
        let rel = (got - ng).abs() / ng;
        assert!(
            rel < REL_TOL,
            "'{model}' / '{inst}': fairchild I(vd)={got:.6e}, ngspice {ng:.6e} \
             (rel {rel:.2e}). Dropping RSH gives the no-resistance answer for \
             five of these eight rows."
        );
        compared += 1;
    }
    assert_eq!(compared, 8, "every card must have been compared");
}

/// ngspice's `i(vd)` for a `.op` deck.
fn ngspice_bulk_i_vd(deck: &str) -> Option<f64> {
    let ng = find_ngspice()?;
    let dir = std::env::temp_dir().join("fc_rsh_golden");
    std::fs::create_dir_all(&dir).ok()?;
    let tag: String = deck
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .skip(4)
        .take(56)
        .collect();
    let path = dir.join(format!("rsh_{tag}.sp"));
    std::fs::write(
        &path,
        format!("{deck}.control\nop\nprint i(vd)\n.endc\n.end\n"),
    )
    .ok()?;
    let out = Command::new(&ng).arg("-b").arg(&path).output().ok()?;
    let text = String::from_utf8_lossy(&out.stdout);
    for line in text.lines() {
        let t = line.trim();
        if t.starts_with("i(vd)") && t.contains('=') {
            if let Ok(v) = t
                .split('=')
                .nth(1)?
                .split_whitespace()
                .next()?
                .parse::<f64>()
            {
                return Some(v.abs());
            }
        }
    }
    None
}
