//! Every reactance a device stamps in transient must also reach `.ac`.
//!
//! # The bug this exists to prevent
//!
//! A device can carry its capacitance two ways. It can report
//! `Device::reactive_branches` and let the integrator build the companion, or it
//! can stamp the companion itself in `Device::load_jacobian_tran`. The second way
//! is invisible to `.ac` and `.noise`, which read
//! `Device::small_signal_reactances` and `Device::load_reactive_jacobian` instead.
//! A device that stamps its own companion and forgets the AC hooks is a transistor
//! with no capacitance in every frequency-domain analysis.
//!
//! That is not hypothetical. The BJT did exactly this. Measured before the fix: a
//! 1 kOhm resistor into the base of a device with `CJE = CJC = CJS = 100p` read
//! `|V(b)| = 1.000000` at 1 kHz, 1 MHz, 10 MHz and 100 MHz, where the corner is
//! 1.04 MHz. Flat, silently. Every transient test passed throughout, because
//! transient takes the other path, and no AC test asked.
//!
//! # Why this shape
//!
//! This is an absence rather than a weakness, so a test per device cannot find it.
//! The missing test looks like no test. What is needed is a structural invariant,
//! measured on every device rather than declared by it. Every matrix cell whose
//! value `load_jacobian_tran` changes must also be written by the AC hooks.
//!
//! Three weaker forms were tried and each let a real fault through.
//! *Reports something* is too weak, because a device that reports four of its five
//! capacitances passes it: dropping only the substrate capacitance from the BJT's
//! AC list walked straight through. *Cells the transient stamp adds and the DC
//! stamp does not* is blind for a two-terminal device, because the diode's junction
//! capacitance sits between the same node pair as its conductance and the cell set
//! never grows. *A count of reported branches* says nothing about which ones.
//! Per-cell values catch all three. `load_jacobian_tran` calls `load_jacobian`
//! first and then adds, so a cell with no reactive part is bit-identical and the
//! comparison needs no tolerance.
//!
//! # What this does not check
//!
//! The cell **values**. A device that reports the right node pairs with the wrong
//! capacitance passes this. The per-device ngspice goldens catch that instead:
//! `cjs_follows_the_depletion_law_in_reverse_and_a_straight_line_forward` for the
//! substrate junction, and the switching-time goldens for the rest. This gate
//! answers one question. Does the reactance reach the frequency domain at all.

use std::collections::BTreeSet;

use fairchild_core::device::{Device, EvalFlags, ReactiveKind};
use fairchild_core::mna::CircuitTopology;
use fairchild_core::newton::build_devices;
use fairchild_core::options::SimOptions;
use fairchild_core::{dc_op_nr_with_registry_opts, DeviceRegistry};
use fairchild_parser::parse_spice;

/// One deck per native device that can carry a reactance, each biased so the
/// reactance is non-zero at the operating point.
///
/// Hand-maintained, because Rust cannot enumerate a trait's implementors at
/// runtime. The count assertion at the end is what makes adding a device a
/// deliberate edit here rather than a silent gap.
const EVERY_REACTIVE_DEVICE: &[(&str, &str)] = &[
    (
        "diode",
        "* diode Cj\n\
         .model dm D (IS=1e-14 CJO=10p VJ=0.75 M=0.5 TT=1n)\n\
         V1 a 0 DC 0.6\n\
         D1 a 0 dm\n\
         .op\n",
    ),
    (
        "mosfet",
        "* mosfet Meyer and junction caps\n\
         .model nm NMOS (VTO=0.7 KP=200u CGSO=1n CGDO=1n CGBO=1n CJ=1m CJSW=1n)\n\
         VG g 0 DC 2\n\
         VD d 0 DC 2\n\
         M1 d g 0 0 nm W=10u L=1u AS=1p AD=1p PS=4u PD=4u\n\
         .op\n",
    ),
    (
        "bjt",
        "* bjt depletion, diffusion and substrate caps\n\
         .model qn NPN (IS=1e-16 BF=100 TF=1n TR=10n CJE=10p CJC=5p \
         CJS=2p VJS=0.75 MJS=0.33)\n\
         VC c 0 DC 5\n\
         VB b 0 DC 0.7\n\
         VS s 0 DC -2\n\
         Q1 c b 0 s qn\n\
         .op\n",
    ),
];

type Cells = BTreeSet<(usize, usize)>;

/// The cells whose value the transient stamp changes: this device's reactance,
/// measured rather than declared.
fn reactive_cells(dev: &dyn Device, n: usize, alpha: f64) -> Cells {
    let mut dc = fairchild_core::mna::MnaMatrix::zeros(n);
    dev.load_jacobian(&mut dc);
    let mut tran = fairchild_core::mna::MnaMatrix::zeros(n);
    dev.load_jacobian_tran(&mut tran, alpha);

    let mut dc_vals = std::collections::BTreeMap::new();
    for (r, row) in dc.a.iter().enumerate() {
        for (c, v) in row.iter() {
            dc_vals.insert((r, c), v);
        }
    }
    let mut out = Cells::new();
    for (r, row) in tran.a.iter().enumerate() {
        for (c, v) in row.iter() {
            if v != dc_vals.get(&(r, c)).copied().unwrap_or(0.0) {
                out.insert((r, c));
            }
        }
    }
    out
}

/// The cells the AC assembly would write for this device, from both hooks.
fn ac_cells(dev: &dyn Device, n: usize) -> Cells {
    let mut out = Cells::new();
    for r in dev.small_signal_reactances() {
        if r.value == 0.0 && matches!(r.kind, ReactiveKind::Inductor) {
            continue;
        }
        // `.ac` stamps a two-port between (pos, neg); a `None` terminal is ground
        // and contributes no row or column, exactly as `stamp_2port_by_id` does.
        for (a, b) in [
            (r.pos, r.pos),
            (r.pos, r.neg),
            (r.neg, r.pos),
            (r.neg, r.neg),
        ] {
            if let (Some(i), Some(j)) = (a, b) {
                out.insert((i, j));
            }
        }
    }
    let mut c_mat = vec![fairchild_core::mna::SparseRow::default(); n];
    dev.load_reactive_jacobian(&mut c_mat);
    for (r, row) in c_mat.iter().enumerate() {
        for (c, v) in row.iter() {
            if v != 0.0 {
                out.insert((r, c));
            }
        }
    }
    out
}

/// The gate.
#[test]
fn every_reactive_cell_a_device_stamps_in_transient_also_reaches_ac() {
    for (name, deck) in EVERY_REACTIVE_DEVICE {
        let net = parse_spice(deck).unwrap_or_else(|e| panic!("{name}: parse: {e:?}"));
        let mut registry = DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let opts = SimOptions::from_netlist(&net);
        let ctx = opts.sim_context();
        let mut topo = CircuitTopology::build(&net);
        let mut devices = build_devices(&net, &mut topo, &ctx, &registry)
            .unwrap_or_else(|e| panic!("{name}: build: {e:?}"));

        // Solve first, so the devices are interrogated at a real operating point
        // rather than at zero, where several capacitances legitimately vanish.
        let sol = dc_op_nr_with_registry_opts(&net, &registry, &opts)
            .unwrap_or_else(|e| panic!("{name}: solve: {e:?}"));
        let n = topo.size;
        let mut x = vec![0.0; n];
        for (node, &row) in topo.node_index.iter() {
            if let Ok(v) = sol.node_voltage(node) {
                x[row] = v;
            }
        }

        let mut saw_reactive = false;
        for dev in devices.iter_mut() {
            // `EvalFlags::tran` on purpose: the capacitance caches are
            // transient-gated in every model here, and `.ac` queries them after
            // exactly this eval.
            dev.eval(&x, EvalFlags::tran(), &ctx);

            let reactive = reactive_cells(dev.as_ref(), n, 1.0 / 1e-9);
            if reactive.is_empty() {
                continue;
            }
            saw_reactive = true;
            let ac = ac_cells(dev.as_ref(), n);
            let missing: Vec<_> = reactive.difference(&ac).copied().collect();
            assert!(
                missing.is_empty(),
                "{name}: `load_jacobian_tran` changes the value of cells \
                 {missing:?}, so they carry a reactance - and the AC hooks write \
                 none of them. `.ac` and `.noise` will run this device with \
                 that reactance missing and report a plausible wrong answer. \
                 Override `small_signal_reactances` (or `load_reactive_jacobian`) \
                 with the *same* reactances `load_jacobian_tran` stamps. This is \
                 how the BJT shipped with no capacitance in any AC sweep."
            );
        }
        assert!(
            saw_reactive,
            "{name}: no device in this deck stamped a reactive companion, so the \
             deck does not exercise what this test checks. Bias it so its \
             capacitance is non-zero, or remove the row."
        );
    }
    assert_eq!(
        EVERY_REACTIVE_DEVICE.len(),
        3,
        "adding a native device with a capacitance means adding a row above. Bump \
         this count in the same commit, so a new device cannot arrive without \
         someone deciding whether its reactance reaches `.ac`."
    );
}
