//! `$abstime` must read the transient clock, not zero.
//!
//! `OsdiSimInfo.abstime` was hardcoded to 0.0, so any Verilog-A model reading
//! `$abstime` — an arbitrary-waveform source, a drift or ageing term, anything
//! time-parametrised — silently saw t = 0 for the whole run.
//!
//! The old fixture recorded the value it was handed into a spare field of its
//! model struct, and the test read it back. That proved the number crossed the
//! ABI and nothing more. `abstime_ramp.va` puts it in the *answer* instead: it
//! drives `k · $abstime` amps into a resistor, so V(out) must be k·R·t at every
//! timepoint, and a dead clock is a flat line.

use fairchild_core::{tran_nr_with_registry, DeviceRegistry};
use fairchild_parser::parse_spice;

use crate::common;

#[test]
fn the_simulation_clock_reaches_the_model() {
    let Some(path) = common::compiled("abstime_ramp") else {
        return;
    };

    // k = 1 A/s into 1 kΩ: V(out) = 1000·t, which is 5 mV at the 5 µs stop.
    let deck = format!(
        "* $abstime through a resistor\n\
         .osdi {}\n\
         Xs out 0 abstime_ramp\n\
         Rl out 0 1k\n\
         .tran 1u 5u\n\
         .end\n",
        path.display()
    );
    let netlist = parse_spice(&deck).expect("parse");
    let mut registry = DeviceRegistry::new();
    registry.register_builtin_models(&netlist.models);
    fairchild_osdi::load_libraries(
        &netlist.osdi_paths,
        &netlist.va_sources,
        None,
        &Default::default(),
        &mut registry,
    )
    .expect("load");

    let r = tran_nr_with_registry(&netlist, 1e-6, 5e-6, &registry).expect("transient failed");

    let times = &r.time;
    let out = r.node_voltages.get("out").expect("node out");
    assert!(times.len() >= 5, "too few timepoints: {}", times.len());
    for (t, v) in times.iter().zip(out) {
        let want = 1000.0 * t;
        assert!(
            (v - want).abs() < 1e-9,
            "t = {t}: V(out) = {v}, want {want} — the model did not see the clock"
        );
    }
    // And it has to have gone somewhere: a clock stuck at zero would satisfy
    // the loop above only if every timepoint were zero too.
    assert!(
        out.last().copied().unwrap_or(0.0) > 1e-3,
        "no ramp at all: {out:?}"
    );
}
