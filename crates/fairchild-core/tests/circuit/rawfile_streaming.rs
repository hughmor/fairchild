//! Streaming a transient, and the binary rawfile spelling.
//!
//! Two claims, and each has an obvious way to pass while being wrong:
//!
//! * **Streaming and collecting agree.** They must produce the same bytes for
//!   the same run, or a long run's output is a second and subtly different
//!   opinion about what the results are.
//! * **Binary and ASCII agree.** Same header, same values. The binary path
//!   exists to stop paying for six formatted digits, not to answer differently.
//!
//! Both are agreement invariants, so on their own they could hide a fault
//! common to both sides. The absolute anchor is at the other end: `tests/
//! ngspice/ngspice_rawfile_golden.rs` compares our binary rawfile with the one
//! ngspice writes for the same run, byte layout included. What is here is the
//! internal consistency that anchor cannot see.

use std::io::Cursor;

use fairchild_core::nutmeg::{self, Encoding};
use fairchild_core::options::SimOptions;
use fairchild_core::tran::{tran_nr_with_registry_opts_into, tran_nr_with_registry_var_opts_into};
use fairchild_core::{
    tran_nr_configured, CollectSink, CsvSink, DeviceRegistry, RawSink, SelectSink, TranSink,
};
use fairchild_parser::{parse_spice, Netlist};

/// An RC step with a source current to read back, so the layout carries a node
/// voltage, a branch current and a `time` column.
const DECK: &str = "\
* rc step
V1 in 0 PULSE(0 1 0 1n 1n 1u 2u)
R1 in out 1k
C1 out 0 1p
.tran 100p 3n
.end
";

fn netlist(extra: &str) -> Netlist {
    parse_spice(&DECK.replace(".end", &format!("{extra}\n.end"))).unwrap()
}

fn registry(net: &Netlist) -> DeviceRegistry {
    let mut r = DeviceRegistry::new();
    r.register_builtin_models(&net.models);
    r
}

/// Run the deck into `sink`, on whichever integrator `opts` selects.
fn run(net: &Netlist, opts: &SimOptions, sink: &mut dyn TranSink) {
    let reg = registry(net);
    let (step, stop) = (100e-12, 3e-9);
    if opts.variable_step {
        tran_nr_with_registry_var_opts_into(net, step, stop, &reg, opts, sink).unwrap();
    } else {
        tran_nr_with_registry_opts_into(net, step, stop, &reg, opts, sink).unwrap();
    }
}

fn streamed(net: &Netlist, opts: &SimOptions, enc: Encoding) -> Vec<u8> {
    let mut sink = RawSink::new(Cursor::new(Vec::new()), "t", enc);
    run(net, opts, &mut sink);
    sink.into_inner().unwrap().into_inner()
}

fn collected(net: &Netlist, opts: &SimOptions, enc: Encoding) -> Vec<u8> {
    let reg = registry(net);
    let r = tran_nr_configured(net, 100e-12, 3e-9, &reg, opts).unwrap();
    let mut buf = Vec::new();
    r.write_raw(&mut buf, "t", enc).unwrap();
    buf
}

/// Only the reserved point-count field's padding may differ: a streamed run
/// does not know its count when it writes the header.
fn same_but_for_padding(a: &[u8], b: &[u8]) -> Result<(), String> {
    let norm = |v: &[u8]| {
        let f = nutmeg::read(v).map_err(|e| format!("unreadable: {e}"))?;
        Ok::<_, String>((f.plotname.clone(), f.vars.len(), f.points.clone()))
    };
    let (pa, va, xa) = norm(a)?;
    let (pb, vb, xb) = norm(b)?;
    if pa != pb || va != vb {
        return Err(format!("headers differ: {pa}/{va} against {pb}/{vb}"));
    }
    if xa.len() != xb.len() {
        return Err(format!("{} points against {}", xa.len(), xb.len()));
    }
    for (i, (ra, rb)) in xa.iter().zip(&xb).enumerate() {
        if ra != rb {
            return Err(format!("point {i}: {ra:?} against {rb:?}"));
        }
    }
    Ok(())
}

/// **Both integrators**, streamed against collected, in both spellings.
///
/// The variable-step path is the one that would be left behind: its point count
/// is not known when the header is written, so it is the only caller that needs
/// the reserved field, and the only one where a wrong count would go unnoticed
/// until a reader tried to use the file.
#[test]
fn streaming_and_collecting_agree_on_both_integrators() {
    for variable in [false, true] {
        let net = netlist("");
        let mut opts = SimOptions::from_netlist(&net);
        opts.variable_step = variable;
        let mode = if variable {
            "variable-step"
        } else {
            "fixed-step"
        };
        for enc in [Encoding::Ascii, Encoding::Binary] {
            let s = streamed(&net, &opts, enc);
            let c = collected(&net, &opts, enc);
            same_but_for_padding(&s, &c)
                .unwrap_or_else(|e| panic!("{mode} {enc:?}: streamed against collected: {e}"));
        }
    }
}

/// The header must state the number of points the body actually carries.
///
/// This is the whole risk of a reserved field: a run that never fills it in
/// leaves a file every reader will mis-slice. `nutmeg::read` checks the count
/// against the body rather than trusting it, so an unfilled field is an error
/// here rather than a silent short read — which is what makes this test able to
/// fail at all.
#[test]
fn a_streamed_rawfile_states_its_real_point_count() {
    for variable in [false, true] {
        let net = netlist("");
        let mut opts = SimOptions::from_netlist(&net);
        opts.variable_step = variable;
        for enc in [Encoding::Ascii, Encoding::Binary] {
            let bytes = streamed(&net, &opts, enc);
            let f = nutmeg::read(&bytes).expect("streamed rawfile must read back");
            assert!(
                f.points.len() > 10,
                "{enc:?}: only {} points",
                f.points.len()
            );
            let stated: usize = String::from_utf8_lossy(&bytes[..400.min(bytes.len())])
                .lines()
                .find_map(|l| l.strip_prefix("No. Points:"))
                .and_then(|v| v.trim().parse().ok())
                .expect("header must carry a parsable point count");
            assert_eq!(stated, f.points.len());
        }
    }
}

/// Binary must be exact where ASCII rounds.
///
/// Agreement to six figures is what the two formats share; the extra ten are
/// the reason to have the binary one at all. Asserting only agreement would
/// pass on a binary writer that quietly formatted and re-parsed.
#[test]
fn binary_carries_more_than_six_figures_and_ascii_does_not() {
    let net = netlist("");
    let opts = SimOptions::from_netlist(&net);
    let a = nutmeg::read(&streamed(&net, &opts, Encoding::Ascii)).unwrap();
    let b = nutmeg::read(&streamed(&net, &opts, Encoding::Binary)).unwrap();

    let va = a.series("v(out)").unwrap();
    let vb = b.series("v(out)").unwrap();
    assert_eq!(va.len(), vb.len());

    let mut ascii_lost = 0;
    for (x, y) in va.iter().zip(&vb) {
        assert!(
            (x - y).abs() <= 1e-6 * y.abs().max(1e-30),
            "spellings disagree beyond ASCII's six figures: {x} against {y}"
        );
        if x.to_bits() != y.to_bits() {
            ascii_lost += 1;
        }
    }
    assert!(
        ascii_lost > va.len() / 2,
        "ASCII kept {}/{} values bit-exact, so this deck cannot tell the two \
         spellings apart and the test proves nothing",
        va.len() - ascii_lost,
        va.len()
    );
}

/// A selected layout must narrow what is *written*, not what is displayed.
///
/// The bug this replaces filtered the rendered text, so asking for one signal
/// of two hundred formatted all two hundred and then discarded 199 — slower
/// than not filtering, since the filter re-parsed what it had just written.
#[test]
fn probe_selection_narrows_the_file_itself() {
    let net = netlist("");
    let opts = SimOptions::from_netlist(&net);

    let mut all = CsvSink::new(Cursor::new(Vec::new()));
    run(&net, &opts, &mut all);

    let mut one = CollectSink::new();
    {
        let mut sel = SelectSink::new(&mut one, &["V(out)".to_string()]);
        run(&net, &opts, &mut sel);
    }
    let r = one.take().unwrap();
    assert_eq!(
        r.node_voltages.len(),
        1,
        "only V(out) should have been kept"
    );
    assert!(r.node_voltages.contains_key("out"));
    assert!(
        r.vsrc_currents.is_empty(),
        "I(v1) was not asked for and must not be carried"
    );
    assert!(!r.time.is_empty());
}

/// A probe that names nothing must stop the run, not narrow it silently (#72).
#[test]
fn an_unmatched_probe_is_an_error() {
    let net = netlist("");
    let opts = SimOptions::from_netlist(&net);
    let reg = registry(&net);
    let mut sink = CollectSink::new();
    let mut sel = SelectSink::new(&mut sink, &["V(nowhere)".to_string()]);
    let err = tran_nr_with_registry_opts_into(&net, 100e-12, 3e-9, &reg, &opts, &mut sel)
        .expect_err("an unmatched probe must fail the run");
    let msg = err.to_string();
    assert!(msg.contains("V(nowhere)"), "{msg}");
}

/// tstart selects what is saved, on both integrators and through a writer.
///
/// The streaming filter has to hold a point it has already decided to drop, in
/// case it turns out to be the last one — an empty waveform is a worse answer
/// than a short one, and that is the case a naive filter gets wrong.
#[test]
fn tstart_reaches_a_streamed_rawfile() {
    for variable in [false, true] {
        let net = netlist(".options tstart=2n");
        let mut opts = SimOptions::from_netlist(&net);
        opts.variable_step = variable;
        let f = nutmeg::read(&streamed(&net, &opts, Encoding::Binary)).unwrap();
        let t = f.series("time").unwrap();
        assert!(!t.is_empty(), "tstart must not empty the run");
        assert!(
            t[0] >= 2e-9 - 1e-15,
            "first saved point is {}, before tstart",
            t[0]
        );

        // Past the end of the run: one point survives.
        let net = netlist(".options tstart=100");
        let mut opts = SimOptions::from_netlist(&net);
        opts.variable_step = variable;
        let f = nutmeg::read(&streamed(&net, &opts, Encoding::Binary)).unwrap();
        assert_eq!(
            f.points.len(),
            1,
            "a tstart past the end must leave exactly the last point"
        );
    }
}
