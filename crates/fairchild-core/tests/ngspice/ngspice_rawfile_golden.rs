//! The binary rawfile's layout, against the simulator that defines it.
//!
//! Everything else about the binary spelling is checked by comparing it with
//! our own ASCII spelling, which cannot see a fault the two share — a byte
//! order, a value stride, a missing separator, or a header line a real reader
//! needs and neither of ours writes. So the anchor is ngspice: it writes the
//! same file for the same run, and it reads ours back.
//!
//! The format, as measured rather than as remembered: header lines identical
//! between the two spellings, `Binary:` in place of `Values:`, then native-
//! endian `f64` immediately after that newline, point-major, no padding and no
//! separators. Checked arithmetically on a 62-point, 4-variable run from
//! ngspice 46: 62 × 4 × 8 = 1984 bytes of body, 2218-byte file, 234-byte
//! header.
//!
//! Skips when ngspice is not installed, like every other golden here.

use std::io::Write;
use std::process::Command;

use fairchild_core::nutmeg::{self, Encoding};
use fairchild_core::options::SimOptions;
use fairchild_core::{tran_nr_configured, DeviceRegistry};
use fairchild_parser::parse_spice;

use super::ngspice_golden::find_ngspice;

const DECK: &str = "\
* rc step
V1 in 0 PULSE(0 1 0 1n 1n 1u 2u)
R1 in out 1k
C1 out 0 1p
.tran 100p 5n
";

/// Have ngspice run `DECK` and write it as a binary rawfile.
fn ngspice_binary_raw(dir: &std::path::Path) -> Option<Vec<u8>> {
    let bin = find_ngspice()?;
    let out = dir.join("ng.raw");
    let deck = dir.join("ng.sp");
    let mut f = std::fs::File::create(&deck).ok()?;
    write!(
        f,
        "{DECK}.control\nset filetype=binary\nrun\nwrite {}\nquit\n.endc\n.end\n",
        out.display()
    )
    .ok()?;
    drop(f);
    let st = Command::new(bin).arg("-b").arg(&deck).output().ok()?;
    if !out.exists() {
        eprintln!(
            "ngspice wrote no rawfile: {}",
            String::from_utf8_lossy(&st.stderr)
        );
        return None;
    }
    std::fs::read(&out).ok()
}

fn ours(enc: Encoding) -> Vec<u8> {
    let net = parse_spice(&format!("{DECK}.end\n")).unwrap();
    let mut reg = DeviceRegistry::new();
    reg.register_builtin_models(&net.models);
    let opts = SimOptions::from_netlist(&net);
    let r = tran_nr_configured(&net, 100e-12, 5e-9, &reg, &opts).unwrap();
    let mut buf = Vec::new();
    r.write_raw(&mut buf, "rc step", enc).unwrap();
    buf
}

/// The number on ngspice's `print v(out)[n] = …` line.
///
/// Anchored on the `v(out)` prefix rather than on "the first line with an `=`":
/// ngspice's banner and plot listing carry `=` too, and reading one of those
/// made this test fail against a rawfile it had in fact read correctly.
fn printed(text: &str) -> Option<f64> {
    text.lines()
        .filter(|l| l.trim_start().starts_with("v(out)"))
        .find_map(|l| l.split('=').nth(1))
        .and_then(|v| v.trim().parse().ok())
}

/// The layout, measured on ngspice's own file rather than assumed.
///
/// This is the test that would catch a wrong stride, a text separator between
/// values, or little-endian written where the reader expects native order: it
/// derives the body length from the header ngspice wrote and requires the file
/// to be exactly that long.
#[test]
fn ngspice_binary_layout_is_header_then_packed_f64() {
    let dir = tempfile::tempdir().unwrap();
    let Some(bytes) = ngspice_binary_raw(dir.path()) else {
        eprintln!("skipping: ngspice not installed");
        return;
    };

    let tag = b"Binary:\n";
    let at = bytes
        .windows(tag.len())
        .position(|w| w == tag)
        .expect("ngspice binary rawfile must carry a `Binary:` line");
    let head = String::from_utf8_lossy(&bytes[..at]);
    let field = |k: &str| -> usize {
        head.lines()
            .find_map(|l| l.strip_prefix(k))
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or_else(|| panic!("no `{k}` in ngspice header:\n{head}"))
    };
    let n_vars = field("No. Variables:");
    let n_pts = field("No. Points:");
    let body = bytes.len() - at - tag.len();
    assert_eq!(
        body,
        n_pts * n_vars * 8,
        "ngspice body is {body} bytes for {n_pts} points × {n_vars} variables; \
         the format is packed f64 with no separators and no padding"
    );

    // And our reader agrees with ngspice's own file, which is the other half:
    // writing the layout correctly is no use if we cannot read one back.
    let f = nutmeg::read(&bytes).expect("our reader must read ngspice's binary rawfile");
    assert_eq!(f.points.len(), n_pts);
    assert_eq!(f.vars.len(), n_vars);
    let v = f
        .series("v(out)")
        .expect("ngspice's rawfile must carry v(out)");
    assert!(
        v.iter().any(|x| *x > 0.5),
        "v(out) never rises; the values were not decoded"
    );
}

/// ngspice reads *our* binary rawfile, and gets our numbers out of it.
///
/// The direction that matters for a user: the file is worth writing because
/// their existing tools open it. This runs ngspice's own `load` on our output
/// and compares a value it prints back with the one we wrote — so a byte order
/// or stride we got wrong shows up as a wrong number rather than as a file that
/// happens to parse.
#[test]
fn ngspice_reads_our_binary_rawfile() {
    let Some(bin) = find_ngspice() else {
        eprintln!("skipping: ngspice not installed");
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let raw = dir.path().join("ours.raw");
    std::fs::write(&raw, ours(Encoding::Binary)).unwrap();

    let mine = nutmeg::read(&std::fs::read(&raw).unwrap()).unwrap();
    let v_out = mine.series("v(out)").unwrap();
    let probe = v_out.len() / 2;
    let expect = v_out[probe];

    // `-b` needs a circuit of its own before a `.control` block will run, so
    // the deck is a throwaway and the rawfile is loaded over the top of it.
    let deck = dir.path().join("read.sp");
    std::fs::write(
        &deck,
        format!(
            "* reader\nR1 a 0 1k\nV1 a 0 1\n.op\n.control\nrun\nload {}\n\
             setplot tran1\nprint v(out)[{probe}]\nquit\n.endc\n.end\n",
            raw.display()
        ),
    )
    .unwrap();
    let out = Command::new(bin).arg("-b").arg(&deck).output().unwrap();
    let text = String::from_utf8_lossy(&out.stdout);
    let got: f64 = printed(&text).unwrap_or_else(|| {
        panic!(
            "ngspice printed no value from our rawfile:\n{text}\n{}",
            String::from_utf8_lossy(&out.stderr)
        )
    });
    assert!(
        (got - expect).abs() <= 1e-6 * expect.abs().max(1e-12),
        "ngspice read {got} where we wrote {expect}"
    );
}

/// A padded point count reads the same as an unpadded one.
///
/// A streamed transient does not know its count when it writes the header, so
/// it reserves the field and fills it in afterwards, right-aligned. That is a
/// deviation from what ngspice writes, and the only thing that makes it safe is
/// that ngspice's reader does not care — which is checked here rather than
/// assumed, because the whole streaming path rests on it.
#[test]
fn ngspice_reads_a_padded_point_count() {
    let Some(bin) = find_ngspice() else {
        eprintln!("skipping: ngspice not installed");
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let bytes = ours(Encoding::Binary);
    let f = nutmeg::read(&bytes).unwrap();
    let probe = f.points.len() / 2;
    let expect = f.series("v(out)").unwrap()[probe];

    // Re-space the count field both ways, as `Writer::finish` and a hand-edit
    // might leave it.
    let text_end = bytes.windows(8).position(|w| w == b"Binary:\n").unwrap();
    let head = String::from_utf8_lossy(&bytes[..text_end]).to_string();
    let line = head
        .lines()
        .find(|l| l.starts_with("No. Points:"))
        .unwrap()
        .to_string();
    let n = line.trim_start_matches("No. Points:").trim().to_string();

    for (label, spaced) in [
        ("right-aligned", format!("No. Points: {n:>20}")),
        ("left-aligned", format!("No. Points: {n:<20}")),
    ] {
        let patched = head.replace(&line, &spaced);
        let mut out = patched.into_bytes();
        out.extend_from_slice(&bytes[text_end..]);

        let raw = dir.path().join(format!("{label}.raw"));
        std::fs::write(&raw, &out).unwrap();
        assert_eq!(
            nutmeg::read(&out).unwrap().points.len(),
            f.points.len(),
            "{label}: our own reader"
        );

        let deck = dir.path().join(format!("{label}.sp"));
        std::fs::write(
            &deck,
            format!(
                "* reader\nR1 a 0 1k\nV1 a 0 1\n.op\n.control\nrun\nload {}\n\
                 setplot tran1\nprint v(out)[{probe}]\nquit\n.endc\n.end\n",
                raw.display()
            ),
        )
        .unwrap();
        let res = Command::new(&bin).arg("-b").arg(&deck).output().unwrap();
        let text = String::from_utf8_lossy(&res.stdout);
        let got: f64 =
            printed(&text).unwrap_or_else(|| panic!("{label}: ngspice printed nothing:\n{text}"));
        assert!(
            (got - expect).abs() <= 1e-6 * expect.abs().max(1e-12),
            "{label}: ngspice read {got} where we wrote {expect}"
        );
    }
}
