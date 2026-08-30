//! Nutmeg rawfiles, in both spellings the format has.
//!
//! Every analysis writes the same header and then either an ASCII or a binary
//! body. Before this module each analysis wrote its own copy of that header —
//! eight copies of the same seven lines — and adding the binary spelling would
//! have made sixteen. One place decides what a rawfile is, so ASCII and binary
//! cannot come to disagree about it, and so a new analysis gets both spellings
//! by writing a `Plot` rather than by copying a writer.
//!
//! The format, as ngspice writes it: seven header lines, a `Variables:` block
//! of one tab-separated line per variable, then `Values:` (ASCII) or `Binary:`
//! (binary). ASCII writes one value per line, with the point index on the first
//! variable's line only. Binary writes native-endian `f64` immediately after
//! the newline, point-major, with no padding and no separators. A complex plot
//! writes two numbers per value in both spellings.
//!
//! ## The point count
//!
//! `No. Points:` sits in the header, above data that has not been produced yet.
//! An analysis that knows its count in advance writes it directly. A streamed
//! transient does not — the run decides how many points there are, and tstart
//! then drops some of them — so [`Writer::start_counting`] reserves a
//! fixed-width field and [`Writer::finish`] seeks back and fills it in,
//! right-aligned so the line carries no trailing whitespace.
//!
//! ngspice reads the padded field exactly as it reads an unpadded one, which
//! was checked against ngspice 46 rather than assumed: a rawfile padded on
//! either side of the number still reports `62 long` and the same value at the
//! same index.
//!
//! Only that one path needs [`std::io::Seek`]. A caller with nowhere to seek
//! to (a pipe, standard output) can wrap a `Cursor<Vec<u8>>` and write the
//! bytes on afterwards.

use std::io::{Seek, SeekFrom, Write};

/// How the values are spelled on disk.
///
/// The header is identical either way — a reader distinguishes the two by the
/// `Values:` / `Binary:` line that ends it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Encoding {
    /// One value per line, six significant figures.
    Ascii,
    /// Native-endian `f64`, point-major, no separators.
    Binary,
}

impl Encoding {
    fn tag(self) -> &'static str {
        match self {
            Encoding::Ascii => "Values:",
            Encoding::Binary => "Binary:",
        }
    }
}

/// One rawfile variable: the name a reader asks for it by, and its ngspice
/// type tag (`time`, `frequency`, `voltage`, `current`, `notype`).
#[derive(Debug, Clone)]
pub struct Var {
    pub name: String,
    pub kind: &'static str,
}

impl Var {
    pub fn new(name: impl Into<String>, kind: &'static str) -> Self {
        Var {
            name: name.into(),
            kind,
        }
    }
    pub fn voltage(name: impl Into<String>) -> Self {
        Var::new(name, "voltage")
    }
    pub fn current(name: impl Into<String>) -> Self {
        Var::new(name, "current")
    }
    pub fn notype(name: impl Into<String>) -> Self {
        Var::new(name, "notype")
    }
}

/// Everything the header states about a plot.
#[derive(Debug, Clone)]
pub struct Plot<'a> {
    pub title: &'a str,
    pub plotname: &'a str,
    /// `Flags: complex` rather than `Flags: real`. A complex plot must be
    /// written with [`Writer::point_complex`].
    pub complex: bool,
    pub vars: Vec<Var>,
    /// The number of points, when the analysis knows it before it runs.
    /// `None` reserves the field and fills it in at [`Writer::finish`].
    pub n_points: Option<usize>,
}

/// Width of the reserved `No. Points:` field, in characters.
///
/// A `usize` cannot print wider than 20 digits, so the field can always hold
/// the count that replaces it.
const COUNT_FIELD: usize = 20;

/// A rawfile part-written: header on disk, points appended one at a time.
pub struct Writer<W> {
    w: W,
    enc: Encoding,
    complex: bool,
    n_vars: usize,
    n_points: usize,
    /// Byte offset of the reserved count field, when the count was not known
    /// up front. `None` means the header already states the right number.
    count_at: Option<u64>,
}

impl<W: Write> Writer<W> {
    /// Write the header of a plot whose point count is already known. Every
    /// later call appends one point.
    ///
    /// An analysis whose count is decided by running — the variable-step
    /// transient, and only it — wants [`Writer::start_counting`] instead.
    pub fn start(w: W, enc: Encoding, plot: &Plot) -> std::io::Result<Self> {
        let n = plot
            .n_points
            .expect("Writer::start needs a known point count; use start_counting");
        Self::header(w, enc, plot, Some(n), None)
    }

    /// Shared between the two constructors: `reserve_at` is the offset of the
    /// count field when one was reserved, which only the seekable path knows.
    fn header(
        mut w: W,
        enc: Encoding,
        plot: &Plot,
        count: Option<usize>,
        reserve_at: Option<u64>,
    ) -> std::io::Result<Self> {
        writeln!(w, "Title: {}", plot.title)?;
        writeln!(w, "Plotname: {}", plot.plotname)?;
        writeln!(
            w,
            "Flags: {}",
            if plot.complex { "complex" } else { "real" }
        )?;
        writeln!(w, "No. Variables: {}", plot.vars.len())?;
        match count {
            Some(n) => writeln!(w, "No. Points: {n}")?,
            None => writeln!(w, "No. Points: {:>width$}", "", width = COUNT_FIELD)?,
        }
        writeln!(w, "Variables:")?;
        for (i, v) in plot.vars.iter().enumerate() {
            writeln!(w, "\t{i}\t{}\t{}", v.name, v.kind)?;
        }
        writeln!(w, "{}", enc.tag())?;

        Ok(Writer {
            w,
            enc,
            complex: plot.complex,
            n_vars: plot.vars.len(),
            n_points: 0,
            count_at: reserve_at,
        })
    }

    /// Append one point of a real plot.
    ///
    /// `values` must carry one entry per variable, in the order the header
    /// declared them. A short or long row is a caller bug that would silently
    /// shift every later value against its name, so it panics rather than
    /// writing a file whose columns have quietly moved.
    pub fn point(&mut self, values: &[f64]) -> std::io::Result<()> {
        assert!(
            !self.complex,
            "point() on a complex plot: use point_complex()"
        );
        assert_eq!(
            values.len(),
            self.n_vars,
            "rawfile point has {} values against {} variables",
            values.len(),
            self.n_vars
        );
        match self.enc {
            Encoding::Ascii => {
                for (k, v) in values.iter().enumerate() {
                    if k == 0 {
                        writeln!(self.w, " {}\t{v:.6e}", self.n_points)?;
                    } else {
                        writeln!(self.w, "\t{v:.6e}")?;
                    }
                }
            }
            Encoding::Binary => {
                let mut buf = Vec::with_capacity(values.len() * 8);
                for v in values {
                    buf.extend_from_slice(&v.to_ne_bytes());
                }
                self.w.write_all(&buf)?;
            }
        }
        self.n_points += 1;
        Ok(())
    }

    /// Append one point of a complex plot, as `(re, im)` per variable.
    pub fn point_complex(&mut self, values: &[(f64, f64)]) -> std::io::Result<()> {
        assert!(self.complex, "point_complex() on a real plot: use point()");
        assert_eq!(
            values.len(),
            self.n_vars,
            "rawfile point has {} values against {} variables",
            values.len(),
            self.n_vars
        );
        match self.enc {
            Encoding::Ascii => {
                for (k, (re, im)) in values.iter().enumerate() {
                    if k == 0 {
                        writeln!(self.w, " {}\t{re:.6e},{im:.6e}", self.n_points)?;
                    } else {
                        writeln!(self.w, "\t{re:.6e},{im:.6e}")?;
                    }
                }
            }
            Encoding::Binary => {
                let mut buf = Vec::with_capacity(values.len() * 16);
                for (re, im) in values {
                    buf.extend_from_slice(&re.to_ne_bytes());
                    buf.extend_from_slice(&im.to_ne_bytes());
                }
                self.w.write_all(&buf)?;
            }
        }
        self.n_points += 1;
        Ok(())
    }

    /// Flush and hand the writer back.
    ///
    /// A writer built by [`Writer::start_counting`] must go through
    /// [`Writer::finish`] instead — this one cannot seek back to the count.
    pub fn done(mut self) -> std::io::Result<W> {
        assert!(
            self.count_at.is_none(),
            "a reserved point count needs finish(), not done()"
        );
        self.w.flush()?;
        Ok(self.w)
    }

    /// How many points have been written so far.
    pub fn points_written(&self) -> usize {
        self.n_points
    }
}

impl<W: Write + Seek> Writer<W> {
    /// Write the header of a plot whose point count is not known yet,
    /// reserving the field for [`Writer::finish`] to fill in.
    ///
    /// `plot.n_points` is ignored. A run that dies before `finish` leaves a
    /// blank count, which is what makes the reader's count check the
    /// difference between a truncated file and a whole one.
    pub fn start_counting(mut w: W, enc: Encoding, plot: &Plot) -> std::io::Result<Self> {
        // Four header lines precede the count; measure rather than assume,
        // since a caller may hand over a writer already partway through a file.
        w.flush()?;
        let before = w.stream_position()?;
        let prefix = format!(
            "Title: {}\nPlotname: {}\nFlags: {}\nNo. Variables: {}\nNo. Points: ",
            plot.title,
            plot.plotname,
            if plot.complex { "complex" } else { "real" },
            plot.vars.len(),
        );
        let at = before + prefix.len() as u64;
        Self::header(w, enc, plot, None, Some(at))
    }

    /// Fill in a reserved point count and hand the writer back.
    pub fn finish(mut self) -> std::io::Result<W> {
        if let Some(at) = self.count_at {
            self.w.flush()?;
            let end = self.w.stream_position()?;
            self.w.seek(SeekFrom::Start(at))?;
            write!(self.w, "{:>width$}", self.n_points, width = COUNT_FIELD)?;
            self.w.flush()?;
            self.w.seek(SeekFrom::Start(end))?;
        }
        self.w.flush()?;
        Ok(self.w)
    }
}

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

/// A rawfile read back: the header, and the values as one row per point.
#[derive(Debug, Clone)]
pub struct RawFile {
    pub title: String,
    pub plotname: String,
    pub complex: bool,
    pub vars: Vec<Var>,
    /// `[point][variable]`, as `(re, im)`. A real plot carries `im == 0.0`.
    pub points: Vec<Vec<(f64, f64)>>,
}

impl RawFile {
    /// The series for one variable by name, real part only.
    pub fn series(&self, name: &str) -> Option<Vec<f64>> {
        let i = self.vars.iter().position(|v| v.name == name)?;
        Some(self.points.iter().map(|p| p[i].0).collect())
    }
}

/// Read a rawfile in either spelling.
///
/// The count in the header is checked against the number of points actually
/// present rather than trusted: a truncated run and a complete one differ in
/// exactly that, and a reader that trusts the header reports a short file as a
/// whole one. This is the reader the round-trip test compares the two
/// spellings with, so it must not paper over a difference between them.
pub fn read(bytes: &[u8]) -> Result<RawFile, String> {
    // The header is ASCII whatever the body is, so it can be split off by
    // scanning for the line that ends it.
    let split = find_body(bytes).ok_or("rawfile has no Values: or Binary: line")?;
    let (head, body, binary) = split;
    let head =
        std::str::from_utf8(head).map_err(|e| format!("rawfile header is not UTF-8: {e}"))?;

    let mut title = String::new();
    let mut plotname = String::new();
    let mut complex = false;
    let mut n_vars = 0usize;
    let mut n_points = 0usize;
    let mut vars: Vec<Var> = Vec::new();
    let mut in_vars = false;

    for line in head.lines() {
        if let Some(v) = line.strip_prefix("Title: ") {
            title = v.to_string();
        } else if let Some(v) = line.strip_prefix("Plotname: ") {
            plotname = v.to_string();
        } else if let Some(v) = line.strip_prefix("Flags: ") {
            complex = v.trim() == "complex";
        } else if let Some(v) = line.strip_prefix("No. Variables: ") {
            n_vars = v.trim().parse().map_err(|_| "bad No. Variables")?;
        } else if let Some(v) = line.strip_prefix("No. Points: ") {
            n_points = v.trim().parse().map_err(|_| "bad No. Points")?;
        } else if line.starts_with("Variables:") {
            in_vars = true;
        } else if in_vars {
            let f: Vec<&str> = line.split('\t').filter(|s| !s.is_empty()).collect();
            if f.len() >= 2 {
                // The type tag is a fixed vocabulary; anything else reads as
                // `notype` rather than failing, since a rawfile from another
                // tool may carry a tag we do not use.
                let kind = match f.get(2).map(|s| s.trim()) {
                    Some("time") => "time",
                    Some("frequency") => "frequency",
                    Some("voltage") => "voltage",
                    Some("current") => "current",
                    Some("voltage-density") => "voltage-density",
                    _ => "notype",
                };
                vars.push(Var::new(f[1].trim(), kind));
            }
        }
    }
    if vars.len() != n_vars {
        return Err(format!(
            "rawfile declares {n_vars} variables and lists {}",
            vars.len()
        ));
    }

    let per_value = if complex { 2 } else { 1 };
    let points = if binary {
        let stride = n_vars * per_value * 8;
        if stride == 0 {
            return Err("rawfile has no variables".into());
        }
        if body.len() % stride != 0 {
            return Err(format!(
                "binary rawfile body is {} bytes, not a whole number of {stride}-byte points",
                body.len()
            ));
        }
        body.chunks_exact(stride)
            .map(|pt| {
                pt.chunks_exact(8 * per_value)
                    .map(|v| {
                        let re = f64::from_ne_bytes(v[..8].try_into().unwrap());
                        let im = if complex {
                            f64::from_ne_bytes(v[8..16].try_into().unwrap())
                        } else {
                            0.0
                        };
                        (re, im)
                    })
                    .collect()
            })
            .collect()
    } else {
        let text =
            std::str::from_utf8(body).map_err(|e| format!("rawfile body is not UTF-8: {e}"))?;
        let mut points: Vec<Vec<(f64, f64)>> = Vec::new();
        let mut cur: Vec<(f64, f64)> = Vec::new();
        for line in text.lines() {
            let t = line.trim();
            if t.is_empty() {
                continue;
            }
            // The point index leads the first variable's line; the value is
            // whatever follows the tab either way.
            let val = t.rsplit('\t').next().unwrap_or(t).trim();
            let (re, im) = match val.split_once(',') {
                Some((a, b)) => (parse(a)?, parse(b)?),
                None => (parse(val)?, 0.0),
            };
            cur.push((re, im));
            if cur.len() == n_vars {
                points.push(std::mem::take(&mut cur));
            }
        }
        if !cur.is_empty() {
            return Err(format!(
                "ASCII rawfile ends {} values into a {n_vars}-variable point",
                cur.len()
            ));
        }
        points
    };

    if points.len() != n_points {
        return Err(format!(
            "rawfile header says {n_points} points and the body carries {}",
            points.len()
        ));
    }

    Ok(RawFile {
        title,
        plotname,
        complex,
        vars,
        points,
    })
}

fn parse(s: &str) -> Result<f64, String> {
    s.trim()
        .parse::<f64>()
        .map_err(|_| format!("rawfile value {s:?} is not a number"))
}

/// Split a rawfile at the line that ends its header.
///
/// Returns `(header, body, is_binary)`. The search is over bytes rather than
/// text because a binary body is not UTF-8 and would fail to decode.
fn find_body(bytes: &[u8]) -> Option<(&[u8], &[u8], bool)> {
    for (tag, binary) in [(&b"Binary:\n"[..], true), (&b"Values:\n"[..], false)] {
        if let Some(at) = find(bytes, tag) {
            return Some((&bytes[..at], &bytes[at + tag.len()..], binary));
        }
    }
    None
}

fn find(hay: &[u8], needle: &[u8]) -> Option<usize> {
    hay.windows(needle.len()).position(|w| w == needle)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn plot(n: Option<usize>) -> Plot<'static> {
        Plot {
            title: "t",
            plotname: "Transient Analysis",
            complex: false,
            vars: vec![Var::new("time", "time"), Var::voltage("v(out)")],
            n_points: n,
        }
    }

    fn write(enc: Encoding, n: Option<usize>, rows: &[[f64; 2]]) -> Vec<u8> {
        let c = Cursor::new(Vec::new());
        let mut wr = match n {
            Some(_) => Writer::start(c, enc, &plot(n)).unwrap(),
            None => Writer::start_counting(c, enc, &plot(n)).unwrap(),
        };
        for r in rows {
            wr.point(r).unwrap();
        }
        match n {
            Some(_) => wr.done().unwrap().into_inner(),
            None => wr.finish().unwrap().into_inner(),
        }
    }

    const ROWS: [[f64; 2]; 3] = [[0.0, 1.0], [1e-9, 0.5], [2e-9, 0.25]];

    /// The two spellings must carry the same numbers. Not "close": the binary
    /// path exists to avoid the ASCII path's rounding, so a mismatch beyond
    /// six figures is the ASCII writer's format and anything larger is a bug.
    #[test]
    fn ascii_and_binary_agree() {
        let a = read(&write(Encoding::Ascii, Some(3), &ROWS)).unwrap();
        let b = read(&write(Encoding::Binary, Some(3), &ROWS)).unwrap();
        assert_eq!(a.vars.len(), b.vars.len());
        assert_eq!(a.points.len(), b.points.len());
        for (pa, pb) in a.points.iter().zip(&b.points) {
            for (va, vb) in pa.iter().zip(pb) {
                assert!(
                    (va.0 - vb.0).abs() <= 1e-6 * vb.0.abs(),
                    "{} against {}",
                    va.0,
                    vb.0
                );
            }
        }
    }

    /// Binary is exact where ASCII rounds. This is most of the point of it, so
    /// it is asserted rather than assumed: a value with more than six figures
    /// must survive the binary round trip bit for bit and must *not* survive
    /// the ASCII one.
    #[test]
    fn binary_is_exact_and_ascii_is_not() {
        let v = std::f64::consts::PI * 1e-7;
        let rows = [[0.0, v]];
        let b = read(&write(Encoding::Binary, Some(1), &rows)).unwrap();
        assert_eq!(b.points[0][1].0.to_bits(), v.to_bits());
        let a = read(&write(Encoding::Ascii, Some(1), &rows)).unwrap();
        assert_ne!(a.points[0][1].0.to_bits(), v.to_bits());
    }

    /// A reserved count must end up stating the real number of points, in both
    /// spellings — the variable-step transient has no other way to write one.
    #[test]
    fn reserved_count_is_filled_in() {
        for enc in [Encoding::Ascii, Encoding::Binary] {
            let bytes = write(enc, None, &ROWS);
            let head = String::from_utf8_lossy(&bytes[..200.min(bytes.len())]).to_string();
            let counted = head
                .lines()
                .find_map(|l| l.strip_prefix("No. Points:"))
                .map(str::trim);
            assert_eq!(counted, Some("3"), "{enc:?}: {head}");
            let r = read(&bytes).unwrap();
            assert_eq!(r.points.len(), 3, "{enc:?}");
        }
    }

    /// A body that stops short of the header's promise is a truncated run.
    /// Reporting it as a whole one is the failure this reader exists to avoid.
    #[test]
    fn short_body_is_an_error() {
        let bytes = write(Encoding::Binary, Some(3), &ROWS);
        let cut = bytes.len() - 16;
        let err = read(&bytes[..cut]).unwrap_err();
        assert!(err.contains("points"), "{err}");
    }

    #[test]
    fn complex_round_trips() {
        let p = Plot {
            title: "t",
            plotname: "AC Analysis",
            complex: true,
            vars: vec![Var::new("frequency", "frequency"), Var::voltage("v(out)")],
            n_points: Some(2),
        };
        let mut wr = Writer::start(Cursor::new(Vec::new()), Encoding::Binary, &p).unwrap();
        wr.point_complex(&[(1.0, 0.0), (0.5, -0.25)]).unwrap();
        wr.point_complex(&[(2.0, 0.0), (0.25, -0.5)]).unwrap();
        let bytes = wr.done().unwrap().into_inner();
        let r = read(&bytes).unwrap();
        assert!(r.complex);
        assert_eq!(r.points[1][1], (0.25, -0.5));
    }
}
