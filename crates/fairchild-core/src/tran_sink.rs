//! Where a transient's timepoints go.
//!
//! Both integrators — fixed step and variable step — hand every accepted point
//! to a [`TranSink`], so a caller chooses between keeping the run in memory and
//! writing it straight out without either integrator knowing which.
//!
//! ## Why this exists
//!
//! [`crate::tran::TranResult`] is the whole run: one `Vec<f64>` per signal, all
//! of it resident until the analysis returns. A 200-node ladder over 200,001
//! timepoints is 323 MB before a byte is written, and the 10⁶ × 10⁴ run that
//! motivated this is 80 GB. No output format fixes that — the cost is in
//! holding the table, not in spelling it.
//!
//! A writing sink holds one row. Memory becomes O(signals) instead of
//! O(signals × timepoints), and the format question ([`RawSink`] against
//! [`CsvSink`], ASCII against binary) becomes a separate and much smaller one.
//!
//! ## What a sink sees
//!
//! [`TranSink::begin`] carries the [`TranLayout`] — the columns, in output
//! order, each knowing where in the solution vector its value lives. Every
//! later [`TranSink::point`] carries a time and the raw solution vector, and
//! the sink reads what it wants through the layout. This keeps the integrators
//! free of any opinion about output: they push `(t, x)` and nothing else.

use indexmap::IndexMap;

use crate::error::SimError;
use crate::mna::CircuitTopology;
use crate::nutmeg::{Encoding, Plot, Var, Writer};
use crate::tran::TranResult;

/// One output column, and where its value comes from.
#[derive(Debug, Clone)]
pub enum TranColumn {
    /// The independent variable.
    Time,
    /// A node voltage: row `index` of the solution vector.
    Voltage { name: String, index: usize },
    /// A voltage-source branch current: row `index` of the solution vector.
    Current { name: String, index: usize },
    /// A λ label. Resolved before the solve and constant across the run, so it
    /// is a value rather than a row — see `TranResult::lambda`.
    Lambda { name: String, value: f64 },
}

impl TranColumn {
    /// The column's value at a solution vector.
    pub fn value(&self, t: f64, x: &[f64]) -> f64 {
        match self {
            TranColumn::Time => t,
            TranColumn::Voltage { index, .. } | TranColumn::Current { index, .. } => x[*index],
            TranColumn::Lambda { value, .. } => *value,
        }
    }

    /// The column's name as a CSV header or `--probe` argument spells it.
    pub fn label(&self) -> String {
        match self {
            TranColumn::Time => "time".to_string(),
            TranColumn::Voltage { name, .. } | TranColumn::Lambda { name, .. } => {
                format!("V({name})")
            }
            TranColumn::Current { name, .. } => format!("I({name})"),
        }
    }

    pub(crate) fn raw_var(&self) -> Var {
        match self {
            TranColumn::Time => Var::new("time", "time"),
            TranColumn::Voltage { name, .. } | TranColumn::Lambda { name, .. } => {
                Var::voltage(format!("v({name})"))
            }
            TranColumn::Current { name, .. } => Var::current(format!("i({name})")),
        }
    }
}

/// The columns a transient writes, settled once the topology is resolved.
///
/// Order is time, node voltages, λ labels, branch currents — the order both
/// existing writers already used, so a reader of either sees no change.
#[derive(Debug, Clone)]
pub struct TranLayout {
    pub columns: Vec<TranColumn>,
}

impl TranLayout {
    pub fn from_topology(topo: &CircuitTopology) -> Self {
        let mut columns = vec![TranColumn::Time];
        for (name, &index) in &topo.node_index {
            columns.push(TranColumn::Voltage {
                name: name.clone(),
                index,
            });
        }
        for (name, wl) in topo.lambda_signals() {
            columns.push(TranColumn::Lambda {
                name: name.to_string(),
                value: wl,
            });
        }
        let n = topo.n_nodes();
        for (name, &index) in &topo.vsrc_index {
            columns.push(TranColumn::Current {
                name: name.clone(),
                index: n + index,
            });
        }
        TranLayout { columns }
    }

    /// Every column's label, in order.
    pub fn labels(&self) -> Vec<String> {
        self.columns.iter().map(|c| c.label()).collect()
    }

    /// Keep only the columns `probes` names, plus time.
    ///
    /// An empty `probes` keeps everything. A probe that names no column comes
    /// back as an error listing it, rather than being dropped — see
    /// [`crate::probe`].
    ///
    /// Selecting here rather than on rendered text is most of why `--probe`
    /// is worth using: the columns nobody asked for are never formatted, so
    /// asking for one signal of two hundred costs a two-hundredth of the work
    /// instead of slightly more than all of it.
    pub fn select(&self, probes: &[String]) -> Result<TranLayout, Vec<String>> {
        if probes.is_empty() {
            return Ok(self.clone());
        }
        let labels = self.labels();
        let missing = crate::probe::unmatched(probes, &labels);
        if !missing.is_empty() {
            return Err(missing);
        }
        let columns = self
            .columns
            .iter()
            .enumerate()
            .filter(|(i, _)| {
                *i == 0 || probes.iter().any(|p| crate::probe::matches(p, &labels[*i]))
            })
            .map(|(_, c)| c.clone())
            .collect();
        Ok(TranLayout { columns })
    }
}

/// Somewhere a transient's timepoints go.
pub trait TranSink {
    /// Called once, before the first point.
    fn begin(&mut self, layout: &TranLayout) -> Result<(), SimError>;
    /// One accepted timepoint, with the full solution vector.
    fn point(&mut self, t: f64, x: &[f64]) -> Result<(), SimError>;
    /// Called once, after the last point.
    fn end(&mut self) -> Result<(), SimError> {
        Ok(())
    }
    /// How many points this sink has taken.
    ///
    /// A streaming caller has no result to count afterwards, and counting the
    /// steps the integrator took would be a different number — tstart drops
    /// points on the way through. This is what was written.
    fn points_written(&self) -> usize;
}

fn io(e: std::io::Error) -> SimError {
    SimError::ParameterError(format!("writing transient output: {e}"))
}

// ---------------------------------------------------------------------------
// Collecting
// ---------------------------------------------------------------------------

/// Builds a [`TranResult`] — the whole run in memory, as before this module.
///
/// This is what every caller that wants a result back uses, and it is why the
/// sink could be introduced without changing a single one of them.
#[derive(Default)]
pub struct CollectSink {
    result: Option<TranResult>,
    /// Where each column of the layout lands in the result, so `point` is an
    /// indexed push rather than a map lookup per signal per timestep.
    plan: Vec<Sub>,
}

enum Sub {
    Time,
    Voltage {
        slot: usize,
        index: usize,
    },
    Current {
        slot: usize,
        index: usize,
    },
    /// λ is stored once in the result, not as a series.
    Skip,
}

impl CollectSink {
    pub fn new() -> Self {
        Self::default()
    }

    /// The finished result. `None` if the sink never saw [`TranSink::begin`].
    pub fn take(&mut self) -> Option<TranResult> {
        self.result.take()
    }
}

impl TranSink for CollectSink {
    fn begin(&mut self, layout: &TranLayout) -> Result<(), SimError> {
        let mut node_voltages: IndexMap<String, Vec<f64>> = IndexMap::new();
        let mut vsrc_currents: IndexMap<String, Vec<f64>> = IndexMap::new();
        let mut lambda: IndexMap<String, f64> = IndexMap::new();
        let mut plan = Vec::with_capacity(layout.columns.len());
        for c in &layout.columns {
            match c {
                TranColumn::Time => plan.push(Sub::Time),
                TranColumn::Voltage { name, index } => {
                    plan.push(Sub::Voltage {
                        slot: node_voltages.len(),
                        index: *index,
                    });
                    node_voltages.insert(name.clone(), Vec::new());
                }
                TranColumn::Current { name, index } => {
                    plan.push(Sub::Current {
                        slot: vsrc_currents.len(),
                        index: *index,
                    });
                    vsrc_currents.insert(name.clone(), Vec::new());
                }
                TranColumn::Lambda { name, value } => {
                    plan.push(Sub::Skip);
                    lambda.insert(name.clone(), *value);
                }
            }
        }
        self.plan = plan;
        self.result = Some(TranResult {
            time: Vec::new(),
            node_voltages,
            vsrc_currents,
            lambda,
        });
        Ok(())
    }

    fn points_written(&self) -> usize {
        self.result.as_ref().map_or(0, |r| r.time.len())
    }

    fn point(&mut self, t: f64, x: &[f64]) -> Result<(), SimError> {
        let r = self
            .result
            .as_mut()
            .expect("CollectSink::point before begin");
        for sub in &self.plan {
            match sub {
                Sub::Time => r.time.push(t),
                Sub::Voltage { slot, index } => r.node_voltages[*slot].push(x[*index]),
                Sub::Current { slot, index } => r.vsrc_currents[*slot].push(x[*index]),
                Sub::Skip => {}
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

/// Streams a Nutmeg rawfile, ASCII or binary.
///
/// The point count is not known when the header is written — a variable-step
/// run decides it by running — so the writer reserves the field and fills it
/// in at [`TranSink::end`]. That is why `W` must seek; a caller with nowhere
/// to seek to can hand over a `Cursor<Vec<u8>>`.
pub struct RawSink<W: std::io::Write + std::io::Seek> {
    w: Option<W>,
    writer: Option<Writer<W>>,
    title: String,
    enc: Encoding,
    row: Vec<f64>,
    cols: Vec<TranColumn>,
    /// Carried past `end`, when the inner writer has been handed back.
    n_points: usize,
}

impl<W: std::io::Write + std::io::Seek> RawSink<W> {
    pub fn new(w: W, title: impl Into<String>, enc: Encoding) -> Self {
        RawSink {
            w: Some(w),
            writer: None,
            title: title.into(),
            enc,
            row: Vec::new(),
            cols: Vec::new(),
            n_points: 0,
        }
    }
}

impl<W: std::io::Write + std::io::Seek> RawSink<W> {
    /// The writer, after [`TranSink::end`].
    ///
    /// A caller that handed over a `Cursor` because it had nowhere to seek to
    /// gets the bytes back this way. `None` before `end`.
    pub fn into_inner(self) -> Option<W> {
        self.w
    }
}

impl<W: std::io::Write + std::io::Seek> TranSink for RawSink<W> {
    fn begin(&mut self, layout: &TranLayout) -> Result<(), SimError> {
        let plot = Plot {
            title: &self.title,
            plotname: "Transient Analysis",
            complex: false,
            vars: layout.columns.iter().map(|c| c.raw_var()).collect(),
            n_points: None,
        };
        let w = self.w.take().expect("RawSink::begin twice");
        self.writer = Some(Writer::start_counting(w, self.enc, &plot).map_err(io)?);
        self.row = vec![0.0; layout.columns.len()];
        self.cols = layout.columns.clone();
        Ok(())
    }

    fn point(&mut self, t: f64, x: &[f64]) -> Result<(), SimError> {
        for (slot, c) in self.cols.iter().enumerate() {
            self.row[slot] = c.value(t, x);
        }
        self.writer
            .as_mut()
            .expect("RawSink::point before begin")
            .point(&self.row)
            .map_err(io)
    }

    fn end(&mut self) -> Result<(), SimError> {
        if let Some(wr) = self.writer.take() {
            self.n_points = wr.points_written();
            self.w = Some(wr.finish().map_err(io)?);
        }
        Ok(())
    }

    fn points_written(&self) -> usize {
        self.writer
            .as_ref()
            .map_or(self.n_points, |wr| wr.points_written())
    }
}

/// Streams CSV.
///
/// The header names the columns the layout carries, so a selected layout emits
/// a narrow file without the full one ever existing.
pub struct CsvSink<W: std::io::Write> {
    w: W,
    cols: Vec<TranColumn>,
    n_points: usize,
}

impl<W: std::io::Write> CsvSink<W> {
    pub fn new(w: W) -> Self {
        CsvSink {
            w,
            cols: Vec::new(),
            n_points: 0,
        }
    }

    /// How many points have been written.
    pub fn points_written(&self) -> usize {
        self.n_points
    }
}

impl<W: std::io::Write> TranSink for CsvSink<W> {
    fn begin(&mut self, layout: &TranLayout) -> Result<(), SimError> {
        self.cols = layout.columns.clone();
        let header = layout.labels().join(",");
        writeln!(self.w, "{header}").map_err(io)
    }

    fn point(&mut self, t: f64, x: &[f64]) -> Result<(), SimError> {
        for (i, c) in self.cols.iter().enumerate() {
            if i > 0 {
                write!(self.w, ",").map_err(io)?;
            }
            write!(self.w, "{:.6e}", c.value(t, x)).map_err(io)?;
        }
        self.n_points += 1;
        writeln!(self.w).map_err(io)
    }

    fn end(&mut self) -> Result<(), SimError> {
        self.w.flush().map_err(io)
    }

    fn points_written(&self) -> usize {
        self.n_points
    }
}

// ---------------------------------------------------------------------------
// tstart
// ---------------------------------------------------------------------------

/// Drops points before `tstart` on the way through.
///
/// `.tran`'s third argument selects what is *saved*, not where integration
/// begins, so it is a filter on the output and not a change to the solve.
/// Applying it here rather than to a finished result is what lets a streaming
/// sink honour it at all, and it means the fixed-step and variable-step paths
/// cannot come to disagree about what tstart means.
///
/// A tstart past the end of the run leaves the last point, because an empty
/// waveform is a worse answer than a short one. That needs one point of
/// lookahead: the most recent dropped point is held, and released at
/// [`TranSink::end`] if nothing else was ever emitted.
pub struct TstartSink<'a> {
    inner: &'a mut dyn TranSink,
    tstart: f64,
    emitted: bool,
    pending: Option<(f64, Vec<f64>)>,
}

impl<'a> TstartSink<'a> {
    pub fn new(inner: &'a mut dyn TranSink, tstart: f64) -> Self {
        TstartSink {
            inner,
            tstart,
            emitted: false,
            pending: None,
        }
    }
}

impl TranSink for TstartSink<'_> {
    fn begin(&mut self, layout: &TranLayout) -> Result<(), SimError> {
        self.inner.begin(layout)
    }

    fn point(&mut self, t: f64, x: &[f64]) -> Result<(), SimError> {
        if self.tstart > 0.0 && t < self.tstart {
            self.pending = Some((t, x.to_vec()));
            return Ok(());
        }
        self.emitted = true;
        self.pending = None;
        self.inner.point(t, x)
    }

    fn end(&mut self) -> Result<(), SimError> {
        if !self.emitted {
            if let Some((t, x)) = self.pending.take() {
                self.inner.point(t, &x)?;
            }
        }
        self.inner.end()
    }

    fn points_written(&self) -> usize {
        self.inner.points_written()
    }
}

/// Narrows the layout to the columns `--probe` names, before anything is
/// written.
///
/// This is what makes `--probe` worth asking for. Filtering a rendered result
/// gives the user the columns they wanted and the simulator all the work: a
/// 200-signal transient asked for one signal formatted all two hundred, then
/// threw away 199 — which cost *more* than not filtering, because the filter
/// re-parsed the text it had just written. Narrowing the layout means the other
/// 199 are never read out of the solution vector at all.
///
/// A probe that names no column is an error rather than a silent drop (#72).
pub struct SelectSink<'a> {
    inner: &'a mut dyn TranSink,
    probes: Vec<String>,
}

impl<'a> SelectSink<'a> {
    pub fn new(inner: &'a mut dyn TranSink, probes: &[String]) -> Self {
        SelectSink {
            inner,
            probes: probes.to_vec(),
        }
    }
}

impl TranSink for SelectSink<'_> {
    fn begin(&mut self, layout: &TranLayout) -> Result<(), SimError> {
        let narrowed = layout.select(&self.probes).map_err(|missing| {
            let list = missing
                .iter()
                .map(|m| format!("'{m}'"))
                .collect::<Vec<_>>()
                .join(", ");
            SimError::ParameterError(format!(
                "--probe {list} matched no signal of the transient output. \
                 Signals are spelled V(<node>) and I(<vsource>), case-insensitively \
                 (run without --probe to see the full list; --list-nodes prints the \
                 deck's nets)"
            ))
        })?;
        self.inner.begin(&narrowed)
    }

    fn point(&mut self, t: f64, x: &[f64]) -> Result<(), SimError> {
        self.inner.point(t, x)
    }

    fn end(&mut self) -> Result<(), SimError> {
        self.inner.end()
    }

    fn points_written(&self) -> usize {
        self.inner.points_written()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn layout() -> TranLayout {
        TranLayout {
            columns: vec![
                TranColumn::Time,
                TranColumn::Voltage {
                    name: "out".into(),
                    index: 0,
                },
                TranColumn::Voltage {
                    name: "mid".into(),
                    index: 1,
                },
                TranColumn::Current {
                    name: "v1".into(),
                    index: 2,
                },
            ],
        }
    }

    /// Selecting a column must drop the others and keep time, and must report
    /// a name that selects nothing rather than quietly narrowing the output.
    #[test]
    fn select_keeps_time_and_reports_misses() {
        let l = layout();
        let got = l.select(&["V(mid)".into()]).unwrap();
        assert_eq!(got.labels(), vec!["time", "V(mid)"]);
        let err = l.select(&["V(nope)".into()]).unwrap_err();
        assert_eq!(err, vec!["V(nope)".to_string()]);
        assert_eq!(l.select(&[]).unwrap().labels().len(), 4);
    }

    /// A selected layout must not merely hide columns downstream — the point
    /// of selecting is that the dropped columns are never written.
    #[test]
    fn csv_writes_only_selected_columns() {
        let l = layout().select(&["V(out)".into()]).unwrap();
        let mut s = CsvSink::new(Cursor::new(Vec::new()));
        s.begin(&l).unwrap();
        s.point(1e-9, &[0.5, 0.25, -1e-3]).unwrap();
        s.end().unwrap();
        let text = String::from_utf8(s.w.into_inner()).unwrap();
        assert_eq!(text, "time,V(out)\n1.000000e-9,5.000000e-1\n");
    }

    /// tstart selects what is saved. The last point survives a tstart past the
    /// end of the run, which is the case a streaming filter is most likely to
    /// get wrong — it has to hold a point it has already decided to drop.
    #[test]
    fn tstart_drops_early_points_but_never_all_of_them() {
        const TIMES: [f64; 4] = [0.0, 1e-9, 2e-9, 3e-9];
        let mut collect = CollectSink::new();
        {
            let mut s = TstartSink::new(&mut collect, 2e-9);
            s.begin(&layout()).unwrap();
            for (k, &t) in TIMES.iter().enumerate() {
                s.point(t, &[k as f64, 0.0, 0.0]).unwrap();
            }
            s.end().unwrap();
        }
        let r = collect.take().unwrap();
        assert_eq!(r.time, vec![2e-9, 3e-9]);

        let mut collect = CollectSink::new();
        {
            let mut s = TstartSink::new(&mut collect, 100.0);
            s.begin(&layout()).unwrap();
            for (k, &t) in TIMES.iter().enumerate() {
                s.point(t, &[k as f64, 0.0, 0.0]).unwrap();
            }
            s.end().unwrap();
        }
        let r = collect.take().unwrap();
        assert_eq!(r.time.len(), 1);
        assert_eq!(r.node_voltages["out"], vec![3.0]);
    }

    /// Selection has to happen at `begin`, before a single value is written —
    /// a wrapper that filtered on the way out would have done the formatting
    /// this exists to avoid.
    #[test]
    fn select_sink_narrows_the_writer_and_refuses_a_bad_probe() {
        let mut collect = CollectSink::new();
        {
            let mut s = SelectSink::new(&mut collect, &["I(v1)".into()]);
            s.begin(&layout()).unwrap();
            s.point(0.0, &[1.0, 2.0, 3.0]).unwrap();
            s.end().unwrap();
        }
        let r = collect.take().unwrap();
        assert!(r.node_voltages.is_empty());
        assert_eq!(r.vsrc_currents["v1"], vec![3.0]);

        let mut collect = CollectSink::new();
        let mut s = SelectSink::new(&mut collect, &["V(nope)".into()]);
        let err = s.begin(&layout()).unwrap_err().to_string();
        assert!(err.contains("V(nope)"), "{err}");
    }

    /// A λ column is a label, not a series: it must reach the rawfile at every
    /// point and the result exactly once.
    #[test]
    fn lambda_is_a_label_in_the_result_and_a_column_in_the_file() {
        let l = TranLayout {
            columns: vec![
                TranColumn::Time,
                TranColumn::Lambda {
                    name: "bus_wl_0".into(),
                    value: 1.31e-6,
                },
            ],
        };
        let mut collect = CollectSink::new();
        collect.begin(&l).unwrap();
        collect.point(0.0, &[]).unwrap();
        collect.point(1e-9, &[]).unwrap();
        let r = collect.take().unwrap();
        assert_eq!(r.lambda["bus_wl_0"], 1.31e-6);
        assert!(r.node_voltages.is_empty());

        let mut raw = RawSink::new(Cursor::new(Vec::new()), "t", Encoding::Binary);
        raw.begin(&l).unwrap();
        raw.point(0.0, &[]).unwrap();
        raw.point(1e-9, &[]).unwrap();
        raw.end().unwrap();
        let bytes = raw.w.take().unwrap().into_inner();
        let f = crate::nutmeg::read(&bytes).unwrap();
        assert_eq!(f.series("v(bus_wl_0)").unwrap(), vec![1.31e-6, 1.31e-6]);
    }
}
