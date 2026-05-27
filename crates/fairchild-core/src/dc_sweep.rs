//! DC sweep analysis (`.dc` directive).
//!
//! For each value in `[start, stop]` (step `step`), the named source's DC
//! amplitude is overridden and a fresh DC operating-point is solved.  The
//! sequence of operating points is returned as a `DcSweepResult` that can
//! be written as CSV or Nutmeg.

use indexmap::IndexMap;
use rayon::prelude::*;

use fairchild_parser::{Element, Netlist, Waveform};

use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::newton::dc_op_nr_with_registry_opts;
use crate::options::SimOptions;

/// One axis of a sweep: name of the swept source plus the linear point grid.
#[derive(Debug, Clone)]
pub struct SweepAxis {
    pub src: String,
    pub values: Vec<f64>,
}

/// Result of a DC sweep.
///
/// For a 1-D sweep `inner` is `None` and each `node_voltages[name][i]` is the
/// operating-point voltage at sweep-point `i`.  For a 2-D nested sweep,
/// timepoints are laid out outer-major: index `i*inner_len + j` corresponds
/// to outer point `i`, inner point `j`.
pub struct DcSweepResult {
    pub outer: SweepAxis,
    pub inner: Option<SweepAxis>,
    pub node_voltages: IndexMap<String, Vec<f64>>,
    pub vsrc_currents: IndexMap<String, Vec<f64>>,
}

impl DcSweepResult {
    /// Total number of sweep points.
    pub fn n_points(&self) -> usize {
        self.outer.values.len() * self.inner.as_ref().map_or(1, |i| i.values.len())
    }

    /// Write the sweep as CSV.  Columns: `<outer_src>`, optional `<inner_src>`,
    /// then `V(node)…` and `I(vsrc)…`.
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        write!(w, "{}", self.outer.src)?;
        if let Some(inner) = &self.inner {
            write!(w, ",{}", inner.src)?;
        }
        for name in self.node_voltages.keys() {
            write!(w, ",V({name})")?;
        }
        for name in self.vsrc_currents.keys() {
            write!(w, ",I({name})")?;
        }
        writeln!(w)?;

        let inner_len = self.inner.as_ref().map_or(1, |i| i.values.len());
        for (i, &outer_v) in self.outer.values.iter().enumerate() {
            for j in 0..inner_len {
                let idx = i * inner_len + j;
                write!(w, "{outer_v:.6e}")?;
                if let Some(inner) = &self.inner {
                    write!(w, ",{:.6e}", inner.values[j])?;
                }
                for series in self.node_voltages.values() {
                    write!(w, ",{:.6e}", series[idx])?;
                }
                for series in self.vsrc_currents.values() {
                    write!(w, ",{:.6e}", series[idx])?;
                }
                writeln!(w)?;
            }
        }
        Ok(())
    }

    /// Write the sweep as an ngspice-compatible Nutmeg ASCII rawfile.
    pub fn write_nutmeg<W: std::io::Write>(&self, mut w: W, title: &str) -> std::io::Result<()> {
        let n_vars =
            1 + self.inner.is_some() as usize + self.node_voltages.len() + self.vsrc_currents.len();
        let n_pts = self.n_points();
        writeln!(w, "Title: {title}")?;
        writeln!(w, "Plotname: DC Sweep")?;
        writeln!(w, "Flags: real")?;
        writeln!(w, "No. Variables: {n_vars}")?;
        writeln!(w, "No. Points: {n_pts}")?;
        writeln!(w, "Variables:")?;
        let mut idx = 0;
        writeln!(w, "\t{idx}\t{}\tvoltage", self.outer.src)?;
        idx += 1;
        if let Some(inner) = &self.inner {
            writeln!(w, "\t{idx}\t{}\tvoltage", inner.src)?;
            idx += 1;
        }
        for name in self.node_voltages.keys() {
            writeln!(w, "\t{idx}\tv({name})\tvoltage")?;
            idx += 1;
        }
        for name in self.vsrc_currents.keys() {
            writeln!(w, "\t{idx}\ti({name})\tcurrent")?;
            idx += 1;
        }
        writeln!(w, "Values:")?;

        let inner_len = self.inner.as_ref().map_or(1, |i| i.values.len());
        let mut point = 0usize;
        for (i, &outer_v) in self.outer.values.iter().enumerate() {
            for j in 0..inner_len {
                let k = i * inner_len + j;
                writeln!(w, " {point}\t{outer_v:.6e}")?;
                if let Some(inner) = &self.inner {
                    writeln!(w, "\t{:.6e}", inner.values[j])?;
                }
                for series in self.node_voltages.values() {
                    writeln!(w, "\t{:.6e}", series[k])?;
                }
                for series in self.vsrc_currents.values() {
                    writeln!(w, "\t{:.6e}", series[k])?;
                }
                point += 1;
            }
        }
        Ok(())
    }
}

/// Compute the linear grid of values for a DC sweep axis.
///
/// Inclusive of both endpoints when the step lands on `stop`; otherwise stops
/// strictly before `stop + 0.5·step`.
fn linspace(start: f64, stop: f64, step: f64) -> Vec<f64> {
    if step == 0.0 {
        return vec![start];
    }
    let sign = if stop >= start { 1.0 } else { -1.0 };
    let s = step.abs() * sign;
    let n = ((stop - start) / s).floor() as i64;
    let n = n.max(0) as usize;
    let mut v = Vec::with_capacity(n + 1);
    for i in 0..=n {
        v.push(start + i as f64 * s);
    }
    // If the final point exactly matches stop (within rounding), keep it.
    if let Some(last) = v.last() {
        if (last - stop).abs() > 1e-12 * stop.abs().max(step.abs()) && (stop - start).abs() > 0.0 {
            v.push(stop);
        }
    }
    v
}

/// Override a named voltage- or current-source's DC value in a netlist clone.
fn override_source_dc(netlist: &mut Netlist, src_lc: &str, value: f64) {
    for el in &mut netlist.elements {
        match el {
            Element::VoltageSource { name, waveform, .. } if name.eq_ignore_ascii_case(src_lc) => {
                *waveform = Waveform::Dc(value);
            }
            Element::CurrentSource { name, waveform, .. } if name.eq_ignore_ascii_case(src_lc) => {
                *waveform = Waveform::Dc(value);
            }
            _ => {}
        }
    }
}

/// Run a `.dc` sweep.
///
/// `src` is the swept source name (case-insensitive); `start`, `stop`, `step`
/// define the linear grid.  `nested` adds an inner-loop sweep run at every
/// outer-point value.
pub fn dc_sweep_with_registry_opts(
    netlist: &Netlist,
    src: &str,
    start: f64,
    stop: f64,
    step: f64,
    nested: Option<(&str, f64, f64, f64)>,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<DcSweepResult, SimError> {
    let outer_vals = linspace(start, stop, step);
    let inner = nested.map(|(name, a, b, s)| SweepAxis {
        src: name.to_lowercase(),
        values: linspace(a, b, s),
    });
    let inner_vals: Vec<f64> = inner.as_ref().map_or(vec![0.0], |i| i.values.clone());

    let total = outer_vals.len() * inner_vals.len();

    // We need the topology to allocate result vectors; just build it once from
    // the unmodified netlist (topology only depends on connectivity, not values).
    let probe = dc_op_nr_with_registry_opts(netlist, registry, opts)?;
    let topo = probe.topo;

    // Each sweep point is independent: clone the netlist, override the
    // swept source(s), run a fresh DC OP.  Solve all points in parallel
    // (rayon), then write into pre-sized result vectors by linear index
    // so the outer-major layout that `write_csv` / `write_nutmeg` expect
    // is preserved.  Under `--verbose` we fall back to serial — the
    // per-point NR diagnostics would otherwise interleave across threads.
    let inner_len = inner_vals.len();
    let mut points: Vec<(usize, f64, f64)> = Vec::with_capacity(total);
    for (i, &ov) in outer_vals.iter().enumerate() {
        for (j, &iv) in inner_vals.iter().enumerate() {
            points.push((i * inner_len + j, ov, iv));
        }
    }

    let solve_one = |outer_v: f64, inner_v: f64| -> Result<Vec<f64>, SimError> {
        let mut nl = netlist.clone();
        override_source_dc(&mut nl, src, outer_v);
        if let Some(axis) = &inner {
            override_source_dc(&mut nl, &axis.src, inner_v);
        }
        let r = dc_op_nr_with_registry_opts(&nl, registry, opts)?;
        Ok(r.x)
    };

    let solved: Vec<Result<(usize, Vec<f64>), SimError>> = if opts.verbose {
        points
            .iter()
            .map(|&(idx, ov, iv)| solve_one(ov, iv).map(|x| (idx, x)))
            .collect()
    } else {
        points
            .par_iter()
            .map(|&(idx, ov, iv)| solve_one(ov, iv).map(|x| (idx, x)))
            .collect()
    };
    let solved: Vec<(usize, Vec<f64>)> = solved.into_iter().collect::<Result<Vec<_>, _>>()?;

    let mut node_voltages: IndexMap<String, Vec<f64>> = topo
        .node_index
        .keys()
        .map(|k| (k.clone(), vec![0.0; total]))
        .collect();
    let mut vsrc_currents: IndexMap<String, Vec<f64>> = topo
        .vsrc_index
        .keys()
        .map(|k| (k.clone(), vec![0.0; total]))
        .collect();
    let n_nodes = topo.n_nodes();
    for (idx, x) in solved {
        for (name, &node_idx) in &topo.node_index {
            node_voltages.get_mut(name).unwrap()[idx] = x[node_idx];
        }
        for (name, &vsrc_idx) in &topo.vsrc_index {
            vsrc_currents.get_mut(name).unwrap()[idx] = x[n_nodes + vsrc_idx];
        }
    }

    Ok(DcSweepResult {
        outer: SweepAxis {
            src: src.to_lowercase(),
            values: outer_vals,
        },
        inner,
        node_voltages,
        vsrc_currents,
    })
}

/// Convenience: run a `.dc` sweep with default `SimOptions` (honoured `.options`).
pub fn dc_sweep_with_registry(
    netlist: &Netlist,
    src: &str,
    start: f64,
    stop: f64,
    step: f64,
    nested: Option<(&str, f64, f64, f64)>,
    registry: &DeviceRegistry,
) -> Result<DcSweepResult, SimError> {
    dc_sweep_with_registry_opts(
        netlist,
        src,
        start,
        stop,
        step,
        nested,
        registry,
        &SimOptions::from_netlist(netlist),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn linspace_basic() {
        let v = linspace(0.0, 1.0, 0.25);
        assert_eq!(v.len(), 5);
        assert!((v[0] - 0.0).abs() < 1e-12);
        assert!((v[4] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn linspace_does_not_overshoot() {
        let v = linspace(0.0, 1.0, 0.3);
        assert!(v.iter().all(|&x| x <= 1.0 + 1e-9));
        assert!((v.last().unwrap() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn dc_sweep_resistor_divider_linear() {
        // V1 → 1k → out → 1k → 0.  V(out) should be V1/2 at each sweep point.
        let net =
            parse_spice("* divider\nV1 in 0 DC 0\nR1 in out 1k\nR2 out 0 1k\n.dc V1 0 5 1\n.end\n")
                .unwrap();
        let reg = DeviceRegistry::new();
        let r = dc_sweep_with_registry(&net, "v1", 0.0, 5.0, 1.0, None, &reg).unwrap();
        assert_eq!(r.outer.values.len(), 6); // 0,1,2,3,4,5
        let v_out = r.node_voltages.get("out").unwrap();
        for (v_in, &v) in r.outer.values.iter().zip(v_out.iter()) {
            assert!((v - 0.5 * v_in).abs() < 1e-6, "V(in)={v_in} V(out)={v}");
        }
    }

    #[test]
    fn dc_sweep_diode_iv_monotonic() {
        // Diode I-V curve: V(b) > 0 and monotonic in V1.
        let net = parse_spice(
            "* diode IV\nV1 a 0 DC 0\nR1 a b 1k\nD1 b 0 myd\n\
             .model myd D (Is=1e-14 N=1)\n.dc V1 0 1 0.1\n.end\n",
        )
        .unwrap();
        let reg = {
            let mut r = DeviceRegistry::new();
            r.register_builtin_diodes(&net.models);
            r
        };
        let r = dc_sweep_with_registry(&net, "v1", 0.0, 1.0, 0.1, None, &reg).unwrap();
        let v_b = r.node_voltages.get("b").unwrap();
        for i in 1..v_b.len() {
            assert!(
                v_b[i] >= v_b[i - 1] - 1e-9,
                "V(b) not monotonic at step {i}: {} → {}",
                v_b[i - 1],
                v_b[i]
            );
        }
    }

    #[test]
    fn dc_sweep_csv_header_and_rows() {
        let net =
            parse_spice("* divider\nV1 in 0 DC 0\nR1 in out 1k\nR2 out 0 1k\n.end\n").unwrap();
        let reg = DeviceRegistry::new();
        let r = dc_sweep_with_registry(&net, "v1", 0.0, 1.0, 0.5, None, &reg).unwrap();
        let mut buf = Vec::new();
        r.write_csv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("v1,"), "header: {s}");
        assert!(s.contains("V(out)"), "{s}");
        // 3 sweep points → 3 data rows + 1 header.
        assert_eq!(s.lines().count(), 4, "{s}");
    }
}
