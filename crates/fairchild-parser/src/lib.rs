mod error;
pub mod expr;
mod spice;

pub use error::{DisciplineError, ParseError};
pub use expr::{Expr, EvalContext, ExprError};
pub use spice::{parse_spice, parse_spice_file};

/// A node name. "0" and "gnd" and "GND" all refer to ground.
pub type NodeName = String;

/// A parsed circuit netlist.
#[derive(Debug, Default, Clone)]
pub struct Netlist {
    pub title: String,
    pub elements: Vec<Element>,
    pub analyses: Vec<Analysis>,
    pub models: Vec<ModelCard>,
    /// Paths from `.osdi <path>` directives — OSDI shared libraries to load.
    pub osdi_paths: Vec<String>,
    /// Net names declared as optical via `.optical <net> ...` directive.
    pub optical_nets: Vec<String>,
    /// Raw `KEY=VALUE` pairs from every `.options` directive in source order.
    /// Values are stored as strings so the consumer (typically `SimOptions::set`)
    /// can parse SPICE suffixes and method names appropriately.
    pub options: Vec<(String, String)>,
    /// `.measure` directives — evaluated by the consumer (CLI / Python) on a
    /// completed TranResult / SimResult.
    pub measurements: Vec<Measurement>,
    /// `.ic V(node)=value …`: initial conditions for transient analysis.
    /// Applied at t=0 *only when* `UIC` is enabled (via `.tran … UIC` or
    /// `.options uic=1`).  Otherwise the DC operating point is used as t=0.
    pub ic: Vec<(NodeName, f64)>,
    /// `.nodeset V(node)=value …`: initial guess for the DC operating point
    /// Newton-Raphson.  Always applied if present.
    pub nodeset: Vec<(NodeName, f64)>,
    /// `.temp <T1> [<T2> …]`: simulation temperatures in Kelvin (converted
    /// from Celsius).  Empty ⇒ use `SimOptions::temp_k` default (300.15 K).
    /// More than one entry asks the CLI/Python driver to repeat every
    /// analysis once per temperature.
    pub temps: Vec<f64>,
    /// `.alter <label>` blocks: each is a set of element / model overrides
    /// applied on top of the base netlist for a re-run pass.  The base run
    /// uses the original netlist; each block runs once after applying.
    pub alters: Vec<AlterBlock>,
    /// `.optical_port NAME [N]` declarations.  Each entry is a single
    /// user-visible port that expands to 3·N underlying wires (re, im, λ
    /// per channel).  Used by the parser's X-line preprocessor.
    pub optical_ports: Vec<OpticalPort>,
}

/// A bundle optical port declared via `.optical_port NAME [N]`.
///
/// A reference to `NAME` in an X-element net list expands to N copies of the
/// 3-wire `(NAME_re_i, NAME_im_i, NAME_wl_i)` tuple; the device instance is
/// replicated once per channel when `channels > 1`.
#[derive(Debug, Clone)]
pub struct OpticalPort {
    pub name: String,
    pub channels: usize,
}

impl OpticalPort {
    /// Underlying wire names for channel `ch` of this port.
    pub fn wires_for_channel(&self, ch: usize) -> [String; 3] {
        [
            format!("{}_re_{}", self.name, ch),
            format!("{}_im_{}", self.name, ch),
            format!("{}_wl_{}", self.name, ch),
        ]
    }

    /// Every underlying wire across every channel.  Used to register the
    /// port's wires as optical nets for discipline-check purposes.
    pub fn all_wires(&self) -> Vec<String> {
        let mut out = Vec::with_capacity(self.channels * 3);
        for ch in 0..self.channels {
            out.extend(self.wires_for_channel(ch));
        }
        out
    }
}

/// One `.alter <label>` block: name-keyed patches applied to the base netlist
/// before a re-run pass.  `element_overrides` replace by name (or append when
/// the name doesn't exist); `model_overrides` follow the same rule.
#[derive(Debug, Clone, Default)]
pub struct AlterBlock {
    pub label: String,
    pub element_overrides: Vec<Element>,
    pub model_overrides:  Vec<ModelCard>,
}

impl Netlist {
    /// Apply an `.alter` block in place: replace any same-named Element /
    /// ModelCard, or append a new entry if the name is unseen.
    pub fn apply_alter(&mut self, block: &AlterBlock) {
        for new_el in &block.element_overrides {
            if let Some(slot) = self.elements.iter_mut().find(|e| element_name(e) == element_name(new_el)) {
                *slot = new_el.clone();
            } else {
                self.elements.push(new_el.clone());
            }
        }
        for new_m in &block.model_overrides {
            if let Some(slot) = self.models.iter_mut().find(|m| m.name == new_m.name) {
                *slot = new_m.clone();
            } else {
                self.models.push(new_m.clone());
            }
        }
    }
}

/// Extract the name (case-folded) from any Element variant, for alter matching.
fn element_name(el: &Element) -> String {
    match el {
        Element::Resistor      { name, .. } |
        Element::Capacitor     { name, .. } |
        Element::Inductor      { name, .. } |
        Element::VoltageSource { name, .. } |
        Element::CurrentSource { name, .. } |
        Element::Diode         { name, .. } |
        Element::Mosfet        { name, .. } |
        Element::XOsdi         { name, .. } |
        Element::Behavioral    { name, .. } => name.to_lowercase(),
    }
}

/// Waveform specification for independent sources.
#[derive(Debug, Clone)]
pub enum Waveform {
    Dc(f64),
    /// PULSE(v0 v1 td tr tf pw per)
    Pulse {
        v0: f64,
        v1: f64,
        td: f64,
        tr: f64,
        tf: f64,
        pw: f64,
        per: f64,
    },
    /// PWL(t0 v0 t1 v1 ...) — piecewise-linear; points must be sorted by time.
    Pwl {
        points: Vec<(f64, f64)>,
    },
    /// SIN(vo va freq td theta phase) — damped sinusoid.
    ///   v(t) = vo + va·exp(−(t−td)·theta)·sin(2π·freq·(t−td) + phase)  for t ≥ td
    ///   v(t) = vo                                                       for t <  td
    /// `phase` is in radians.
    Sin {
        vo:    f64,
        va:    f64,
        freq:  f64,
        td:    f64,
        theta: f64,
        phase: f64,
    },
    /// EXP(v1 v2 td1 tau1 td2 tau2) — rising/falling exponentials.
    ///   v(t) = v1 for t < td1
    ///        = v1 + (v2−v1)·(1−exp(−(t−td1)/tau1))                                          for td1 ≤ t < td2
    ///        = v1 + (v2−v1)·(1−exp(−(t−td1)/tau1)) + (v1−v2)·(1−exp(−(t−td2)/tau2))         for t ≥ td2
    Exp {
        v1:   f64,
        v2:   f64,
        td1:  f64,
        tau1: f64,
        td2:  f64,
        tau2: f64,
    },
    /// SFFM(vo va fc mdi fs) — single-frequency FM.
    ///   v(t) = vo + va·sin(2π·fc·t + mdi·sin(2π·fs·t))
    Sffm {
        vo:  f64,
        va:  f64,
        fc:  f64,
        mdi: f64,
        fs:  f64,
    },
    /// AM(va vo mf fc td) — amplitude-modulated sinusoid (ngspice form).
    ///   v(t) = vo                                                         for t < td
    ///   v(t) = va·sin(2π·mf·(t−td))·sin(2π·fc·(t−td)) + vo                for t ≥ td
    Am {
        va: f64,
        vo: f64,
        mf: f64,
        fc: f64,
        td: f64,
    },
}

impl Waveform {
    /// Value used for DC operating-point (t = 0).
    pub fn dc_value(&self) -> f64 {
        match self {
            Waveform::Dc(v) => *v,
            Waveform::Pulse { v0, .. } => *v0,
            Waveform::Pwl { points } => points.first().map(|(_, v)| *v).unwrap_or(0.0),
            // All continuous shapes use their value at t=0 as the DC point.
            // For SIN with td>0, this is just vo; for EXP, v1; etc.
            Waveform::Sin  { vo, .. } => *vo,
            Waveform::Exp  { v1, .. } => *v1,
            Waveform::Sffm { vo, .. } => *vo,  // sin(0)=0 → vo dominates
            Waveform::Am   { vo, .. } => *vo,
        }
    }

    /// Next time strictly after `t` at which this waveform has a slope discontinuity.
    ///
    /// Returns `None` for smooth (DC) waveforms or when all breakpoints are in the past.
    pub fn next_breakpoint(&self, t: f64) -> Option<f64> {
        match self {
            Waveform::Dc(_) => None,
            Waveform::Pulse { td, tr, tf, pw, per, .. } => {
                if t < *td {
                    return Some(*td);
                }
                // Offsets from the start of a period where slope changes.
                let offsets: [f64; 4] = [0.0, *tr, tr + pw, tr + pw + tf];
                if *per <= 0.0 {
                    return offsets.iter()
                        .map(|b| td + b)
                        .filter(|&bp| bp > t)
                        .reduce(f64::min);
                }
                let phase = (t - td) % per;
                let base = t - phase; // start of current period
                offsets.iter()
                    .map(|b| base + b)
                    .chain(std::iter::once(base + per))
                    .filter(|&bp| bp > t)
                    .reduce(f64::min)
            }
            Waveform::Pwl { points } => {
                points.iter().map(|(pt, _)| *pt).find(|&pt| pt > t)
            }
            Waveform::Sin { td, .. } | Waveform::Am { td, .. } => {
                // Slope discontinuity at td (constant → oscillatory).
                if t < *td { Some(*td) } else { None }
            }
            Waveform::Exp { td1, td2, .. } => {
                if t < *td1 { Some(*td1) }
                else if t < *td2 { Some(*td2) }
                else { None }
            }
            // SFFM is smooth everywhere; let the LTE controller pick steps.
            Waveform::Sffm { .. } => None,
        }
    }

    /// Value at time t (seconds).
    pub fn at(&self, t: f64) -> f64 {
        match self {
            Waveform::Dc(v) => *v,
            Waveform::Pulse { v0, v1, td, tr, tf, pw, per } => {
                if t < *td {
                    return *v0;
                }
                // Time within the current period.
                let tp = if *per > 0.0 { (t - td) % per } else { t - td };
                if tp < *tr {
                    v0 + (v1 - v0) * tp / tr
                } else if tp < tr + pw {
                    *v1
                } else if tp < tr + pw + tf {
                    v1 + (v0 - v1) * (tp - tr - pw) / tf
                } else {
                    *v0
                }
            }
            Waveform::Pwl { points } => {
                if points.is_empty() {
                    return 0.0;
                }
                if t <= points[0].0 {
                    return points[0].1;
                }
                if t >= points[points.len() - 1].0 {
                    return points[points.len() - 1].1;
                }
                // Binary search for the segment containing t.
                let idx = points.partition_point(|(pt, _)| *pt <= t);
                let (t0, v0) = points[idx - 1];
                let (t1, v1) = points[idx];
                let frac = (t - t0) / (t1 - t0);
                v0 + (v1 - v0) * frac
            }
            Waveform::Sin { vo, va, freq, td, theta, phase } => {
                if t < *td {
                    *vo
                } else {
                    let dt = t - td;
                    let damp = if *theta > 0.0 { (-theta * dt).exp() } else { 1.0 };
                    vo + va * damp * (2.0 * std::f64::consts::PI * freq * dt + phase).sin()
                }
            }
            Waveform::Exp { v1, v2, td1, tau1, td2, tau2 } => {
                if t < *td1 {
                    *v1
                } else {
                    let rise = (v2 - v1) * (1.0 - (-(t - td1) / tau1).exp());
                    if t < *td2 {
                        v1 + rise
                    } else {
                        let fall = (v1 - v2) * (1.0 - (-(t - td2) / tau2).exp());
                        v1 + rise + fall
                    }
                }
            }
            Waveform::Sffm { vo, va, fc, mdi, fs } => {
                let mod_arg = 2.0 * std::f64::consts::PI * fs * t;
                vo + va * (2.0 * std::f64::consts::PI * fc * t + mdi * mod_arg.sin()).sin()
            }
            Waveform::Am { va, vo, mf, fc, td } => {
                if t < *td {
                    *vo
                } else {
                    let dt = t - td;
                    let m_lo = (2.0 * std::f64::consts::PI * mf * dt).sin();
                    let m_hi = (2.0 * std::f64::consts::PI * fc * dt).sin();
                    va * m_lo * m_hi + vo
                }
            }
        }
    }
}

/// A single circuit element.
#[derive(Debug, Clone)]
pub enum Element {
    Resistor {
        name: String,
        pos: NodeName,
        neg: NodeName,
        resistance: f64,
    },
    Capacitor {
        name: String,
        pos: NodeName,
        neg: NodeName,
        capacitance: f64,
    },
    Inductor {
        name: String,
        pos: NodeName,
        neg: NodeName,
        inductance: f64,
    },
    VoltageSource {
        name: String,
        pos: NodeName,
        neg: NodeName,
        waveform: Waveform,
    },
    CurrentSource {
        name: String,
        pos: NodeName,
        neg: NodeName,
        waveform: Waveform,
    },
    Diode {
        name: String,
        anode: NodeName,
        cathode: NodeName,
        model_name: String,
    },
    Mosfet {
        name: String,
        drain: NodeName,
        gate: NodeName,
        source: NodeName,
        bulk: NodeName,
        model_name: String,
        params: Vec<(String, f64)>,
    },
    /// Generic OSDI instance: `X<name> <net0> <net1> ... <model_name> [param=value ...]`
    /// Port order matches terminal order in the OSDI descriptor.
    XOsdi {
        name: String,
        nets: Vec<NodeName>,
        model_name: String,
        params: Vec<(String, f64)>,
    },
    /// Behavioural source: `B<name> n+ n- V=<expr>` or `I=<expr>`.
    ///
    /// `V=` form contributes a voltage between (pos, neg) equal to `expr(x)`
    /// (an MNA aux row, like an ordinary voltage source).  `I=` form drives
    /// a current of `expr(x)` from pos→neg.
    Behavioral {
        name: String,
        pos:  NodeName,
        neg:  NodeName,
        kind: BehavioralKind,
        expr: expr::Expr,
    },
}

/// Whether a B-element is a voltage or current source.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BehavioralKind { Voltage, Current }

/// A `.measure tran|dc|ac NAME …` directive.
#[derive(Debug, Clone)]
pub struct Measurement {
    /// User-visible name (e.g. "tpd", "vmax").
    pub name: String,
    /// Analysis context — currently only `Tran` is honoured by the post-processor.
    pub analysis: MeasAnalysis,
    pub kind: MeasKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeasAnalysis { Tran, Dc, Ac }

/// What to compute.
#[derive(Debug, Clone)]
pub enum MeasKind {
    /// FIND expr AT=t   — evaluate `expr` at time `t`.
    FindAt    { expr: expr::Expr, at: f64 },
    /// FIND expr WHEN cond — first time `cond` is true; report `expr` there.
    FindWhen  { expr: expr::Expr, cond: expr::Expr, cross: usize },
    /// Aggregate operation over [from, to].
    Aggregate { op: MeasOp, expr: expr::Expr, from: Option<f64>, to: Option<f64> },
    /// DERIV expr AT=t — numerical derivative at t.
    DerivAt   { expr: expr::Expr, at: f64 },
    /// TRIG cond1 [val=v1] [cross=n] TARG cond2 [val=v2] [cross=n] — delay measurement.
    TrigTarg  {
        trig_expr: expr::Expr, trig_val: f64, trig_cross: usize,
        targ_expr: expr::Expr, targ_val: f64, targ_cross: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeasOp { Max, Min, Avg, Rms, Pp, Integ }

/// A model card parsed from `.model <name> <kind> [param=value ...]`.
#[derive(Debug, Clone)]
pub struct ModelCard {
    pub name: String,
    pub kind: String,
    pub params: Vec<(String, f64)>,
}

/// Frequency spacing for AC sweep.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AcVariation {
    /// Points per decade (logarithmic).
    Dec,
    /// Points per octave (logarithmic).
    Oct,
    /// Total points (linear).
    Lin,
}

/// Check for optical↔electrical discipline mismatches in the netlist.
///
/// Any net declared via `.optical` is flagged if it is also connected to a
/// purely-electrical element (R, L, C, V, I, D, M).  Mixed-domain elements
/// (XOsdi) may legitimately connect optical nets to electrical ones and are
/// not checked here.
///
/// Returns the first mismatch found, or `Ok(())` if the netlist is clean.
pub fn check_disciplines(netlist: &Netlist) -> Result<(), DisciplineError> {
    use std::collections::HashSet;

    let optical: HashSet<&str> = netlist.optical_nets.iter().map(|s| s.as_str()).collect();
    if optical.is_empty() {
        return Ok(());
    }

    let check = |element_name: &str, net: &str| -> Result<(), DisciplineError> {
        if optical.contains(net) {
            Err(DisciplineError {
                element: element_name.to_string(),
                net: net.to_string(),
            })
        } else {
            Ok(())
        }
    };

    for el in &netlist.elements {
        match el {
            Element::Resistor  { name, pos, neg, .. } => { check(name, pos)?; check(name, neg)?; }
            Element::Capacitor { name, pos, neg, .. } => { check(name, pos)?; check(name, neg)?; }
            Element::Inductor  { name, pos, neg, .. } => { check(name, pos)?; check(name, neg)?; }
            Element::VoltageSource { name, pos, neg, .. } => { check(name, pos)?; check(name, neg)?; }
            Element::CurrentSource { name, pos, neg, .. } => { check(name, pos)?; check(name, neg)?; }
            Element::Diode  { name, anode, cathode, .. } => { check(name, anode)?; check(name, cathode)?; }
            Element::Mosfet { name, drain, gate, source, bulk, .. } => {
                check(name, drain)?; check(name, gate)?;
                check(name, source)?; check(name, bulk)?;
            }
            // XOsdi is intentionally not checked: mixed-domain connections are valid.
            Element::XOsdi { .. } => {}
            Element::Behavioral { name, pos, neg, .. } => { check(name, pos)?; check(name, neg)?; }
        }
    }
    Ok(())
}

/// A requested simulation analysis.
#[derive(Debug, Clone)]
pub enum Analysis {
    Op,
    Tran { step: f64, stop: f64 },
    /// `.ac DEC|OCT|LIN <points> <fstart> <fstop>`
    Ac { variation: AcVariation, points: usize, fstart: f64, fstop: f64 },
    /// `.dc SRC START STOP STEP [SRC2 START2 STOP2 STEP2]`
    ///
    /// First sweep is the outer loop; optional nested sweep is the inner loop.
    /// Source names are stored lowercased to match element naming.
    Dc {
        src: String,
        start: f64,
        stop: f64,
        step: f64,
        nested: Option<DcSweepSpec>,
    },
    /// `.noise V(<out_node>[,<ref_node>]) <input_src> DEC|OCT|LIN <pts> <fstart> <fstop>`
    ///
    /// Small-signal noise sweep.  Reports output noise PSD at the observation
    /// node and input-referred PSD through the named excitation source.
    Noise {
        out_pos: String,
        out_neg: String,           // "0" if user writes `V(out)`
        input_src: String,
        variation: AcVariation,
        points: usize,
        fstart: f64,
        fstop: f64,
    },
}

/// One leg of a `.dc` sweep (the nested form's inner specification).
#[derive(Debug, Clone)]
pub struct DcSweepSpec {
    pub src: String,
    pub start: f64,
    pub stop: f64,
    pub step: f64,
}
