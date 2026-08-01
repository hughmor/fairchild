//! Small-signal noise analysis.
//!
//! For each swept frequency:
//!   1. Build the same complex linearized system A(jω) = G + jB used by AC.
//!   2. Solve the adjoint problem  A^T · λ = e_out  to obtain the transfer
//!      impedance from any internal current injection to the observation
//!      node V(out_pos) − V(out_neg).
//!   3. For every uncorrelated noise source k at nodes (p,n) with one-sided
//!      current PSD S_ik(f) [A²/Hz]:
//!      S_V_out_k(f) = |λ[p] − λ[n]|² · S_ik(f).
//!   4. Sum the contributions and divide by the squared signal-path gain to
//!      get the input-referred PSD.
//!
//! This commit ships the full plumbing (parser → CLI/Python) plus resistor
//! thermal noise (4kT/R).  Diode shot noise and MOSFET channel noise land in
//! the next commit once the Device trait grows a `noise_sources()` hook.

use indexmap::IndexMap;

use fairchild_parser::{Element, Netlist};

use crate::device::{Device, EvalFlags, ReactiveKind};
use crate::device_registry::DeviceRegistry;
use crate::error::SimError;
use crate::mna::{stamp_2port_by_id, stamp_netlist_scaled, stamp_passive_2port, CircuitTopology};
use crate::newton::build_devices;
use crate::options::SimOptions;
use crate::solver::LinearSolver;

/// Boltzmann's constant in J/K.
const KB: f64 = 1.380649e-23;

/// Result of a `.noise` sweep.  All PSDs are one-sided (V²/Hz).
pub struct NoiseResult {
    pub freq: Vec<f64>,
    /// Total output-referred voltage noise PSD at each frequency (V²/Hz).
    pub onoise_psd: Vec<f64>,
    /// Input-referred PSD = onoise / |H(f)|² where H is the small-signal gain
    /// from `input_src` to (out_pos, out_neg).  NaN at frequencies where the
    /// transfer function is too small to invert reliably.
    pub inoise_psd: Vec<f64>,
}

impl NoiseResult {
    /// Write a CSV with columns `freq_hz, onoise_v2hz, onoise_vrthz, inoise_v2hz, inoise_vrthz`.
    pub fn write_csv<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        writeln!(
            w,
            "freq_hz,onoise_v2hz,onoise_vrthz,inoise_v2hz,inoise_vrthz"
        )?;
        for i in 0..self.freq.len() {
            let on = self.onoise_psd[i];
            let in_ = self.inoise_psd[i];
            writeln!(
                w,
                "{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
                self.freq[i],
                on,
                on.max(0.0).sqrt(),
                in_,
                in_.max(0.0).sqrt()
            )?;
        }
        Ok(())
    }
}

/// Run a `.noise` analysis.
///
/// `freqs` — already-expanded sweep points (use the same `freq_decade` /
///   `freq_linear` / `freq_oct` helpers as `ac::ac_analysis`).
/// `out_pos`, `out_neg` — observation node pair (`"0"` for ground).
/// `input_src` — name of the voltage source that defines the input port; used
///   only for input-referred PSD via the signal-path gain.
pub fn noise_analysis(
    netlist: &Netlist,
    freqs: &[f64],
    out_pos: &str,
    out_neg: &str,
    input_src: &str,
    registry: &DeviceRegistry,
    opts: &SimOptions,
) -> Result<NoiseResult, SimError> {
    crate::connectivity::check_connectivity(netlist)?;
    let ctx = opts.sim_context();
    let mut topo = CircuitTopology::build(netlist);
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();

    // DC operating point.
    let mut devices = build_devices(netlist, &mut topo, &ctx, registry)?;
    let n_nodes = topo.n_nodes();
    let size = topo.size;
    let dc_solver = opts.linear_solver(size);
    let noise_solver = opts.linear_solver(2 * size);
    let x0 = run_dc_op(&topo, netlist, &mut devices, &ctx, opts, &*dc_solver)?;

    // Real (G) and imaginary-coefficient (C, L⁻¹) parts of Y(jω) = G + j(ωC − L⁻¹/ω).
    let mut g_mat = stamp_netlist_scaled(&topo, netlist, 1.0, &empty, &empty).a;
    for dev in devices.iter_mut() {
        // tran() flags so device small-signal cap caches populate (see ac.rs).
        dev.eval(&x0, EvalFlags::tran(), &ctx);
        let mut tmp = crate::mna::MnaMatrix::zeros(size);
        dev.load_jacobian(&mut tmp);
        // `tmp.a` is sparse now, so this walks only the cells the device
        // actually stamped instead of the whole row.
        for (g_row, t_row) in g_mat.iter_mut().zip(tmp.a.iter()) {
            for (j, t) in t_row.iter() {
                g_row[j] += t;
            }
        }
    }
    topo.stamp_gmin(&mut g_mat, opts.gmin);
    // ponytail: dense G/C/L and a dense 2n×2n adjoint system, same trade-off
    // and same upgrade path as `ac.rs`. Tracked as task #12.
    let mut c_mat = vec![vec![0.0f64; size]; size];
    let mut l_mat = vec![vec![0.0f64; size]; size];
    for el in &netlist.elements {
        match el {
            Element::Capacitor {
                pos,
                neg,
                capacitance,
                ..
            } => stamp_passive_2port(&mut c_mat, &topo.node_index, pos, neg, *capacitance),
            Element::Inductor {
                pos,
                neg,
                inductance,
                ..
            } => stamp_passive_2port(&mut l_mat, &topo.node_index, pos, neg, 1.0 / *inductance),
            _ => {}
        }
    }
    // Device-internal small-signal reactances (diode Cj, MOSFET caps, photonic
    // parasitics) — previously absent from noise; see ac.rs for the rationale.
    for dev in devices.iter() {
        for r in dev.small_signal_reactances() {
            match r.kind {
                ReactiveKind::Capacitor => stamp_2port_by_id(&mut c_mat, r.pos, r.neg, r.value),
                ReactiveKind::Inductor if r.value != 0.0 => {
                    stamp_2port_by_id(&mut l_mat, r.pos, r.neg, 1.0 / r.value)
                }
                ReactiveKind::Inductor => {}
            }
        }
        // Devices whose reactance is a general ∂q/∂x matrix rather than a set
        // of two-terminal branches (OSDI/Verilog-A) stamp it themselves.
        dev.load_reactive_jacobian(&mut c_mat);
    }

    // Locate the named input source so we can compute its signal-path gain.
    let input_vsrc_idx = topo.vsrc_index.get(input_src).copied().ok_or_else(|| {
        SimError::ParameterError(format!("noise: input source '{input_src}' not found"))
    })?;

    let out_pos_idx = if out_pos == "0" {
        None
    } else {
        topo.node_index.get(out_pos).copied()
    };
    let out_neg_idx = if out_neg == "0" {
        None
    } else {
        topo.node_index.get(out_neg).copied()
    };
    if out_pos == "0" && out_neg == "0" {
        return Err(SimError::ParameterError(
            "noise: output node cannot be ground".into(),
        ));
    }

    let temp_k = opts.temp_k;
    let four_kt = 4.0 * KB * temp_k;

    let mut result = NoiseResult {
        freq: freqs.to_vec(),
        onoise_psd: Vec::with_capacity(freqs.len()),
        inoise_psd: Vec::with_capacity(freqs.len()),
    };

    for &f in freqs {
        let omega = 2.0 * std::f64::consts::PI * f;

        // B = ωC − L⁻¹/ω (purely real coefficient matrix; imaginary part of Y).
        let mut b_mat = vec![vec![0.0f64; size]; size];
        for i in 0..size {
            for j in 0..size {
                b_mat[i][j] = omega * c_mat[i][j] - l_mat[i][j] / omega;
            }
        }

        // Forward 2N×2N block system for the signal-path gain H(f):
        //   [ G  -B ] [V_re]   [b_re]
        //   [ B   G ] [V_im] = [b_im]
        // where the RHS injects 1 V at the input source's branch row.
        let n2 = 2 * size;
        let mut a_fwd = vec![vec![0.0f64; n2]; n2];
        for i in 0..size {
            for j in 0..size {
                a_fwd[i][j] = g_mat[i][j];
                a_fwd[i][size + j] = -b_mat[i][j];
                a_fwd[size + i][j] = b_mat[i][j];
                a_fwd[size + i][size + j] = g_mat[i][j];
            }
        }
        let mut rhs_fwd = vec![0.0f64; n2];
        rhs_fwd[n_nodes + input_vsrc_idx] = 1.0; // unit AC amplitude on V source
        let v_fwd = noise_solver.solve(&CircuitTopology::sparse_from_dense(&a_fwd), &rhs_fwd)?;
        let v_re_fwd = &v_fwd[..size];
        let v_im_fwd = &v_fwd[size..];
        let h_re = pick(v_re_fwd, out_pos_idx) - pick(v_re_fwd, out_neg_idx);
        let h_im = pick(v_im_fwd, out_pos_idx) - pick(v_im_fwd, out_neg_idx);
        let h_mag_sq = h_re * h_re + h_im * h_im;

        // Adjoint solve.  Real-block transpose of [G -B; B G] is
        //   [ G^T   B^T ]
        //   [-B^T   G^T ]
        // which, by the advisor's note, represents M^H rather than M^T.  For
        // |λ[p]−λ[n]|² in the noise sum this is fine — conjugation preserves
        // magnitude — but we MUST NOT compare λ values directly to an A^T λ
        // solution from a different convention.
        let mut a_adj = vec![vec![0.0f64; n2]; n2];
        for i in 0..size {
            for j in 0..size {
                a_adj[i][j] = g_mat[j][i];
                a_adj[i][size + j] = b_mat[j][i];
                a_adj[size + i][j] = -b_mat[j][i];
                a_adj[size + i][size + j] = g_mat[j][i];
            }
        }
        let mut rhs_adj = vec![0.0f64; n2];
        // e_out: +1 at out_pos, -1 at out_neg, both in the real block (imag
        // RHS is zero — we observe a real-valued node voltage).
        if let Some(i) = out_pos_idx {
            rhs_adj[i] += 1.0;
        }
        if let Some(i) = out_neg_idx {
            rhs_adj[i] -= 1.0;
        }
        let lam = noise_solver.solve(&CircuitTopology::sparse_from_dense(&a_adj), &rhs_adj)?;
        let lam_re = &lam[..size];
        let lam_im = &lam[size..];

        // Resistor thermal noise: 4kT/R between (pos, neg) for every linear R.
        let mut s_v_out = 0.0_f64;
        for el in &netlist.elements {
            if let Element::Resistor {
                pos,
                neg,
                resistance,
                ..
            } = el
            {
                if *resistance <= 0.0 {
                    continue;
                }
                let s_i = four_kt / resistance; // 4kT/R [A²/Hz]
                let p_idx = topo.node_index.get(pos).copied();
                let n_idx = topo.node_index.get(neg).copied();
                let z_re = pick(lam_re, p_idx) - pick(lam_re, n_idx);
                let z_im = pick(lam_im, p_idx) - pick(lam_im, n_idx);
                let z_mag_sq = z_re * z_re + z_im * z_im;
                s_v_out += s_i * z_mag_sq;
            }
        }
        // Device-internal noise (diode shot, MOSFET channel thermal, …).
        // Each device contributes one or more uncorrelated current sources
        // between specific terminal indices; magnitude squared of the
        // transfer impedance picks up the output PSD contribution.
        for dev in devices.iter() {
            for (p_idx, n_idx, s_i) in dev.noise_sources(&ctx) {
                let z_re = pick(lam_re, p_idx) - pick(lam_re, n_idx);
                let z_im = pick(lam_im, p_idx) - pick(lam_im, n_idx);
                let z_mag_sq = z_re * z_re + z_im * z_im;
                s_v_out += s_i * z_mag_sq;
            }
        }

        let s_v_in = if h_mag_sq > 1e-30 {
            s_v_out / h_mag_sq
        } else {
            f64::NAN
        };
        result.onoise_psd.push(s_v_out);
        result.inoise_psd.push(s_v_in);
    }

    Ok(result)
}

/// Pick a node value from a vector indexed by `NodeId`-style Option<usize>.
fn pick(v: &[f64], idx: Option<usize>) -> f64 {
    match idx {
        Some(i) => v[i],
        None => 0.0,
    }
}

/// DC OP shared with ac.rs (kept private here to avoid coupling the modules).
fn run_dc_op(
    topo: &CircuitTopology,
    netlist: &Netlist,
    devices: &mut [Box<dyn Device>],
    ctx: &crate::device::SimContext,
    opts: &SimOptions,
    solver: &dyn LinearSolver,
) -> Result<Vec<f64>, SimError> {
    let empty: IndexMap<String, (f64, f64)> = IndexMap::new();
    let n_nodes = topo.n_nodes();
    // Not every unknown is a volt — see `crate::tolerance`.
    let tol = crate::tolerance::Tolerances::build(netlist, topo, opts);
    let mut x = vec![0.0f64; topo.size];

    for _ in 0..opts.itl1 {
        let mut mat = stamp_netlist_scaled(topo, netlist, 1.0, &empty, &empty);
        for dev in devices.iter_mut() {
            dev.eval(&x, EvalFlags::dc(), ctx);
            dev.load_residual(&mut mat.b);
            dev.load_jacobian(&mut mat);
        }
        topo.stamp_gmin(&mut mat.a, opts.gmin);
        let x_new = solver.solve(&mat.a, &mat.b)?;
        let max_dv = x_new
            .iter()
            .zip(x.iter())
            .take(n_nodes)
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);
        let x_next: Vec<f64> = if max_dv > opts.vmax {
            let scale = opts.vmax / max_dv;
            x.iter()
                .zip(x_new.iter())
                .map(|(o, n)| o + scale * (n - o))
                .collect()
        } else {
            x_new
        };
        let converged = tol.converged(&x_next, &x);
        x = x_next;
        if converged {
            return Ok(x);
        }
    }
    Err(SimError::NoConvergence { iters: opts.itl1 })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    /// Resistor-divider thermal noise.  At any frequency:
    ///   S_V_out = 4kT·(R1 ‖ R2)
    /// for a 1-V AC input driving R1→out and R2→0.  Independent of frequency
    /// in this purely resistive circuit.
    #[test]
    fn resistor_divider_thermal_noise() {
        // The current parser doesn't have an `AC` source tag; AC analysis
        // treats every voltage source with a nonzero DC value as the
        // excitation.  noise() routes its own RHS injection by name, so the
        // DC value here is irrelevant for the test.
        let net = parse_spice(
            "* thermal\nV1 in 0 DC 1\nR1 in out 1k\nR2 out 0 1k\n\
             .noise V(out) V1 DEC 1 1k 1k\n.end\n",
        )
        .unwrap();
        let mut registry = crate::device_registry::DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let opts = SimOptions::default();
        let r = noise_analysis(&net, &[1e3], "out", "0", "v1", &registry, &opts).unwrap();

        let kt = KB * opts.temp_k;
        let expected = 4.0 * kt * 500.0; // 4kT·(R1||R2)
        let s = r.onoise_psd[0];
        let rel = (s - expected).abs() / expected;
        assert!(
            rel < 0.01,
            "S_V_out={s:.3e} expected={expected:.3e} rel={rel:.3e}"
        );
    }

    /// Diode shot noise.  A current source biases a diode to ~1 mA so its
    /// internal Id is set by the source, not by an exponential.  At low f
    /// (cap negligible) the diode's small-signal resistance r_d = V_T / Id
    /// terminates the noise current into the output node.  The expected
    /// output PSD is:
    ///     S_V_out = (R_load || r_d)² · (2qId + 4kT/R_load · (r_d/(R_load+r_d))²)
    /// We just verify shot dominates: S_V_out > 4kT/R_load · z² and the
    /// magnitude is roughly 2qId · (R_load || r_d)² within 30%.
    #[test]
    fn diode_shot_noise_dominates_at_high_bias() {
        let src = "* diode shot\n\
                   Vbias bias 0 DC 1\n\
                   Ib   0 b  1m\n\
                   D1   b 0  myd\n\
                   .model myd D (Is=1e-14 N=1)\n\
                   .noise V(b) Vbias DEC 1 1k 1k\n.end\n";
        let net = parse_spice(src).unwrap();
        let mut registry = crate::device_registry::DeviceRegistry::new();
        registry.register_builtin_models(&net.models);
        let opts = SimOptions::default();
        let r = noise_analysis(&net, &[1e3], "b", "0", "vbias", &registry, &opts).unwrap();

        // The output is the diode anode; its small-signal resistance is the
        // only impedance to ground at that node aside from gmin.  λ at the
        // diode terminals therefore ≈ r_d (= V_T / Id ≈ 25.85 Ω at 27 °C,
        // Id=1mA), so output PSD ≈ 2qId · r_d².
        const Q: f64 = 1.602176634e-19;
        let id = 1e-3;
        let vt = KB * opts.temp_k / Q;
        let r_d = vt / id;
        let expected = 2.0 * Q * id * r_d * r_d;
        let s = r.onoise_psd[0];
        assert!(
            (s - expected).abs() / expected < 0.1,
            "diode shot S_V_out={s:.3e} expected≈{expected:.3e}"
        );
    }
}
