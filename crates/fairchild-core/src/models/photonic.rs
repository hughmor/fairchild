//! Native Rust photonic passive devices (B3).
//!
//! Each device implements `Device` directly — no Verilog-A round-trip, no
//! OSDI shared library, no Norton hack.  Outputs are stamped via the
//! direct-potential pattern (one auxiliary MNA row per output potential
//! equation, requested through `Device::num_extra_nodes` and bound via
//! `bind_extra_nodes`).
//!
//! Port convention (matches B1 discipline scheme, single-channel
//! forward-only): each optical port is a 3-wire bundle (re, im, λ).  The
//! .optical_port bundle directive (B2) maps a user-visible port name to
//! the underlying wires.
//!
//! All models in this module are equivalent-physics implementations from
//! public textbook formulas (Saleh-Teich, Heuck-Englund, Pozar) and carry
//! no PDK-specific calibration.

use crate::device::{Device, EvalFlags, NodeId, SimContext};
use crate::mna::MnaMatrix;

// ────────────────────────────────────────────────────────────────────────
// Native straight waveguide
// ────────────────────────────────────────────────────────────────────────

/// Straight optical waveguide — propagation loss + accumulated phase.
///
/// Physics: `A_out = A_in · exp(-α·L/2) · exp(-j·β·L)` with `β = 2π·n_g/λ`.
///
/// Implementation: 3 direct-potential equations enforce
///   V(out_re) = T·(cos(φ)·V(in_re) + sin(φ)·V(in_im))
///   V(out_im) = T·(-sin(φ)·V(in_re) + cos(φ)·V(in_im))
///   V(out_λ)  = V(in_λ)
/// stamped via three auxiliary branch rows reserved at construction time.
pub struct NativeWaveguide {
    // Parameters (SI internally; SPICE-style "_um"/"_nm" entry via set_real_param).
    length_m:        f64,
    n_g:             f64,
    alpha_neper_m:   f64,
    wavelength_m:    f64,
    // Terminal node indices: [in_re, in_im, in_lambda, out_re, out_im, out_lambda].
    nodes: [NodeId; 6],
    // Internal branch rows (one per potential equation).  Populated by
    // `bind_extra_nodes`.
    branches: [Option<usize>; 3],
    // Cached rotation coefficients (computed by `eval`, read by `load_jacobian`).
    c_cached: f64,
    s_cached: f64,
}

impl NativeWaveguide {
    pub fn new() -> Self {
        NativeWaveguide {
            length_m:      100e-6,
            n_g:           4.2,
            alpha_neper_m: dB_per_cm_to_neper_per_m(2.0),
            wavelength_m:  1550e-9,
            nodes:    [None; 6],
            branches: [None; 3],
            c_cached: 1.0,
            s_cached: 0.0,
        }
    }
}

impl Device for NativeWaveguide {
    fn num_terminals(&self) -> usize { 6 }

    fn setup_model(&mut self, _ctx: &SimContext) {}

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert_eq!(terminals.len(), 6,
            "NativeWaveguide: expected 6 terminals [in_re, in_im, in_λ, out_re, out_im, out_λ]");
        for i in 0..6 {
            self.nodes[i] = terminals[i];
        }
    }

    fn num_extra_nodes(&self) -> usize { 3 }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        self.branches[0] = Some(first_idx);
        self.branches[1] = Some(first_idx + 1);
        self.branches[2] = Some(first_idx + 2);
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "l_um"          => { self.length_m       = value * 1e-6;                       true }
            "l_m" | "length"=> { self.length_m       = value;                              true }
            "n_g"           => { self.n_g            = value;                              true }
            "alpha_db_cm"   => { self.alpha_neper_m  = dB_per_cm_to_neper_per_m(value);    true }
            "wavelength_nm" => { self.wavelength_m   = value * 1e-9;                       true }
            "wavelength_m"  => { self.wavelength_m   = value;                              true }
            _ => false,
        }
    }

    fn eval(&mut self, x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {
        // Read λ from the input wavelength node.  Bootstrap with the design
        // wavelength when the wire is near zero (initial NR iterate at x=0).
        let lambda = match self.nodes[2] {
            Some(i) => {
                let v = x[i];
                if v.abs() > self.wavelength_m * 0.5 { v } else { self.wavelength_m }
            }
            None => self.wavelength_m,
        };
        let beta  = 2.0 * std::f64::consts::PI * self.n_g / lambda;
        let phi   = beta * self.length_m;
        let t_amp = (-self.alpha_neper_m * self.length_m / 2.0).exp();
        self.c_cached = t_amp * phi.cos();
        self.s_cached = t_amp * phi.sin();
    }

    fn load_residual(&self, b: &mut [f64]) {
        // All three branch equations are homogeneous (target − Σ k_i·V_i = 0),
        // so b stays zero.  The branch-flow contribution to terminal-row
        // KCL is also zero for an ideal passive: net current at any optical
        // port is by convention zero in this discipline scheme.
        let _ = b;
    }

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        self.stamp(mat);
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) {
        self.load_residual(b);
    }

    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) {
        self.load_jacobian(mat);
    }
}

impl NativeWaveguide {
    /// Stamp the three branch equations into the MNA Jacobian using
    /// coefficients cached by the most recent `eval`.
    fn stamp(&self, mat: &mut MnaMatrix) {
        let c = self.c_cached;
        let s = self.s_cached;
        // Equation 1: V(out_re) − c·V(in_re) − s·V(in_im) = 0
        self.stamp_potential_eq(mat, 0, self.nodes[3], &[
            (self.nodes[0], -c),
            (self.nodes[1], -s),
        ]);
        // Equation 2: V(out_im) + s·V(in_re) − c·V(in_im) = 0
        self.stamp_potential_eq(mat, 1, self.nodes[4], &[
            (self.nodes[0],  s),
            (self.nodes[1], -c),
        ]);
        // Equation 3: V(out_λ) − V(in_λ) = 0
        self.stamp_potential_eq(mat, 2, self.nodes[5], &[
            (self.nodes[2], -1.0),
        ]);
    }

    /// Stamp one direct-potential equation `V(out) = Σ k_i · V(in_i)`.
    ///
    /// Allocates one auxiliary branch row at `self.branches[branch_idx]`.
    /// Stamp pattern:
    ///   branch row: +1 at out, k_i at each in_i, RHS = 0.
    ///   KCL at out: +1 at branch column (branch current leaves through out).
    fn stamp_potential_eq(
        &self,
        mat: &mut MnaMatrix,
        branch_idx: usize,
        out_node: NodeId,
        ins: &[(NodeId, f64)],
    ) {
        let (Some(out), Some(j)) = (out_node, self.branches[branch_idx]) else {
            // Output is ground or branch wasn't bound — skip.  Either is a
            // misconfiguration; let downstream singularity surface it.
            return;
        };
        // Branch row.
        mat.a[j][out] += 1.0;
        for &(in_node, k) in ins {
            if let Some(in_i) = in_node {
                mat.a[j][in_i] += k;
            }
        }
        // KCL at output: branch column carries the current that enforces V(out).
        mat.a[out][j] += 1.0;
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native directional coupler (2×2)
// ────────────────────────────────────────────────────────────────────────

/// 2×2 directional coupler with length-coupled cross-coefficient.
///
/// Lossless coupling matrix:
///   [c]   [ t   k] [a]    with t = cos(κL), k = -j sin(κL),
///   [d] = [ k   t] [b]    so |t|² + |k|² = 1.
///
/// In SVEA real/imag form:
///   c_re = t·a_re + s·b_im     d_re = t·b_re + s·a_im
///   c_im = t·a_im − s·b_re     d_im = t·b_im − s·a_re
/// (with t = cos(κL), s = sin(κL)).  Wavelength passes through unchanged
/// to both outputs from the corresponding input.
pub struct NativeDirectionalCoupler {
    kappa_per_m: f64,
    length_m:    f64,
    // Terminals: [a_re, a_im, a_λ, b_re, b_im, b_λ, c_re, c_im, c_λ, d_re, d_im, d_λ]
    nodes: [NodeId; 12],
    // 6 direct potential equations (c_re, c_im, c_λ, d_re, d_im, d_λ).
    branches: [Option<usize>; 6],
}

impl NativeDirectionalCoupler {
    pub fn new() -> Self {
        Self {
            kappa_per_m: 100.0,    // 100 rad/m — gives κL = 0.5 at L = 5 mm
            length_m:    5e-3,
            nodes:    [None; 12],
            branches: [None; 6],
        }
    }
}

impl Device for NativeDirectionalCoupler {
    fn num_terminals(&self) -> usize { 12 }

    fn setup_model(&mut self, _ctx: &SimContext) {}

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert_eq!(terminals.len(), 12);
        for i in 0..12 { self.nodes[i] = terminals[i]; }
    }

    fn num_extra_nodes(&self) -> usize { 6 }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..6 { self.branches[i] = Some(first_idx + i); }
    }

    fn set_real_param(&mut self, name: &str, value: f64) -> bool {
        match name.to_lowercase().as_str() {
            "kappa_per_m" | "kappa" => { self.kappa_per_m = value; true }
            "l_um"   => { self.length_m = value * 1e-6; true }
            "l_m" | "length" => { self.length_m = value; true }
            "kappa_l" | "kappal" => {
                // Direct coupling-angle override — keeps L fixed, scales κ.
                self.kappa_per_m = if self.length_m > 0.0 { value / self.length_m } else { 0.0 };
                true
            }
            _ => false,
        }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let kl = self.kappa_per_m * self.length_m;
        let t = kl.cos();
        let s = kl.sin();
        // c_re = t·a_re + s·b_im   (branch 0)
        self.stamp_potential_eq(mat, 0, self.nodes[6], &[
            (self.nodes[0], -t), (self.nodes[4], -s),
        ]);
        // c_im = t·a_im − s·b_re   (branch 1)
        self.stamp_potential_eq(mat, 1, self.nodes[7], &[
            (self.nodes[1], -t), (self.nodes[3],  s),
        ]);
        // c_λ  = a_λ              (branch 2)
        self.stamp_potential_eq(mat, 2, self.nodes[8], &[
            (self.nodes[2], -1.0),
        ]);
        // d_re = t·b_re + s·a_im   (branch 3)
        self.stamp_potential_eq(mat, 3, self.nodes[9], &[
            (self.nodes[3], -t), (self.nodes[1], -s),
        ]);
        // d_im = t·b_im − s·a_re   (branch 4)
        self.stamp_potential_eq(mat, 4, self.nodes[10], &[
            (self.nodes[4], -t), (self.nodes[0],  s),
        ]);
        // d_λ  = b_λ              (branch 5)
        self.stamp_potential_eq(mat, 5, self.nodes[11], &[
            (self.nodes[5], -1.0),
        ]);
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

impl NativeDirectionalCoupler {
    fn stamp_potential_eq(
        &self,
        mat: &mut MnaMatrix,
        branch_idx: usize,
        out_node: NodeId,
        ins: &[(NodeId, f64)],
    ) {
        let (Some(out), Some(j)) = (out_node, self.branches[branch_idx]) else { return };
        mat.a[j][out] += 1.0;
        for &(in_node, k) in ins {
            if let Some(in_i) = in_node { mat.a[j][in_i] += k; }
        }
        mat.a[out][j] += 1.0;
    }
}

// ────────────────────────────────────────────────────────────────────────
// Native 1×2 Y-junction splitter (3 dB lossless)
// ────────────────────────────────────────────────────────────────────────

/// Equal-power lossless splitter.  c = d = a / √2.  Wavelength duplicated.
pub struct NativeSplitter {
    // Terminals: [a_re, a_im, a_λ, c_re, c_im, c_λ, d_re, d_im, d_λ]
    nodes: [NodeId; 9],
    branches: [Option<usize>; 6],
}

impl NativeSplitter {
    pub fn new() -> Self {
        Self { nodes: [None; 9], branches: [None; 6] }
    }
}

impl Device for NativeSplitter {
    fn num_terminals(&self) -> usize { 9 }

    fn setup_model(&mut self, _ctx: &SimContext) {}

    fn setup_instance(&mut self, terminals: &[NodeId], _ctx: &SimContext) {
        debug_assert_eq!(terminals.len(), 9);
        for i in 0..9 { self.nodes[i] = terminals[i]; }
    }

    fn num_extra_nodes(&self) -> usize { 6 }

    fn bind_extra_nodes(&mut self, first_idx: usize) {
        for i in 0..6 { self.branches[i] = Some(first_idx + i); }
    }

    fn eval(&mut self, _x: &[f64], _flags: EvalFlags, _ctx: &SimContext) {}

    fn load_residual(&self, _b: &mut [f64]) {}

    fn load_jacobian(&self, mat: &mut MnaMatrix) {
        let k = 1.0 / 2.0_f64.sqrt();
        // c_re = k · a_re
        self.stamp_potential_eq(mat, 0, self.nodes[3], &[(self.nodes[0], -k)]);
        // c_im = k · a_im
        self.stamp_potential_eq(mat, 1, self.nodes[4], &[(self.nodes[1], -k)]);
        // c_λ  = a_λ
        self.stamp_potential_eq(mat, 2, self.nodes[5], &[(self.nodes[2], -1.0)]);
        // d_re = k · a_re
        self.stamp_potential_eq(mat, 3, self.nodes[6], &[(self.nodes[0], -k)]);
        // d_im = k · a_im
        self.stamp_potential_eq(mat, 4, self.nodes[7], &[(self.nodes[1], -k)]);
        // d_λ  = a_λ
        self.stamp_potential_eq(mat, 5, self.nodes[8], &[(self.nodes[2], -1.0)]);
    }

    fn load_residual_tran(&self, b: &mut [f64], _alpha: f64) { self.load_residual(b); }
    fn load_jacobian_tran(&self, mat: &mut MnaMatrix, _alpha: f64) { self.load_jacobian(mat); }
}

impl NativeSplitter {
    fn stamp_potential_eq(
        &self,
        mat: &mut MnaMatrix,
        branch_idx: usize,
        out_node: NodeId,
        ins: &[(NodeId, f64)],
    ) {
        let (Some(out), Some(j)) = (out_node, self.branches[branch_idx]) else { return };
        mat.a[j][out] += 1.0;
        for &(in_node, k) in ins {
            if let Some(in_i) = in_node { mat.a[j][in_i] += k; }
        }
        mat.a[out][j] += 1.0;
    }
}

// ────────────────────────────────────────────────────────────────────────
// Shared utilities
// ────────────────────────────────────────────────────────────────────────

#[allow(non_snake_case)]
fn dB_per_cm_to_neper_per_m(alpha_db_cm: f64) -> f64 {
    // 1 dB = ln(10)/20 Np ≈ 0.1151 Np; 1 cm = 0.01 m → multiply by 100/cm.
    alpha_db_cm * 100.0 * std::f64::consts::LN_10 / 20.0
}

#[cfg(test)]
mod tests {
    use crate::newton::dc_op_nr_with_registry;
    use crate::device_registry::DeviceRegistry;
    use fairchild_parser::parse_spice;

    /// Drive a native waveguide directly through voltage sources on its
    /// underlying re/im/λ wires and verify the output amplitude matches the
    /// closed-form `exp(-α·L/2) · exp(-jβL)` formula.
    ///
    /// `.optical` declarations are skipped here — we use XOsdi (which the
    /// discipline check exempts) for the native device and treat the wires
    /// as plain electrical nets for the voltage sources.
    #[test]
    fn native_waveguide_amplitude_matches_closed_form() {
        // L = 100 µm, n_g = 4.2, α = 2 dB/cm, λ = 1.55 µm.
        // T = exp(-α_Np/m · L / 2), where α_Np/m = 2·100·ln(10)/20 ≈ 23.026 Np/m.
        //   → T = exp(-23.026·100e-6/2) = exp(-1.1513e-3) ≈ 0.998850
        // φ = 2π·n_g·L/λ ≈ 2π·4.2·100/1550 wraps; |A_out| = T regardless of φ.
        let netlist = parse_spice(
            "* native waveguide test\n\
             V_re in_re 0 DC 1.0\n\
             V_im in_im 0 DC 0.0\n\
             V_wl in_wl 0 DC 1.55e-6\n\
             X1 in_re in_im in_wl out_re out_im out_wl fc_waveguide \
                L_um=100 n_g=4.2 alpha_dB_cm=2.0 wavelength_nm=1550\n\
             .op\n.end\n"
        ).unwrap();
        let registry = DeviceRegistry::new();
        let r = dc_op_nr_with_registry(&netlist, &registry)
            .expect("DC OP should converge");
        let v_re = r.node_voltage("out_re").unwrap();
        let v_im = r.node_voltage("out_im").unwrap();
        let amp  = (v_re * v_re + v_im * v_im).sqrt();
        let expected = (-23.0258509_f64 * 100e-6 / 2.0).exp(); // 0.99885
        assert!((amp - expected).abs() < 1e-5,
            "|A_out|={amp:.6} expected={expected:.6}");
        // Output wavelength must equal input wavelength.
        let v_wl = r.node_voltage("out_wl").unwrap();
        assert!((v_wl - 1.55e-6).abs() < 1e-15);
    }

    /// 3 dB splitter — equal-power split, in-phase.  |c|² + |d|² = |a|²,
    /// and c = d (both halves of the input).
    #[test]
    fn native_splitter_equal_power_split() {
        let netlist = parse_spice(
            "* splitter test\n\
             V_re a_re 0 DC 1.0\n\
             V_im a_im 0 DC 0.0\n\
             V_wl a_wl 0 DC 1.55e-6\n\
             X1 a_re a_im a_wl c_re c_im c_wl d_re d_im d_wl fc_splitter\n\
             .op\n.end\n"
        ).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new())
            .expect("DC OP should converge");
        let c_re = r.node_voltage("c_re").unwrap();
        let c_im = r.node_voltage("c_im").unwrap();
        let d_re = r.node_voltage("d_re").unwrap();
        let d_im = r.node_voltage("d_im").unwrap();
        // Each output: a / √2 ≈ 0.7071
        let expected_re = 1.0 / 2.0_f64.sqrt();
        assert!((c_re - expected_re).abs() < 1e-9);
        assert!((d_re - expected_re).abs() < 1e-9);
        assert!(c_im.abs() < 1e-9);
        assert!(d_im.abs() < 1e-9);
        // Power conservation: c² + d² ≈ 1 (input was 1.0)
        let p_total = c_re * c_re + c_im * c_im + d_re * d_re + d_im * d_im;
        assert!((p_total - 1.0).abs() < 1e-9);
    }

    /// Directional coupler: at κL = π/4, transmission and coupling are equal
    /// (cos²(π/4) = sin²(π/4) = 0.5).  With input only on port a, the cross
    /// port d gets power transferred via the imaginary cross-coupling.
    #[test]
    fn native_dcoupler_half_power_at_kl_pi_over_4() {
        let netlist = parse_spice(
            "* dcoupler test\n\
             V_ar a_re 0 DC 1.0\n\
             V_ai a_im 0 DC 0.0\n\
             V_awl a_wl 0 DC 1.55e-6\n\
             V_br b_re 0 DC 0.0\n\
             V_bi b_im 0 DC 0.0\n\
             V_bwl b_wl 0 DC 1.55e-6\n\
             X1 a_re a_im a_wl b_re b_im b_wl \
                c_re c_im c_wl d_re d_im d_wl fc_dcoupler kappa_L=0.7853981633974483\n\
             .op\n.end\n"
        ).unwrap();
        let r = dc_op_nr_with_registry(&netlist, &DeviceRegistry::new())
            .expect("DC OP should converge");
        // With a=(1,0), b=(0,0), κL=π/4 → t=s=1/√2:
        //   c_re = t·a_re = 1/√2,    c_im = -s·b_re = 0
        //   d_re = t·b_re + s·a_im = 0,   d_im = -s·a_re = -1/√2
        let half = 1.0 / 2.0_f64.sqrt();
        let c_re = r.node_voltage("c_re").unwrap();
        let c_im = r.node_voltage("c_im").unwrap();
        let d_re = r.node_voltage("d_re").unwrap();
        let d_im = r.node_voltage("d_im").unwrap();
        assert!((c_re - half).abs() < 1e-9, "c_re={c_re} expected {half}");
        assert!(c_im.abs() < 1e-9);
        assert!(d_re.abs() < 1e-9);
        assert!((d_im + half).abs() < 1e-9, "d_im={d_im} expected {}", -half);
        // Power: |c|² + |d|² = 1 (lossless).
        let p = c_re * c_re + c_im * c_im + d_re * d_re + d_im * d_im;
        assert!((p - 1.0).abs() < 1e-9, "p={p}");
    }
}
