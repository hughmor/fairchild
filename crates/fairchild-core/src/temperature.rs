//! How device parameters move with temperature.
//!
//! `.temp` used to scale `kT/q` and nothing else. So a `.temp 125` run gave the
//! thermal voltage at 125 °C and every device parameter at nominal, and
//! `TNOM`/`EG`/`XTI` were on the accepted-but-not-modelled list for all three
//! native models (#77 §5). The warnings fired, so it was not silent — but
//! "temperature sweeps work" was a much stronger claim than what was true.
//!
//! # This file owns the laws
//!
//! Every temperature dependence in the native models comes from here, because
//! the same bandgap expression appears in the diode's saturation current, the
//! BJT's, and the MOSFET's threshold, and three copies are three chances to
//! disagree. The devices hold a *factor* derived once per solve rather than
//! calling in from `eval`.
//!
//! # Everything here was fitted from ngspice, not read from the SPICE source
//!
//! The SPICE literature carries several variants of each of these, and a PDK
//! relies on the one the reference simulator actually ships. Each law below was
//! back-solved from measured currents and then checked against explicit parameter
//! settings. Residuals, from `ngspice_temperature_golden.rs`:
//!
//! | law | how it was measured | residual |
//! |---|---|---|
//! | diode `IS(T)` | reverse saturation at −1 V *is* `−IS(T)` | fitted `XTI` 2.9999, `EG` 1.1100 |
//! | BJT `IS(T)` | `IC` at fixed `VBE`, dividing out `exp(VBE/vt(T))` | 1.6e−6 |
//! | BJT `BF(T)` | `IC/IB` | 1.5278 measured against 1.5278 |
//! | MOSFET `KP(T)` | `sqrt(Id)` against two `Vgs`, which separates it from `VTO` | exact |
//! | MOSFET `VTO(T)` | the same fit's intercept | 0.04 mV |

/// SPICE3's Boltzmann constant and electron charge, to the digits it uses.
///
/// Not `crate::device`'s: the bandgap expressions below are curve fits whose
/// coefficients were published against these values, so using a more accurate
/// `k` moves `PHI(T)` in the fourth decimal and stops matching ngspice.
const BOLTZ: f64 = 1.3806226e-23;
const CHARGE: f64 = 1.6021918e-19;

/// The temperature SPICE's built-in coefficients are referenced to, 27 °C.
///
/// Distinct from a card's `TNOM`: `TNOM` says what temperature the card's
/// parameters were extracted at, and this is a constant inside the fits.
pub const REFTEMP_K: f64 = 300.15;

/// Default `TNOM` when a card does not give one — also 27 °C.
pub const TNOM_DEFAULT_K: f64 = 300.15;

/// Silicon's activation energy, `EG`'s default, in eV.
pub const EG_DEFAULT: f64 = 1.11;

/// `XTI`'s default: the saturation current's temperature exponent.
pub const XTI_DEFAULT: f64 = 3.0;

/// `XTB`'s default: no beta temperature dependence.
pub const XTB_DEFAULT: f64 = 0.0;

/// Silicon's bandgap at `t` kelvin, the SPICE3 fit.
///
/// `1.16 − 7.02e-4·T²/(T + 1108)`. Used by the MOSFET threshold shift and by
/// [`mos_junction_is_factor`]. The diode's and the BJT's saturation currents take
/// `EG` from the card instead, because a card may describe a material this fit
/// does not. The MOSFET card has no `EG`, which is why it uses the fit.
pub fn bandgap_ev(t: f64) -> f64 {
    1.16 - (7.02e-4 * t * t) / (t + 1108.0)
}

/// Multiplier on a **MOSFET's** bulk junction saturation current at `t`.
///
/// ```text
/// Isat(T) = Isat · exp(Eg(TNOM)/vt(TNOM) − Eg(T)/vt(T))
/// ```
///
/// A third law, not either of the two above, and measured to be so. The diode's
/// and the BJT's use a **constant** `EG` from the card in `exp(EG·(T/TNOM −
/// 1)/vt(T))` and multiply by `(T/TNOM)^XTI`. A MOSFET card carries neither `EG`
/// nor `XTI`, and SPICE puts the temperature-dependent bandgap in the exponent
/// instead.
///
/// The difference is not small. Against ngspice at `TNOM = 27 °C`, reading `Isat`
/// off a bulk held one volt reverse with `gmin = 0`:
///
/// | T | ngspice | this law | the diode law at `EG = 1.11` |
/// |---|---|---|---|
/// | −40 °C | 1.835405e-20 | 1.834720e-20 | 4.404638e-20 |
/// | 0 °C | 1.047853e-16 | 1.047717e-16 | 1.437655e-16 |
/// | 75 °C | 5.986135e-12 | 5.987214e-12 | 3.712831e-12 |
/// | 125 °C | 9.446510e-10 | 9.449566e-10 | 3.869226e-10 |
///
/// The last column is out by up to 2.4×. This one agrees to 3.7e-4 worst case,
/// and the residual grows with `|T − TNOM|` in the exponent, which is where a
/// small difference in `k/q` would show. Solving the measurements for `k/q` gives
/// 8.61638e-5 against SPICE3's 8.617087e-5, so the remaining gap is a constant
/// this tree deliberately does not copy, on a leakage current between 1e-20 and
/// 1e-9 A.
pub fn mos_junction_is_factor(t: f64, tnom: f64) -> f64 {
    if t <= 0.0 || tnom <= 0.0 {
        return 1.0;
    }
    let vt = t * BOLTZ / CHARGE;
    let vtnom = tnom * BOLTZ / CHARGE;
    (bandgap_ev(tnom) / vtnom - bandgap_ev(t) / vt).exp()
}

/// The built-in-potential temperature term SPICE calls `pbfact`.
///
/// Isolated here because it appears twice in the MOSFET threshold — once at
/// `TNOM` to back out an unshifted `PHI`, once at `T` to shift it — and the two
/// uses must be the same function.
pub fn pbfact(t: f64) -> f64 {
    let vt = t * BOLTZ / CHARGE;
    let arg = -bandgap_ev(t) / (2.0 * BOLTZ * t) + 1.1150877 / (BOLTZ * 2.0 * REFTEMP_K);
    -2.0 * vt * (1.5 * (t / REFTEMP_K).ln() + CHARGE * arg)
}

/// Multiplier on a **diode's** `IS` at temperature `t`.
///
/// ```text
/// IS(T) = IS · exp((T/TNOM − 1)·EG/(N·vt(T))) · (T/TNOM)^(XTI/N)
/// ```
///
/// The emission coefficient divides both terms, which is what distinguishes this
/// from the BJT form below. Returns a factor rather than a current so the caller
/// can apply `AREA` afterwards and stay idempotent however many times setup runs.
pub fn diode_is_factor(t: f64, tnom: f64, eg: f64, xti: f64, n: f64) -> f64 {
    if t <= 0.0 || tnom <= 0.0 || n <= 0.0 {
        return 1.0;
    }
    let ratio = t / tnom;
    let vt = t * BOLTZ / CHARGE;
    ((ratio - 1.0) * eg / (n * vt)).exp() * ratio.powf(xti / n)
}

/// Multiplier on a **BJT's** `IS` at temperature `t`.
///
/// ```text
/// IS(T) = IS · (T/TNOM)^XTI · exp(EG·(T/TNOM − 1)/vt(T))
/// ```
///
/// The same shape as the diode's with no emission coefficient in it. Written as
/// its own function rather than `diode_is_factor(.., n = 1.0)` because they are
/// two model definitions that happen to coincide at `N = 1`, and collapsing them
/// would make a future divergence in either a silent change to the other.
pub fn bjt_is_factor(t: f64, tnom: f64, eg: f64, xti: f64) -> f64 {
    if t <= 0.0 || tnom <= 0.0 {
        return 1.0;
    }
    let ratio = t / tnom;
    let vt = t * BOLTZ / CHARGE;
    ratio.powf(xti) * (eg * (ratio - 1.0) / vt).exp()
}

/// Multiplier on `BF`/`BR` at temperature `t`: `(T/TNOM)^XTB`.
pub fn beta_factor(t: f64, tnom: f64, xtb: f64) -> f64 {
    if t <= 0.0 || tnom <= 0.0 {
        return 1.0;
    }
    (t / tnom).powf(xtb)
}

/// Multiplier on a MOSFET's `KP` at temperature `t`: `(T/TNOM)^-1.5`.
///
/// This is the mobility law. The exponent is SPICE's and is not a card
/// parameter — `UO`/`THETA`-style mobility degradation is a different thing and
/// is not modelled.
pub fn mobility_factor(t: f64, tnom: f64) -> f64 {
    if t <= 0.0 || tnom <= 0.0 {
        return 1.0;
    }
    (t / tnom).powf(-1.5)
}

/// A MOSFET's surface potential at `t`, from its nominal `PHI`.
pub fn scaled_phi(phi_nom: f64, t: f64, tnom: f64) -> f64 {
    if t <= 0.0 || tnom <= 0.0 {
        return phi_nom;
    }
    // Back out the temperature-free part at TNOM, then put it back at T. The
    // division by `fact1` is why `pbfact` has to be one function: a `pbfact` that
    // differed between the two calls would leave a residual shift at `T = TNOM`.
    let fact1 = tnom / REFTEMP_K;
    let phio = (phi_nom - pbfact(tnom)) / fact1;
    (t / REFTEMP_K) * phio + pbfact(t)
}

/// A junction's built-in potential at `t`, from its nominal value.
///
/// The same law as [`scaled_phi`], because a MOSFET's surface potential and a
/// diode's junction potential are the same quantity in SPICE and move by the same
/// expression. Named separately so a reader of `diode.rs` finds it, and delegating
/// so the two cannot drift.
pub fn scaled_junction_potential(vj_nom: f64, t: f64, tnom: f64) -> f64 {
    scaled_phi(vj_nom, t, tnom)
}

/// Multiplier on a zero-bias junction capacitance (`CJO`, `CJE`, `CJC`, `CJ`) at
/// temperature `t`.
///
/// ```text
/// pbo     = (VJ − pbfact(TNOM)) / (TNOM/300.15)
/// gmaold  = (VJ    − pbo)/pbo      cjfact  = 1 / (1 + M·(4e-4·(TNOM−300.15) − gmaold))
/// gmanew  = (VJ(T) − pbo)/pbo      cjfact1 =      1 + M·(4e-4·(T   −300.15) − gmanew)
/// CJO(T)  = CJO · cjfact · cjfact1
/// ```
///
/// The two halves are not redundant: `cjfact` un-references the card from its own
/// `TNOM` and `cjfact1` re-references it to `T`, so a card extracted at `TNOM`
/// and run at `TNOM` comes back unshifted. That identity is the one property a
/// cross-simulator comparison cannot check — both would share an offset — and it
/// is asserted directly in the tests.
///
/// Measured against ngspice to 1.2e-6…8.2e-5 across −40/27/75/125 °C at two
/// values of `M`. `M` appears three times in the law, which is why it is swept.
pub fn junction_cap_factor(vj_nom: f64, m: f64, t: f64, tnom: f64) -> f64 {
    if t <= 0.0 || tnom <= 0.0 || vj_nom <= 0.0 {
        return 1.0;
    }
    let fact1 = tnom / REFTEMP_K;
    let pbo = (vj_nom - pbfact(tnom)) / fact1;
    if pbo == 0.0 {
        return 1.0;
    }
    let vj_t = scaled_junction_potential(vj_nom, t, tnom);
    let gmaold = (vj_nom - pbo) / pbo;
    let denom = 1.0 + m * (4e-4 * (tnom - REFTEMP_K) - gmaold);
    if denom == 0.0 {
        return 1.0;
    }
    let gmanew = (vj_t - pbo) / pbo;
    let cjfact1 = 1.0 + m * (4e-4 * (t - REFTEMP_K) - gmanew);
    cjfact1 / denom
}

/// A MOSFET's threshold at `t`.
///
/// ```text
/// VTO(T) = VTO − type·Γ·√PHI + ½(EG(TNOM) − EG(T)) + type·½(PHI(T) − PHI)
///                + type·Γ·√PHI(T)
/// ```
///
/// `type` is +1 for NMOS and −1 for PMOS. With `GAMMA = 0` the two body-effect
/// terms drop and the shift is the bandgap and the surface potential; that is the
/// case the 0.04 mV residual was measured in, and `GAMMA` is covered separately
/// because `√PHI(T)` enters twice and a sign error there cancels at `Γ = 0`.
pub fn scaled_vto(
    vto: f64,
    gamma: f64,
    phi_nom: f64,
    phi_t: f64,
    t: f64,
    tnom: f64,
    is_pmos: bool,
) -> f64 {
    let ty = if is_pmos { -1.0 } else { 1.0 };
    let vbi = vto - ty * gamma * phi_nom.max(0.0).sqrt()
        + 0.5 * (bandgap_ev(tnom) - bandgap_ev(t))
        + ty * 0.5 * (phi_t - phi_nom);
    vbi + ty * gamma * phi_t.max(0.0).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every law must be the identity at `T = TNOM`. This is the one property
    /// that cannot be checked against another simulator — both would agree about
    /// a shared offset — and getting it wrong shifts every deck that never asked
    /// for a temperature at all.
    #[test]
    fn every_law_is_the_identity_at_nominal() {
        for tnom in [300.15, 273.15, 400.0] {
            let t = tnom;
            assert!(
                (diode_is_factor(t, tnom, EG_DEFAULT, XTI_DEFAULT, 1.0) - 1.0).abs() < 1e-15,
                "diode IS factor at T = TNOM = {tnom}"
            );
            assert!(
                (bjt_is_factor(t, tnom, EG_DEFAULT, XTI_DEFAULT) - 1.0).abs() < 1e-15,
                "BJT IS factor at T = TNOM = {tnom}"
            );
            assert!((beta_factor(t, tnom, 1.5) - 1.0).abs() < 1e-15);
            assert!((mobility_factor(t, tnom) - 1.0).abs() < 1e-15);
            let phi = scaled_phi(0.6, t, tnom);
            assert!(
                (phi - 0.6).abs() < 1e-12,
                "PHI at T = TNOM = {tnom} is {phi}, not 0.6"
            );
            let vto = scaled_vto(0.7, 0.5, 0.6, phi, t, tnom, false);
            assert!(
                (vto - 0.7).abs() < 1e-12,
                "VTO at T = TNOM = {tnom} is {vto}, not 0.7"
            );
        }
    }

    /// The measured MOSFET threshold, to the precision it was measured at.
    ///
    /// ngspice at `VTO=0.7 PHI=0.6 GAMMA=0`, fitted from `sqrt(Id)` against two
    /// `Vgs` so the fit separates the threshold from the mobility.
    #[test]
    fn vto_matches_the_measured_shift() {
        for (t, want) in [(348.15, 0.6521), (398.15, 0.6014)] {
            let phi = scaled_phi(0.6, t, 300.15);
            let got = scaled_vto(0.7, 0.0, 0.6, phi, t, 300.15, false);
            assert!(
                (got - want).abs() < 1e-4,
                "VTO({t} K) is {got:.4} and ngspice fits {want:.4}"
            );
        }
    }

    /// The junction-capacitance law is the identity at `T = TNOM`, for every
    /// `TNOM` and every `M`.
    ///
    /// This is the property no cross-simulator comparison can check: both would
    /// share an offset. It is also the one that breaks every existing transient
    /// golden at once if the two halves of the law disagree — `cjfact`
    /// un-references the card and `cjfact1` re-references it, and they have to
    /// cancel exactly.
    #[test]
    fn the_junction_cap_law_is_the_identity_at_nominal() {
        for tnom in [300.15, 273.15, 398.15] {
            for m in [0.0, 0.33, 0.5, 1.0] {
                for vj in [0.6, 0.75, 1.0] {
                    let f = junction_cap_factor(vj, m, tnom, tnom);
                    assert!(
                        (f - 1.0).abs() < 1e-12,
                        "CJO factor at T = TNOM = {tnom}, M = {m}, VJ = {vj} is \
                         {f}, not 1"
                    );
                    let p = scaled_junction_potential(vj, tnom, tnom);
                    assert!(
                        (p - vj).abs() < 1e-12,
                        "VJ at T = TNOM = {tnom} is {p}, not {vj}"
                    );
                }
            }
        }
    }

    /// A junction gets *wider* capacitance as it heats, because the built-in
    /// potential narrows with the bandgap. A sign error in either half of the
    /// law reverses this while still passing the identity test above.
    #[test]
    fn a_hotter_junction_has_more_capacitance_and_less_potential() {
        let (cold, hot, tnom) = (233.15, 398.15, 300.15);
        let vj = 0.75;
        assert!(
            scaled_junction_potential(vj, hot, tnom) < vj,
            "the built-in potential falls as the bandgap narrows"
        );
        assert!(
            scaled_junction_potential(vj, cold, tnom) > vj,
            "and rises when cold"
        );
        assert!(
            junction_cap_factor(vj, 0.5, hot, tnom) > 1.0,
            "less potential means a thinner depletion layer, so more capacitance"
        );
        assert!(
            junction_cap_factor(vj, 0.5, cold, tnom) < 1.0,
            "and less when cold"
        );
        // M scales the whole correction, so M = 0 leaves only the 4e-4 term.
        let m0 = junction_cap_factor(vj, 0.0, hot, tnom);
        assert!(
            (m0 - 1.0).abs() < 1e-12,
            "M = 0 removes the grading correction entirely: got {m0}"
        );
    }

    /// The direction of each law, which a transposed exponent would break while
    /// still landing on 1.0 at nominal.
    #[test]
    fn each_law_moves_the_right_way() {
        let (t, tnom) = (398.15, 300.15);
        assert!(
            diode_is_factor(t, tnom, EG_DEFAULT, XTI_DEFAULT, 1.0) > 1e4,
            "leakage rises steeply with temperature — it is the reason a hot \
             diode leaks"
        );
        assert!(
            mobility_factor(t, tnom) < 1.0,
            "mobility falls with temperature, so KP falls and drive current with it"
        );
        assert!(
            scaled_phi(0.6, t, tnom) < 0.6,
            "the surface potential falls as the bandgap narrows"
        );
        assert!(
            scaled_vto(0.7, 0.0, 0.6, scaled_phi(0.6, t, tnom), t, tnom, false) < 0.7,
            "and the threshold falls with it"
        );
        assert!(beta_factor(t, tnom, 1.5) > 1.0, "XTB > 0 raises beta");
    }
}
