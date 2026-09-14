//! Is a lumped device short enough to *be* lumped?
//!
//! A lumped phase shifter holds its whole electrode at one voltage. That is
//! true while the device is short against the RF wavelength and false exactly
//! where a travelling-wave modulator earns its name — and the lumped answer is
//! plausible either way. The eye closes smoothly, the bandwidth comes out
//! finite, nothing looks wrong. `fc_tw_ps` models the electrode; nothing told a
//! user which one their deck needed (#122).
//!
//! # Why this warns rather than switching
//!
//! The obvious version is to auto-segment when the deck "needs" it. Two things
//! stop that being a good idea.
//!
//! The failure condition is electrical **length** — `L·n_m·f/c` — which belongs
//! to the device and the drive, not to the timestep. A 3 mm electrode at
//! `n_m = 4` is a tenth of a wavelength at 2.5 GHz whether it is stepped at 1 ps
//! or 100 ps, and fine at 100 MHz either way. And the drive's bandwidth is not
//! knowable before solving: a `PULSE` has content to `1/t_r`, a `PWL` has
//! whatever its corners imply, and a compiled driver has whatever it decides to
//! do at run time. Inferring a bandwidth from that and silently changing the
//! device is how a confidently wrong answer gets made.
//!
//! So: say what was noticed, name the element and the number, and name the
//! parameter that fixes it. `crate::unmodelled` is the precedent — the codebase
//! already has a place for "this deck names an effect it is not getting", and
//! each entry says what the deck *loses*.
//!
//! # Why the transient and not `.ac`
//!
//! `.ac` is where a bandwidth gets measured, so it looks like the place this
//! matters most. It is the wrong place to *ask*: a sweep's top frequency is
//! chosen to bracket a response, not to describe a drive. Sweeping to 1 THz to
//! find a 17 GHz pole is normal, and reporting a 3 mm arm as "42 wavelengths
//! long at 1000 GHz" is true and useless. Measured on
//! `examples/photonic/noisy_eye_and_ber.py`, which emitted 28 of them.
//!
//! A transient's sources are different: their edges are a property of the deck
//! the user wrote, and the fastest of them is a real statement about what the
//! circuit is being asked to carry. So the check runs there, once per element.

use std::collections::HashSet;
use std::sync::{Mutex, OnceLock};

use fairchild_parser::{Element, Netlist, Waveform};

use crate::device::Device;
use crate::warn_user;

/// Fraction of an RF wavelength past which a lumped electrode is worth
/// mentioning.
///
/// A tenth is the usual engineering line for "lumped is fine", and it is where
/// the walk-off `sinc` has fallen by about 1.6 %. Below it the two models agree
/// to better than a component tolerance; above it they start to part.
const LUMPED_FRACTION: f64 = 0.1;

/// The highest frequency the deck's own sources put into the circuit.
///
/// Edges, not fundamentals: a 1 GHz square wave with 10 ps edges is a 35 GHz
/// signal, and it is the edge that finds a long electrode. The knee `0.35/t_r`
/// is the standard rule of thumb for where a linear ramp's spectrum rolls off.
///
/// `None` when every source is DC, which is the honest answer — a deck with no
/// time-varying drive has no bandwidth for this to be measured against.
/// Compiled drivers (OSDI/Verilog-A) are invisible here for the same reason
/// they are invisible everywhere before a solve: they have not run yet.
pub fn source_knee_hz(netlist: &Netlist) -> Option<f64> {
    let knee_of_rise = |t: f64| (t > 0.0).then(|| 0.35 / t);
    let mut hi: Option<f64> = None;
    let mut note = |f: Option<f64>| {
        if let Some(f) = f.filter(|f| f.is_finite() && *f > 0.0) {
            hi = Some(hi.map_or(f, |h: f64| h.max(f)));
        }
    };
    for el in &netlist.elements {
        let w = match el {
            Element::VoltageSource { waveform, .. } | Element::CurrentSource { waveform, .. } => {
                waveform
            }
            _ => continue,
        };
        match w {
            Waveform::Dc(_) => {}
            Waveform::Pulse { tr, tf, .. } => {
                note(knee_of_rise(tr.min(*tf).max(f64::MIN_POSITIVE)));
            }
            Waveform::Pwl { points } => {
                // The shortest segment is the fastest edge in the list.
                let shortest = points
                    .windows(2)
                    .map(|w| w[1].0 - w[0].0)
                    .filter(|dt| *dt > 0.0)
                    .fold(f64::INFINITY, f64::min);
                note(knee_of_rise(shortest));
            }
            Waveform::Sin { freq, .. } => note(Some(*freq)),
            Waveform::Exp { tau1, tau2, .. } => note(knee_of_rise(tau1.min(*tau2))),
            // The highest component of an FM/AM pair is the carrier plus its
            // sidebands; the modulation index sets how far they reach.
            Waveform::Sffm { fc, mdi, fs, .. } => note(Some(fc + (1.0 + mdi.abs()) * fs)),
            Waveform::Am { fc, mf, .. } => note(Some(fc + mf)),
        }
    }
    hi
}

/// What every lumped device in a circuit assumes, collected once.
///
/// Devices are not always still in scope when the frequency is known — `.ac`
/// assembles its matrices before it is handed a sweep — so the limits are
/// gathered at build time and the comparison happens later.
#[derive(Clone, Debug, Default)]
pub struct LumpedLimits {
    /// `(element name, highest frequency the lumped assumption holds to)`.
    entries: Vec<(String, f64)>,
}

impl LumpedLimits {
    /// Ask every device, keeping the ones that make a lumped assumption.
    ///
    /// `names` is parallel to `devices` — `crate::newton::build_device_names`
    /// builds it — so a message can say which element rather than which index.
    pub fn collect(devices: &[Box<dyn Device>], names: &[String]) -> Self {
        LumpedLimits {
            entries: devices
                .iter()
                .enumerate()
                .filter_map(|(i, d)| {
                    let f = d.lumped_valid_to_hz()?;
                    let name = names
                        .get(i)
                        .cloned()
                        .unwrap_or_else(|| "a phase shifter".to_string());
                    Some((name, f))
                })
                .collect(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Warn about every entry that `f_hz` is past, once per element for the
    /// life of the process.
    ///
    /// A script that runs twenty transients over one netlist — a sweep, a
    /// fit, a Monte Carlo — would otherwise repeat itself twenty times, and a
    /// warning nobody finishes reading is a warning nobody reads.
    pub fn warn_above(&self, f_hz: f64) {
        if !(f_hz.is_finite() && f_hz > 0.0) {
            return;
        }
        static SAID: OnceLock<Mutex<HashSet<String>>> = OnceLock::new();
        for (name, f_lumped) in &self.entries {
            if f_hz <= *f_lumped {
                continue;
            }
            let fresh = SAID
                .get_or_init(|| Mutex::new(HashSet::new()))
                .lock()
                .map(|mut said| said.insert(name.clone()))
                .unwrap_or(true);
            if !fresh {
                continue;
            }
            warn_user!(
                "{name} is {:.2} RF wavelengths long at {:.1} GHz, and a lumped \
                 phase shifter holds its whole electrode at one voltage. Its \
                 lumped assumption holds to about {:.1} GHz. `fc_tw_ps` models \
                 the electrode as a transmission line — see \
                 docs/photonic-models.md",
                LUMPED_FRACTION * f_hz / f_lumped,
                f_hz / 1e9,
                f_lumped / 1e9,
            );
        }
    }
}

/// [`LumpedLimits::warn_above`] against the deck's own sources.
pub fn warn_from_sources(netlist: &Netlist, devices: &[Box<dyn Device>], names: &[String]) {
    if let Some(f) = source_knee_hz(netlist) {
        LumpedLimits::collect(devices, names).warn_above(f);
    }
}

/// The frequency at which a device of length `l_m` on an electrode of index
/// `n_m` reaches [`LUMPED_FRACTION`] of a wavelength.
///
/// Shared so a device implementing
/// [`Device::lumped_valid_to_hz`](crate::device::Device::lumped_valid_to_hz)
/// does not have to restate the criterion, and so there is one place to change
/// it.
pub fn lumped_limit_hz(l_m: f64, n_m: f64) -> Option<f64> {
    (l_m > 0.0 && n_m > 0.0).then(|| LUMPED_FRACTION * super::models::photonic::C0 / (l_m * n_m))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fairchild_parser::parse_spice;

    #[test]
    fn the_knee_is_the_fastest_edge_not_the_fundamental() {
        // 1 GHz square wave with 10 ps edges: a 35 GHz signal, not a 1 GHz one.
        let nl = parse_spice("* k\nV1 a 0 PULSE(0 1 0 10p 10p 500p 1n)\nR1 a 0 50\n").unwrap();
        let f = source_knee_hz(&nl).expect("a pulse has edges");
        assert!(
            (f - 35e9).abs() < 1e9,
            "0.35/10 ps is 35 GHz, got {:.1} GHz",
            f / 1e9
        );
    }

    #[test]
    fn a_dc_only_deck_has_no_bandwidth_to_report() {
        let nl = parse_spice("* dc\nV1 a 0 DC 1\nR1 a 0 50\n").unwrap();
        assert_eq!(source_knee_hz(&nl), None);
    }

    #[test]
    fn a_sine_reports_its_own_frequency() {
        let nl = parse_spice("* sin\nV1 a 0 SIN(0 1 2.5e9 0 0 0)\nR1 a 0 50\n").unwrap();
        let f = source_knee_hz(&nl).expect("a sine has a frequency");
        assert!((f - 2.5e9).abs() < 1.0, "got {f}");
    }

    #[test]
    fn the_fastest_source_in_the_deck_wins() {
        let nl = parse_spice(
            "* two\nV1 a 0 SIN(0 1 1e9 0 0 0)\nV2 b 0 PULSE(0 1 0 5p 5p 1n 2n)\n\
             R1 a 0 50\nR2 b 0 50\n",
        )
        .unwrap();
        let f = source_knee_hz(&nl).expect("both have content");
        assert!((f - 70e9).abs() < 1e9, "0.35/5 ps is 70 GHz, got {f:e}");
    }

    /// Only a device that *makes* a lumped assumption reports one.
    ///
    /// The on/off pair is the point. A long lumped shifter must report a limit,
    /// and `fc_tw_ps` — which models its electrode as a transmission line —
    /// must report nothing, or the check would be warning about the very
    /// device it recommends. A passive waveguide has no electrode at all.
    #[test]
    fn only_a_lumped_electrode_reports_a_limit() {
        use crate::models::{pn_phase_shifter_cap, NativeTwPhaseShifter, NativeWaveguide};

        let arm = |l_um: f64| {
            let mut d = pn_phase_shifter_cap();
            d.set_real_param("l_um", l_um);
            d.lumped_valid_to_hz()
        };
        // 3 mm at the assumed n_m = 4.2: a tenth of a wavelength at 2.4 GHz.
        let long = arm(3000.0).expect("a driven shifter has an electrode");
        assert!(
            (long - 2.379e9).abs() < 1e7,
            "3 mm should be lumped to about 2.4 GHz, got {:.3} GHz",
            long / 1e9
        );
        // 10 um is lumped past any frequency anyone sweeps.
        let short = arm(10.0).expect("still an electrode");
        assert!(
            short > 500e9,
            "10 um should be lumped well past any sweep, got {:.1} GHz",
            short / 1e9
        );

        assert_eq!(
            NativeTwPhaseShifter::new().lumped_valid_to_hz(),
            None,
            "fc_tw_ps models its electrode, so it makes no lumped assumption — \
             warning about it would be recommending the device to itself"
        );
        assert_eq!(
            NativeWaveguide::new().lumped_valid_to_hz(),
            None,
            "a passive waveguide has no electrode to be long"
        );
    }

    /// And the limits collect against their element names, which is what the
    /// message needs to be useful.
    #[test]
    fn the_limits_carry_the_element_name() {
        use crate::models::pn_phase_shifter_cap;

        let mut d = pn_phase_shifter_cap();
        d.set_real_param("l_um", 3000.0);
        let devices: Vec<Box<dyn Device>> = vec![Box::new(d)];
        let names = vec!["xarm1".to_string()];
        let limits = LumpedLimits::collect(&devices, &names);
        assert!(!limits.is_empty());
        assert_eq!(limits.entries[0].0, "xarm1");
        // Nothing to say when the drive is slower than the limit.
        assert!(LumpedLimits::default().is_empty());
    }

    #[test]
    fn the_limit_is_a_tenth_of_a_wavelength() {
        // 3 mm at n_m = 4: a wavelength is c/(4·f), and a tenth of it at
        // f = 0.1·c/(3e-3·4) = 2.5 GHz.
        let f = lumped_limit_hz(3e-3, 4.0).expect("positive geometry");
        assert!((f - 2.498e9).abs() < 1e7, "got {:.3} GHz", f / 1e9);
        assert_eq!(lumped_limit_hz(0.0, 4.0), None);
    }
}
