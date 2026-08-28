//! Binned `.model` cards: one model name, several parameter sets, geometry
//! picks which one.
//!
//! Essentially every foundry PDK bins its transistors. The card name carries the
//! bin index and the card carries the window it covers:
//!
//! ```spice
//! .model nch.1 nmos (LMIN=0.18u LMAX=0.30u WMIN=0.22u WMAX=1u  VTO=0.42 …)
//! .model nch.2 nmos (LMIN=0.30u LMAX=1.00u WMIN=0.22u WMAX=1u  VTO=0.45 …)
//! M1 d g s b nch W=0.5u L=0.25u
//! ```
//!
//! `M1` gets bin 1. Before this, `.model nch.1` registered a model *called*
//! `nch.1`, and the element asking for `nch` failed as `unknown model` — which is
//! the whole of why a PDK would not load.
//!
//! # This file owns the concept
//!
//! "What is a bin", "which bin does this geometry get", and "what happens when
//! none of them match" are answered here and nowhere else. The registry stores a
//! [`BinGroup`] per model name and asks it; the Spectre front end rewrites its
//! braced form into the `.model name.N` form so there is one representation to
//! reason about rather than two.
//!
//! # The selection rule, and why this one
//!
//! A bin matches when `LMIN ≤ L ≤ LMAX` and `WMIN ≤ W ≤ WMAX`, closed on both
//! ends. Real PDKs write contiguous windows (`…LMAX=0.30u` then `LMIN=0.30u…`),
//! so half-open intervals would drop `L = 0.30u` into a gap and refuse a deck
//! that is not wrong. Closed intervals never gap.
//!
//! Closed intervals do let two bins share a boundary, and there the tightest
//! window wins, then the lowest bin index. That is arbitrary but *deterministic*,
//! and it is safe because a PDK's parameters are continuous across a boundary it
//! chose to write twice.
//!
//! Bins whose **interiors** overlap are a different thing — a malformed card set
//! where the choice really does change the answer. Those are named on stderr once
//! at registration, not silently resolved.
//!
//! Geometry outside every bin is a **hard error** naming the geometry and every
//! window, per `docs/spice_support.md`. Picking the nearest bin would be a wrong
//! answer with nothing to read.

use crate::error::SimError;

/// The four card parameters that choose a bin rather than describing a device.
pub const SELECTORS: [&str; 4] = ["lmin", "lmax", "wmin", "wmax"];

/// True for a parameter name that selects a bin. Registration filters these out
/// of its unknown-parameter warning, because they are not model parameters.
#[allow(non_snake_case)]
pub fn IS_SELECTOR(key: &str) -> bool {
    SELECTORS.iter().any(|s| key.eq_ignore_ascii_case(s))
}

/// The geometry window one bin covers, in metres.
///
/// An unbinned card gets the unbounded window, so a plain `.model` and a binned
/// one take the same path through [`BinGroup::select`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Window {
    pub lmin: f64,
    pub lmax: f64,
    pub wmin: f64,
    pub wmax: f64,
}

impl Window {
    /// Covers every geometry — what an unbinned card gets.
    pub const UNBOUNDED: Window = Window {
        lmin: f64::NEG_INFINITY,
        lmax: f64::INFINITY,
        wmin: f64::NEG_INFINITY,
        wmax: f64::INFINITY,
    };

    fn is_unbounded(&self) -> bool {
        *self == Window::UNBOUNDED
    }

    fn contains(&self, l: f64, w: f64) -> bool {
        l >= self.lmin && l <= self.lmax && w >= self.wmin && w <= self.wmax
    }

    /// Area of the window, for "tightest bin wins". Infinite for the unbounded
    /// one, which is what makes an explicit bin beat a catch-all card.
    fn extent(&self) -> f64 {
        (self.lmax - self.lmin) * (self.wmax - self.wmin)
    }

    /// True when the two windows share interior, as opposed to touching at an
    /// edge. Touching gives `a.lmax == b.lmin`, so one of the strict tests fails.
    fn overlaps_interior(&self, other: &Window) -> bool {
        let l = self.lmin < other.lmax && other.lmin < self.lmax;
        let w = self.wmin < other.wmax && other.wmin < self.wmax;
        l && w
    }

    fn describe(&self) -> String {
        let f = |v: f64, unbounded: &str| {
            if v.is_infinite() {
                unbounded.to_string()
            } else {
                format!("{v:.4e}")
            }
        };
        format!(
            "L {}..{}, W {}..{}",
            f(self.lmin, "-inf"),
            f(self.lmax, "+inf"),
            f(self.wmin, "-inf"),
            f(self.wmax, "+inf")
        )
    }
}

/// One card of a binned set: the window it covers and its parameters.
#[derive(Debug, Clone)]
pub struct ModelBin {
    /// The `N` from `name.N`. `None` for an unbinned card.
    pub index: Option<u32>,
    pub window: Window,
    /// The card's parameters, with the four geometry selectors removed — they
    /// select a bin and are not model parameters, so leaving them in would have
    /// every binned card warn about four unknown parameters.
    pub params: Vec<(String, f64)>,
}

/// Every card registered under one model name.
#[derive(Debug, Clone, Default)]
pub struct BinGroup {
    bins: Vec<ModelBin>,
}

impl BinGroup {
    pub fn push(&mut self, bin: ModelBin) {
        self.bins.push(bin);
    }

    /// True when this name carries more than one card, or one card with a
    /// window. Used only for diagnostics.
    pub fn is_binned(&self) -> bool {
        self.bins.len() > 1 || self.bins.iter().any(|b| !b.window.is_unbounded())
    }

    pub fn bins(&self) -> &[ModelBin] {
        &self.bins
    }

    /// Any card's parameters, for a caller that has no geometry to select with.
    ///
    /// Only correct for an unbinned name, and every caller of it is in that
    /// position by construction (a device family with no `W`/`L`). A binned name
    /// reached here would silently take bin one, so it does not: the caller gets
    /// `None` and has to go through [`Self::select`].
    pub fn unbinned(&self) -> Option<&[(String, f64)]> {
        match self.bins.as_slice() {
            [only] if only.window.is_unbounded() => Some(&only.params),
            _ => None,
        }
    }

    /// The parameters for an instance of geometry `l` × `w`.
    ///
    /// See the module docs for the rule. `name` is only for the error text.
    pub fn select(&self, name: &str, l: f64, w: f64) -> Result<&[(String, f64)], SimError> {
        let mut matching: Vec<&ModelBin> = self
            .bins
            .iter()
            .filter(|b| b.window.contains(l, w))
            .collect();
        if matching.is_empty() {
            let mut windows: Vec<String> = self
                .bins
                .iter()
                .map(|b| {
                    let which = b
                        .index
                        .map_or("(unbinned)".into(), |i| format!("{name}.{i}"));
                    format!("  {which}: {}", b.window.describe())
                })
                .collect();
            windows.sort();
            return Err(SimError::NoMatchingBin {
                model: name.to_string(),
                l,
                w,
                windows: windows.join("\n"),
            });
        }
        // Tightest window first, then lowest bin index. Deterministic at a shared
        // boundary; see the module docs for why arbitrary is acceptable there and
        // why a real interior overlap is reported instead.
        matching.sort_by(|a, b| {
            a.window
                .extent()
                .partial_cmp(&b.window.extent())
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.index.cmp(&b.index))
        });
        Ok(&matching[0].params)
    }

    /// Interior overlaps, as `(index_a, index_b)` pairs. Empty for a well-formed
    /// set. Called once per model name at registration.
    pub fn interior_overlaps(&self) -> Vec<(String, String)> {
        let label = |b: &ModelBin| b.index.map_or("(unbinned)".into(), |i| format!(".{i}"));
        let mut out = Vec::new();
        for (i, a) in self.bins.iter().enumerate() {
            for b in &self.bins[i + 1..] {
                if a.window.overlaps_interior(&b.window) {
                    out.push((label(a), label(b)));
                }
            }
        }
        out
    }
}

/// Read a card as a bin, returning the name it should be registered under.
///
/// `nch.1` with a geometry selector becomes bin 1 of `nch`. A dotted name
/// *without* a selector stays whole: dots are legal in model names, and
/// reinterpreting `my.model` as a bin of `my` would break a deck that never
/// asked for binning.
pub fn classify(name: &str, params: &[(String, f64)]) -> (String, ModelBin) {
    let get = |key: &str| {
        params
            .iter()
            .find(|(k, _)| k.eq_ignore_ascii_case(key))
            .map(|(_, v)| *v)
    };
    let has_selector = SELECTORS.iter().any(|s| get(s).is_some());

    let window = if has_selector {
        Window {
            lmin: get("lmin").unwrap_or(f64::NEG_INFINITY),
            lmax: get("lmax").unwrap_or(f64::INFINITY),
            wmin: get("wmin").unwrap_or(f64::NEG_INFINITY),
            wmax: get("wmax").unwrap_or(f64::INFINITY),
        }
    } else {
        Window::UNBOUNDED
    };

    // Only a selector-carrying card may claim a `.N` suffix as a bin index.
    let (base, index) = match name.rsplit_once('.') {
        Some((base, tail)) if has_selector && !base.is_empty() => match tail.parse::<u32>() {
            Ok(i) => (base.to_string(), Some(i)),
            Err(_) => (name.to_string(), None),
        },
        _ => (name.to_string(), None),
    };

    let params = params
        .iter()
        .filter(|(k, _)| !IS_SELECTOR(k))
        .cloned()
        .collect();

    (
        base,
        ModelBin {
            index,
            window,
            params,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(pairs: &[(&str, f64)]) -> Vec<(String, f64)> {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    fn bin(index: u32, lmin: f64, lmax: f64, wmin: f64, wmax: f64) -> ModelBin {
        ModelBin {
            index: Some(index),
            window: Window {
                lmin,
                lmax,
                wmin,
                wmax,
            },
            params: p(&[("vto", index as f64)]),
        }
    }

    /// A dotted name is a bin only when the card says it is. Dots are legal in
    /// model names, so the selector is what makes the difference.
    #[test]
    fn a_dotted_name_without_a_selector_is_not_a_bin() {
        let (base, b) = classify("my.model", &p(&[("vto", 0.5)]));
        assert_eq!(base, "my.model");
        assert_eq!(b.index, None);
        assert_eq!(b.window, Window::UNBOUNDED);

        let (base, b) = classify("nch.1", &p(&[("lmin", 1e-7), ("lmax", 2e-7)]));
        assert_eq!(base, "nch");
        assert_eq!(b.index, Some(1));
    }

    /// The selectors pick a bin, so they must not survive into the parameters —
    /// every binned card would otherwise warn about four unknown parameters.
    #[test]
    fn the_selectors_do_not_reach_the_model() {
        let (_, b) = classify(
            "nch.1",
            &p(&[
                ("lmin", 1e-7),
                ("LMAX", 2e-7),
                ("wmin", 1e-7),
                ("wmax", 1e-6),
                ("vto", 0.4),
            ]),
        );
        assert_eq!(b.params, p(&[("vto", 0.4)]));
        assert_eq!(b.window.lmax, 2e-7, "a selector is read case-insensitively");
    }

    /// The case the rule exists for: PDKs write contiguous windows, so the
    /// shared endpoint must land in a bin rather than in a gap.
    #[test]
    fn a_geometry_on_a_shared_boundary_selects_deterministically() {
        let mut g = BinGroup::default();
        g.push(bin(1, 0.18e-6, 0.30e-6, 0.22e-6, 1e-6));
        g.push(bin(2, 0.30e-6, 1.00e-6, 0.22e-6, 1e-6));
        // Interior of bin 1, interior of bin 2, and the boundary they share.
        for (l, want) in [(0.25e-6, 1.0), (0.5e-6, 2.0), (0.30e-6, 1.0)] {
            let got = g.select("nch", l, 0.5e-6).expect("a bin covers this");
            assert_eq!(got, &p(&[("vto", want)]), "L={l:e} should take bin {want}");
        }
        assert!(
            g.interior_overlaps().is_empty(),
            "touching at an edge is not an overlap"
        );
    }

    /// Nested bins: an explicit window beats a catch-all card of the same name,
    /// which is what "tightest wins" is for.
    #[test]
    fn a_tight_bin_beats_a_catch_all() {
        let mut g = BinGroup::default();
        g.push(ModelBin {
            index: None,
            window: Window::UNBOUNDED,
            params: p(&[("vto", 99.0)]),
        });
        g.push(bin(1, 0.18e-6, 0.30e-6, 0.22e-6, 1e-6));
        let got = g.select("nch", 0.25e-6, 0.5e-6).unwrap();
        assert_eq!(
            got,
            &p(&[("vto", 1.0)]),
            "the bounded bin is the tighter one"
        );
    }

    /// Outside every window is a hard error naming the geometry and the windows.
    /// The alternative is picking the nearest bin, which is a wrong answer.
    #[test]
    fn geometry_outside_every_bin_is_an_error_that_names_the_windows() {
        let mut g = BinGroup::default();
        g.push(bin(1, 0.18e-6, 0.30e-6, 0.22e-6, 1e-6));
        g.push(bin(2, 0.30e-6, 1.00e-6, 0.22e-6, 1e-6));
        let err = g
            .select("nch", 5e-6, 0.5e-6)
            .expect_err("5 um is past LMAX");
        let msg = err.to_string();
        for needle in ["nch", "5.0000e-6", "nch.1", "nch.2", "1.8000e-7"] {
            assert!(msg.contains(needle), "error should name {needle}: {msg}");
        }
    }

    /// A genuine interior overlap is a malformed card set, and the choice there
    /// does change the answer, so it is reported rather than resolved silently.
    #[test]
    fn interiors_that_overlap_are_reported() {
        let mut g = BinGroup::default();
        g.push(bin(1, 0.18e-6, 0.50e-6, 0.22e-6, 1e-6));
        g.push(bin(2, 0.30e-6, 1.00e-6, 0.22e-6, 1e-6));
        assert_eq!(
            g.interior_overlaps(),
            vec![(".1".to_string(), ".2".to_string())]
        );
    }

    /// `unbinned` is the escape hatch for a device family with no geometry, and
    /// it must refuse a binned name rather than hand back bin one.
    #[test]
    fn unbinned_refuses_a_binned_group() {
        let mut g = BinGroup::default();
        g.push(bin(1, 0.18e-6, 0.30e-6, 0.22e-6, 1e-6));
        assert!(g.unbinned().is_none());

        let mut plain = BinGroup::default();
        plain.push(ModelBin {
            index: None,
            window: Window::UNBOUNDED,
            params: p(&[("is", 1e-14)]),
        });
        assert_eq!(plain.unbinned(), Some(p(&[("is", 1e-14)]).as_slice()));
    }
}
