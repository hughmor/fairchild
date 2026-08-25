//! The `.alter` × `.temp` grid: what "every corner this deck declares" means.
//!
//! One definition, because two frontends now need it. The CLI has always
//! expanded this grid; `Circuit.run_all()` expands the same one, and a second
//! copy would be a second opinion about whether a deck with one `.temp` card is
//! one corner or two — the sort of disagreement #33 exists to stop.
//!
//! ## Why a single `.temp` is not a sweep
//!
//! `SimOptions::from_netlist` already folds `temps.first()` into `temp_k`, so a
//! deck with one `.temp 75` is *already* at 75 °C before any corner expansion
//! happens. Treating it as a one-point sweep here would be harmless but
//! misleading — it would report a corner axis to a caller that has none. Only
//! two or more temperatures make a sweep.
//!
//! ## What a corner is not
//!
//! It is not a *run*. A corner is a resolved `(netlist, options)` pair — what to
//! simulate, not which analyses to simulate. Which analyses run is
//! `netlist.analyses`, and the same list runs at every corner.

use fairchild_parser::Netlist;

use crate::options::SimOptions;

/// One leaf of the `.alter` × `.temp` grid: a fully-resolved netlist plus a
/// `SimOptions` carrying that corner's temperature.
#[derive(Clone, Debug)]
pub struct Corner {
    /// Position in the alter axis; `0` is the base netlist.
    pub alter_idx: usize,
    /// Position in the temperature axis.
    pub temp_idx: usize,
    /// `"base"`, or the `.alter` block's label.
    pub alter_label: String,
    pub temp_k: f64,
    /// The netlist with this corner's `.alter` overrides applied.
    pub netlist: Netlist,
    /// `opts` with this corner's temperature.
    pub opts: SimOptions,
}

impl Corner {
    /// Temperature in °C, which is the unit every user-facing surface quotes.
    pub fn temp_c(&self) -> f64 {
        self.temp_k - 273.15
    }
}

/// The expanded grid, with the axis lengths that decide whether a corner needs
/// naming at all — a single-corner run should not grow a `.alter_base.temp_27c`
/// suffix.
#[derive(Clone, Debug)]
pub struct CornerGrid {
    pub corners: Vec<Corner>,
    /// Base run plus one per `.alter` block, so always ≥ 1.
    pub n_alters: usize,
    /// `≥ 2` only when the deck declares a real temperature sweep.
    pub n_temps: usize,
}

/// Expand `netlist`'s `.alter` blocks and `.temp` sweep into the full grid.
///
/// Always yields at least one corner (`alter_label = "base"`), so a caller
/// never has to special-case a deck that declares no corners at all.
pub fn expand_corners(netlist: &Netlist, opts: &SimOptions) -> CornerGrid {
    // A single `.temp` is already in `opts.temp_k`; see the module note.
    let temp_sweep: Vec<f64> = if netlist.temps.len() > 1 {
        netlist.temps.clone()
    } else {
        vec![opts.temp_k]
    };

    let mut alter_runs: Vec<(String, Netlist)> = vec![("base".into(), netlist.clone())];
    for block in &netlist.alters {
        let mut patched = netlist.clone();
        patched.apply_alter(block);
        alter_runs.push((block.label.clone(), patched));
    }

    let mut corners = Vec::with_capacity(alter_runs.len() * temp_sweep.len());
    for (ai, (alter_label, alter_netlist)) in alter_runs.iter().enumerate() {
        for (ti, &temp_k) in temp_sweep.iter().enumerate() {
            let mut corner_opts = opts.clone();
            corner_opts.temp_k = temp_k;
            corners.push(Corner {
                alter_idx: ai,
                temp_idx: ti,
                alter_label: alter_label.clone(),
                temp_k,
                netlist: alter_netlist.clone(),
                opts: corner_opts,
            });
        }
    }

    CornerGrid {
        n_alters: alter_runs.len(),
        n_temps: temp_sweep.len(),
        corners,
    }
}
