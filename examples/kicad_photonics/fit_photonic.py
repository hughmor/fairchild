#!/usr/bin/env python3
"""
fit_photonic.py — add-drop ring resonator parameter fitting for fairchild.

Fits PN + thermal phase-shifter parameters to experimental transmission
spectra from a lightlab NdSweeper dataset.

Topology: add-drop MRM
---------------------------------------------------------------------------
  IN ──► WG1 ──► CPL2(bus) ──► WG3 ──► THRU
                  │                     │
               CPL2_a2 ◄─── PS1 ◄─── CPL1_b1
               CPL2_b2 ──► PS2  ──► CPL1_a1
                                     CPL1_b2 ──► WG2 ──► DROP
---------------------------------------------------------------------------

Supported models (both have PN junction + heater):
  fc_pn_th_ps       — linear EO (L1)
  fc_pn_th_ps_cap   — depletion-mode with C_j(V) + da/dV (L2)

Staged fitting strategy:
  Stage 1 — passive  : fit n_eff, alpha_dB_cm, kappa_l from 0 V / 0 mA spectrum
  Stage 2 — thermal  : fit p_pi_th from heater current sweeps
  Stage 3 — EO       : fit v_pi_l from junction voltage sweeps
  Pre-fit            : r_heater from heater I/V; g_pn from junction I/V

Setup
-----
  cd /path/to/fairchild
  maturin develop --release -m crates/fairchild-py/Cargo.toml
  source .venv/bin/activate
"""

#%% ── imports ────────────────────────────────────────────────────────────────

import sys
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.optimize import differential_evolution, minimize

try:
    import matplotlib
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

try:
    import fairchild as fc
except ImportError:
    sys.exit(
        "fairchild Python extension not found.\n"
        "Build it with:  maturin develop --release -m crates/fairchild-py/Cargo.toml"
    )

import lightlab as ll
from lightlab.util.data import Spectrum
from lightlab.util.sweep import NdSweeper


#%% ── dataset path & window ──────────────────────────────────────────────────

DATA_PATH = str(
    Path(__file__).parents[2] / "data" / "giona_neuron2_mod_joint_IV_spec"
)

# Restrict to the resonance window of our target ring (nm).
WL_LO, WL_HI = 1545.8, 1547.5
# Number of wavelength points for simulation (downsampled from measured).
N_SIM_PTS = 100


#%% ── data loading helpers ───────────────────────────────────────────────────

def load_sweep(path: str = DATA_PATH) -> NdSweeper:
    """Load the NdSweeper dataset from disk."""
    sweep = NdSweeper()
    sweep.load(path)
    return sweep


def _trim_spectrum(spec: Spectrum, lo: float = WL_LO, hi: float = WL_HI):
    """Return (wl_nm, T_dB) trimmed to [lo, hi] nm window."""
    wl = np.asarray(spec.absc, dtype=float)
    T  = np.asarray(spec.ordi, dtype=float)
    mask = (wl >= lo) & (wl <= hi)
    return wl[mask], T[mask]


def _downsample(wl: np.ndarray, T: np.ndarray, n: int = N_SIM_PTS):
    """Uniformly downsample (wl, T) to n points via linear interpolation."""
    wl_new = np.linspace(wl[0], wl[-1], n)
    T_new  = np.interp(wl_new, wl, T)
    return wl_new, T_new


def _sweep_axes(sweep: NdSweeper):
    """
    Return (hc_mA, jv_V) 1-D axis arrays from a loaded NdSweeper.

    NdSweeper.save() stores each actuated variable as an N-D grid in
    sweep.data, so axis values are recovered by slicing one column/row.
    """
    hc_grid = np.asarray(sweep.data["Heater Current (mA)"], dtype=float)
    jv_grid = np.asarray(sweep.data["Junction Voltage (V)"], dtype=float)
    # HC varies along axis-0 (rows), JV along axis-1 (cols).
    hc_mA = hc_grid[:, 0]
    jv_V  = jv_grid[0, :]
    return hc_mA, jv_V


@dataclass
class SweepData:
    """Extracted and pre-processed measurements from the NdSweeper dataset."""
    # 1-D axis arrays
    hc_mA: np.ndarray          # heater current (mA), shape (n_hc,)
    jv_V:  np.ndarray          # junction voltage (V), shape (n_jv,)
    # Trimmed + downsampled spectra grids — shape (n_hc, n_jv, N_SIM_PTS)
    wl_nm:   np.ndarray        # shared wavelength axis (nm), shape (N_SIM_PTS,)
    T_dB:    np.ndarray        # transmission (dB), shape (n_hc, n_jv, N_SIM_PTS)
    # Electrical measurements — shape (n_hc, n_jv)
    v_heat_V:  np.ndarray      # heater voltage (V)
    i_junc_mA: np.ndarray      # junction current (mA)


def extract_data(sweep: NdSweeper) -> SweepData:
    """Pre-process the NdSweeper into arrays ready for fitting."""
    hc_mA, jv_V = _sweep_axes(sweep)
    spectra   = sweep.data["Spectrum"]
    v_heat_V  = np.asarray(sweep.data["Heater Voltage (V)"],    dtype=float)
    i_junc_mA = np.asarray(sweep.data["Junction Current (mA)"], dtype=float)

    n_hc, n_jv = len(hc_mA), len(jv_V)
    # Use the first valid spectrum to set the shared wavelength axis.
    wl0, _ = _trim_spectrum(spectra[0, 0])
    wl_sim, _ = _downsample(wl0, wl0)   # just the axis

    T_dB = np.empty((n_hc, n_jv, N_SIM_PTS), dtype=float)
    for i in range(n_hc):
        for j in range(n_jv):
            wl, T = _trim_spectrum(spectra[i, j])
            _, T_ds = _downsample(wl, T)
            T_dB[i, j] = T_ds

    return SweepData(
        hc_mA=hc_mA,
        jv_V=jv_V,
        wl_nm=wl_sim,
        T_dB=T_dB,
        v_heat_V=v_heat_V,
        i_junc_mA=i_junc_mA,
    )


#%% ── parameter specification ────────────────────────────────────────────────

@dataclass
class ParamSpec:
    """One tunable or fixed simulation parameter."""
    name: str
    value: float
    fixed: bool = False
    bounds: tuple = (0.0, 1.0)
    description: str = ""

    def copy(self, **overrides):
        p = deepcopy(self)
        for k, v in overrides.items():
            setattr(p, k, v)
        return p


_WAVEGUIDE = [
    ParamSpec("n_g",         4.2,       fixed=True,  bounds=(3.5, 5.5),
              description="group index (from waveguide simulation)"),
    # n_eff=2.403 puts the resonance near 1546.5 nm for this ring geometry
    ParamSpec("n_eff",       2.403,     fixed=False, bounds=(2.0, 3.5),
              description="effective index at lambda_ref (sets resonance position)"),
    # l_m = pi * r = pi * 8.078e-6 m per arm (half ring circumference)
    ParamSpec("l_m",         25.378e-6, fixed=True,  bounds=(1e-6, 1e-3),
              description="PS arm length (m) — half ring circumference"),
    ParamSpec("alpha_db_cm", 10.0,      fixed=False, bounds=(0.5, 80.0),
              description="propagation loss (dB/cm)"),
    # MUST be 0 for ring resonators — pin_at_ref=1 suppresses resonances.
    ParamSpec("pin_at_ref",  0.0,       fixed=True,  bounds=(0.0, 1.0),
              description="phase reference mode (0=absolute, 1=pin-at-ref)"),
]

_PN_EO = [
    ParamSpec("v_pi_l",  0.02,  fixed=False, bounds=(1e-3, 0.5),
              description="Vpi*L (V*m) — EO tuning efficiency"),
    ParamSpec("g_pn",    1e-3,  fixed=True,  bounds=(1e-9, 1e3),
              description="PN junction ohmic conductance (S)"),
]

_PN_CAP = [
    ParamSpec("c_j0",  20e-15, fixed=False, bounds=(1e-16, 1e-12),
              description="zero-bias junction capacitance (F)"),
    ParamSpec("v_bi",  0.7,    fixed=False, bounds=(0.3, 1.4),
              description="built-in voltage (V)"),
    ParamSpec("m_j",   0.5,    fixed=False, bounds=(0.2, 0.95),
              description="junction grading coefficient"),
    ParamSpec("da_dv", 0.0,    fixed=True,  bounds=(0.0, 50.0),
              description="bias-dependent loss slope (Np/m/V)"),
]

_HEATER = [
    # r_heater measured from heater I/V: R_total=368.2 Ohm → R_arm=184.1 Ohm
    ParamSpec("r_heater", 184.1,  fixed=True,  bounds=(50.0, 2000.0),
              description="heater resistance per arm (Ohm)"),
    # p_pi_th estimated from thermal tuning: ~52 mW/arm at 2 mA
    ParamSpec("p_pi_th",  52e-3,  fixed=False, bounds=(1e-3, 0.5),
              description="heater pi-power (W) — thermal tuning efficiency"),
]

_COUPLER = [
    # kappa_L=0.15 gives ~17 dB extinction at the measured ring loss level
    ParamSpec("kappa_l", 0.15, fixed=False, bounds=(0.001, 0.99),
              description="coupler kappa*L (shared for both CPL1 and CPL2)"),
]

MODELS: dict[str, dict] = {
    "fc_pn_th_ps": {
        "description": "PN + thermal — combined, linear EO (L1)",
        "ps_params":   _WAVEGUIDE + _PN_EO + _HEATER,
        "has_pn":      True,
        "has_heater":  True,
        "has_cap":     False,
    },
    "fc_pn_th_ps_cap": {
        "description": "PN + thermal — depletion-mode with C_j(V) + da/dV (L2)",
        "ps_params":   _WAVEGUIDE + _PN_EO + _PN_CAP + _HEATER,
        "has_pn":      True,
        "has_heater":  True,
        "has_cap":     True,
    },
}


#%% ── netlist builder ─────────────────────────────────────────────────────────

# fc_pn_th_ps terminal order (from phase_shifters.rs):
#   in, out, anode, cathode, heat_p, heat_n
# anode=GND, cathode=PN_BIAS → forward convention: V_pn = PN_BIAS (positive = reverse bias)
# heat_p / heat_n filled per-arm for series heater wiring.
_PS_LINE_TEMPLATES = {
    "fc_pn_th_ps":     "X{name} {in_} {out} GND PN_BIAS {heat_p} {heat_n} fc_pn_th_ps {params}",
    "fc_pn_th_ps_cap": "X{name} {in_} {out} GND PN_BIAS {heat_p} {heat_n} fc_pn_th_ps_cap {params}",
}

# Ring arm optical connections: CPL1_b1 → PS1 → CPL2_a2, CPL2_b2 → PS2 → CPL1_a1
_ARM_CONNECTIONS = [
    ("CPL1_b1", "CPL2_a2"),  # arm 1 (upper half-ring)
    ("CPL2_b2", "CPL1_a1"),  # arm 2 (lower half-ring)
]

# Series heater wiring: current source drives HEAT_BIAS → PS1 → HEAT_MID → PS2 → GND
_HEATER_CONNECTIONS = [
    ("HEAT_BIAS", "HEAT_MID"),  # arm 1
    ("HEAT_MID",  "GND"),       # arm 2
]


def _params_to_spice(params: dict[str, float]) -> str:
    return " ".join(f"{k}={v:.8g}" for k, v in params.items())


def build_netlist(
    model: str,
    ps_params: dict[str, float],
    coupler_kappa_l: float,
    wavelength_nm: float,
    v_pn: float = 0.0,
    i_heat: float = 0.0,
) -> str:
    """
    Build a complete SPICE netlist for the add-drop ring resonator.

    Optical power is measured at CWL_OUT (laser) and GC_THRU_OUT (after
    output GC) so T_dB = 10*log10(P_GC_THRU_OUT / P_CWL_OUT) matches the
    experimental fibre-to-fibre transmission convention.

    Parameters
    ----------
    model          : one of the keys in MODELS
    ps_params      : fairchild SPICE parameter name → value (both arms share)
    coupler_kappa_l: kappa*L for both directional couplers
    wavelength_nm  : CW laser wavelength (nm)
    v_pn           : PN reverse-bias voltage (V)
    i_heat         : heater drive current (A); series through both arms
    """
    ps_p = {k: v for k, v in ps_params.items() if k != "kappa_l"}
    ps_spice = _params_to_spice(ps_p)

    tmpl = _PS_LINE_TEMPLATES[model]
    ps_lines = []
    for idx, ((in_, out), (heat_p, heat_n)) in enumerate(
        zip(_ARM_CONNECTIONS, _HEATER_CONNECTIONS), start=1
    ):
        ps_lines.append(tmpl.format(
            name=f"PS{idx}", in_=in_, out=out, params=ps_spice,
            heat_p=heat_p, heat_n=heat_n,
        ))

    optical_ports = [
        "CWL_OUT", "GC_THRU_OUT", "GC_DROP_OUT",
        "IN_OPT", "THRU_OPT", "DROP_OPT",
        "CPL1_a1", "CPL1_a2", "CPL1_b1", "CPL1_b2",
        "CPL2_a1", "CPL2_a2", "CPL2_b1", "CPL2_b2",
    ]

    return "\n".join([
        *[f".optical_port {p}" for p in optical_ports],
        ".op",
        f"XCWL1 CWL_OUT fc_cw_laser power_mW=1.0 wavelength_nm={wavelength_nm:.6f}",
        "V1 VDD GND DC 2",
        f"V_PN  PN_BIAS  GND DC {v_pn:.6g}",
        # Series heater: conventional current flows from HEAT_BIAS through
        # arm-1 heater, then arm-2 heater, to GND.
        f"IHEAT HEAT_BIAS GND DC {i_heat:.6g}",
        # Grating couplers — loss only; we read optical power on both sides
        "XGC_IN   CWL_OUT  IN_OPT      fc_grating_coupler alpha_dB=9.0",
        "XGC_THRU THRU_OPT GC_THRU_OUT fc_grating_coupler alpha_dB=9.0",
        "XGC_DROP DROP_OPT GC_DROP_OUT fc_grating_coupler alpha_dB=9.0",
        # Bus waveguides
        "XWG1 IN_OPT   CPL2_a1  fc_waveguide l_m=304.5e-6 n_g=4.2 alpha_dB_cm=1.0",
        "XWG2 CPL1_b2  DROP_OPT fc_waveguide l_m=304.5e-6 n_g=4.2 alpha_dB_cm=1.0",
        "XWG3 CPL2_b1  THRU_OPT fc_waveguide l_m=304.5e-6 n_g=4.2 alpha_dB_cm=1.0",
        # Directional couplers (CPL2=bus, CPL1=drop)
        f"XCPL1 CPL1_a1 CPL1_a2 CPL1_b1 CPL1_b2 fc_dcoupler kappa_L={coupler_kappa_l:.8g}",
        f"XCPL2 CPL2_a1 CPL2_a2 CPL2_b1 CPL2_b2 fc_dcoupler kappa_L={coupler_kappa_l:.8g}",
        *ps_lines,
        ".end",
    ]) + "\n"


#%% ── simulation helpers ─────────────────────────────────────────────────────

def _optical_power(r, net: str) -> float:
    """Compute optical power (amplitude²) from the DC-OP result at optical bundle net.

    Fairchild expands each .optical_port NET into nodes NET_re_0, NET_im_0,
    NET_wl_0, accessed via r["V(net_re_0)"] (lowercase net name).
    """
    n = net.lower()
    try:
        re = float(r[f"V({n}_re_0)"][0])
        im = float(r[f"V({n}_im_0)"][0])
    except Exception:
        re, im = 0.0, 0.0
    return re * re + im * im


def wavelength_sweep(
    model: str,
    ps_params: dict[str, float],
    coupler_kappa_l: float,
    wavelengths_nm: np.ndarray,
    v_pn: float = 0.0,
    i_heat: float = 0.0,
) -> np.ndarray:
    """
    DC-OP sweep over wavelengths; return T_dB = 10·log10(P_out / P_in).

    Uses optical power at CWL_OUT (before input GC) and GC_THRU_OUT (after
    output GC) to match the fibre-to-fibre experimental convention.
    """
    netlist0 = build_netlist(model, ps_params, coupler_kappa_l,
                             wavelengths_nm[0], v_pn, i_heat)
    ckt = fc.Circuit()
    ckt.load_str(netlist0)

    t_dB = np.empty(len(wavelengths_nm))
    for i, wl in enumerate(wavelengths_nm):
        ckt.set_param("XCWL1", "wavelength_nm", wl)
        ckt.set_param("V_PN",  "dc", v_pn)
        ckt.set_param("IHEAT", "dc", i_heat)
        try:
            r = ckt.run("op")
            p_in  = max(_optical_power(r, "CWL_OUT"),      1e-30)
            p_out = max(_optical_power(r, "GC_THRU_OUT"),  1e-30)
            t_dB[i] = 10.0 * np.log10(p_out / p_in)
        except Exception:
            t_dB[i] = -100.0
    return t_dB


def _find_resonance(wl: np.ndarray, t_dB: np.ndarray) -> float:
    """Wavelength of the deepest dip via quadratic interpolation."""
    idx = int(np.argmin(t_dB))
    if 1 <= idx < len(t_dB) - 1:
        y0, y1, y2 = t_dB[idx - 1], t_dB[idx], t_dB[idx + 1]
        d = 2 * y1 - y0 - y2
        if abs(d) > 1e-12:
            frac = (y0 - y2) / (2 * d)
            dw = wl[idx] - wl[max(idx - 1, 0)]
            return float(wl[idx] + frac * dw)
    return float(wl[idx])


def _normalise_spectrum(t_dB: np.ndarray, percentile: float = 95.0) -> np.ndarray:
    """Subtract the background level (high percentile) so peak → 0 dB."""
    bg = np.percentile(t_dB, percentile)
    return t_dB - bg


def _spectrum_loss(sim_t: np.ndarray, meas_t: np.ndarray) -> float:
    """
    Mean-squared error on background-normalised spectra.

    Both are shifted so their background (95th percentile) is 0 dB, which
    removes sensitivity to absolute insertion loss (GC coupling efficiency).
    """
    s = _normalise_spectrum(sim_t)
    m = _normalise_spectrum(meas_t)
    return float(np.mean((s - m) ** 2))


#%% ── parameter vector helpers ───────────────────────────────────────────────

def _free_specs(specs: list[ParamSpec]) -> list[ParamSpec]:
    return [s for s in specs if not s.fixed]


def _pack(free: list[ParamSpec]) -> np.ndarray:
    return np.array([s.value for s in free])


def _unpack(x: np.ndarray, specs: list[ParamSpec]) -> dict[str, float]:
    result, xi = {}, 0
    for s in specs:
        if s.fixed:
            result[s.name] = s.value
        else:
            result[s.name] = float(x[xi])
            xi += 1
    return result


def _bounds(specs: list[ParamSpec]) -> list[tuple]:
    return [s.bounds for s in specs if not s.fixed]


#%% ── pre-fitting from electrical I/V ────────────────────────────────────────

def prefit_r_heater(sd: SweepData) -> float:
    """
    Fit heater resistance from heater I/V.

    V_heat = 2 * R_arm * I_heat (two arms in series), so
    R_arm = slope(V_heat vs I_heat) / 2.
    """
    hc_A = sd.hc_mA * 1e-3                          # convert to amperes
    v_heat_mean = sd.v_heat_V.mean(axis=1)           # average over JV axis

    # Only use points where |I| > 0.1 mA to avoid near-zero noise.
    mask = np.abs(hc_A) > 1e-4
    if mask.sum() < 2:
        return 184.1                                 # fallback to known value

    slope = float(np.polyfit(hc_A[mask], v_heat_mean[mask], 1)[0])
    r_arm = slope / 2.0
    print(f"  Pre-fit R_heater: {slope:.1f} Ω total → {r_arm:.1f} Ω/arm")
    return max(r_arm, 10.0)


def prefit_g_pn(sd: SweepData) -> float:
    """
    Fit PN junction conductance from the forward-bias I/V slope.

    Uses junction current / junction voltage in a linear regime near 0 V.
    """
    i_junc_A = sd.i_junc_mA.mean(axis=0) * 1e-3    # average over HC axis
    jv_V = sd.jv_V

    # Linear region: 0 V to +0.4 V (before strong injection)
    mask = (jv_V >= 0.0) & (jv_V <= 0.4)
    if mask.sum() < 2:
        return 1e-3

    slope = float(np.polyfit(jv_V[mask], i_junc_A[mask], 1)[0])
    g_pn = max(slope, 1e-9)
    print(f"  Pre-fit g_pn: {g_pn:.3e} S (from forward I/V slope)")
    return g_pn


#%% ── staged fitting ─────────────────────────────────────────────────────────

def _run_de(objective, specs, x0, label="", maxiter=200, popsize=12, seed=0, verbose=True):
    """Thin wrapper around differential_evolution; always returns a param dict."""
    bounds = _bounds(specs)
    free = _free_specs(specs)
    if not free:
        if verbose:
            print(f"  {label}: no free parameters — skipping.")
        return _unpack(np.array([]), specs)

    if verbose:
        print(f"  {label}: fitting {[s.name for s in free]}")
        print(f"    differential_evolution (maxiter={maxiter}, popsize={popsize})")

    count = [0]
    def wrapped(x):
        count[0] += 1
        return objective(x)

    result = differential_evolution(
        wrapped, bounds,
        x0=x0,
        maxiter=maxiter, popsize=popsize, seed=seed,
        tol=1e-6, mutation=(0.5, 1.5), recombination=0.7,
        polish=True, disp=False,
    )
    if verbose:
        print(f"    converged={result.success}  loss={result.fun:.4e}"
              f"  ({count[0]} evals)")
    return _unpack(result.x, specs)


def stage1_passive(
    model: str,
    sd: SweepData,
    base_specs: Optional[list[ParamSpec]] = None,
    coupler_specs: Optional[list[ParamSpec]] = None,
    maxiter: int = 300,
    popsize: int = 15,
    verbose: bool = True,
) -> dict[str, float]:
    """
    Stage 1: fit n_eff, alpha_dB_cm, kappa_l from the passive spectrum
    (heater current ≈ 0 mA, junction voltage ≈ 0 V).

    EO and thermal parameters are held fixed.
    """
    if base_specs is None:
        base_specs = deepcopy(MODELS[model]["ps_params"])
    if coupler_specs is None:
        coupler_specs = deepcopy(_COUPLER)

    # Fix EO and thermal; free only passive optical params + kappa_l
    freeze = {"v_pi_l", "g_pn", "c_j0", "v_bi", "m_j", "da_dv",
              "r_heater", "p_pi_th", "tau_th"}
    for s in base_specs:
        if s.name in freeze:
            s.fixed = True
    for s in coupler_specs:
        s.fixed = False   # kappa_l is free in stage 1

    # Passive spectrum: closest to 0 mA and 0 V
    i0 = int(np.argmin(np.abs(sd.hc_mA)))
    j0 = int(np.argmin(np.abs(sd.jv_V)))
    meas_T = sd.T_dB[i0, j0]
    wls    = sd.wl_nm

    combined = base_specs + coupler_specs

    def objective(x):
        params = _unpack(x, combined)
        kl     = params.pop("kappa_l")
        sim_T  = wavelength_sweep(model, params, kl, wls)
        return _spectrum_loss(sim_T, meas_T)

    x0   = _pack(_free_specs(combined))
    best = _run_de(objective, combined, x0,
                   label="Stage 1 (passive)", maxiter=maxiter,
                   popsize=popsize, verbose=verbose)

    if verbose:
        ps_p = {k: v for k, v in best.items() if k != "kappa_l"}
        res  = _find_resonance(wls, wavelength_sweep(
            model, ps_p, best["kappa_l"], wls))
        print(f"    Fitted resonance: {res:.4f} nm")
    return best


def stage2_thermal(
    model: str,
    sd: SweepData,
    stage1_result: dict[str, float],
    base_specs: Optional[list[ParamSpec]] = None,
    coupler_specs: Optional[list[ParamSpec]] = None,
    n_hc_points: int = 8,
    maxiter: int = 200,
    popsize: int = 12,
    verbose: bool = True,
) -> dict[str, float]:
    """
    Stage 2: fit p_pi_th from heater current sweep (JV ≈ 0 V).

    n_eff, alpha_dB_cm, kappa_l are fixed at Stage-1 values.
    """
    if base_specs is None:
        base_specs = deepcopy(MODELS[model]["ps_params"])
    if coupler_specs is None:
        coupler_specs = deepcopy(_COUPLER)

    # Fix everything except p_pi_th (and tau_th if present)
    for s in base_specs + coupler_specs:
        s.fixed = True
        if s.name in stage1_result:
            s.value = stage1_result[s.name]

    thermal_free = {"p_pi_th", "tau_th"}
    for s in base_specs:
        if s.name in thermal_free:
            s.fixed = False

    # Select a subset of heater currents (skip 0-crossing, use positive side)
    hc_all = sd.hc_mA
    hc_pos = hc_all[hc_all > 0]
    if len(hc_pos) == 0:
        hc_pos = hc_all[hc_all != 0]
    n_pts = min(n_hc_points, len(hc_pos))
    hc_sel = hc_pos[np.round(np.linspace(0, len(hc_pos) - 1, n_pts)).astype(int)]

    j0 = int(np.argmin(np.abs(sd.jv_V)))
    wls = sd.wl_nm

    combined = base_specs + coupler_specs

    def objective(x):
        params = _unpack(x, combined)
        kl     = params.pop("kappa_l")
        loss   = 0.0
        for hc in hc_sel:
            i_hc   = int(np.argmin(np.abs(sd.hc_mA - hc)))
            meas_T = sd.T_dB[i_hc, j0]
            sim_T  = wavelength_sweep(model, params, kl, wls, i_heat=hc * 1e-3)
            loss  += _spectrum_loss(sim_T, meas_T)
        return loss / len(hc_sel)

    x0   = _pack(_free_specs(combined))
    best = {**stage1_result,
            **_run_de(objective, combined, x0,
                      label="Stage 2 (thermal)", maxiter=maxiter,
                      popsize=popsize, verbose=verbose)}
    return best


def stage3_eo(
    model: str,
    sd: SweepData,
    stage2_result: dict[str, float],
    base_specs: Optional[list[ParamSpec]] = None,
    coupler_specs: Optional[list[ParamSpec]] = None,
    n_jv_points: int = 8,
    maxiter: int = 200,
    popsize: int = 12,
    verbose: bool = True,
) -> dict[str, float]:
    """
    Stage 3: fit v_pi_l (and optionally g_pn) from junction voltage sweep
    (heater current ≈ 0 mA).

    Uses only reverse-bias (JV ≤ 0) to avoid injection-regime complications.
    n_eff, alpha_dB_cm, kappa_l, p_pi_th are fixed at prior-stage values.
    """
    if base_specs is None:
        base_specs = deepcopy(MODELS[model]["ps_params"])
    if coupler_specs is None:
        coupler_specs = deepcopy(_COUPLER)

    for s in base_specs + coupler_specs:
        s.fixed = True
        if s.name in stage2_result:
            s.value = stage2_result[s.name]

    eo_free = {"v_pi_l"}
    for s in base_specs:
        if s.name in eo_free:
            s.fixed = False

    # Reverse-bias subset of junction voltages
    jv_rev = sd.jv_V[sd.jv_V <= 0]
    if len(jv_rev) == 0:
        print("  Stage 3: no reverse-bias JV points — skipping.")
        return stage2_result
    n_pts = min(n_jv_points, len(jv_rev))
    jv_sel = jv_rev[np.round(np.linspace(0, len(jv_rev) - 1, n_pts)).astype(int)]

    i0 = int(np.argmin(np.abs(sd.hc_mA)))
    wls = sd.wl_nm

    combined = base_specs + coupler_specs

    def objective(x):
        params = _unpack(x, combined)
        kl     = params.pop("kappa_l")
        loss   = 0.0
        for jv in jv_sel:
            j_jv   = int(np.argmin(np.abs(sd.jv_V - jv)))
            meas_T = sd.T_dB[i0, j_jv]
            sim_T  = wavelength_sweep(model, params, kl, wls, v_pn=jv)
            loss  += _spectrum_loss(sim_T, meas_T)
        return loss / len(jv_sel)

    x0   = _pack(_free_specs(combined))
    best = {**stage2_result,
            **_run_de(objective, combined, x0,
                      label="Stage 3 (EO)", maxiter=maxiter,
                      popsize=popsize, verbose=verbose)}
    return best


def fit_staged(
    model: str = "fc_pn_th_ps",
    path: str = DATA_PATH,
    verbose: bool = True,
) -> dict[str, float]:
    """
    Full staged fitting pipeline.

    Loads data, pre-fits electrical parameters, then runs Stages 1–3.
    Returns the final best-fit parameter dict.
    """
    print(f"\n{'='*60}")
    print(f"Staged fitting: model={model}")
    print(f"Data: {path}")
    print(f"{'='*60}")

    print("\nLoading dataset …")
    sweep = load_sweep(path)
    sd = extract_data(sweep)
    print(f"  Grid: {len(sd.hc_mA)} HC × {len(sd.jv_V)} JV")
    print(f"  HC range: {sd.hc_mA.min():.2f} … {sd.hc_mA.max():.2f} mA")
    print(f"  JV range: {sd.jv_V.min():.2f} … {sd.jv_V.max():.2f} V")
    print(f"  Wavelength window: {sd.wl_nm[0]:.3f} … {sd.wl_nm[-1]:.3f} nm"
          f" ({N_SIM_PTS} pts)")

    # ── pre-fits from electrical data ──────────────────────────────────────
    print("\nPre-fitting electrical parameters …")
    r_heater = prefit_r_heater(sd)
    g_pn     = prefit_g_pn(sd)

    # Apply pre-fit values to all relevant ParamSpec lists
    def _apply_prefit(specs):
        for s in specs:
            if s.name == "r_heater":
                s.value = r_heater
                s.fixed = True
            elif s.name == "g_pn":
                s.value = g_pn
                s.fixed = True
        return specs

    base1 = _apply_prefit(deepcopy(MODELS[model]["ps_params"]))
    cpl1  = deepcopy(_COUPLER)

    # ── Stage 1: passive ──────────────────────────────────────────────────
    print("\nStage 1 — passive (n_eff, alpha_dB_cm, kappa_l) …")
    s1 = stage1_passive(model, sd, base_specs=base1, coupler_specs=cpl1,
                        verbose=verbose)

    # ── Stage 2: thermal ──────────────────────────────────────────────────
    base2 = _apply_prefit(deepcopy(MODELS[model]["ps_params"]))
    cpl2  = deepcopy(_COUPLER)
    print("\nStage 2 — thermal (p_pi_th) …")
    s2 = stage2_thermal(model, sd, s1, base_specs=base2, coupler_specs=cpl2,
                        verbose=verbose)

    # ── Stage 3: EO ───────────────────────────────────────────────────────
    base3 = _apply_prefit(deepcopy(MODELS[model]["ps_params"]))
    cpl3  = deepcopy(_COUPLER)
    print("\nStage 3 — EO (v_pi_l) …")
    s3 = stage3_eo(model, sd, s2, base_specs=base3, coupler_specs=cpl3,
                   verbose=verbose)

    print_results(s3)
    return s3


#%% ── result reporting & plotting ────────────────────────────────────────────

def print_results(best: dict[str, float]):
    print("\n── Best-fit parameters ──────────────────────────────────")
    for k, v in best.items():
        print(f"  {k:<18s} = {v:>14.6g}")
    print()


def plot_ring_fit(
    model: str,
    best: dict[str, float],
    sd: SweepData,
    out_path: str = "ring_fit.png",
    n_panels: int = 4,
):
    """
    Plot measured vs fitted spectra for a selection of bias conditions.

    Panel layout:
      [0] passive (0 mA, 0 V)
      [1] max heater current
      [2] max reverse bias
      [3] max forward bias
    """
    if not _HAS_MPL:
        print("matplotlib not available — skipping plot.")
        return

    kl = best.get("kappa_l", _COUPLER[0].value)
    wls = sd.wl_nm

    i0 = int(np.argmin(np.abs(sd.hc_mA)))
    j0 = int(np.argmin(np.abs(sd.jv_V)))
    i_max = int(np.argmax(sd.hc_mA))
    j_rev = int(np.argmin(sd.jv_V))    # most negative (max reverse bias)
    j_fwd = int(np.argmax(sd.jv_V))    # most positive (max forward bias)

    cases = [
        ("passive (0 mA, 0 V)",
         dict(v_pn=0.0, i_heat=0.0), i0, j0),
        (f"thermal ({sd.hc_mA[i_max]:.1f} mA, 0 V)",
         dict(v_pn=0.0, i_heat=sd.hc_mA[i_max] * 1e-3), i_max, j0),
        (f"EO reverse ({sd.jv_V[j_rev]:.2f} V, 0 mA)",
         dict(v_pn=sd.jv_V[j_rev], i_heat=0.0), i0, j_rev),
        (f"EO forward ({sd.jv_V[j_fwd]:.2f} V, 0 mA)",
         dict(v_pn=sd.jv_V[j_fwd], i_heat=0.0), i0, j_fwd),
    ]

    fig, axes = plt.subplots(1, len(cases), figsize=(5 * len(cases), 4),
                             sharey=False)
    for ax, (title, kwargs, i, j) in zip(axes, cases):
        meas_T = sd.T_dB[i, j]
        sim_T  = wavelength_sweep(model, best, kl, wls, **kwargs)
        ax.plot(wls, _normalise_spectrum(sim_T),  lw=2, label="sim")
        ax.plot(wls, _normalise_spectrum(meas_T), "--", lw=1.5, label="meas")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Transmission (dB, norm.)")
        ax.legend(fontsize=8)

    fig.suptitle(f"Ring fit — {model}", y=1.01)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved: {out_path}")
    plt.close(fig)


def plot_2d_sweep(
    model: str,
    best: dict[str, float],
    sd: SweepData,
    out_path: str = "ring_sweep.png",
):
    """
    2-D colour map of resonance shift vs heater current and junction voltage.
    Left panel: measured; right panel: simulated.
    """
    if not _HAS_MPL:
        return

    kl  = best.get("kappa_l", _COUPLER[0].value)
    wls = sd.wl_nm

    i0 = int(np.argmin(np.abs(sd.hc_mA)))
    j0 = int(np.argmin(np.abs(sd.jv_V)))
    res0_meas = _find_resonance(wls, sd.T_dB[i0, j0])

    n_hc, n_jv = len(sd.hc_mA), len(sd.jv_V)
    delta_meas = np.empty((n_hc, n_jv))
    delta_sim  = np.empty((n_hc, n_jv))

    for i, hc in enumerate(sd.hc_mA):
        for j, jv in enumerate(sd.jv_V):
            delta_meas[i, j] = _find_resonance(wls, sd.T_dB[i, j]) - res0_meas
            sim_T = wavelength_sweep(model, best, kl, wls,
                                     v_pn=jv, i_heat=hc * 1e-3)
            delta_sim[i, j]  = _find_resonance(wls, sim_T) - res0_meas

    vmax = max(np.abs(delta_meas).max(), np.abs(delta_sim).max())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, data, title in [
        (ax1, delta_meas, "Measured"),
        (ax2, delta_sim,  "Simulated"),
    ]:
        im = ax.pcolormesh(sd.jv_V, sd.hc_mA, data,
                           cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
        ax.set_xlabel("Junction Voltage (V)")
        ax.set_ylabel("Heater Current (mA)")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label="Δλ_res (nm)")
    fig.suptitle(f"Resonance shift map — {model}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved: {out_path}")
    plt.close(fig)


#%% ── quick diagnostics ──────────────────────────────────────────────────────

def quick_sim(
    model: str = "fc_pn_th_ps",
    v_pn: float = 0.0,
    i_heat_mA: float = 0.0,
) -> None:
    """
    Simulate and plot a single spectrum with default parameters.
    Useful for sanity-checking the topology before fitting.
    """
    ps_p = {s.name: s.value for s in MODELS[model]["ps_params"]}
    kl   = _COUPLER[0].value
    wls  = np.linspace(WL_LO, WL_HI, N_SIM_PTS)

    print(f"Quick sim: model={model}  v_pn={v_pn} V  i_heat={i_heat_mA} mA")
    t_dB = wavelength_sweep(model, ps_p, kl, wls,
                            v_pn=v_pn, i_heat=i_heat_mA * 1e-3)
    res = _find_resonance(wls, t_dB)
    print(f"  Resonance: {res:.4f} nm")
    print(f"  Min T_dB: {t_dB.min():.2f} dB  at {wls[np.argmin(t_dB)]:.4f} nm")
    print(f"  T range: {t_dB.min():.2f} … {t_dB.max():.2f} dB")

    if _HAS_MPL:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(wls, t_dB)
        ax.axvline(res, color="r", ls="--", alpha=0.5, label=f"res={res:.3f} nm")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("T_dB (fibre-to-fibre)")
        ax.set_title(f"Quick sim — {model}")
        ax.legend()
        plt.tight_layout()
        plt.show()


#%% ── entry point ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default="fc_pn_th_ps",
                    choices=list(MODELS.keys()))
    ap.add_argument("--data", default=DATA_PATH,
                    help="path to NdSweeper pickle (without extension)")
    ap.add_argument("--quick-sim", action="store_true",
                    help="run a quick simulation with default params and exit")
    ap.add_argument("--plot-result", metavar="PKL",
                    help="load a previously-saved result dict and plot it")
    ap.add_argument("--maxiter", type=int, default=200)
    ap.add_argument("--popsize", type=int, default=12)
    args = ap.parse_args()

    if args.quick_sim:
        quick_sim(args.model)
    else:
        best = fit_staged(model=args.model, path=args.data, verbose=True)

        sweep = load_sweep(args.data)
        sd    = extract_data(sweep)
        plot_ring_fit(args.model, best, sd)
        plot_2d_sweep(args.model, best, sd)
