#!/usr/bin/env python3
"""
ringfit.py — add-drop ring machinery + staged CW fit (giona chip).

Shared library for every fit in this directory (dataset -> observables,
netlist build, ring wavelength sweep, staged fitter, plots) AND the CLI
driver for the May sparse-sweep staged fit.

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

Supported models (all have PN junction + heater):
  fc_pn_th_ps       — linear EO (L1)
  fc_pn_th_ps_cap   — depletion-mode with C_j(V) + da/dV (L2)
  fc_pn_th_ps_full  — piecewise PN: depletion + injection + TPA/self-heat (L3)

Staged fitting strategy:
  Stage 1 — passive      : fit n_eff, alpha_dB_cm, kappa_l from 0 V / 0 mA spectrum
  Stage 2 — thermal      : fit p_pi_th from heater current sweeps
  Stage 3 — EO           : fit v_pi_l (L1/L2) or dn_dv_rev + da_dv_rev (L3) from JV≤0
  Stage 4 — injection    : fit dn_dv_inj + da_dv_inj (L3 only) from 0 < JV ≤ 0.9 V
  Pre-fit                : r_heater from heater I/V; g_pn or (i_sat, n_diode) from I/V

V=0 self-consistency note (L3): the model guarantees continuity at V=0 by construction —
  both phi_eo_rev and phi_eo_inj are identically zero at v_pn=0, so separate depletion
  and injection fits always "meet at V=0" without any explicit constraint.

Setup
-----
  cd /path/to/fairchild
  maturin develop --release -m crates/fairchild-py/Cargo.toml
  source .venv/bin/activate
"""

#%% ── imports ────────────────────────────────────────────────────────────────

import json
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

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"        # raw lightlab captures (gitignored, large)
RESULTS = HERE / "results"     # fitted params (.json) + plots (.png)

DATA_PATH = str(DATA_DIR / "giona_neuron2_mod_joint_IV_spec")

# Restrict to the resonance window of our target ring (nm).
WL_LO, WL_HI = 1545.8, 1547.5
# Number of wavelength points for simulation (downsampled from measured).
N_SIM_PTS = 250
# Discard forward-bias data beyond this voltage (dataset-specific bad-data limit).
JV_HI = 1.0
# Parallel shunt resistor present on the physical device under test (Ω).
R_SHUNT_OHM = 2000.0


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

    # Discard bad data beyond JV_HI (dataset-specific cutoff).
    jv_mask  = jv_V <= JV_HI
    jv_V     = jv_V[jv_mask]
    jv_cols  = np.where(jv_mask)[0]     # original column indices kept
    v_heat_V  = v_heat_V[:, jv_mask]
    i_junc_mA = i_junc_mA[:, jv_mask]

    n_hc, n_jv = len(hc_mA), len(jv_V)
    # Use the first valid spectrum to set the shared wavelength axis.
    wl0, _ = _trim_spectrum(spectra[0, jv_cols[0]])
    wl_sim, _ = _downsample(wl0, wl0)   # just the axis

    T_dB = np.empty((n_hc, n_jv, N_SIM_PTS), dtype=float)
    for i in range(n_hc):
        for j, orig_j in enumerate(jv_cols):
            wl, T = _trim_spectrum(spectra[i, orig_j])
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

# ── Full-model (fc_pn_th_ps_full) parameter groups ──────────────────────────

_PN_REV_FULL = [
    # phi_eo_rev = 2π·L·dn_dv_rev·v_pn/λ.  With swapped terminal wiring v_pn = JV,
    # so for a typical Si PN ring (reverse bias → red shift) dn_dv_rev < 0.
    # Allow both signs so DE can find the correct direction.
    ParamSpec("dn_dv_rev", -5.024e-5, fixed=False, bounds=(-5e-3, 5e-3),
              description="reverse-bias Δn_eff/dV (depletion, 1/V); negative → red shift"),
    ParamSpec("da_dv_rev", 7.83,      fixed=False, bounds=(0.0, 200.0),
              description="reverse-bias FCA loss slope (Np/m/V)"),
    # Junction cap: unfittable without dynamic (AC) data — held fixed
    ParamSpec("c_j0",      1.375e-13, fixed=True, bounds=(1e-15, 1e-11),
              description="zero-bias junction capacitance (F) [fixed: no AC data]"),
    ParamSpec("v_bi",      0.917,     fixed=True, bounds=(0.3, 1.4),
              description="built-in voltage (V) [fixed: no AC data]"),
    ParamSpec("m_j",       0.5,       fixed=True, bounds=(0.2, 0.95),
              description="junction grading coefficient [fixed: no AC data]"),
]

_PN_FWD_FULL = [
    # i_sat and n_diode pre-fit from log-linear region of forward I/V
    ParamSpec("i_sat",      1e-12,    fixed=True,  bounds=(1e-20, 1e-6),
              description="Shockley saturation current (A) [pre-fit from I/V]"),
    ParamSpec("n_diode",    1.05,     fixed=True,  bounds=(0.5, 5.0),
              description="diode ideality factor [pre-fit from I/V]"),
    # Carrier lifetime: only observable from AC / transient data
    ParamSpec("tau_carrier", 10e-9,   fixed=True,  bounds=(0.1e-9, 100e-9),
              description="carrier lifetime (s) [fixed: AC-only observable]"),
    ParamSpec("dn_dv_inj",  1.311e-4, fixed=False, bounds=(1e-7, 1e-2),
              description="injection Δn_eff per (exp(V/Vt)−1) [dimensionless coeff]"),
    ParamSpec("da_dv_inj",  150.0,    fixed=False, bounds=(0.0, 2000.0),
              description="injection FCA loss prefactor (Np/m per carrier unit)"),
    # Series resistance pre-fit from diode I/V; held fixed during EO stages
    ParamSpec("r_series",   0.0,      fixed=True,  bounds=(0.0, 5000.0),
              description="series resistance (Ω) [pre-fit from diode I/V]"),
]

_SELF_HEATING = [
    # TPA and thermal: negligible at ~−10 dBm optical power at the ring.
    # r_th can be fit from optical-power-dependent sweeps, but not this dataset.
    ParamSpec("beta_tpa",  7.9e-12,   fixed=True, bounds=(1e-13, 1e-10),
              description="TPA coefficient (m/W) [fixed: Si literature, low power]"),
    ParamSpec("a_eff_m2",  1.257e-13, fixed=True, bounds=(1e-14, 1e-12),
              description="effective mode area (m²) [fixed: WG simulation]"),
    ParamSpec("r_th",      0.0,       fixed=True, bounds=(0.0, 1e5),
              description="thermal resistance (K/W) [fixed=0: constant optical power]"),
    ParamSpec("dn_dt",     1.86e-4,   fixed=True, bounds=(1e-5, 1e-3),
              description="thermo-optic dn/dT (1/K) [fixed: crystalline Si]"),
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
    "fc_pn_th_ps_full": {
        "description": "PN + thermal — full piecewise model: depletion + injection + TPA (L3)",
        "ps_params":   _WAVEGUIDE + _PN_REV_FULL + _PN_FWD_FULL + _SELF_HEATING + _HEATER,
        "has_pn":      True,
        "has_heater":  True,
        "has_cap":     True,
        "has_injection": True,
    },
}


#%% ── netlist builder ─────────────────────────────────────────────────────────

# fc_pn_th_ps / fc_pn_th_ps_cap terminal order: [in, out, anode, cathode, heat_p, heat_n]
#   anode=GND, cathode=PN_BIAS → v_pn_physics = PN_BIAS; positive PN_BIAS = reverse bias ✓
#
# fc_pn_th_ps_full terminal order: [in, out, anode, cathode, heat_p, heat_n]
#   SWAPPED: anode=PN_BIAS, cathode=GND → v_pn_physics = PN_BIAS (same dataset convention).
#   With this wiring: PN_BIAS < 0 → depletion branch activates; PN_BIAS > 0 → injection. ✓
_PS_LINE_TEMPLATES = {
    "fc_pn_th_ps":      "X{name} {in_} {out} GND PN_BIAS {heat_p} {heat_n} fc_pn_th_ps {params}",
    "fc_pn_th_ps_cap":  "X{name} {in_} {out} GND PN_BIAS {heat_p} {heat_n} fc_pn_th_ps_cap {params}",
    "fc_pn_th_ps_full": "X{name} {in_} {out} PN_BIAS GND {heat_p} {heat_n} fc_pn_th_ps_full {params}",
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
    # kappa_l is a coupler param (separate element), not a phase-shifter param.
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
        # COMMENTED OUT BECAUSE DATA DOESN'T FULLY DETERMINE THIS FIT MODEL # Physical 2 kΩ shunt resistor present on the measured device (pads to GND).
        # *(["R_SHUNT  PN_BIAS  GND  {:.0f}".format(R_SHUNT_OHM)]
        #   if model == "fc_pn_th_ps_full" else []),
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


#%% ── diode model helpers ────────────────────────────────────────────────────

def _diode_current_rs(
    v_arr: np.ndarray,
    i_sat: float,
    n: float,
    vt: float,
    r_s: float,
) -> np.ndarray:
    """
    Vectorized Newton-Raphson for the Shockley-with-series-resistance model.

    Solves V = V_j + R_s · I_d(V_j) for each element of v_arr, where
    I_d(V_j) = i_sat · (exp(V_j / (n·Vt)) − 1).

    For R_s ≤ 0 falls back to direct Shockley (no iteration needed).
    """
    clamp = 40.0
    if r_s <= 0.0:
        return i_sat * (np.exp(np.clip(v_arr / (n * vt), -clamp, clamp)) - 1.0)

    vj = np.array(v_arr, dtype=float)          # initial guess: V_j = V
    for _ in range(60):
        e   = np.exp(np.clip(vj / (n * vt), -clamp, clamp))
        i_d = i_sat * (e - 1.0)
        F   = vj + r_s * i_d - v_arr           # residual
        dF  = 1.0 + r_s * i_sat * e / (n * vt) # Jacobian
        delta = F / dF
        vj   -= delta
        if np.max(np.abs(delta)) < 1e-14:
            break
    e_final = np.exp(np.clip(vj / (n * vt), -clamp, clamp))
    return i_sat * (e_final - 1.0)


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


def prefit_diode_iv(sd: SweepData) -> tuple[float, float, float]:
    """
    Fit Shockley+series-resistance diode (i_sat, n_diode, r_series) from forward I/V.

    First bootstraps (i_sat, n_diode) from a log-linear fit in 0.3–0.6 V,
    then jointly fits all three parameters over the full 0.1 V–JV_HI range
    in log-current space using L-BFGS-B.

    Returns (i_sat [A], n_diode, r_series [Ω]).
    """
    Vt = 0.025852   # kT/q at 300 K
    i_meas_A = sd.i_junc_mA.mean(axis=0) * 1e-3    # average over HC axis, A
    # Subtract the physical 2 kΩ shunt so we fit the diode branch alone.
    i_junc_A = i_meas_A - sd.jv_V / R_SHUNT_OHM

    # Full fitting range: forward bias with positive diode current.
    fit_mask = (sd.jv_V >= 0.1) & (i_junc_A > 1e-10)
    if fit_mask.sum() < 3:
        print("  Pre-fit diode I/V: insufficient points — using defaults.")
        return 1e-12, 1.05, 0.0

    v_fit = sd.jv_V[fit_mask]
    i_fit = i_junc_A[fit_mask]

    # Bootstrap: log-linear fit in 0.3–0.6 V (exponential regime, R_s negligible).
    # After shunt subtraction, the 0.3–0.6 V region may have near-zero diode residual
    # (shunt current ~0.15–0.3 mA dominates), so fall back to safe defaults if needed.
    narrow = (sd.jv_V >= 0.3) & (sd.jv_V <= 0.6) & (i_junc_A > 0)
    if narrow.sum() >= 2:
        slope0, ic0 = np.polyfit(sd.jv_V[narrow], np.log(i_junc_A[narrow]), 1)
        i_sat0 = float(np.exp(ic0))
        n0     = float(np.clip(1.0 / (slope0 * Vt), 0.5, 5.0))
        if not (1e-20 < i_sat0 < 1e-3) or not (0.5 < n0 < 5.0):
            i_sat0, n0 = 1e-12, 1.05
    else:
        i_sat0, n0 = 1e-12, 1.05

    # Bootstrap R_s: at the highest V point, V_j_ideal = n·Vt·log(I/i_sat+1),
    # so R_s ≈ (V - V_j_ideal) / I.
    v_hi, i_hi = v_fit[-1], i_fit[-1]
    v_j_ideal  = n0 * Vt * np.log(max(i_hi / i_sat0 + 1.0, 1.0))
    r_s0       = max((v_hi - v_j_ideal) / i_hi, 1.0) if i_hi > 0 else 1.0

    # Joint L-BFGS-B fit in log-current space: x = [log(i_sat), n, r_s]
    def objective(x):
        i_s = float(np.exp(np.clip(x[0], -50, 0)))
        n_d = float(x[1])
        r_s = float(x[2])
        i_pred = _diode_current_rs(v_fit, i_s, n_d, Vt, r_s)
        log_pred = np.log(np.maximum(i_pred, 1e-30))
        return float(np.mean((log_pred - np.log(i_fit)) ** 2))

    x0 = [np.log(max(i_sat0, 1e-30)), n0, r_s0]
    res = minimize(
        objective, x0, method="L-BFGS-B",
        bounds=[(-45.0, -5.0), (0.5, 5.0), (0.0, 5000.0)],
        options={"ftol": 1e-12, "gtol": 1e-9, "maxiter": 2000},
    )
    i_sat   = float(np.exp(res.x[0]))
    n_diode = float(np.clip(res.x[1], 0.5, 5.0))
    r_series = float(max(res.x[2], 0.0))
    print(f"  Pre-fit diode: i_sat={i_sat:.3e} A, n_diode={n_diode:.3f}, "
          f"r_series={r_series:.1f} Ω  (loss={res.fun:.4e})")
    return i_sat, n_diode, r_series


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
    freeze = {
        "v_pi_l", "g_pn", "c_j0", "v_bi", "m_j", "da_dv",
        "r_heater", "p_pi_th", "tau_th",
        # Full-model EO/injection/self-heating params (frozen in Stage 1)
        "dn_dv_rev", "da_dv_rev", "dn_dv_inj", "da_dv_inj",
        "i_sat", "n_diode", "tau_carrier",
        "beta_tpa", "a_eff_m2", "r_th", "dn_dt",
    }
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


def stage3_eo_full(
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
    Stage 3 (full model): fit dn_dv_rev and da_dv_rev from reverse-bias spectra (JV ≤ 0).

    V=0 self-consistency is guaranteed by model construction: both phi_eo_rev and
    phi_eo_inj are exactly zero at v_pn=0, so depletion and injection fits always
    meet at V=0 without any explicit constraint.  The slope discontinuity at V=0 is
    physically intentional (depletion ↔ injection are different mechanisms).
    """
    if base_specs is None:
        base_specs = deepcopy(MODELS[model]["ps_params"])
    if coupler_specs is None:
        coupler_specs = deepcopy(_COUPLER)

    for s in base_specs + coupler_specs:
        s.fixed = True
        if s.name in stage2_result:
            s.value = stage2_result[s.name]

    for s in base_specs:
        if s.name in {"dn_dv_rev", "da_dv_rev"}:
            s.fixed = False

    jv_rev = sd.jv_V[sd.jv_V <= 0]
    if len(jv_rev) == 0:
        print("  Stage 3 (depletion): no reverse-bias JV points — skipping.")
        return stage2_result
    n_pts  = min(n_jv_points, len(jv_rev))
    jv_sel = jv_rev[np.round(np.linspace(0, len(jv_rev) - 1, n_pts)).astype(int)]

    i0  = int(np.argmin(np.abs(sd.hc_mA)))
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
                      label="Stage 3 (EO depletion, full)", maxiter=maxiter,
                      popsize=popsize, verbose=verbose)}
    return best


def stage4_eo_injection(
    model: str,
    sd: SweepData,
    stage3_result: dict[str, float],
    base_specs: Optional[list[ParamSpec]] = None,
    coupler_specs: Optional[list[ParamSpec]] = None,
    n_jv_points: int = 8,
    maxiter: int = 200,
    popsize: int = 12,
    verbose: bool = True,
) -> dict[str, float]:
    """
    Stage 4 (full model): fit dn_dv_inj and da_dv_inj from forward-bias spectra.

    Restricted to 0 < JV ≤ 0.9 V — the Shockley exponential clamps at 40·Vt ≈ 1.04 V,
    so fitting beyond ~0.9 V is unreliable.
    """
    if base_specs is None:
        base_specs = deepcopy(MODELS[model]["ps_params"])
    if coupler_specs is None:
        coupler_specs = deepcopy(_COUPLER)

    for s in base_specs + coupler_specs:
        s.fixed = True
        if s.name in stage3_result:
            s.value = stage3_result[s.name]

    for s in base_specs:
        if s.name in {"dn_dv_inj", "da_dv_inj"}:
            s.fixed = False

    # Restrict to 0 < JV ≤ 0.9 V (below Shockley clamp at 40·Vt ≈ 1.04 V)
    jv_fwd = sd.jv_V[(sd.jv_V > 0) & (sd.jv_V <= 0.9)]
    if len(jv_fwd) == 0:
        print("  Stage 4 (injection): no forward-bias JV points (0–0.9 V) — skipping.")
        return stage3_result
    n_pts  = min(n_jv_points, len(jv_fwd))
    jv_sel = jv_fwd[np.round(np.linspace(0, len(jv_fwd) - 1, n_pts)).astype(int)]

    i0  = int(np.argmin(np.abs(sd.hc_mA)))
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
    best = {**stage3_result,
            **_run_de(objective, combined, x0,
                      label="Stage 4 (EO injection, full)", maxiter=maxiter,
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

    is_full = (model == "fc_pn_th_ps_full")

    # ── pre-fits from electrical data ──────────────────────────────────────
    print("\nPre-fitting electrical parameters …")
    r_heater = prefit_r_heater(sd)

    if is_full:
        i_sat, n_diode, r_series = prefit_diode_iv(sd)
    else:
        g_pn = prefit_g_pn(sd)

    # Apply pre-fit values to all relevant ParamSpec lists
    def _apply_prefit(specs):
        for s in specs:
            if s.name == "r_heater":
                s.value = r_heater
                s.fixed = True
            elif not is_full and s.name == "g_pn":
                s.value = g_pn
                s.fixed = True
            elif is_full and s.name == "i_sat":
                s.value = i_sat
                s.fixed = True
            elif is_full and s.name == "n_diode":
                s.value = n_diode
                s.fixed = True
            elif is_full and s.name == "r_series":
                s.value = r_series
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
    if is_full:
        print("\nStage 3 — EO depletion (dn_dv_rev, da_dv_rev) …")
        s3 = stage3_eo_full(model, sd, s2, base_specs=base3, coupler_specs=cpl3,
                            verbose=verbose)
    else:
        print("\nStage 3 — EO (v_pi_l) …")
        s3 = stage3_eo(model, sd, s2, base_specs=base3, coupler_specs=cpl3,
                       verbose=verbose)

    # ── Stage 4: injection EO (full model only) ───────────────────────────
    if is_full:
        base4 = _apply_prefit(deepcopy(MODELS[model]["ps_params"]))
        cpl4  = deepcopy(_COUPLER)
        print("\nStage 4 — EO injection (dn_dv_inj, da_dv_inj) …")
        s4 = stage4_eo_injection(model, sd, s3, base_specs=base4, coupler_specs=cpl4,
                                 verbose=verbose)
    else:
        s4 = s3

    print_results(s4)
    return s4


#%% ── result save / load ─────────────────────────────────────────────────────

def save_params(best: dict[str, float], model: str, path: str) -> None:
    """Save best-fit params and model name to a JSON file."""
    with open(path, "w") as f:
        json.dump({"model": model, "params": best}, f, indent=2)
    print(f"Parameters saved: {path}")


def load_params(path: str) -> tuple[str, dict[str, float]]:
    """Load model name and params dict from a JSON file saved by save_params()."""
    with open(path) as f:
        data = json.load(f)
    return data["model"], {k: float(v) for k, v in data["params"].items()}


#%% ── result reporting & plotting ────────────────────────────────────────────

def print_results(best: dict[str, float]):
    print("\n── Best-fit parameters ──────────────────────────────────")
    for k, v in best.items():
        print(f"  {k:<18s} = {v:>14.6g}")
    print()


def _sel_indices(arr: np.ndarray, n: int, lo: float = -np.inf, hi: float = np.inf) -> np.ndarray:
    """Return n evenly-spaced indices into arr restricted to [lo, hi]."""
    mask = (arr >= lo) & (arr <= hi)
    full_idx = np.where(mask)[0]
    if len(full_idx) == 0:
        return np.array([], dtype=int)
    chosen = np.round(np.linspace(0, len(full_idx) - 1, min(n, len(full_idx)))).astype(int)
    return full_idx[chosen]


def _style_legend(ax, loc="lower right"):
    """Add a solid/dashed sim/data legend entry."""
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color="k", lw=1.8, label="sim"),
               Line2D([0], [0], color="k", lw=1.2, ls="--", label="data")]
    ax.legend(handles=handles, fontsize=8, loc=loc)


def plot_ring_fit(
    model: str,
    best: dict[str, float],
    sd: SweepData,
    out_path: str = "ring_fit.png",
):
    """Four-panel overview: passive, thermal, EO reverse, EO forward."""
    if not _HAS_MPL:
        print("matplotlib not available — skipping plot.")
        return

    kl  = best.get("kappa_l", _COUPLER[0].value)
    wls = sd.wl_nm

    i0    = int(np.argmin(np.abs(sd.hc_mA)))
    j0    = int(np.argmin(np.abs(sd.jv_V)))
    i_max = int(np.argmax(sd.hc_mA))
    j_rev = int(np.argmin(sd.jv_V))
    j_fwd = int(np.argmax(sd.jv_V))

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

    fig, axes = plt.subplots(1, 4, figsize=(20, 4), sharey=False)
    for ax, (title, kwargs, i, j) in zip(axes, cases):
        meas_T = sd.T_dB[i, j]
        sim_T  = wavelength_sweep(model, best, kl, wls, **kwargs)
        ax.plot(wls, _normalise_spectrum(sim_T),  lw=1.8)
        ax.plot(wls, _normalise_spectrum(meas_T), ls="--", lw=1.2)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Transmission (dB, norm.)")
    _style_legend(axes[0])

    fig.suptitle(f"Ring fit — {model}")
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
    """2-D colour map of resonance shift vs heater current and junction voltage."""
    if not _HAS_MPL:
        return

    kl  = best.get("kappa_l", _COUPLER[0].value)
    wls = sd.wl_nm

    i0        = int(np.argmin(np.abs(sd.hc_mA)))
    j0        = int(np.argmin(np.abs(sd.jv_V)))
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


def plot_iv_curves(
    model: str,
    best: dict[str, float],
    sd: SweepData,
    out_path: str = "iv_curves.png",
):
    """
    Heater I-V and diode I-V: model vs data.

    Left panel  — heater: V = 2·R_arm·I (model) vs measured heater voltage.
    Right panel — diode: Shockley or linear conductance (model) vs measured junction current.
                  Lower inset uses log scale on the forward-bias half to show
                  Shockley fit quality.
    """
    if not _HAS_MPL:
        return

    # ── heater I-V ────────────────────────────────────────────────────────
    hc_A         = sd.hc_mA * 1e-3
    v_heat_meas  = sd.v_heat_V.mean(axis=1)      # average over JV axis
    r_heater     = best.get("r_heater", 184.1)
    v_heat_model = 2.0 * r_heater * hc_A

    # ── diode I-V ─────────────────────────────────────────────────────────
    jv_V          = sd.jv_V
    i_junc_mA_meas = sd.i_junc_mA.mean(axis=0)   # average over HC axis

    Vt = 0.025852   # kT/q at 300 K
    if model == "fc_pn_th_ps_full":
        n_d    = best.get("n_diode",  1.05)
        i_sat  = best.get("i_sat",    1e-12)
        r_s    = best.get("r_series", 0.0)
        # Total current = diode branch + physical 2 kΩ shunt branch.
        i_model_A = _diode_current_rs(jv_V, i_sat, n_d, Vt, r_s) + jv_V / R_SHUNT_OHM
    else:
        g_pn      = best.get("g_pn", 1e-3)
        i_model_A = g_pn * jv_V
    i_model_mA = i_model_A * 1e3

    fig = plt.figure(figsize=(12, 4))
    gs  = fig.add_gridspec(1, 2, wspace=0.35)

    # Left: heater I-V
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(sd.hc_mA, v_heat_model, lw=1.8, label=f"model (2·{r_heater:.0f} Ω)")
    ax1.plot(sd.hc_mA, v_heat_meas,  ls="--", lw=1.2, label="data")
    ax1.set_xlabel("Heater Current (mA)")
    ax1.set_ylabel("Heater Voltage (V)")
    ax1.set_title("Heater I-V")
    ax1.legend(fontsize=8)

    # Right: diode I-V — two stacked sub-axes sharing X
    gs_r = gs[1].subgridspec(2, 1, hspace=0.08, height_ratios=[1, 1])
    ax2t = fig.add_subplot(gs_r[0])   # linear, full range
    ax2b = fig.add_subplot(gs_r[1], sharex=ax2t)   # log, forward only

    for ax, yscale, fwd_only, ylabel in [
        (ax2t, "linear", False, "I_junc (mA)"),
        (ax2b, "log",    True,  "I_junc (mA, log)"),
    ]:
        if fwd_only:
            mask = jv_V > 0
            x    = jv_V[mask]
            ym   = i_model_mA[mask]
            yd   = i_junc_mA_meas[mask]
            # only keep positive values for log scale
            pos  = (ym > 0) & (yd > 0)
            ax.plot(x[pos], ym[pos], lw=1.8)
            ax.plot(x[pos], yd[pos], ls="--", lw=1.2)
            ax.set_yscale("log")
        else:
            ax.plot(jv_V, i_model_mA,     lw=1.8,         label="model")
            ax.plot(jv_V, i_junc_mA_meas, ls="--", lw=1.2, label="data")
            ax.legend(fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=8)

    ax2b.set_xlabel("Junction Voltage (V)")
    ax2t.set_title("Diode I-V")
    plt.setp(ax2t.get_xticklabels(), visible=False)

    fig.suptitle(f"I-V curves — {model}")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved: {out_path}")
    plt.close(fig)


def plot_spectra_vs_jv(
    model: str,
    best: dict[str, float],
    sd: SweepData,
    out_path: str = "spectra_vs_jv.png",
    n_curves: int = 6,
):
    """
    Spectra at two fixed heater currents (0 mA and 2 mA), sweeping junction voltage.

    n_curves uniformly-spaced JV points are shown per panel, coloured by JV using
    the coolwarm map (blue = most negative, red = most positive).
    Solid lines = simulation; dashed lines = data.
    """
    if not _HAS_MPL:
        return

    kl  = best.get("kappa_l", _COUPLER[0].value)
    wls = sd.wl_nm

    jv_idx = _sel_indices(sd.jv_V, n_curves)
    cmap   = plt.get_cmap("coolwarm")
    norm   = plt.Normalize(vmin=sd.jv_V.min(), vmax=sd.jv_V.max())

    # Fixed heater currents: 0 mA and 2 mA (closest available)
    hc_targets = [0.0, 2.0]
    hc_indices = [int(np.argmin(np.abs(sd.hc_mA - t))) for t in hc_targets]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, i_hc in zip(axes, hc_indices):
        hc_mA = sd.hc_mA[i_hc]
        for j_jv in jv_idx:
            jv    = sd.jv_V[j_jv]
            color = cmap(norm(jv))
            sim_T  = wavelength_sweep(model, best, kl, wls,
                                      v_pn=jv, i_heat=hc_mA * 1e-3)
            meas_T = sd.T_dB[i_hc, j_jv]
            ax.plot(wls, _normalise_spectrum(sim_T),  color=color, lw=1.8)
            ax.plot(wls, _normalise_spectrum(meas_T), color=color, lw=1.0, ls="--")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Transmission (dB, norm.)")
        ax.set_title(f"HC = {hc_mA:.1f} mA")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="Junction Voltage (V)")
    _style_legend(axes[0])

    fig.suptitle(f"Spectra vs junction voltage — {model}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved: {out_path}")
    plt.close(fig)


def plot_spectra_vs_hc(
    model: str,
    best: dict[str, float],
    sd: SweepData,
    out_path: str = "spectra_vs_hc.png",
    n_curves: int = 6,
):
    """
    Spectra at three fixed junction voltages, sweeping heater current.

    Three panels: max reverse bias, moderate forward (~1 V), max forward.
    n_curves uniformly-spaced HC points (≥ 0 mA) per panel, coloured by heater
    current using the plasma map (dark = low, bright = high).
    Solid lines = simulation; dashed lines = data.
    """
    if not _HAS_MPL:
        return

    kl  = best.get("kappa_l", _COUPLER[0].value)
    wls = sd.wl_nm

    hc_idx = _sel_indices(sd.hc_mA, n_curves, lo=0.0)
    hc_sel = sd.hc_mA[hc_idx]
    cmap   = plt.get_cmap("plasma")
    norm   = plt.Normalize(vmin=hc_sel.min(), vmax=hc_sel.max())

    # Three JV slices
    j_rev = int(np.argmin(sd.jv_V))
    j_mod = int(np.argmin(np.abs(sd.jv_V - 1.0)))
    j_fwd = int(np.argmax(sd.jv_V))
    jv_cases = [
        (j_rev, f"JV = {sd.jv_V[j_rev]:.2f} V  (max reverse)"),
        (j_mod, f"JV = {sd.jv_V[j_mod]:.2f} V  (≈+1 V)"),
        (j_fwd, f"JV = {sd.jv_V[j_fwd]:.2f} V  (max forward)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 4), sharey=True)
    for ax, (j_jv, jv_label) in zip(axes, jv_cases):
        jv = sd.jv_V[j_jv]
        for i_hc in hc_idx:
            hc    = sd.hc_mA[i_hc]
            color = cmap(norm(hc))
            sim_T  = wavelength_sweep(model, best, kl, wls,
                                      v_pn=jv, i_heat=hc * 1e-3)
            meas_T = sd.T_dB[i_hc, j_jv]
            ax.plot(wls, _normalise_spectrum(sim_T),  color=color, lw=1.8)
            ax.plot(wls, _normalise_spectrum(meas_T), color=color, lw=1.0, ls="--")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Transmission (dB, norm.)")
        ax.set_title(jv_label, fontsize=9)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="Heater Current (mA)")
    _style_legend(axes[0])

    fig.suptitle(f"Spectra vs heater current — {model}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved: {out_path}")
    plt.close(fig)


def explore_diode_models(
    sd: SweepData,
    best: Optional[dict[str, float]] = None,
    out_path: str = "diode_models.png",
) -> None:
    """
    Compare diode I-V models against measured forward-bias data.

    Top panel  — log I-V: data points, pure Shockley (narrow fit), Shockley+Rs.
    Bottom panel — log10 residuals (model / data) with RMS in legend.

    If best is provided and contains i_sat/n_diode/r_series, those values are used
    for the Shockley+Rs model instead of re-running the fit.
    """
    if not _HAS_MPL:
        return

    Vt = 0.025852
    i_meas_A = sd.i_junc_mA.mean(axis=0) * 1e-3
    # Subtract physical 2 kΩ shunt so we compare diode-only model to diode-only data.
    i_diode_A = i_meas_A - sd.jv_V / R_SHUNT_OHM

    # Model 1: pure Shockley bootstrapped from 0.3–0.6 V of shunt-corrected current.
    narrow = (sd.jv_V >= 0.3) & (sd.jv_V <= 0.6) & (i_diode_A > 0)
    if narrow.sum() >= 2:
        slope0, ic0 = np.polyfit(sd.jv_V[narrow], np.log(i_diode_A[narrow]), 1)
        i_sat0 = float(np.exp(ic0))
        n0     = float(np.clip(1.0 / (slope0 * Vt), 0.5, 5.0))
        if not (1e-20 < i_sat0 < 1e-3) or not (0.5 < n0 < 5.0):
            i_sat0, n0 = 1e-12, 1.05
    else:
        i_sat0, n0 = 1e-12, 1.05
    i_shockley = i_sat0 * (np.exp(np.clip(sd.jv_V / (n0 * Vt), -40, 40)) - 1.0)

    # Model 2: Shockley+Rs — use best dict if available to avoid re-fitting.
    if best is not None and "i_sat" in best and "r_series" in best:
        i_sat_rs = best["i_sat"]
        n_rs     = best.get("n_diode", 1.05)
        r_s      = best["r_series"]
    else:
        i_sat_rs, n_rs, r_s = prefit_diode_iv(sd)
    i_rs = _diode_current_rs(sd.jv_V, i_sat_rs, n_rs, Vt, r_s)

    # Restrict comparison to forward bias with positive shunt-corrected data.
    fwd = (sd.jv_V >= 0.1) & (i_diode_A > 1e-10)
    v_fwd  = sd.jv_V[fwd]
    id_fwd = i_diode_A[fwd]
    is_fwd = np.maximum(i_shockley[fwd], 1e-30)
    ir_fwd = np.maximum(i_rs[fwd],       1e-30)

    rms_s  = float(np.sqrt(np.mean((np.log10(is_fwd / np.maximum(id_fwd, 1e-30))) ** 2)))
    rms_rs = float(np.sqrt(np.mean((np.log10(ir_fwd / np.maximum(id_fwd, 1e-30))) ** 2)))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                   gridspec_kw={"hspace": 0.06})

    # Top: log I-V
    ax1.semilogy(v_fwd, id_fwd * 1e3, "ko", ms=3, zorder=5, label="data (shunt-corrected)")
    ax1.semilogy(v_fwd, is_fwd * 1e3, "b-", lw=1.8,
                 label=f"Shockley  (I_sat={i_sat0:.2e} A, n={n0:.3f})")
    ax1.semilogy(v_fwd, ir_fwd * 1e3, "r-", lw=1.8,
                 label=f"Shockley+Rs  (Rs={r_s:.0f} Ω, n={n_rs:.3f})")
    ax1.set_ylabel("I_diode (mA, shunt-corrected)")
    ax1.legend(fontsize=8)
    ax1.set_title("Diode model comparison (shunt-corrected data)")
    ax1.grid(True, which="both", alpha=0.3)

    # Bottom: log-residuals
    res_s  = np.log10(is_fwd / np.maximum(id_fwd, 1e-30))
    res_rs = np.log10(ir_fwd / np.maximum(id_fwd, 1e-30))
    ax2.plot(v_fwd, res_s,  "b-", lw=1.5, label=f"Shockley   RMS={rms_s:.3f} dec")
    ax2.plot(v_fwd, res_rs, "r-", lw=1.5, label=f"Shockley+Rs  RMS={rms_rs:.3f} dec")
    ax2.axhline(0, color="k", lw=0.6)
    ax2.set_ylabel("log10(I_model / I_data)")
    ax2.set_xlabel("Junction Voltage (V)")
    ax2.legend(fontsize=8)
    ax2.grid(True, which="both", alpha=0.3)

    fig.suptitle("Diode I-V model exploration")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved: {out_path}")
    plt.close(fig)


def plot_all(model: str, best: dict[str, float], sd: SweepData, prefix: str = "") -> None:
    """Run all plot functions with a shared filename prefix."""
    RESULTS.mkdir(exist_ok=True)
    p = str(RESULTS / (prefix + "_")) if prefix else str(RESULTS) + "/"
    plot_ring_fit(model, best, sd,         out_path=f"{p}ring_fit.png")
    plot_2d_sweep(model, best, sd,         out_path=f"{p}ring_sweep.png")
    plot_iv_curves(model, best, sd,        out_path=f"{p}iv_curves.png")
    plot_spectra_vs_jv(model, best, sd,    out_path=f"{p}spectra_vs_jv.png")
    plot_spectra_vs_hc(model, best, sd,    out_path=f"{p}spectra_vs_hc.png")
    if model == "fc_pn_th_ps_full":
        explore_diode_models(sd, best=best, out_path=f"{p}diode_models.png")


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
                    choices=list(MODELS.keys()),
                    help="phase-shifter model (default: fc_pn_th_ps)")
    ap.add_argument("--data", default=DATA_PATH,
                    help="path to NdSweeper pickle (without extension)")
    ap.add_argument("--quick-sim", action="store_true",
                    help="run a quick simulation with default params and exit")
    ap.add_argument("--load-params", metavar="JSON",
                    help="skip fitting; load params from JSON and plot only")
    ap.add_argument("--save-params", metavar="JSON", default=None,
                    help="path to save best-fit params (default: <model>_fit.json)")
    ap.add_argument("--plot-prefix", default="",
                    help="filename prefix for all output plots")
    ap.add_argument("--maxiter", type=int, default=200)
    ap.add_argument("--popsize", type=int, default=12)
    args = ap.parse_args()

    if args.quick_sim:
        quick_sim(args.model)
        sys.exit(0)

    sweep = load_sweep(args.data)
    sd    = extract_data(sweep)

    if args.load_params:
        model, best = load_params(args.load_params)
        print(f"Loaded params from {args.load_params}  (model={model})")
        print_results(best)
    else:
        best  = fit_staged(model=args.model, path=args.data, verbose=True)
        model = args.model
        json_path = args.save_params or str(RESULTS / f"{model}_fit.json")
        save_params(best, model, json_path)

    plot_all(model, best, sd, prefix=args.plot_prefix or model)
