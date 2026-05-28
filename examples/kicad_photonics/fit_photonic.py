#!/usr/bin/env python3
"""
fit_photonic.py — photonic device parameter fitting for fairchild.

Fits PN and/or thermal phase-shifter parameters (and coupler coupling
coefficients) to experimental transmission-vs-wavelength data and
electro-optic / thermal tuning curves.

Topology
--------
MZI modulator with two directional couplers and two phase-shifter arms,
based on the Fairchild-format single_modulator_tb.sp netlist.

    IN ──► WG1 ──► CPL1 ──┬──► PS_arm1 ──┬──► CPL2 ──► WG3 ──► THRU
                           └──► PS_arm2 ──┘

Experimental data formats (CSV, first row is header)
-----------------------------------------------------
  spectrum.csv       :  wavelength_nm,  transmission_dB
  eo_tuning.csv      :  v_pn_V,         resonance_nm
  thermal_tuning.csv :  i_heat_A,        resonance_nm

Usage
-----
  # Smoke-test with synthetic data (no real files needed)
  python fit_photonic.py

  # Fit to your measured spectrum
  python fit_photonic.py --spectrum spectrum.csv

  # Full fit: spectrum + EO + thermal tuning, specific model
  python fit_photonic.py --model fc_pn_th_ps_cap \\
      --spectrum spectrum.csv --eo eo_tuning.csv --thermal thermal_tuning.csv

  # List available models and their parameters
  python fit_photonic.py --list-models

Setup
-----
  cd /path/to/fairchild
  maturin develop --release -m crates/fairchild-py/Cargo.toml
"""

import argparse
import re
import sys
from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.optimize import differential_evolution, minimize

try:
    import matplotlib
    matplotlib.use("Agg")
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



# ── parameter specification ───────────────────────────────────────────────────

@dataclass
class ParamSpec:
    """One tunable or fixed simulation parameter."""
    name: str                       # fairchild SPICE parameter name
    value: float                    # initial / fixed value
    fixed: bool = False             # True → held constant during optimisation
    bounds: tuple = (0.0, 1.0)
    description: str = ""

    def copy(self, **overrides):
        p = deepcopy(self)
        for k, v in overrides.items():
            setattr(p, k, v)
        return p


# ── per-model parameter tables ────────────────────────────────────────────────
# Edit `fixed` and `bounds` here to control the optimisation space.

_WAVEGUIDE = [
    ParamSpec("n_g",         4.2,      fixed=True,  bounds=(3.5, 5.0),
              description="group index (dispersion slope; set from waveguide sim)"),
    ParamSpec("n_eff",       2.40,     fixed=False, bounds=(2.0, 3.5),
              description="effective index at λ_ref (sets fringe spacing / resonance position)"),
    ParamSpec("l_m",         25.13e-6, fixed=True,  bounds=(1e-6, 5e-3),
              description="modulator arm length (m) — usually known from layout"),
    ParamSpec("alpha_db_cm", 10.0,     fixed=False, bounds=(0.5, 60.0),
              description="propagation loss (dB/cm); sets extinction / insertion loss"),
    ParamSpec("pin_at_ref",  1.0,      fixed=True,  bounds=(0.0, 1.0),
              description="pin resonance to λ_ref (1=on-resonance at laser wl, 0=absolute phase)"),
]

_PN_EO = [
    ParamSpec("v_pi_l",      0.02,     fixed=False, bounds=(1e-3, 0.5),
              description="Vπ·L (V·m) — EO tuning efficiency"),
    ParamSpec("g_pn",        1e-3,     fixed=True,  bounds=(1e-9, 1e3),
              description="PN junction ohmic conductance (S)"),
]

_PN_CAP = [
    ParamSpec("c_j0",        20e-15,   fixed=False, bounds=(1e-16, 1e-12),
              description="zero-bias junction capacitance (F)"),
    ParamSpec("v_bi",        0.7,      fixed=False, bounds=(0.3, 1.4),
              description="built-in voltage (V)"),
    ParamSpec("m_j",         0.5,      fixed=False, bounds=(0.2, 0.95),
              description="junction grading coefficient"),
    ParamSpec("da_dv",       0.0,      fixed=True,  bounds=(0.0, 50.0),
              description="bias-dependent loss slope (Np/m/V)"),
]

_HEATER = [
    ParamSpec("r_heater",    200.0,    fixed=True,  bounds=(50.0, 5000.0),
              description="heater resistance (Ω) — usually measured directly"),
    ParamSpec("p_pi_th",     5e-3,     fixed=False, bounds=(1e-4, 0.2),
              description="heater π-power (W); thermal tuning efficiency"),
]

_HEATER_TAU = [
    ParamSpec("tau_th",      10e-6,    fixed=False, bounds=(1e-7, 1e-3),
              description="thermal time constant (s)"),
]

_COUPLER = [
    ParamSpec("kappa_l",     0.0769,   fixed=False, bounds=(0.01, 0.99),
              description="coupler κ·L [same value applied to both CPL1 and CPL2]"),
]

# Aggregate model definitions.
MODELS: dict[str, dict] = {
    "fc_pn_ps": {
        "description": "PN phase shifter — linear EO, no junction cap",
        "ps_params":   _WAVEGUIDE + _PN_EO,
        "has_pn":      True,
        "has_heater":  False,
        "has_cap":     False,
    },
    "fc_pn_ps_cap": {
        "description": "PN phase shifter — depletion-mode with bias-dependent C_j(V)",
        "ps_params":   _WAVEGUIDE + _PN_EO + _PN_CAP,
        "has_pn":      True,
        "has_heater":  False,
        "has_cap":     True,
    },
    "fc_thermal_ps": {
        "description": "Thermal phase shifter — instantaneous Joule heating",
        "ps_params":   _WAVEGUIDE + _HEATER,
        "has_pn":      False,
        "has_heater":  True,
        "has_cap":     False,
    },
    "fc_thermal_ps_rc": {
        "description": "Thermal phase shifter with first-order RC time constant",
        "ps_params":   _WAVEGUIDE + _HEATER + _HEATER_TAU,
        "has_pn":      False,
        "has_heater":  True,
        "has_cap":     False,
    },
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


# ── netlist builder ───────────────────────────────────────────────────────────

# Terminal layout for each model's PS SPICE line (Fairchild subcircuit format).
# Optical bundle ports each expand to 3 wires internally; electrical ports follow.
# anode=GND, cathode=PN_BIAS → V_pn = -PN_BIAS (reverse-bias convention).
# {heat_p}/{heat_n} are filled per-arm for series heater wiring.
_PS_LINE_TEMPLATES = {
    "fc_pn_ps":         "X{name} {in_} {out} GND PN_BIAS fc_pn_ps {params}",
    "fc_pn_ps_cap":     "X{name} {in_} {out} GND PN_BIAS fc_pn_ps_cap {params}",
    "fc_thermal_ps":    "X{name} {in_} {out} {heat_p} {heat_n} fc_thermal_ps {params}",
    "fc_thermal_ps_rc": "X{name} {in_} {out} {heat_p} {heat_n} fc_thermal_ps_rc {params}",
    "fc_pn_th_ps":      "X{name} {in_} {out} GND PN_BIAS {heat_p} {heat_n} fc_pn_th_ps {params}",
    "fc_pn_th_ps_cap":  "X{name} {in_} {out} GND PN_BIAS {heat_p} {heat_n} fc_pn_th_ps_cap {params}",
}

# MZI arm optical connections (cross-coupled topology).
_ARM_CONNECTIONS = [
    ("CPL2_b1", "CPL1_a2"),  # arm 1
    ("CPL1_b2", "CPL2_a1"),  # arm 2
]

# Series heater wiring: arm1 top → shared midpoint → arm2 bottom → GND.
# A single current source drives current through both heaters in series.
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
    Build a complete SPICE netlist string for the MZI modulator.

    Parameters
    ----------
    model          : one of the keys in MODELS
    ps_params      : dict of fairchild parameter name → value (both arms share these)
    coupler_kappa_l: κ·L for both directional couplers
    wavelength_nm  : CW laser wavelength (nm)
    v_pn           : PN bias voltage (V)
    i_heat         : heater drive current (A); flows through both heaters in series
    """
    # Remove coupler param from ps_params if it leaked in
    ps_p = {k: v for k, v in ps_params.items() if k != "kappa_l"}
    ps_spice = _params_to_spice(ps_p)

    tmpl = _PS_LINE_TEMPLATES[model]
    ps_lines = []
    for i, ((in_, out), (heat_p, heat_n)) in enumerate(
        zip(_ARM_CONNECTIONS, _HEATER_CONNECTIONS), start=1
    ):
        ps_lines.append(tmpl.format(
            name=f"PS{i}", in_=in_, out=out, params=ps_spice,
            heat_p=heat_p, heat_n=heat_n,
        ))

    # All optical nets must be declared so the parser expands them to 3-wire bundles.
    optical_ports = [
        "IN_OPT", "THRU_OPT", "DROP_OPT", "ADD_OPT",
        "CWL_OUT", "GC_DROP_OUT", "GC_THRU_OUT",
        "CPL1_a1", "CPL1_a2", "CPL1_b1", "CPL1_b2",
        "CPL2_a1", "CPL2_a2", "CPL2_b1", "CPL2_b2",
    ]
    port_decls = [f".optical_port {p}" for p in optical_ports]

    return "\n".join([
        ".title MZI modulator — parameter fit",
        *port_decls,
        ".op",
        f"XCWL1 CWL_OUT fc_cw_laser power_mW=1.0 wavelength_nm={wavelength_nm:.6f}",
        "V1 VDD GND DC 2",
        f"V_PN  PN_BIAS  GND DC {v_pn:.6g}",
        f"IHEAT HEAT_BIAS GND DC {i_heat:.6g}",
        "R1 V_DROP GND 1k",
        "R2 V_THRU GND 1k",
        "XGC_IN  CWL_OUT IN_OPT fc_grating_coupler alpha_dB=6.5",
        "XGC_DROP DROP_OPT GC_DROP_OUT fc_grating_coupler alpha_dB=6.5",
        "XGC_THRU THRU_OPT GC_THRU_OUT fc_grating_coupler alpha_dB=6.5",
        "XPD1 GC_DROP_OUT VDD V_DROP fc_photodetector"
        " responsivity=0.7 i_dark_a=1e-9 r_shunt=20k",
        "XPD2 GC_THRU_OUT VDD V_THRU fc_photodetector"
        " responsivity=0.7 i_dark_a=1e-9 r_shunt=20k",
        "XWG1 IN_OPT CPL1_a1 fc_waveguide l_m=304.5e-6 n_g=4.2 alpha_dB_cm=1.0",
        "XWG2 CPL2_b2 DROP_OPT fc_waveguide l_m=304.5e-6 n_g=4.2 alpha_dB_cm=1.0",
        "XWG3 CPL1_b1 THRU_OPT fc_waveguide l_m=304.5e-6 n_g=4.2 alpha_dB_cm=1.0",
        "XWG4 ADD_OPT CPL2_a2 fc_waveguide l_m=304.5e-6 n_g=4.2 alpha_dB_cm=1.0",
        f"XCPL1 CPL1_a1 CPL1_a2 CPL1_b1 CPL1_b2"
        f" fc_dcoupler kappa_L={coupler_kappa_l:.8g}",
        f"XCPL2 CPL2_a1 CPL2_a2 CPL2_b1 CPL2_b2"
        f" fc_dcoupler kappa_L={coupler_kappa_l:.8g}",
        *ps_lines,
        ".end",
    ]) + "\n"


# ── simulation helpers ────────────────────────────────────────────────────────

def wavelength_sweep(
    model: str,
    ps_params: dict[str, float],
    coupler_kappa_l: float,
    wavelengths_nm: np.ndarray,
    v_pn: float = 0.0,
    i_heat: float = 0.0,
) -> np.ndarray:
    """
    Run a DC OP at each wavelength and return V_THRU (V).

    The result is proportional to transmitted optical power:
        P_thru = V_THRU / (responsivity · R_load)
    """
    netlist0 = build_netlist(model, ps_params, coupler_kappa_l,
                             wavelengths_nm[0], v_pn, i_heat)
    ckt = fc.Circuit()
    ckt.load_str(netlist0)

    v_thru = np.empty(len(wavelengths_nm))
    for i, wl in enumerate(wavelengths_nm):
        ckt.set_param("XCWL1", "wavelength_nm", wl)
        ckt.set_param("V_PN", "dc", v_pn)
        ckt.set_param("IHEAT", "dc", i_heat)
        try:
            r = ckt.run("op")
            v_thru[i] = max(float(r["V_THRU"][0]), 1e-30)
        except Exception:
            v_thru[i] = 1e-30
    return v_thru


def find_extremum_wavelength(
    wavelengths_nm: np.ndarray,
    v_thru: np.ndarray,
    kind: str = "min",
) -> float:
    """
    Return the wavelength of the first transmission minimum (ring dip) or
    maximum (MZI fringe peak) using quadratic interpolation around the
    coarse-grid extremum.
    """
    idx = int(np.argmin(v_thru) if kind == "min" else np.argmax(v_thru))
    # Quadratic refinement using neighbours.
    if 1 <= idx < len(v_thru) - 1:
        y0, y1, y2 = v_thru[idx - 1], v_thru[idx], v_thru[idx + 1]
        d = 2 * y1 - y0 - y2
        if abs(d) > 1e-30:
            frac = (y0 - y2) / (2 * d)
            dw = wavelengths_nm[idx] - wavelengths_nm[max(idx - 1, 0)]
            return wavelengths_nm[idx] + frac * dw
    return wavelengths_nm[idx]


def to_dB(v: np.ndarray) -> np.ndarray:
    """Convert linear voltage to dB, normalised to the max value."""
    v = np.clip(v, 1e-30, None)
    return 10.0 * np.log10(v / v.max())


# ── parameter vector helpers ──────────────────────────────────────────────────

def _free_specs(all_specs: list[ParamSpec]) -> list[ParamSpec]:
    return [s for s in all_specs if not s.fixed]


def _pack(free_specs: list[ParamSpec]) -> np.ndarray:
    return np.array([s.value for s in free_specs])


def _unpack(x: np.ndarray, all_specs: list[ParamSpec]) -> dict[str, float]:
    """Merge optimiser values back into the full named parameter dict."""
    result = {}
    xi = 0
    for s in all_specs:
        if s.fixed:
            result[s.name] = s.value
        else:
            result[s.name] = float(x[xi])
            xi += 1
    return result


def _all_bounds(all_specs: list[ParamSpec]) -> list[tuple[float, float]]:
    return [s.bounds for s in all_specs if not s.fixed]


# ── datasets ──────────────────────────────────────────────────────────────────

@dataclass
class Dataset:
    """One experimental dataset to contribute to the loss function."""
    kind: str                          # "spectrum" | "eo_tuning" | "thermal_tuning"
    x: np.ndarray                      # wavelength_nm (spectrum) or voltage (tuning)
    y: np.ndarray                      # transmission_dB (spectrum) or resonance_nm (tuning)
    weight: float = 1.0
    extremum: str = "min"              # "min" (ring dip) or "max" (MZI peak) for tuning curves


def load_spectrum(path: str) -> Dataset:
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return Dataset(kind="spectrum", x=d[:, 0], y=d[:, 1])


def load_eo_tuning(path: str) -> Dataset:
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return Dataset(kind="eo_tuning", x=d[:, 0], y=d[:, 1])


def load_thermal_tuning(path: str) -> Dataset:
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return Dataset(kind="thermal_tuning", x=d[:, 0], y=d[:, 1])


# ── synthetic placeholder data ────────────────────────────────────────────────

def _synthetic_spectrum(
    model: str,
    all_specs: list[ParamSpec],
    coupler_kappa_l: float,
    center_nm: float = 1550.0,
    span_nm: float = 5.0,
    n_points: int = 61,
    noise_dB: float = 0.2,
    rng_seed: int = 42,
) -> Dataset:
    """
    Generate a plausible synthetic MZI spectrum using the simulator itself.
    Useful for smoke-testing the fitting pipeline without real measurements.
    """
    rng = np.random.default_rng(rng_seed)
    wls = np.linspace(center_nm - span_nm / 2, center_nm + span_nm / 2, n_points)
    ps_p = {s.name: s.value for s in all_specs}
    v = wavelength_sweep(model, ps_p, coupler_kappa_l, wls)
    t_dB = to_dB(v) + rng.normal(0.0, noise_dB, size=len(wls))
    return Dataset(kind="spectrum", x=wls, y=t_dB)


def _synthetic_eo_tuning(
    model: str,
    all_specs: list[ParamSpec],
    coupler_kappa_l: float,
    v_pn_range: tuple = (-4.0, 0.0),
    n_points: int = 9,
    center_nm: float = 1550.0,
    scan_span_nm: float = 2.0,
    noise_nm: float = 0.005,
    rng_seed: int = 43,
) -> Optional[Dataset]:
    """Synthetic resonance-shift vs PN-bias dataset."""
    info = MODELS[model]
    if not info["has_pn"]:
        return None
    rng = np.random.default_rng(rng_seed)
    wls = np.linspace(center_nm - scan_span_nm / 2,
                      center_nm + scan_span_nm / 2, 41)
    voltages = np.linspace(v_pn_range[0], v_pn_range[1], n_points)
    ps_p = {s.name: s.value for s in all_specs}
    resonances = []
    for v in voltages:
        vt = wavelength_sweep(model, ps_p, coupler_kappa_l, wls, v_pn=v)
        res = find_extremum_wavelength(wls, vt, kind="min")
        resonances.append(res)
    resonances = np.array(resonances) + rng.normal(0.0, noise_nm, size=len(voltages))
    return Dataset(kind="eo_tuning", x=voltages, y=resonances)


def _synthetic_thermal_tuning(
    model: str,
    all_specs: list[ParamSpec],
    coupler_kappa_l: float,
    i_heat_range: tuple = (0.0, 20e-3),
    n_points: int = 9,
    center_nm: float = 1550.0,
    scan_span_nm: float = 4.0,
    noise_nm: float = 0.02,
    rng_seed: int = 44,
) -> Optional[Dataset]:
    """Synthetic resonance-shift vs heater-current dataset."""
    info = MODELS[model]
    if not info["has_heater"]:
        return None
    rng = np.random.default_rng(rng_seed)
    wls = np.linspace(center_nm - scan_span_nm / 2,
                      center_nm + scan_span_nm / 2, 41)
    currents = np.linspace(i_heat_range[0], i_heat_range[1], n_points)
    ps_p = {s.name: s.value for s in all_specs}
    resonances = []
    for i in currents:
        vt = wavelength_sweep(model, ps_p, coupler_kappa_l, wls, i_heat=i)
        res = find_extremum_wavelength(wls, vt, kind="min")
        resonances.append(res)
    resonances = np.array(resonances) + rng.normal(0.0, noise_nm, size=len(currents))
    return Dataset(kind="thermal_tuning", x=currents, y=resonances)


# ── loss function ─────────────────────────────────────────────────────────────

def _loss(
    x: np.ndarray,
    all_specs: list[ParamSpec],
    coupler_specs: list[ParamSpec],
    model: str,
    datasets: list[Dataset],
    scan_wls_nm: np.ndarray,
    center_nm: float,
) -> float:
    ps_params = _unpack(x, all_specs)
    kappa_l = _unpack(x, coupler_specs)["kappa_l"]

    total = 0.0

    for ds in datasets:
        if ds.kind == "spectrum":
            wls = ds.x
            vt = wavelength_sweep(model, ps_params, kappa_l, wls)
            sim_dB = to_dB(vt)
            # Align baselines: remove DC offset (overall insertion-loss
            # calibration) by subtracting means.
            diff = (ds.y - ds.y.mean()) - (sim_dB - sim_dB.mean())
            total += ds.weight * np.mean(diff ** 2)

        elif ds.kind == "eo_tuning":
            residuals = []
            for v_pn, res_meas in zip(ds.x, ds.y):
                vt = wavelength_sweep(model, ps_params, kappa_l, scan_wls_nm,
                                      v_pn=v_pn)
                res_sim = find_extremum_wavelength(scan_wls_nm, vt, ds.extremum)
                residuals.append(res_sim - res_meas)
            # Normalise residuals in nm; remove a common offset (absolute
            # resonance position is set by n_eff which may be fit separately).
            r = np.array(residuals)
            total += ds.weight * np.mean((r - r.mean()) ** 2) * 1e4  # nm²→scale

        elif ds.kind == "thermal_tuning":
            residuals = []
            for i_heat, res_meas in zip(ds.x, ds.y):
                vt = wavelength_sweep(model, ps_params, kappa_l, scan_wls_nm,
                                      i_heat=i_heat)
                res_sim = find_extremum_wavelength(scan_wls_nm, vt, ds.extremum)
                residuals.append(res_sim - res_meas)
            r = np.array(residuals)
            total += ds.weight * np.mean((r - r.mean()) ** 2) * 1e4

    return total


# ── optimiser ─────────────────────────────────────────────────────────────────

def fit(
    model: str,
    datasets: list[Dataset],
    *,
    all_ps_specs: Optional[list[ParamSpec]] = None,
    coupler_specs: Optional[list[ParamSpec]] = None,
    center_nm: float = 1550.0,
    scan_span_nm: float = 3.0,
    scan_points: int = 51,
    maxiter: int = 200,
    popsize: int = 12,
    seed: int = 0,
    verbose: bool = True,
) -> dict[str, float]:
    """
    Fit PS and coupler parameters to the provided datasets.

    Returns the best-fit parameter dict.
    """
    if all_ps_specs is None:
        all_ps_specs = deepcopy(MODELS[model]["ps_params"])
    if coupler_specs is None:
        coupler_specs = deepcopy(_COUPLER)

    # Combine into one vector for the optimiser.
    combined_specs = all_ps_specs + coupler_specs
    free = _free_specs(combined_specs)
    if not free:
        print("No free parameters — nothing to fit.")
        return _unpack(np.array([]), combined_specs)

    x0 = _pack(free)
    bounds = _all_bounds(combined_specs)
    scan_wls = np.linspace(center_nm - scan_span_nm / 2,
                           center_nm + scan_span_nm / 2, scan_points)

    call_count = [0]

    def objective(x):
        call_count[0] += 1
        return _loss(x, all_ps_specs, coupler_specs, model, datasets,
                     scan_wls, center_nm)

    if verbose:
        print(f"\nFitting model: {model}")
        print(f"  Free params ({len(free)}): {[s.name for s in free]}")
        print(f"  Datasets: {[d.kind for d in datasets]}")
        print(f"  Running differential_evolution"
              f" (maxiter={maxiter}, popsize={popsize}) …")

    result = differential_evolution(
        objective,
        bounds,
        x0=x0,
        maxiter=maxiter,
        popsize=popsize,
        seed=seed,
        tol=1e-5,
        mutation=(0.5, 1.5),
        recombination=0.7,
        polish=True,
        disp=False,
    )

    if verbose:
        print(f"  Optimisation converged: {result.success}  "
              f"(loss={result.fun:.4e}, {call_count[0]} evals)")

    best = _unpack(result.x, combined_specs)
    return best


# ── result reporting & plotting ───────────────────────────────────────────────

def print_results(best: dict, all_specs: list[ParamSpec], coupler_specs: list[ParamSpec]):
    all_combined = all_specs + coupler_specs
    print("\n── Best-fit parameters ──────────────────────────────────")
    for s in all_combined:
        val = best.get(s.name, s.value)
        flag = "  [fixed]" if s.fixed else ""
        print(f"  {s.name:<18s} = {val:>14.6g}   # {s.description}{flag}")
    print()


def plot_fit(
    model: str,
    best: dict,
    datasets: list[Dataset],
    center_nm: float = 1550.0,
    span_nm: float = 5.0,
    out_path: str = "fit_result.png",
):
    """Save a comparison plot of measured vs fitted simulation."""
    if not _HAS_MPL:
        print("matplotlib not available — skipping plot.")
        return

    kappa_l = best.get("kappa_l", _COUPLER[0].value)
    wls = np.linspace(center_nm - span_nm / 2, center_nm + span_nm / 2, 201)
    vt_fit = wavelength_sweep(model, best, kappa_l, wls)
    sim_dB = to_dB(vt_fit)

    fig, axes = plt.subplots(1, max(1, sum(
        1 for d in datasets if d.kind == "spectrum"
    ) + sum(
        1 for d in datasets if d.kind in ("eo_tuning", "thermal_tuning")
    )), figsize=(5 * max(1, len(datasets)), 4), squeeze=False)
    axes = axes.ravel()

    panel = 0
    for ds in datasets:
        ax = axes[panel]
        if ds.kind == "spectrum":
            ax.plot(wls, sim_dB, label="simulation (fit)", lw=2)
            ax.plot(ds.x, ds.y - ds.y.mean() + sim_dB.mean(),
                    "o", ms=4, label="measured (aligned)")
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Transmission (dB, normalised)")
            ax.set_title(f"Spectrum — {model}")
            ax.legend(fontsize=8)
            panel += 1
        elif ds.kind in ("eo_tuning", "thermal_tuning"):
            scan_wls = np.linspace(center_nm - span_nm / 2,
                                   center_nm + span_nm / 2, 61)
            sim_res = []
            for v in ds.x:
                v_pn = v if ds.kind == "eo_tuning" else 0.0
                i_heat = v if ds.kind == "thermal_tuning" else 0.0
                vt = wavelength_sweep(model, best, kappa_l, scan_wls,
                                      v_pn=v_pn, i_heat=i_heat)
                sim_res.append(find_extremum_wavelength(scan_wls, vt, ds.extremum))
            xlabel = "V_PN (V)" if ds.kind == "eo_tuning" else "I_heat (A)"
            ax.plot(ds.x, sim_res, label="simulation (fit)", lw=2)
            ax.plot(ds.x, ds.y, "o", ms=4, label="measured")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Resonance / fringe wavelength (nm)")
            ax.set_title(f"Tuning — {ds.kind}")
            ax.legend(fontsize=8)
            panel += 1

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Plot saved: {out_path}")
    plt.close(fig)


# ── CLI ───────────────────────────────────────────────────────────────────────

def _list_models():
    print("\nAvailable phase-shifter models:\n")
    for name, info in MODELS.items():
        print(f"  {name}")
        print(f"    {info['description']}")
        free_p = [s.name for s in info["ps_params"] if not s.fixed]
        fixed_p = [s.name for s in info["ps_params"] if s.fixed]
        print(f"    Free  : {free_p}")
        print(f"    Fixed : {fixed_p}")
        print()


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--model", default="fc_pn_th_ps",
                    choices=list(MODELS.keys()),
                    help="phase-shifter model to fit (default: fc_pn_th_ps)")
    ap.add_argument("--spectrum", metavar="CSV",
                    help="path to transmission spectrum CSV (wavelength_nm, transmission_dB)")
    ap.add_argument("--eo", metavar="CSV",
                    help="path to EO tuning CSV (v_pn_V, resonance_nm)")
    ap.add_argument("--thermal", metavar="CSV",
                    help="path to thermal tuning CSV (i_heat_A, resonance_nm)")
    ap.add_argument("--center-nm", type=float, default=1550.0,
                    help="centre wavelength for sweeps (default: 1550 nm)")
    ap.add_argument("--span-nm", type=float, default=5.0,
                    help="wavelength span for sweeps (default: 5 nm)")
    ap.add_argument("--maxiter", type=int, default=200,
                    help="differential_evolution max iterations (default: 200)")
    ap.add_argument("--popsize", type=int, default=12,
                    help="population size per free parameter (default: 12)")
    ap.add_argument("--plot", metavar="PNG", default="fit_result.png",
                    help="output plot filename (default: fit_result.png)")
    ap.add_argument("--no-plot", action="store_true",
                    help="suppress plotting")
    ap.add_argument("--list-models", action="store_true",
                    help="print available models and parameters, then exit")
    args = ap.parse_args()

    if args.list_models:
        _list_models()
        return

    model = args.model
    all_ps_specs = deepcopy(MODELS[model]["ps_params"])
    coupler_specs = deepcopy(_COUPLER)

    # ── load or synthesise datasets ──────────────────────────────────────────
    datasets: list[Dataset] = []

    if args.spectrum:
        print(f"Loading spectrum from {args.spectrum}")
        datasets.append(load_spectrum(args.spectrum))
    if args.eo:
        print(f"Loading EO tuning from {args.eo}")
        datasets.append(load_eo_tuning(args.eo))
    if args.thermal:
        print(f"Loading thermal tuning from {args.thermal}")
        datasets.append(load_thermal_tuning(args.thermal))

    if not datasets:
        print("No data files provided — running smoke-test with synthetic data.")
        print("(Pass --spectrum / --eo / --thermal to use real measurements.)\n")

        nominal_ps = {s.name: s.value for s in all_ps_specs}
        nominal_kl = coupler_specs[0].value

        # Perturb defaults slightly to give the optimiser something to recover.
        perturbed = {**nominal_ps,
                     "n_eff": nominal_ps.get("n_eff", 2.40) * 1.005,
                     "alpha_db_cm": nominal_ps.get("alpha_db_cm", 10.0) * 0.85}

        syn_ps_specs = deepcopy(all_ps_specs)
        for s in syn_ps_specs:
            if s.name in perturbed:
                s.value = perturbed[s.name]

        spec_ds = _synthetic_spectrum(model, syn_ps_specs, nominal_kl,
                                      center_nm=args.center_nm,
                                      span_nm=args.span_nm)
        datasets.append(spec_ds)

        eo_ds = _synthetic_eo_tuning(model, syn_ps_specs, nominal_kl,
                                     center_nm=args.center_nm,
                                     scan_span_nm=args.span_nm)
        if eo_ds:
            datasets.append(eo_ds)

        th_ds = _synthetic_thermal_tuning(model, syn_ps_specs, nominal_kl,
                                          center_nm=args.center_nm,
                                          scan_span_nm=args.span_nm)
        if th_ds:
            datasets.append(th_ds)

    # ── fit ──────────────────────────────────────────────────────────────────
    best = fit(
        model,
        datasets,
        all_ps_specs=all_ps_specs,
        coupler_specs=coupler_specs,
        center_nm=args.center_nm,
        scan_span_nm=args.span_nm,
        scan_points=51,
        maxiter=args.maxiter,
        popsize=args.popsize,
        verbose=True,
    )

    print_results(best, all_ps_specs, coupler_specs)

    if not args.no_plot:
        plot_fit(model, best, datasets,
                 center_nm=args.center_nm,
                 span_nm=args.span_nm,
                 out_path=args.plot)


if __name__ == "__main__":
    main()
