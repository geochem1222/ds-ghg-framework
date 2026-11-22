"""Two-layer implementation of the depth-state-greenhouse-gas model (DS-GHG).

The model is an executable hypothesis generator, not a calibrated forecasting
tool. Concentrations are in micromolar (numerically equal to mmol m-3), areal
process rates are in mmol m-2 d-1, depth is in m, and time is in days.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent

STATE_NAMES = [
    "Te_C",
    "Th_C",
    "Ts_C",
    "O2e_uM",
    "O2h_uM",
    "NO3h_uM",
    "SO4h_uM",
    "FeIII_mmol_m2",
    "CH4e_uM",
    "CH4h_uM",
    "CH4s_uM",
    "CO2e_uM",
    "CO2h_uM",
    "CO2s_uM",
    "N2Oe_uM",
    "N2Oh_uM",
    "N2Os_uM",
    "bubble_store_mmol_m2",
    "stratification_index",
    "cum_CO2_flux_mmol_m2",
    "cum_CH4_diff_flux_mmol_m2",
    "cum_CH4_ebu_flux_mmol_m2",
    "cum_N2O_flux_mmol_m2",
    "cum_methanogenesis_mmol_m2",
    "cum_AeMO_mmol_m2",
    "cum_AOM_NO3_mmol_m2",
    "cum_AOM_SO4_mmol_m2",
    "cum_AOM_Fe_mmol_m2",
    "cum_nitrification_N2O_mmol_m2",
    "cum_denitrification_N2O_mmol_m2",
    "cum_N2O_reduction_mmol_m2",
]
IDX = {name: i for i, name in enumerate(STATE_NAMES)}
NONNEGATIVE_STATE_NAMES = [
    "O2e_uM",
    "O2h_uM",
    "NO3h_uM",
    "SO4h_uM",
    "FeIII_mmol_m2",
    "CH4e_uM",
    "CH4h_uM",
    "CH4s_uM",
    "CO2e_uM",
    "CO2h_uM",
    "CO2s_uM",
    "N2Oe_uM",
    "N2Oh_uM",
    "N2Os_uM",
    "bubble_store_mmol_m2",
]


@dataclass(frozen=True)
class Scenario:
    scenario_id: str
    forcing: str
    depth_m: float
    duration_d: float
    Kd_m_1: float
    mixing_factor: float
    productivity_mmolO2_m2_d: float
    NO3_boundary_uM: float
    SO4_boundary_uM: float
    hydrodynamic_regime: str = "mixed"
    notes: str = ""


def read_parameters(path: Path = ROOT / "parameters.csv") -> Dict[str, float]:
    table = pd.read_csv(path)
    return {str(row["name"]): float(row["value"]) for _, row in table.iterrows()}


def read_scenarios(path: Path = ROOT / "scenarios.csv") -> list[Scenario]:
    rows: list[Scenario] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                Scenario(
                    scenario_id=row["scenario_id"],
                    forcing=row["forcing"],
                    depth_m=float(row["depth_m"]),
                    duration_d=float(row["duration_d"]),
                    Kd_m_1=float(row["Kd_m_1"]),
                    mixing_factor=float(row["mixing_factor"]),
                    productivity_mmolO2_m2_d=float(row["productivity_mmolO2_m2_d"]),
                    NO3_boundary_uM=float(row["NO3_boundary_uM"]),
                    SO4_boundary_uM=float(row["SO4_boundary_uM"]),
                    hydrodynamic_regime=(row.get("hydrodynamic_regime", "mixed").strip() or "mixed"),
                    notes=row.get("notes", ""),
                )
            )
    return rows


def smooth_box(t: float, start: float, stop: float, edge: float = 0.35) -> float:
    return 0.5 * (np.tanh((t - start) / edge) - np.tanh((t - stop) / edge))


def surface_temperature(t: float, forcing: str, p: Dict[str, float]) -> float:
    base = p["surface_temperature_base_C"]
    diel = p["surface_temperature_diel_amplitude_C"] * np.sin(2 * np.pi * (t - 0.25))
    if forcing in {"diurnal", "drawdown", "rise_fall"}:
        return base + diel
    if forcing == "heatwave":
        return base + diel + p["heatwave_delta_C"] * smooth_box(t, 6.0, 14.0)
    if forcing == "cooling":
        return base + diel - p["cooling_delta_C"] * smooth_box(t, 6.0, 14.0)
    if forcing == "seasonal":
        return base + p["seasonal_amplitude_C"] * np.sin(2 * np.pi * t / 90.0)
    return base


def depth_trajectory(t: float, s: Scenario) -> tuple[float, float]:
    """Return water depth and its derivative in m and m d-1."""
    d0 = s.depth_m
    if s.forcing == "drawdown":
        low = max(0.05, 0.2 * d0)
        if 8.0 <= t < 10.0:
            return d0 + (low - d0) * (t - 8.0) / 2.0, (low - d0) / 2.0
        if t >= 10.0:
            return low, 0.0
    if s.forcing == "rise_fall":
        low = max(0.05, 0.2 * d0)
        if t < 5.0:
            return low, 0.0
        if t < 8.0:
            return low + (d0 - low) * (t - 5.0) / 3.0, (d0 - low) / 3.0
        if t < 12.0:
            return d0, 0.0
        if t < 14.0:
            return d0 + (low - d0) * (t - 12.0) / 2.0, (low - d0) / 2.0
        if t < 18.0:
            return low, 0.0
        if t < 20.0:
            target = 0.75 * d0
            return low + (target - low) * (t - 18.0) / 2.0, (target - low) / 2.0
        return 0.75 * d0, 0.0
    return d0, 0.0


def layer_depths(depth: float, p: Dict[str, float]) -> tuple[float, float]:
    """Return epilimnetic and hypolimnetic layer thicknesses."""
    depth = max(depth, 0.02)
    fraction = np.clip(p["epilimnion_fraction"], 0.1, 0.9)
    minimum = min(p["minimum_layer_depth_m"], 0.45 * depth)
    epi = float(np.clip(fraction * depth, minimum, depth - minimum))
    return epi, depth - epi


def daylight_fraction(t: float) -> float:
    return max(0.0, float(np.sin(2.0 * np.pi * (t - 0.25))))


def ventilation_window(t: float, s: Scenario) -> float:
    if s.hydrodynamic_regime != "ventilation":
        return 0.0
    return smooth_box(t, 0.45 * s.duration_d, 0.55 * s.duration_d, edge=0.12)


def mixing_multiplier(t: float, s: Scenario, p: Dict[str, float]) -> float:
    regime = s.hydrodynamic_regime
    if regime == "mixed":
        return s.mixing_factor * p["mixed_regime_multiplier"]
    if regime == "diel":
        night = 1.0 - daylight_fraction(t)
        return s.mixing_factor * (1.0 + p["diel_night_mixing_multiplier"] * night)
    if regime == "ventilation":
        event = ventilation_window(t, s)
        return s.mixing_factor * (1.0 + (p["ventilation_mixing_multiplier"] - 1.0) * event)
    if regime != "persistent":
        raise ValueError(f"Unknown hydrodynamic regime: {regime}")
    return s.mixing_factor


def hydrodynamic_diagnostics(
    t: float,
    y: np.ndarray,
    s: Scenario,
    p: Dict[str, float],
) -> Dict[str, float | str]:
    depth, depth_dot = depth_trajectory(t, s)
    epi_depth, hypo_depth = layer_depths(depth, p)
    strength = float(np.clip(y[IDX["stratification_index"]], 0.0, 1.0))
    mixing = max(mixing_multiplier(t, s, p), 0.05)
    mixed_exchange = p["interlayer_exchange_mixed_m_d"]
    stratified_exchange = p["interlayer_exchange_stratified_m_d"]
    exchange = stratified_exchange + (mixed_exchange - stratified_exchange) * (1.0 - strength) ** 2
    exchange *= max(1.0, mixing / max(s.mixing_factor, 0.05))

    event = ventilation_window(t, s)
    if event > 0.1:
        state = "ventilating"
    elif strength >= p["stratification_state_threshold"]:
        state = "stratified"
    else:
        state = "mixed"
    return {
        "depth": depth,
        "depth_dot": depth_dot,
        "epi_depth": epi_depth,
        "hypo_depth": hypo_depth,
        "mixing_multiplier": mixing,
        "interlayer_exchange": exchange,
        "ventilation_intensity": event,
        "hydrodynamic_state": state,
    }


def stratification_tendency(
    t: float,
    y: np.ndarray,
    s: Scenario,
    p: Dict[str, float],
    hydro: Dict[str, float | str],
) -> float:
    strength = float(np.clip(y[IDX["stratification_index"]], 0.0, 1.0))
    mixing = float(hydro["mixing_multiplier"])
    if s.hydrodynamic_regime == "mixed":
        return -p["stratification_decay_d_1"] * mixing * strength

    depth = float(hydro["depth"])
    surface = surface_temperature(t, s.forcing, p)
    bottom = float(y[IDX["Th_C"]])
    thermal_drive = max(surface - bottom, 0.0) / p["stratification_temperature_scale_C"]
    depth_factor = depth / (depth + p["stratification_depth_half_m"])
    build = (
        p["stratification_build_d_1"]
        * depth_factor
        * daylight_fraction(t)
        * thermal_drive
        * (1.0 - strength)
        / (1.0 + mixing)
    )
    convective_drive = max(bottom - surface, 0.0) / p["stratification_temperature_scale_C"]
    decay = p["stratification_decay_d_1"] * mixing * strength
    decay += p["convective_mixing_decay_d_1"] * convective_drive * strength
    decay += p["ventilation_stratification_decay_d_1"] * float(hydro["ventilation_intensity"]) * strength
    return build - decay


def oxygen_saturation_uM(temp_c: float) -> float:
    mg_l = 14.621 - 0.41022 * temp_c + 0.007991 * temp_c**2 - 0.000077774 * temp_c**3
    return max(0.0, mg_l) * 1000.0 / 32.0


def q10(temp_c: float, value: float, reference_c: float = 20.0) -> float:
    return value ** ((temp_c - reference_c) / 10.0)


def monod(value: float, half_saturation: float) -> float:
    value = max(0.0, value)
    return value / (half_saturation + value + 1e-12)


def process_rates(
    t: float,
    y: np.ndarray,
    s: Scenario,
    p: Dict[str, float],
) -> Dict[str, float | str]:
    v = {name: max(0.0, float(y[i])) for name, i in IDX.items()}
    hydro = hydrodynamic_diagnostics(t, y, s, p)
    depth = float(hydro["depth"])
    depth_dot = float(hydro["depth_dot"])
    epi_depth = float(hydro["epi_depth"])
    hypo_depth = float(hydro["hypo_depth"])
    exchange = float(hydro["interlayer_exchange"])
    temp = v["Ts_C"]
    o2e = v["O2e_uM"]
    o2h = v["O2h_uM"]
    no3 = v["NO3h_uM"]
    so4 = v["SO4h_uM"]
    fe = v["FeIII_mmol_m2"]
    ch4e = v["CH4e_uM"]
    ch4h = v["CH4h_uM"]
    ch4s = v["CH4s_uM"]
    n2os = v["N2Os_uM"]

    oxic_e = monod(o2e, p["K_O2_uM"])
    oxic_h = monod(o2h, p["K_O2_uM"])
    anoxic = p["K_O2_inhibition_uM"] / (p["K_O2_inhibition_uM"] + o2h + 1e-12)
    labile_carbon = p["labile_carbon_fraction"]

    mineralization = p["V_OM_mmol_m2_d"] * q10(temp, p["Q10_OM"]) * labile_carbon
    meth_inhibition = (
        anoxic
        * p["K_NO3_inhibition_uM"] / (p["K_NO3_inhibition_uM"] + no3 + 1e-12)
        * p["K_SO4_inhibition_uM"] / (p["K_SO4_inhibition_uM"] + so4 + 1e-12)
        * p["K_Fe_inhibition_mmol_m2"] / (p["K_Fe_inhibition_mmol_m2"] + fe + 1e-12)
    )
    methanogenesis = (
        p["V_meth_mmol_m2_d"]
        * q10(temp, p["Q10_meth"])
        * labile_carbon
        * meth_inhibition
    )
    aemo_s = (
        p["V_AeMO_sed_mmol_m2_d"]
        * q10(temp, p["Q10_AeMO"])
        * monod(ch4s, p["K_CH4_AeMO_uM"])
        * oxic_h
    )
    aemo_e_areal = (
        p["V_AeMO_water_uM_d"]
        * q10(v["Te_C"], p["Q10_AeMO"])
        * monod(ch4e, p["K_CH4_AeMO_uM"])
        * oxic_e
        * epi_depth
    )
    aemo_h_areal = (
        p["V_AeMO_water_uM_d"]
        * q10(v["Th_C"], p["Q10_AeMO"])
        * monod(ch4h, p["K_CH4_AeMO_uM"])
        * oxic_h
        * hypo_depth
    )

    ch4_lim = monod(ch4s, p["K_CH4_AOM_uM"])
    aom_no3 = (
        p["V_AOM_NO3_mmol_m2_d"]
        * q10(temp, p["Q10_AOM"])
        * ch4_lim
        * monod(no3, p["K_NO3_AOM_uM"])
        * anoxic
    )
    aom_so4 = (
        p["V_AOM_SO4_mmol_m2_d"]
        * q10(temp, p["Q10_AOM"])
        * ch4_lim
        * monod(so4, p["K_SO4_AOM_uM"])
        * anoxic
    )
    aom_fe = (
        p["V_AOM_Fe_mmol_m2_d"]
        * q10(temp, p["Q10_AOM"])
        * ch4_lim
        * monod(fe, p["K_Fe_AOM_mmol_m2"])
        * anoxic
    )

    nitrification = p["V_nitrification_mmol_m2_d"] * q10(temp, p["Q10_nitrification"]) * oxic_h
    suboxic_peak = 1.0 + p["nitrification_suboxic_multiplier"] * np.exp(
        -((o2h - p["nitrification_O2_optimum_uM"]) / p["nitrification_O2_width_uM"]) ** 2
    )
    nitrification_n2o = nitrification * p["nitrification_N2O_yield"] * suboxic_peak
    denitrification = (
        p["V_denitrification_mmol_m2_d"]
        * q10(temp, p["Q10_denitrification"])
        * monod(no3, p["K_NO3_denitrification_uM"])
        * anoxic
        * labile_carbon
    )
    denitrification_n2o = denitrification * p["denitrification_N2O_yield"]
    n2o_reduction = (
        p["V_N2O_reduction_mmol_m2_d"]
        * q10(temp, p["Q10_N2O_reduction"])
        * monod(n2os, p["K_N2O_reduction_uM"])
        * anoxic
        * labile_carbon
    )

    ch4_swi = p["k_CH4_swi_m_d"] * max(ch4s - ch4h, 0.0)
    co2_swi = p["k_CO2_swi_m_d"] * max(v["CO2s_uM"] - v["CO2h_uM"], 0.0)
    n2o_swi = p["k_N2O_swi_m_d"] * (v["N2Os_uM"] - v["N2Oh_uM"])

    pressure_ratio = 1.0 + p["rho_water_kg_m3"] * p["gravity_m_s2"] * depth / p["atmospheric_pressure_Pa"]
    bubble_threshold = p["bubble_threshold_1atm_uM"] * pressure_ratio + p["bubble_capillary_equivalent_uM"]
    bubble_formation = (
        max(ch4s - bubble_threshold, 0.0)
        * p["sediment_porosity"]
        * p["active_sediment_depth_m"]
        / p["bubble_nucleation_time_d"]
    )
    bubble_capacity = p["bubble_capacity_base_mmol_m2"] + p["bubble_capacity_depth_slope_mmol_m3"] * depth
    pressure_release = p["drawdown_release_factor"] * max(-depth_dot, 0.0) * v["bubble_store_mmol_m2"] / max(depth, 0.02)
    excess_release = max(v["bubble_store_mmol_m2"] - bubble_capacity, 0.0) / p["bubble_release_time_d"]
    bubble_release = pressure_release + excess_release
    bubble_survival = np.exp(-depth / p["bubble_dissolution_length_m"])
    bubble_to_air = bubble_release * bubble_survival
    bubble_to_water = bubble_release - bubble_to_air

    ch4_interlayer = exchange * (ch4h - ch4e)
    co2_interlayer = exchange * (v["CO2h_uM"] - v["CO2e_uM"])
    n2o_interlayer = exchange * (v["N2Oh_uM"] - v["N2Oe_uM"])
    o2_interlayer = exchange * (o2e - o2h)

    co2_flux = p["k_CO2_air_m_d"] * (v["CO2e_uM"] - p["CO2_equilibrium_uM"])
    ch4_diff_flux = p["k_CH4_air_m_d"] * (ch4e - p["CH4_equilibrium_uM"])
    n2o_flux = p["k_N2O_air_m_d"] * (v["N2Oe_uM"] - p["N2O_equilibrium_uM"])

    return {
        "depth": depth,
        "depth_dot": depth_dot,
        "epi_depth": epi_depth,
        "hypo_depth": hypo_depth,
        "mixing_multiplier": float(hydro["mixing_multiplier"]),
        "interlayer_exchange": exchange,
        "ventilation_intensity": float(hydro["ventilation_intensity"]),
        "hydrodynamic_state": str(hydro["hydrodynamic_state"]),
        "surface_temperature": surface_temperature(t, s.forcing, p),
        "mineralization": mineralization,
        "methanogenesis": methanogenesis,
        "aemo_s": aemo_s,
        "aemo_e": aemo_e_areal,
        "aemo_h": aemo_h_areal,
        "aemo_w": aemo_e_areal + aemo_h_areal,
        "aom_no3": aom_no3,
        "aom_so4": aom_so4,
        "aom_fe": aom_fe,
        "nitrification": nitrification,
        "nitrification_n2o": nitrification_n2o,
        "denitrification": denitrification,
        "denitrification_n2o": denitrification_n2o,
        "n2o_reduction": n2o_reduction,
        "ch4_swi": ch4_swi,
        "co2_swi": co2_swi,
        "n2o_swi": n2o_swi,
        "ch4_interlayer": ch4_interlayer,
        "co2_interlayer": co2_interlayer,
        "n2o_interlayer": n2o_interlayer,
        "o2_interlayer": o2_interlayer,
        "bubble_formation": bubble_formation,
        "bubble_release": bubble_release,
        "bubble_to_air": bubble_to_air,
        "bubble_to_water": bubble_to_water,
        "co2_flux": co2_flux,
        "ch4_diff_flux": ch4_diff_flux,
        "n2o_flux": n2o_flux,
    }


def model_rhs(t: float, y: np.ndarray, s: Scenario, p: Dict[str, float]) -> np.ndarray:
    v = {name: max(0.0, float(y[i])) for name, i in IDX.items()}
    r = process_rates(t, y, s, p)
    depth = max(float(r["depth"]), 0.02)
    epi_depth = max(float(r["epi_depth"]), 0.005)
    hypo_depth = max(float(r["hypo_depth"]), 0.005)
    exchange = float(r["interlayer_exchange"])
    depth_dot = float(r["depth_dot"])
    pore_inventory_depth = p["sediment_porosity"] * p["active_sediment_depth_m"]
    dy = np.zeros_like(y)

    tau_epi = (
        p["rho_water_kg_m3"] * p["heat_capacity_water_J_kg_K"] * epi_depth
        / (p["effective_surface_heat_transfer_W_m2_K"] * 86400.0)
    )
    bed_light = np.exp(-s.Kd_m_1 * depth)
    daylight = daylight_fraction(t)
    heat_exchange = exchange
    dy[IDX["Te_C"]] = (
        (float(r["surface_temperature"]) - v["Te_C"]) / tau_epi
        + heat_exchange * (v["Th_C"] - v["Te_C"]) / epi_depth
    )
    dy[IDX["Th_C"]] = (
        heat_exchange * (v["Te_C"] - v["Th_C"]) / hypo_depth
        + (v["Ts_C"] - v["Th_C"]) / p["water_sediment_heat_exchange_d"]
    )
    dy[IDX["Ts_C"]] = (
        (v["Th_C"] - v["Ts_C"]) / p["sediment_thermal_response_d"]
        + p["bed_light_heating_C_d"] * bed_light * daylight
    )
    dy[IDX["stratification_index"]] = stratification_tendency(t, y, s, p, hydrodynamic_diagnostics(t, y, s, p))

    o2_sat = oxygen_saturation_uM(v["Te_C"])
    reaeration = (
        p["k_O2_air_m_d"]
        * float(r["mixing_multiplier"])
        * (o2_sat - v["O2e_uM"])
    )
    oxygen_availability_e = monod(v["O2e_uM"], p["K_O2_consumption_uM"])
    oxygen_availability_h = monod(v["O2h_uM"], p["K_O2_consumption_uM"])
    respiration_epi = (
        p["water_respiration_uM_d"] * q10(v["Te_C"], p["Q10_OM"]) * epi_depth * oxygen_availability_e
    )
    respiration_hypo = (
        p["water_respiration_uM_d"] * q10(v["Th_C"], p["Q10_OM"]) * hypo_depth * oxygen_availability_h
    )
    oxygen_demand = (
        p["sediment_oxygen_demand_mmol_m2_d"]
        * q10(v["Ts_C"], p["Q10_OM"])
        * oxygen_availability_h
        + 2.0 * (float(r["aemo_s"]) + float(r["aemo_h"]))
        + 2.0 * float(r["nitrification"])
    )
    renewal_o2 = (o2_sat - v["O2e_uM"]) * max(depth_dot, 0.0)
    dy[IDX["O2e_uM"]] = (
        reaeration
        + s.productivity_mmolO2_m2_d
        - respiration_epi
        - 2.0 * float(r["aemo_e"])
        - float(r["o2_interlayer"])
        + renewal_o2
    ) / epi_depth
    dy[IDX["O2h_uM"]] = (
        float(r["o2_interlayer"]) - respiration_hypo - oxygen_demand
    ) / hypo_depth

    delivery_ratio = max(exchange / p["interlayer_exchange_mixed_m_d"], 0.01)
    acceptor_exchange = p["acceptor_exchange_d_1"] * delivery_ratio
    dilution_no3 = (s.NO3_boundary_uM - v["NO3h_uM"]) * max(depth_dot, 0.0) / depth
    dilution_so4 = (s.SO4_boundary_uM - v["SO4h_uM"]) * max(depth_dot, 0.0) / depth
    dy[IDX["NO3h_uM"]] = (
        acceptor_exchange * (s.NO3_boundary_uM - v["NO3h_uM"])
        + dilution_no3
        - (float(r["denitrification"]) + p["NO3_per_CH4_AOM"] * float(r["aom_no3"])) / hypo_depth
    )
    dy[IDX["SO4h_uM"]] = (
        acceptor_exchange * (s.SO4_boundary_uM - v["SO4h_uM"])
        + dilution_so4
        - p["SO4_per_CH4_AOM"] * float(r["aom_so4"]) / hypo_depth
    )
    dy[IDX["FeIII_mmol_m2"]] = -p["Fe_per_CH4_AOM"] * float(r["aom_fe"]) + p["Fe_regeneration_mmol_m2_d"]

    dy[IDX["CH4s_uM"]] = (
        float(r["methanogenesis"])
        - float(r["aemo_s"])
        - float(r["aom_no3"])
        - float(r["aom_so4"])
        - float(r["aom_fe"])
        - float(r["ch4_swi"])
        - float(r["bubble_formation"])
    ) / pore_inventory_depth
    dy[IDX["CH4h_uM"]] = (
        float(r["ch4_swi"])
        + float(r["bubble_to_water"])
        - float(r["aemo_h"])
        - float(r["ch4_interlayer"])
    ) / hypo_depth
    dy[IDX["CH4e_uM"]] = (
        float(r["ch4_interlayer"])
        - float(r["aemo_e"])
        - float(r["ch4_diff_flux"])
        + (p["CH4_equilibrium_uM"] - v["CH4e_uM"]) * max(depth_dot, 0.0)
    ) / epi_depth
    dy[IDX["bubble_store_mmol_m2"]] = float(r["bubble_formation"]) - float(r["bubble_release"])

    co2_production = (
        float(r["mineralization"])
        + float(r["aemo_s"])
        + float(r["aom_no3"])
        + float(r["aom_so4"])
        + float(r["aom_fe"])
        + 0.5 * float(r["methanogenesis"])
    )
    dy[IDX["CO2s_uM"]] = (co2_production - float(r["co2_swi"])) / pore_inventory_depth
    water_co2_epi = p["water_CO2_production_uM_d"] * q10(v["Te_C"], p["Q10_OM"]) * epi_depth
    water_co2_hypo = p["water_CO2_production_uM_d"] * q10(v["Th_C"], p["Q10_OM"]) * hypo_depth
    productivity_co2_uptake = s.productivity_mmolO2_m2_d * p["photosynthetic_CO2_per_O2"]
    dy[IDX["CO2h_uM"]] = (
        float(r["co2_swi"])
        + float(r["aemo_h"])
        + water_co2_hypo
        - float(r["co2_interlayer"])
    ) / hypo_depth
    dy[IDX["CO2e_uM"]] = (
        float(r["co2_interlayer"])
        + float(r["aemo_e"])
        + water_co2_epi
        - productivity_co2_uptake
        - float(r["co2_flux"])
        + (p["CO2_equilibrium_uM"] - v["CO2e_uM"]) * max(depth_dot, 0.0)
    ) / epi_depth

    n2o_net_sed = float(r["nitrification_n2o"]) + float(r["denitrification_n2o"]) - float(r["n2o_reduction"])
    dy[IDX["N2Os_uM"]] = (n2o_net_sed - float(r["n2o_swi"])) / pore_inventory_depth
    dy[IDX["N2Oh_uM"]] = (float(r["n2o_swi"]) - float(r["n2o_interlayer"])) / hypo_depth
    dy[IDX["N2Oe_uM"]] = (
        float(r["n2o_interlayer"])
        - float(r["n2o_flux"])
        + (p["N2O_equilibrium_uM"] - v["N2Oe_uM"]) * max(depth_dot, 0.0)
    ) / epi_depth

    dy[IDX["cum_CO2_flux_mmol_m2"]] = float(r["co2_flux"])
    dy[IDX["cum_CH4_diff_flux_mmol_m2"]] = float(r["ch4_diff_flux"])
    dy[IDX["cum_CH4_ebu_flux_mmol_m2"]] = float(r["bubble_to_air"])
    dy[IDX["cum_N2O_flux_mmol_m2"]] = float(r["n2o_flux"])
    dy[IDX["cum_methanogenesis_mmol_m2"]] = float(r["methanogenesis"])
    dy[IDX["cum_AeMO_mmol_m2"]] = float(r["aemo_s"]) + float(r["aemo_w"])
    dy[IDX["cum_AOM_NO3_mmol_m2"]] = float(r["aom_no3"])
    dy[IDX["cum_AOM_SO4_mmol_m2"]] = float(r["aom_so4"])
    dy[IDX["cum_AOM_Fe_mmol_m2"]] = float(r["aom_fe"])
    dy[IDX["cum_nitrification_N2O_mmol_m2"]] = float(r["nitrification_n2o"])
    dy[IDX["cum_denitrification_N2O_mmol_m2"]] = float(r["denitrification_n2o"])
    dy[IDX["cum_N2O_reduction_mmol_m2"]] = float(r["n2o_reduction"])

    for name in NONNEGATIVE_STATE_NAMES:
        i = IDX[name]
        if y[i] <= 0.0 and dy[i] < 0.0:
            dy[i] = 0.0
    stratification_index = IDX["stratification_index"]
    if y[stratification_index] <= 0.0 and dy[stratification_index] < 0.0:
        dy[stratification_index] = 0.0
    if y[stratification_index] >= 1.0 and dy[stratification_index] > 0.0:
        dy[stratification_index] = 0.0
    return dy


def initial_state(s: Scenario, p: Dict[str, float]) -> np.ndarray:
    y0 = np.zeros(len(STATE_NAMES), dtype=float)
    initial = {
        "Te_C": p["initial_temperature_C"],
        "Th_C": p["initial_temperature_C"],
        "Ts_C": p["initial_temperature_C"],
        "O2e_uM": p["initial_O2_uM"],
        "O2h_uM": p["initial_O2_uM"],
        "NO3h_uM": s.NO3_boundary_uM,
        "SO4h_uM": s.SO4_boundary_uM,
        "FeIII_mmol_m2": p["initial_FeIII_mmol_m2"],
        "CH4e_uM": p["initial_CH4w_uM"],
        "CH4h_uM": p["initial_CH4w_uM"],
        "CH4s_uM": p["initial_CH4s_uM"],
        "CO2e_uM": p["initial_CO2w_uM"],
        "CO2h_uM": p["initial_CO2w_uM"],
        "CO2s_uM": p["initial_CO2s_uM"],
        "N2Oe_uM": p["initial_N2Ow_uM"],
        "N2Oh_uM": p["initial_N2Ow_uM"],
        "N2Os_uM": p["initial_N2Os_uM"],
        "bubble_store_mmol_m2": p["initial_bubble_store_mmol_m2"],
        "stratification_index": p["initial_stratification_index"],
    }
    for name, value in initial.items():
        y0[IDX[name]] = value
    return y0


def simulate(s: Scenario, p: Dict[str, float], output_step_h: float = 1.0) -> pd.DataFrame:
    t_eval = np.arange(0.0, s.duration_d + 1e-9, output_step_h / 24.0)
    solution = solve_ivp(
        lambda t, y: model_rhs(t, y, s, p),
        (0.0, s.duration_d),
        initial_state(s, p),
        t_eval=t_eval,
        method="LSODA",
        rtol=1e-6,
        atol=1e-8,
        max_step=1.0 / 24.0,
    )
    if not solution.success:
        raise RuntimeError(f"{s.scenario_id}: {solution.message}")
    frame = pd.DataFrame(solution.y.T, columns=STATE_NAMES)
    frame.insert(0, "time_d", solution.t)
    frame.insert(0, "scenario_id", s.scenario_id)
    rates = [process_rates(t, y, s, p) for t, y in zip(solution.t, solution.y.T)]
    frame["depth_m"] = [r["depth"] for r in rates]
    frame["epilimnion_depth_m"] = [r["epi_depth"] for r in rates]
    frame["hypolimnion_depth_m"] = [r["hypo_depth"] for r in rates]
    frame["surface_temperature_C"] = [r["surface_temperature"] for r in rates]
    frame["mixing_multiplier"] = [r["mixing_multiplier"] for r in rates]
    frame["interlayer_exchange_m_d"] = [r["interlayer_exchange"] for r in rates]
    frame["ventilation_intensity"] = [r["ventilation_intensity"] for r in rates]
    frame["hydrodynamic_state"] = [r["hydrodynamic_state"] for r in rates]
    for key in [
        "methanogenesis",
        "aemo_s",
        "aemo_e",
        "aemo_h",
        "aemo_w",
        "aom_no3",
        "aom_so4",
        "aom_fe",
        "nitrification_n2o",
        "denitrification_n2o",
        "n2o_reduction",
        "bubble_to_air",
        "co2_flux",
        "ch4_diff_flux",
        "n2o_flux",
        "ch4_interlayer",
        "co2_interlayer",
        "n2o_interlayer",
        "o2_interlayer",
    ]:
        frame[key] = [r[key] for r in rates]
    frame[NONNEGATIVE_STATE_NAMES] = frame[NONNEGATIVE_STATE_NAMES].clip(lower=0.0)
    frame["stratification_index"] = frame["stratification_index"].clip(lower=0.0, upper=1.0)
    epi = frame["epilimnion_depth_m"]
    hypo = frame["hypolimnion_depth_m"]
    depth = frame["depth_m"]
    frame["Tw_C"] = (frame["Te_C"] * epi + frame["Th_C"] * hypo) / depth
    frame["O2_uM"] = frame["O2h_uM"]
    frame["CH4w_uM"] = (frame["CH4e_uM"] * epi + frame["CH4h_uM"] * hypo) / depth
    frame["CO2w_uM"] = (frame["CO2e_uM"] * epi + frame["CO2h_uM"] * hypo) / depth
    frame["N2Ow_uM"] = (frame["N2Oe_uM"] * epi + frame["N2Oh_uM"] * hypo) / depth
    return frame


def summarize(frame: pd.DataFrame, s: Scenario) -> Dict[str, float | str]:
    tail = frame[frame["time_d"] >= max(0.0, s.duration_d - 5.0)]
    last = frame.iloc[-1]
    return {
        "scenario_id": s.scenario_id,
        "forcing": s.forcing,
        "hydrodynamic_regime": s.hydrodynamic_regime,
        "nominal_depth_m": s.depth_m,
        "minimum_depth_m": frame["depth_m"].min(),
        "maximum_depth_m": frame["depth_m"].max(),
        "sediment_temperature_amplitude_C": 0.5 * (tail["Ts_C"].max() - tail["Ts_C"].min()),
        "mean_O2_uM": frame["O2h_uM"].mean(),
        "minimum_O2_uM": frame["O2h_uM"].min(),
        "anoxic_hours": 24.0 * np.trapezoid(
            (frame["O2h_uM"] < 10.0).astype(float), frame["time_d"]
        ),
        "mean_stratification_index": frame["stratification_index"].mean(),
        "stratified_hours": 24.0 * np.trapezoid(
            (frame["hydrodynamic_state"] == "stratified").astype(float), frame["time_d"]
        ),
        "ventilating_hours": 24.0 * np.trapezoid(
            (frame["hydrodynamic_state"] == "ventilating").astype(float), frame["time_d"]
        ),
        "minimum_interlayer_exchange_m_d": frame["interlayer_exchange_m_d"].min(),
        "maximum_interlayer_exchange_m_d": frame["interlayer_exchange_m_d"].max(),
        "maximum_hypolimnetic_CH4_uM": frame["CH4h_uM"].max(),
        "final_hypolimnetic_CH4_uM": last["CH4h_uM"],
        "CO2_flux_cumulative_mmol_m2": last["cum_CO2_flux_mmol_m2"],
        "CH4_diffusive_cumulative_mmol_m2": last["cum_CH4_diff_flux_mmol_m2"],
        "CH4_ebullitive_cumulative_mmol_m2": last["cum_CH4_ebu_flux_mmol_m2"],
        "N2O_flux_cumulative_mmol_m2": last["cum_N2O_flux_mmol_m2"],
        "gross_methanogenesis_mmol_m2": last["cum_methanogenesis_mmol_m2"],
        "gross_AeMO_mmol_m2": last["cum_AeMO_mmol_m2"],
        "gross_AOM_mmol_m2": last["cum_AOM_NO3_mmol_m2"] + last["cum_AOM_SO4_mmol_m2"] + last["cum_AOM_Fe_mmol_m2"],
        "AOM_NO3_mmol_m2": last["cum_AOM_NO3_mmol_m2"],
        "AOM_SO4_mmol_m2": last["cum_AOM_SO4_mmol_m2"],
        "AOM_Fe_mmol_m2": last["cum_AOM_Fe_mmol_m2"],
        "gross_N2O_production_mmol_m2": last["cum_nitrification_N2O_mmol_m2"] + last["cum_denitrification_N2O_mmol_m2"],
        "gross_N2O_reduction_mmol_m2": last["cum_N2O_reduction_mmol_m2"],
    }


def make_demo_plot(summary: pd.DataFrame, output_path: Path) -> None:
    static = summary[
        summary["scenario_id"].str.match(r"^D\d+_diurnal$")
    ].sort_values("nominal_depth_m")
    if static.empty:
        return
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.2))
    axes[0].plot(static["nominal_depth_m"], static["sediment_temperature_amplitude_C"], "o-", color="#D55E00")
    axes[0].set(ylabel="Bed-temperature amplitude (C)", xlabel="Water depth (m)")
    axes[1].plot(static["nominal_depth_m"], static["mean_O2_uM"], "o-", color="#0072B2")
    axes[1].set(ylabel="Mean DO (uM)", xlabel="Water depth (m)")
    total_ch4 = static["CH4_diffusive_cumulative_mmol_m2"] + static["CH4_ebullitive_cumulative_mmol_m2"]
    axes[2].plot(static["nominal_depth_m"], total_ch4, "o-", color="#CC79A7")
    axes[2].set(ylabel="Cumulative CH4 flux (mmol m-2)", xlabel="Water depth (m)")
    states = summary[summary["scenario_id"].isin(["D100_diurnal", "D100_persistent_state", "D100_ventilation_state"])]
    states = states.set_index("scenario_id").reindex(
        ["D100_diurnal", "D100_persistent_state", "D100_ventilation_state"]
    )
    axes[3].bar(
        ["Mixed", "Persistent", "Ventilated"],
        states["mean_stratification_index"],
        color=["#56B4E9", "#E69F00", "#009E73"],
    )
    axes[3].set(ylabel="Mean stratification index", xlabel="1-m hydrodynamic regime", ylim=(0, 1))
    axes[3].tick_params(axis="x", rotation=20)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#D9DDE0", lw=0.6)
    fig.suptitle("Illustrative output: parameters require local calibration", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_scenarios(
    scenarios: Iterable[Scenario],
    p: Dict[str, float],
    output_dir: Path,
    make_plot: bool = True,
) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    summaries = []
    for scenario in scenarios:
        frame = simulate(scenario, p)
        frames.append(frame)
        summaries.append(summarize(frame, scenario))
    all_timeseries = pd.concat(frames, ignore_index=True)
    summary = pd.DataFrame(summaries)
    all_timeseries.to_csv(output_dir / "timeseries.csv", index=False)
    summary.to_csv(output_dir / "summary.csv", index=False)
    if make_plot:
        make_demo_plot(summary, output_dir / "illustrative_depth_response.png")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=ROOT / "parameters.csv")
    parser.add_argument("--scenarios", type=Path, default=ROOT / "scenarios.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    summary = run_scenarios(
        read_scenarios(args.scenarios),
        read_parameters(args.parameters),
        args.output_dir,
        make_plot=not args.no_plot,
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
