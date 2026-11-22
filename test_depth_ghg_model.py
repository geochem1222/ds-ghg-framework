import unittest

import numpy as np

from depth_ghg_model import (
    IDX,
    NONNEGATIVE_STATE_NAMES,
    Scenario,
    initial_state,
    process_rates,
    read_parameters,
    simulate,
    summarize,
)


class DepthGHGModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = read_parameters()

    def scenario(
        self,
        scenario_id,
        forcing,
        depth,
        days=20,
        no3=20,
        so4=100,
        regime="mixed",
        mixing=1.0,
    ):
        return Scenario(
            scenario_id,
            forcing,
            depth,
            days,
            2.0,
            mixing,
            0.0,
            no3,
            so4,
            hydrodynamic_regime=regime,
        )

    def test_deeper_water_attenuates_diurnal_bed_temperature(self):
        shallow = self.scenario("shallow", "diurnal", 0.05)
        deep = self.scenario("deep", "diurnal", 1.0)
        shallow_amp = summarize(simulate(shallow, self.p, 3), shallow)["sediment_temperature_amplitude_C"]
        deep_amp = summarize(simulate(deep, self.p, 3), deep)["sediment_temperature_amplitude_C"]
        self.assertGreater(shallow_amp, deep_amp)

    def test_nitrate_and_sulfate_aom_can_operate_under_anoxia(self):
        s = self.scenario("anoxic", "static", 0.3, no3=300, so4=2000)
        y = initial_state(s, self.p)
        y[IDX["O2h_uM"]] = 0.0
        y[IDX["CH4s_uM"]] = 1000.0
        rates = process_rates(0.0, y, s, self.p)
        self.assertGreater(rates["aom_no3"], 0.0)
        self.assertGreater(rates["aom_so4"], 0.0)

    def test_warming_increases_intrinsic_methanogenesis(self):
        s = self.scenario("kinetics", "static", 0.3, no3=0, so4=0)
        cold = initial_state(s, self.p)
        warm = cold.copy()
        cold[IDX["Ts_C"]] = 10.0
        warm[IDX["Ts_C"]] = 25.0
        cold[IDX["O2h_uM"]] = warm[IDX["O2h_uM"]] = 0.0
        self.assertGreater(
            process_rates(0.0, warm, s, self.p)["methanogenesis"],
            process_rates(0.0, cold, s, self.p)["methanogenesis"],
        )

    def test_drawdown_releases_more_bubble_methane_than_static_depth(self):
        static = self.scenario("static", "diurnal", 0.6, 20)
        drawdown = self.scenario("drawdown", "drawdown", 0.6, 20)
        static_total = summarize(simulate(static, self.p, 3), static)["CH4_ebullitive_cumulative_mmol_m2"]
        drawdown_total = summarize(simulate(drawdown, self.p, 3), drawdown)["CH4_ebullitive_cumulative_mmol_m2"]
        self.assertGreater(drawdown_total, static_total)

    def test_states_are_finite_and_nonnegative(self):
        s = self.scenario("rise_fall", "rise_fall", 0.6, 20)
        frame = simulate(s, self.p, 3)
        state = frame[NONNEGATIVE_STATE_NAMES].to_numpy()
        self.assertTrue(np.isfinite(state).all())
        self.assertGreaterEqual(state.min(), 0.0)
        self.assertTrue(frame["stratification_index"].between(0.0, 1.0).all())

    def test_layer_depths_sum_to_total_depth(self):
        s = self.scenario("geometry", "rise_fall", 0.6, 20)
        frame = simulate(s, self.p, 3)
        reconstructed = frame["epilimnion_depth_m"] + frame["hypolimnion_depth_m"]
        self.assertTrue(np.allclose(reconstructed, frame["depth_m"], atol=1e-12))

    def test_persistent_regime_develops_stratification(self):
        mixed = self.scenario("mixed", "diurnal", 1.0, 30, regime="mixed", mixing=0.35)
        persistent = self.scenario(
            "persistent", "diurnal", 1.0, 30, regime="persistent", mixing=0.35
        )
        mixed_summary = summarize(simulate(mixed, self.p, 3), mixed)
        persistent_summary = summarize(simulate(persistent, self.p, 3), persistent)
        self.assertGreater(
            persistent_summary["mean_stratification_index"],
            mixed_summary["mean_stratification_index"],
        )
        self.assertGreater(persistent_summary["stratified_hours"], 0.0)

    def test_ventilation_event_temporarily_restores_exchange(self):
        s = self.scenario(
            "ventilation", "diurnal", 1.0, 30, regime="ventilation", mixing=0.35
        )
        frame = simulate(s, self.p, 1)
        event = frame[frame["ventilation_intensity"] > 0.5]
        pre_event = frame[(frame["time_d"] >= 10.0) & (frame["time_d"] <= 12.0)]
        self.assertFalse(event.empty)
        self.assertTrue((event["hydrodynamic_state"] == "ventilating").all())
        self.assertGreater(
            event["interlayer_exchange_m_d"].mean(),
            pre_event["interlayer_exchange_m_d"].mean(),
        )
        self.assertLess(event["stratification_index"].min(), pre_event["stratification_index"].mean())

    def test_gross_process_integrals_are_nonnegative(self):
        s = self.scenario("budget", "heatwave", 0.3, 20)
        last = simulate(s, self.p, 3).iloc[-1]
        names = [
            "cum_methanogenesis_mmol_m2",
            "cum_AeMO_mmol_m2",
            "cum_AOM_NO3_mmol_m2",
            "cum_AOM_SO4_mmol_m2",
            "cum_AOM_Fe_mmol_m2",
            "cum_nitrification_N2O_mmol_m2",
            "cum_denitrification_N2O_mmol_m2",
            "cum_N2O_reduction_mmol_m2",
        ]
        self.assertTrue((last[names] >= 0.0).all())


if __name__ == "__main__":
    unittest.main()
