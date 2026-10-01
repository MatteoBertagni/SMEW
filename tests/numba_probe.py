"""Exercise numerical functions in a fresh process with the selected JIT setting."""

import importlib.abc
import inspect
import os
from pathlib import Path
import sys

import numpy as np


def collect_outputs(*, helpers_only=False):
    disabled = os.environ.get("NUMBA_DISABLE_JIT", "0") != "0"
    if disabled or helpers_only:
        class BlockNative(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                if fullname.startswith("smew._native"):
                    raise AssertionError(f"Unexpected native import: {fullname}")
        sys.meta_path.insert(0, BlockNative())

    import smew
    from smew import equations, vegetation, soil_pores
    from smew.biogeochem import _biogeochem_balance
    from smew.biogeochem2psd import _biogeochem_balance2psd
    from smew.ic import _conc_to_f_CEC, _total_to_f_CEC_and_conc, _Kelland

    outputs = {}

    def record(name, value):
        if isinstance(value, (tuple, list)):
            for i, item in enumerate(value):
                record(f"{name}/{i}", item)
        else:
            outputs[name] = np.asarray(value)

    def call(fn, *args, **kwargs):
        result = fn(*args, **kwargs)
        record(fn.__name__, result)
        return result

    soil_names = (
        "sand", "loamy sand", "sandy loam", "loam", "silt", "silt loam",
        "sandy clay loam", "silty clay loam", "clay loam", "sandy clay", "silty clay", "clay",
    )
    for soil in soil_names:
        record("soil/" + soil, smew.soil_const(soil))
        for model in ("ding2016", "campbell"):
            record(f"pores/{soil}/{model}", smew.soil_pore_pdf(soil, model, 50))
    call(soil_pores._gradient, np.array([1., 2., 8., 10.]), np.array([0., .2, 1., 2.]))
    for fn in (smew.D_0, smew.Dw_0, smew.plant_nutr_f):
        call(fn)
    for fn in (smew.CO2_atm, smew.MM, smew.K_Al, smew.carb_weath_const):
        call(fn, 1.)
    call(smew.K_C, np.array([273.15, 293.15, 313.15]), 1.)
    call(smew.K_GT_CEC, "loam", 1.)
    for mineral in ("albite", "anorthite", "analcime", "forsterite", "wollastonite", "diopside", "muscovite", "labradorite", "augite", "alkali_feldspar", "Fe_forsterite", "nepheline", "apatite", "leucite"):
        record("mineral/" + mineral, smew.min_const(mineral, 1.))
        record("omega/" + mineral, smew.sil_Omega(mineral, .001, .001, .001, .001, 1e-9, 1e-9, 1e-5, 1e-6, 1., 1., 1.))
    call(smew.sil_Wr, "forsterite", .5, 1e-6, .001, .002, .003, .5, .2, 1., 1.)
    call(smew.carb_W, 1., 1., .5, 1.5, .6, .3, 1., 1., 10., 10.)
    call(smew.temp, .8, 15., 10., 5., .3, 8., 1., 1.)
    x = np.ones(8)
    call(smew.veg, 1., 20., 10., 2., np.zeros(100), .5)
    for i, stages in enumerate(([], [0, 0], [0, 1] + [2]*12 + [3, 4, 0, 0], [0, 1] + [2]*12 + [3, 4], ([0, 1] + [2]*12 + [3, 4, 0, 0])*2)):
        stages = np.array(stages, dtype=np.int8)
        record(f"season/{i}", smew.veg_seasonal(stages, 10., .1))
        record(f"boundaries/{i}/stage", vegetation.get_stage_boundaries(stages, 2))
        record(f"boundaries/{i}/season", vegetation.get_season_boundaries(stages))
    # Demand-limited, diffusion-limited, passive-only and zero-concentration uptake.
    call(smew.up_act, .5, .01, np.array([.0001, .01, 1., 1.]), 1., .001,
         .001, .001, .001, 0., .001, 1., 1., 1., .1)
    temperatures = np.array([-2., -1., 2., 3., 8., 15., 20., 15.])
    for mode in (0, 1):
        for frozen in (False, True):
            record(f"moisture/{mode}/{frozen}", smew.moisture_balance(
                x*.001, .3, "loam", x*.002, x, 1., mode, .6, 8., 1.,
                temp_soil=temperatures if frozen else None,
            ))
    for label, soc, co2, tau, litter in (("tau", 1000., None, 1000., None), ("co2", 1000., .01, None, None), ("estimate_soc", None, None, 1000., .1)):
        record("respiration/" + label, smew.respiration(
            litter, soc, co2, 1., "loam", x*.6, x, 1., .3, x*15., 1., 1., tau,
        ))
    try:
        smew.respiration(None, 1000., None, 1., "loam", x*.6, x, 1., .3, -x, 1., 1., 1000.)
    except ValueError as exc:
        assert "Mean decomposition activity is zero" in str(exc)
    else:
        raise AssertionError("Frozen respiration should reject a steady-state estimate")
    call(smew.mov_avg, np.arange(10.), 2)
    d = np.array([50e-6, 100e-6, 200e-6])
    widths = np.array([50e-6, 50e-6, 100e-6])
    call(smew.psd_evol, d, widths, d, widths, np.ones(3), 3, 100., .35, 3e6)
    number = call(smew.psd_number_from_mass, np.array([.2, .5, .3]), d, 3e6)
    pore_d, pore_pdf = smew.soil_pore_pdf("loam", "ding2016", 50)
    call(smew.normalized_cumulative_area, np.array([4., 1., 2., 2., 3.]), np.array([2., -1., 1., 7., 3.]))
    call(smew.wet_f_Anand, pore_d, pore_pdf, .6, d, number)
    for mode in ("constant", "linear", "nonlinear"):
        record("wetness/" + mode, smew.wetness_SA(.6, mode, pore_d, pore_pdf, d, number, 1., 1e-4))
    # Numba's random generator is independent of Python's; check structural properties.
    for rain in (smew.rain_stoc(.2, .01, 365., 1.), smew.rain_stoc_season(np.full(12, .2), np.full(12, .01), 365., 1.)):
        assert rain.shape == (365,) and np.isfinite(rain).all() and (rain >= 0.).all()
    fractions = np.array([.6, .2, .08, .05, .04, .03])
    conc, _ = call(smew.f_CEC_to_conc, fractions, 6., "loam", 1., 1.)
    call(smew.f_CEC_and_conc_to_K, fractions, np.asarray(conc), 6., "loam", 1., 1.)
    call(smew.Amann, fractions, 6., .001, "loam", 1., 1.)
    # Exercise every equation directly as well as through cminpack.
    for name, fn in inspect.getmembers(equations):
        if name.endswith(("_equations", "_equation", "_residual")):
            size = 16 if name.startswith("biogeochem") else 13 if name.startswith("total_to_cec") else 5 if name.startswith("kelland") else 1
            out = np.zeros(size)
            args = [np.ones(size) if param == "p" else out if param == "out" else 1.
                    for param in inspect.signature(fn).parameters]
            record(name, fn(*args))
            if "out" in inspect.signature(fn).parameters:
                record(name + "/out", out)
    if not helpers_only:
        call(smew.conc_to_f_CEC, conc, 6., "loam", 1., 1.)
        totals = [conc[i] * (.4*.3*.6*1000) + fractions[i]*.1/charge
                  for i, charge in enumerate((2, 2, 1, 1))]
        call(smew.total_to_f_CEC_and_conc, totals, 6., .07, [.6], "loam", .4, .3, .1, 1., 1.)
        totals[0] += .01
        call(smew.Kelland, totals, 6., conc, [.6], "loam", .4, .3, .1, 1., 1.)
        from tests.example_marimo_notebook import app
        original = smew.biogeochem_balance
        captured = {}
        freeze = False

        def balance(**inputs):
            if freeze:
                inputs["temp_soil"] = inputs["temp_soil"].copy()
                inputs["temp_soil"][24:] = -2.
            captured.update(inputs)
            return original(**inputs)

        smew.biogeochem_balance = balance
        for variant in ("default", "no_rock", "nonlinear_frozen"):
            settings = dict(duration_days=2, timestep_minutes=60)
            freeze = variant == "nonlinear_frozen"
            if variant == "no_rock":
                settings.update(rock_mass_g_m2=0., application_day=0.)
            elif freeze:
                settings.update(
                    dissolution_factor=1., mineral_mass_fractions=[1.], minerals=["forsterite"],
                    particle_diameters_um=[50., 100., 200.], particle_mass_fractions=[.2, .5, .3],
                    pore_model="ding2016", pore_particle_mixing=1., rock_density_g_m3=3e6,
                    wet_surface_model="nonlinear", initial_vegetation_g_m2=1000.,
                    vegetation_capacity_g_m2=3000, vegetation_growth_days=100,
                    vegetation_start_day=0, root_area_index=10, root_diameter_m=.4e-3,
                )
            _, definitions = app.run(defs=settings)
            for name, value in definitions["results"].items():
                record(f"simulation/{variant}/{name}", value)
            for name in ("UP_Ca", "UP_Mg", "UP_K", "UP_Si", "wet_f", "frozen"):
                record(f"simulation/{variant}/{name}", definitions["chemistry"][name])
            if variant == "default":
                two_inputs = {name: captured[name] for name in inspect.signature(smew.biogeochem_balance2psd).parameters if name in captured}
        smew.biogeochem_balance = original
        two_inputs.update(M_rock_in=100., t_app=0., M_rock_in2=50., t_app2=1.,
                          mineral2=two_inputs["mineral"], rock_f_in2=two_inputs["rock_f_in"],
                          d_in2=two_inputs["d_in"], psd_perc_in2=two_inputs["psd_perc_in"], SSA_in2=two_inputs["SSA_in"])
        for label, first, second in (("both", 100., 50.), ("none", 0., 0.)):
            two_inputs.update(M_rock_in=first, M_rock_in2=second)
            result = smew.biogeochem_balance2psd(**two_inputs)
            for name in ("pH", "Ca", "Mg", "IC_tot", "M_rock", "M_rock2", "SA", "SA2"):
                record(f"two_rocks/{label}/{name}", result[name])
    # A successful call must produce native signatures, not an object-mode fallback.
    numerical = [getattr(smew, name) for name in (
        "CO2_atm", "D_0", "Dw_0", "MM", "K_Al", "K_C", "K_GT_CEC", "plant_nutr_f",
        "soil_hydraulic_const", "soil_const", "min_const", "carb_weath_const", "temp",
        "rain_stoc", "rain_stoc_season", "moisture_balance", "respiration", "soil_pore_pdf",
        "pore_pdf_ding2016", "pore_pdf_campbell", "veg", "veg_seasonal", "up_act",
        "carb_W", "sil_Omega", "sil_Wr", "psd_evol", "psd_number_from_mass",
        "normalized_cumulative_area", "wet_f_Anand", "wetness_SA", "mov_avg",
        "f_CEC_to_conc", "f_CEC_and_conc_to_K", "Amann",
    )] + [vegetation.get_stage_boundaries, vegetation.get_season_boundaries]
    if not helpers_only:
        numerical += [_biogeochem_balance, _biogeochem_balance2psd, _conc_to_f_CEC, _total_to_f_CEC_and_conc, _Kelland]
    for fn in numerical:
        if disabled:
            assert not hasattr(fn, "py_func"), fn.__name__
        else:
            assert fn.nopython_signatures, fn.__name__
    return outputs


if __name__ == "__main__":
    np.savez(Path(sys.argv[1]), **collect_outputs(helpers_only="--helpers-only" in sys.argv[2:]))
