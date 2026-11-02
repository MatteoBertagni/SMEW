"""Compatibility checks for the deprecated two-application model."""

from pathlib import Path

import numpy as np
import pytest

import smew
from smew.equations import water_alk_equations


@pytest.mark.parametrize("rock_mass", (0., 1.))
def test_two_application_model_warns_and_preserves_rate_flux_contract(rock_mass):
    x = np.ones(3)
    fractions = np.array((.6, .2, .08, .05, .04, .03))
    concentrations, constants = smew.f_CEC_to_conc(
        fractions, 6., "loam", 1., 1.,
    )
    with pytest.warns(DeprecationWarning, match="biogeochem_balance2psd is deprecated") as caught:
        result = smew.biogeochem_balance2psd(
            n=.4, s=x*.6, L=x*0, T=x*0, I=x*0, v=x*0,
            k_v=1., RAI=1., root_d=.001, Zr=.3,
            r_het=x*.001, r_aut=x*.001, D=x*.01, temp_soil=x*20.,
            pH_in=6., conc_in=concentrations, f_CEC_in=fractions,
            K_CEC=constants, CEC_tot=.1, Si_in=.001,
            CaCO3_in=0., MgCO3_in=0.,
            M_rock_in=rock_mass, t_app=0., mineral=("forsterite",),
            rock_f_in=np.array([1.]), d_in=np.array([1e-4]),
            psd_perc_in=np.array([1.]), SSA_in=np.nan,
            M_rock_in2=rock_mass, t_app2=0., mineral2=("forsterite",),
            rock_f_in2=np.array([1.]), d_in2=np.array([1e-4]),
            psd_perc_in2=np.array([1.]), SSA_in2=np.nan,
            diss_f=1., dt=1/24., conv_Al=1., conv_mol=1., keyword_add=0,
        )

    assert Path(caught[0].filename).resolve() == Path(__file__).resolve()
    assert np.isfinite(result["pH"]).all()
    assert (result["H_rain"] > 0.).all()
    for i, hydrogen in enumerate(result["H_rain"]):
        residual = water_alk_equations(
            np.array([hydrogen]), 0., result["k1"][i], result["k2"][i],
            result["CO2_w_rain"][i], result["k_w"][i],
        )
        assert abs(residual) < 1e-12
    # The last timestep also solves the same balance at the soil's trial Alk0.
    assert result["H0_2"] > 0.
    residual = water_alk_equations(
        np.array([result["H0_2"]]), result["Alk0"], result["k1"][-1],
        result["k2"][-1], result["CO2_w0"], result["k_w"][-1],
    )
    assert abs(residual) < 1e-12
    assert result["M_rock"][0] == rock_mass
    assert result["M_rock2"][0] == rock_mass

    if rock_mass > 0:
        for suffix in ("", "2"):
            # Both distributions use the chemical rate directly, including at
            # initialization. Moisture scales the flux, as in biogeochem.
            rate = np.array([
                smew.sil_Wr(
                    "forsterite", result["Omega" + suffix][0, i], result["H"][i],
                    result["k_H_T" + suffix][0, i], result["k_w_T" + suffix][0, i],
                    result["k_OH_T" + suffix][0, i], result["n_H" + suffix][0],
                    result["n_OH" + suffix][0], 1., 1.,
                )
                for i in range(x.size)
            ])
            assert (rate > 0.).all()
            np.testing.assert_allclose(result["Wr" + suffix][0], rate, rtol=1e-14, atol=0)
            wet_area = result["SA" + suffix] * result["s"]
            np.testing.assert_allclose(
                result["EW" + suffix][0],
                rate * wet_area * result["rock_f" + suffix][0],
                rtol=1e-14, atol=0,
            )
