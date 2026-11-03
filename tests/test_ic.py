"""CEC initialization produces scalar fractions usable by the balance model."""

import numpy as np
import pytest

import smew


@pytest.fixture
def balance_inputs():
    x = np.ones(3)
    fractions = np.array((.6, .2, .08, .05, .04, .03))
    concentrations, constants = smew.f_CEC_to_conc(fractions, 6., "loam", 1., 1.)
    return dict(
        n=.4, s=x*.6, L=x*0, T=x*0, I=x*0, v=x*0,
        k_v=1., RAI=1., root_d=.001, Zr=.3,
        r_het=x*.001, r_aut=x*.001, D=x*.01, temp_soil=x*20.,
        pH_in=6., conc_in=concentrations, f_CEC_in=fractions,
        K_CEC=constants, CEC_tot=.1, Si_in=.001,
        CaCO3_in=0., MgCO3_in=0., M_rock_in=0., t_app=0.,
        mineral=("forsterite",), rock_f_in=np.array([1.]),
        d_in=np.array([1e-4]), psd_perc_in=np.array([1.]), SSA_in=np.nan,
        diss_f=1., dt=1/24., conv_Al=1., conv_mol=1., keyword_add=0,
    )


def test_concentration_to_cec_round_trip_initializes_balance(balance_inputs):
    fractions, constants = smew.conc_to_f_CEC(
        balance_inputs["conc_in"], 6., "loam", 1., 1.,
    )
    assert np.asarray(fractions).shape == (6,)
    np.testing.assert_allclose(fractions, balance_inputs["f_CEC_in"], rtol=1e-9)
    np.testing.assert_allclose(constants, balance_inputs["K_CEC"])

    result = smew.biogeochem_balance(
        **{**balance_inputs, "f_CEC_in": fractions, "K_CEC": constants},
    )
    assert np.isfinite(result["pH"]).all()
    for species, fraction in zip(("Ca", "Mg", "K", "Na", "Al", "H"), fractions):
        assert result["f_" + species][0] == pytest.approx(fraction)


@pytest.mark.parametrize("fractions", (
    np.ones((6, 1)), np.ones((1, 6)), np.ones((2, 3)),
    np.ones(5), np.ones(7), np.array(.5), None,
))
def test_balance_rejects_invalid_cec_shapes(balance_inputs, fractions):
    with pytest.raises(ValueError, match="f_CEC_in must be a one-dimensional array with exactly six values"):
        smew.biogeochem_balance(**{**balance_inputs, "f_CEC_in": fractions})
