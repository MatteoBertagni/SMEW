"""CEC-supported soil textures must also work with hydrology and pores."""

import numpy as np
import pytest

import smew


SOIL_GROUPS = (
    ("sand", "sand"),
    ("loamy sand", "sand"),
    ("sandy loam", "sand"),
    ("loam", "loam"),
    ("silt loam", "loam"),
    ("silt", "loam"),
    ("clay loam", "clay"),
    ("silty clay", "clay"),
    ("clay", "clay"),
)


@pytest.mark.parametrize("soil,cec_group", SOIL_GROUPS)
@pytest.mark.parametrize("conv_mol", (1., 1e6))
def test_cec_supported_soil_texture_works_across_model(soil, cec_group, conv_mol):
    assert np.isfinite(smew.soil_const(soil)).all()
    for pore_model in ("ding2016", "campbell"):
        diameters, density = smew.soil_pore_pdf(soil, pore_model, n_points=50)
        assert np.isfinite(diameters).all()
        assert np.isfinite(density).all()

    fractions = np.array((.6, .2, .08, .05, .04, .03))
    concentrations, constants = smew.f_CEC_to_conc(
        fractions, 6., soil, conv_mol, 1.,
    )
    expected_concentrations, expected_constants = smew.f_CEC_to_conc(
        fractions, 6., cec_group, conv_mol, 1.,
    )
    assert np.isfinite(concentrations).all()
    assert (np.asarray(concentrations) > 0.).all()
    np.testing.assert_allclose(concentrations, expected_concentrations, rtol=1e-14)
    np.testing.assert_allclose(constants, expected_constants, rtol=1e-14)


def test_legacy_silty_loam_cec_alias_is_preserved():
    np.testing.assert_array_equal(
        smew.K_GT_CEC("silty loam", 1.), smew.K_GT_CEC("silt loam", 1.),
    )


@pytest.mark.parametrize("soil", ("peat", "sandy clay loam", "silty clay loam", "sandy clay"))
def test_unsupported_cec_soil_reports_input_and_supported_names(soil):
    with pytest.raises(ValueError) as caught:
        smew.K_GT_CEC(soil, 1.)
    message = str(caught.value)
    assert f"Unknown soil type '{soil}'" in message
    assert "K_GT_CEC" in message
    for supported_soil, _ in SOIL_GROUPS:
        assert supported_soil in message

    fractions = np.array((.6, .2, .08, .05, .04, .03))
    with pytest.raises(ValueError, match=f"Unknown soil type '{soil}'"):
        smew.f_CEC_to_conc(fractions, 6., soil, 1., 1.)
