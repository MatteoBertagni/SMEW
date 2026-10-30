# SPDX-License-Identifier: AGPL-3.0-only
# -*- coding: utf-8 -*-

import numpy as np
from numba import njit
from .constants import soil_hydraulic_const

#------------------------------------------------------------------------------

# Constants for the double exponential (DE) model of soil pores
# Ding et al. (2016) http://dx.doi.org/10.1016/j.still.2015.10.007 

ding2016_param = {
    "clay": dict(C=0.288, A1=0.121, h1=5447.5, A2=0.079, h2=269.0),
    "silty clay": dict(C=0.192, A1=0.140, h1=4031.7, A2=0.201, h2=269.9),
    "silty clay loam": dict(C=0.258, A1=0.092, h1=7883.8, A2=0.081, h2=761.6),
    "clay loam": dict(C=0.224, A1=0.154, h1=2679.4, A2=0.209, h2=55.7),
    "sandy clay loam": dict(C=0.205, A1=0.143, h1=24893.3, A2=0.116, h2=248.5),
    "silt loam": dict(C=0.130, A1=0.155, h1=9683.5, A2=0.122, h2=293.2),
    "loam": dict(C=0.170, A1=0.174, h1=5621.3, A2=0.164, h2=145.5),
    "silt": dict(C=0.099, A1=0.156, h1=919.3, A2=0.156, h2=919.0),
    "sandy loam": dict(C=0.089, A1=0.143, h1=5445.5, A2=0.190, h2=240.1),
    "sand": dict(C=0.047, A1=0.050, h1=960.2, A2=0.346, h2=32.4),
}

# Immutable lookup values are usable directly from Numba.
_DING2016_PARAMETERS = tuple(
    (soil, (p["C"], p["A1"], p["h1"], p["A2"], p["h2"]))
    for soil, p in ding2016_param.items()
)

MPA_TO_CM_H2O = 10197.16213

#------------------------------------------------------------------------------

@njit(nogil=True, error_model="numpy")
def _ding2016_parameters(soil):
    if soil == "loamy sand":
        soil = "sand"
    elif soil == "sandy clay":
        soil = "sandy clay loam"
    for name, parameters in _DING2016_PARAMETERS:
        if name == soil:
            return parameters
    raise ValueError("Invalid soil type!")


@njit(nogil=True, error_model="numpy")
def _gradient(values, coordinates):
    """NumPy's default 1D gradient on nonuniform coordinates."""
    if len(values) < 2:
        raise ValueError("At least two points are required for a pore distribution.")
    result = np.empty(len(values))
    result[0] = (values[1] - values[0]) / (coordinates[1] - coordinates[0])
    result[-1] = (values[-1] - values[-2]) / (coordinates[-1] - coordinates[-2])
    for i in range(1, len(values) - 1):
        left = coordinates[i] - coordinates[i - 1]
        right = coordinates[i + 1] - coordinates[i]
        result[i] = (
            -right / (left * (left + right)) * values[i - 1]
            + (right - left) / (left * right) * values[i]
            + left / (right * (left + right)) * values[i + 1]
        )
    return result


@njit(nogil=True, error_model="numpy")
def pore_pdf_ding2016(soil, n_points=500, h_min=10.0, h_max=10*MPA_TO_CM_H2O):
    
    """
    Pore-size distribution from Ding et al. (2016) double-exponential
    water retention parameters.
    
    The water retention curve is mapped to pore diameter with
    d = 2 * 0.149 / h, where h is in cm H2O and d is returned in m.
    
    The density is computed in log-diameter space, dF/dln(d), and then
    normalized over pore diameter. This gives an effective density [1/m]
    for constructing the pore CDF used by wet_f_Anand.
    
    h_min and h_max are suction limits [cm H2O]. Lower h gives larger pores;
    higher h gives smaller pores. The default h_max corresponds to SMEW's
    10 MPa hygroscopic point.
    
    Returns
    -------
    pore_d : ndarray
        Pore diameter [m].
    pore_pdf : ndarray
        Effective pore-size density [1/m].
    """
    
    C, A1, h1, A2, h2 = _ding2016_parameters(soil)

    h = np.logspace(np.log10(h_max), np.log10(h_min), n_points)

    theta = C + A1 * np.exp(-h / h1) + A2 * np.exp(-h / h2)

    pore_d = 2.0 * (0.149 / h) * 1e-2

    F = (theta - theta[0]) / (theta[-1] - theta[0])
        
    # Convention as for Anand implementation: distribution per log-diameter interval.
    pore_pdf = _gradient(F, np.log(pore_d))
    pore_pdf = np.maximum(pore_pdf, 0.0)

    area = np.trapezoid(pore_pdf, pore_d)
    pore_pdf = pore_pdf / area

    return pore_d, pore_pdf
 
#------------------------------------------------------------------------------

@njit(nogil=True, error_model="numpy")
def pore_pdf_campbell(soil, n_points=500, h_max=10*MPA_TO_CM_H2O):
    
    """
    Campbell-consistent pore-size distribution.
    
    The Campbell water retention curve is mapped to pore diameter with
    d = 2 * 0.149 / h, where h is in cm H2O and d is returned in m.
    
    The density is computed in log-diameter space, dF/dln(d), and then
    normalized over pore diameter. This gives an effective density [1/m]
    for constructing the pore CDF used by wet_f_Anand.
    
    Returns
    -------
    pore_d : ndarray
        Pore diameter [m].
    pore_pdf : ndarray
        Effective pore-size density [1/m].
    """

    psi_s_log_cm, b, K_s, n = soil_hydraulic_const(soil)

    pore_d_max = 2.0 * (0.149 / psi_s_log_cm) * 1e-2
    pore_d_min = 2.0 * (0.149 / h_max) * 1e-2

    pore_d = np.logspace(np.log10(pore_d_min), np.log10(pore_d_max), n_points)

    a = 1.0 / b

    F = (pore_d**a - pore_d_min**a) / (
        pore_d_max**a - pore_d_min**a
    )

    # Same convention as Ding/Anand: distribution per log-diameter interval.
    pore_pdf = _gradient(F, np.log(pore_d))

    pore_pdf = np.maximum(pore_pdf, 0.0)
    pore_pdf = pore_pdf / np.trapezoid(pore_pdf, pore_d)

    return pore_d, pore_pdf

#------------------------------------------------------------------------------

@njit(nogil=True, error_model="numpy")
def soil_pore_pdf(soil, pore_model="ding2016", n_points=500, h_min=10.0, h_max=10*MPA_TO_CM_H2O):
        
    """
    Estimate pore diameter PDF from soil texture.

    pore_model='ding2016' uses Ding et al. (2016) DE parameters.
    pore_model='campbell' uses the SMEW Campbell/Clapp-Hornberger closure.

    h_min and h_max define the suction range [cm H2O] used to map water retention to pore diameters.
    Lower h gives larger pores; higher h gives smaller pores, with d = 2 * 0.149 / h.
    Default h_max corresponds to 10 MPa, SMEW's hygroscopic-point suction
    """

    if pore_model == "ding2016":
        return pore_pdf_ding2016(
            soil,
            n_points=n_points,
            h_min=h_min,
            h_max=h_max,
        )

    if pore_model == "campbell":
        return pore_pdf_campbell(
            soil,
            n_points=n_points,
            h_max=h_max,
        )

    raise ValueError("pore_model must be 'ding2016' or 'campbell'.")
