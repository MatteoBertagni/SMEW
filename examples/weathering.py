# SPDX-License-Identifier: AGPL-3.0-only
"""A repeatable SMEW example, adapted from Example.ipynb."""

import marimo

__generated_with = "0.24.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import smew

    return mo, np, smew


@app.cell
def _(mo):
    mo.md(r"""
    # Soil weathering: a repeatable example

    Follow water, organic carbon and mineral weathering in a loamy soil.
    This example follows **Example.ipynb**, retaining its **365-day duration**
    and ten-minute timestep. It uses **fixed weekly rainfall** so that repeated
    runs can be compared.

    Read the explanations and change the inputs below to explore the model.
    Units are shown next to each input. The calculation runs automatically;
    the **Generate plots** button at the end draws the results on demand.

    ## 1. Simulation time and soil
    A ten-minute timestep resolves the coupled water and chemistry calculations.
    Moisture is the fraction of the soil pore volume filled with water.
    """)
    return


@app.cell
def _():
    duration_days = 365
    timestep_minutes = 10
    return duration_days, timestep_minutes


@app.cell
def _():
    # 1: moisture changes with rainfall and losses; 0: maintain constant moisture.
    moisture_mode = 1
    return (moisture_mode,)


@app.cell
def _():
    # Account for snow storage and melting when soil temperature crosses 0 °C.
    frozen_soil_hydrology = True
    snow_melt_rate_m_day_c = 0.005
    return frozen_soil_hydrology, snow_melt_rate_m_day_c


@app.cell
def _():
    soil_type = "loam"
    soil_depth_m = 0.3
    soil_bulk_density_g_m3 = 1.2e6
    initial_moisture = 0.5
    return initial_moisture, soil_bulk_density_g_m3, soil_depth_m, soil_type


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Weather and vegetation
    Seasonal air temperature determines soil temperature and evaporation demand.
    Rain falls on day 1 and every seven days thereafter; each event supplies
    seven days' share of the annual rainfall. These fixed events replace the
    random rainfall used in the original example.

    Vegetation starts at its carrying capacity, as in the original example.
    """)
    return


@app.cell
def _():
    latitude_degrees = 40
    first_day_of_year = 1
    altitude_m = 33
    wind_speed_m_s = 1.0
    albedo = 0.25
    coastal = False
    mean_air_temperature_c = 13
    annual_temperature_amplitude_c = 11
    daily_temperature_amplitude_c = 5
    annual_rainfall_m = 1.2
    rain_interval_days = 7
    first_rain_day = 1
    return (
        albedo,
        altitude_m,
        annual_rainfall_m,
        annual_temperature_amplitude_c,
        coastal,
        daily_temperature_amplitude_c,
        first_day_of_year,
        first_rain_day,
        latitude_degrees,
        mean_air_temperature_c,
        rain_interval_days,
        wind_speed_m_s,
    )


@app.cell
def _():
    vegetation_capacity_g_m2 = 3000
    initial_vegetation_g_m2 = 3000
    vegetation_growth_days = 100
    vegetation_start_day = 0
    root_area_index = 10  # m² of roots / m² of ground
    root_diameter_m = 0.4e-3
    return (
        initial_vegetation_g_m2,
        root_area_index,
        root_diameter_m,
        vegetation_capacity_g_m2,
        vegetation_growth_days,
        vegetation_start_day,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. Carbon and initial soil chemistry
    Organic carbon decomposition and root respiration supply soil CO₂.
    Exchange sites on the soil hold positively charged ions; the initial
    fractions below share these sites and must add up to one.

    The soil initially contains no carbonate minerals. Background ion inputs
    compensate for initial losses when `balance_background_inputs` is 1.
    """)
    return


@app.cell
def _():
    litter_input_g_c_m2_day = 1.0
    organic_carbon_percent = 0.05
    soil_co2_multiple_of_atmosphere = 10
    root_to_microbial_respiration_ratio = 1.0
    initial_ph = 4.0
    exchange_capacity_mmol_per_100g = 10.0
    # Order: calcium, magnesium, potassium, sodium, aluminium, hydrogen.
    initial_exchange_fractions = [0.30, 0.15, 0.10, 0.05, 0.00, 0.40]
    initial_silicon_umol_l = 0.0
    initial_caco3_umol_m2 = 0.0
    initial_mgco3_umol_m2 = 0.0
    balance_background_inputs = 1
    return (
        balance_background_inputs,
        exchange_capacity_mmol_per_100g,
        initial_caco3_umol_m2,
        initial_exchange_fractions,
        initial_mgco3_umol_m2,
        initial_ph,
        initial_silicon_umol_l,
        litter_input_g_c_m2_day,
        organic_carbon_percent,
        root_to_microbial_respiration_ratio,
        soil_co2_multiple_of_atmosphere,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ## 4. Rock application
    Add forsterite powder at the start of the simulation. Set the rock mass
    to zero to explore the same soil without an enhanced-weathering treatment.
    Particle diameters are specified in micrometres (µm).
    """)
    return


@app.cell
def _():
    rock_mass_g_m2 = 1000.0
    application_day = 0.0
    return application_day, rock_mass_g_m2


@app.cell
def _():
    minerals = ["forsterite"]
    mineral_mass_fractions = [1.0]
    particle_diameters_um = [100.0]
    particle_mass_fractions = [1.0]
    dissolution_factor = 1.0
    # Surface wetness can be "constant", "linear", or "nonlinear".
    wet_surface_model = "linear"
    pore_model = "ding2016"  # Used only by the nonlinear option.
    rock_density_g_m3 = 3.0e6
    pore_particle_mixing = 1.0
    return (
        dissolution_factor,
        mineral_mass_fractions,
        minerals,
        particle_diameters_um,
        particle_mass_fractions,
        pore_model,
        pore_particle_mixing,
        rock_density_g_m3,
        wet_surface_model,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ## 5. Run the model
    The following cells convert the inputs into model units, calculate weather
    and water balance, then organic carbon and soil chemistry. The `results`
    dictionary collects numerical time series for plotting and automated tests.
    """)
    return


@app.cell
def _(duration_days, np, timestep_minutes):
    dt = timestep_minutes / (24 * 60)  # days
    time_days = np.arange(0, duration_days, dt)
    conv_mol = 1e6  # mol → µmol
    conv_al = 1e3  # additional scaling for aluminium species
    return conv_al, conv_mol, dt, time_days


@app.cell
def _(
    albedo,
    altitude_m,
    annual_rainfall_m,
    annual_temperature_amplitude_c,
    coastal,
    daily_temperature_amplitude_c,
    dt,
    duration_days,
    first_day_of_year,
    first_rain_day,
    frozen_soil_hydrology,
    initial_moisture,
    initial_vegetation_g_m2,
    latitude_degrees,
    mean_air_temperature_c,
    moisture_mode,
    np,
    rain_interval_days,
    smew,
    snow_melt_rate_m_day_c,
    soil_depth_m,
    soil_type,
    time_days,
    vegetation_capacity_g_m2,
    vegetation_growth_days,
    vegetation_start_day,
    wind_speed_m_s,
):
    latitude = np.deg2rad(latitude_degrees)
    temp_air, temp_soil, temp_min, temp_max = smew.temp(
        latitude, mean_air_temperature_c, annual_temperature_amplitude_c,
        daily_temperature_amplitude_c, soil_depth_m, duration_days, dt,
        first_day_of_year,
    )
    wind = np.full(time_days.shape, wind_speed_m_s)
    et0 = smew.ET0(
        latitude, altitude_m, temp_air, temp_soil, temp_min, temp_max,
        wind, albedo, soil_depth_m, coastal, duration_days, dt, first_day_of_year,
    )
    rain = np.zeros_like(time_days)
    for _day in np.arange(first_rain_day, duration_days, rain_interval_days):
        _index = int(round(_day / dt))
        if _index < len(rain):
            rain[_index] += annual_rainfall_m * rain_interval_days / 365

    vegetation = smew.veg(
        initial_vegetation_g_m2, vegetation_growth_days, vegetation_capacity_g_m2,
        vegetation_start_day, temp_soil, dt,
    )
    moisture, s_w, s_i, infiltration, leakage, transpiration, evaporation, runoff, irrigation, porosity = smew.moisture_balance(
        rain, soil_depth_m, soil_type, et0, vegetation, vegetation_capacity_g_m2,
        moisture_mode, initial_moisture, duration_days, dt,
        temp_soil=temp_soil if frozen_soil_hydrology else None,
        melt_rate=snow_melt_rate_m_day_c,
    )
    return (
        evaporation,
        infiltration,
        leakage,
        moisture,
        porosity,
        rain,
        runoff,
        temp_soil,
        transpiration,
        vegetation,
    )


@app.cell
def _(
    conv_mol,
    dt,
    litter_input_g_c_m2_day,
    moisture,
    organic_carbon_percent,
    root_to_microbial_respiration_ratio,
    smew,
    soil_bulk_density_g_m3,
    soil_co2_multiple_of_atmosphere,
    soil_depth_m,
    soil_type,
    temp_soil,
    vegetation,
    vegetation_capacity_g_m2,
):
    initial_soc = soil_bulk_density_g_m3 * organic_carbon_percent / 100
    initial_soil_co2 = soil_co2_multiple_of_atmosphere * smew.CO2_atm(conv_mol)
    soc, microbial_respiration, root_respiration, diffusivity = smew.respiration(
        litter_input_g_c_m2_day, initial_soc, initial_soil_co2,
        root_to_microbial_respiration_ratio, soil_type, moisture, vegetation,
        vegetation_capacity_g_m2, soil_depth_m, temp_soil, dt, conv_mol,
    )
    return diffusivity, microbial_respiration, root_respiration, soc


@app.cell
def _(
    application_day,
    balance_background_inputs,
    conv_al,
    conv_mol,
    diffusivity,
    dissolution_factor,
    dt,
    exchange_capacity_mmol_per_100g,
    infiltration,
    initial_caco3_umol_m2,
    initial_exchange_fractions,
    initial_mgco3_umol_m2,
    initial_ph,
    initial_silicon_umol_l,
    leakage,
    microbial_respiration,
    mineral_mass_fractions,
    minerals,
    moisture,
    np,
    particle_diameters_um,
    particle_mass_fractions,
    pore_model,
    pore_particle_mixing,
    porosity,
    rock_density_g_m3,
    rock_mass_g_m2,
    root_area_index,
    root_diameter_m,
    root_respiration,
    smew,
    soil_bulk_density_g_m3,
    soil_depth_m,
    soil_type,
    temp_soil,
    transpiration,
    vegetation,
    vegetation_capacity_g_m2,
    wet_surface_model,
):
    cec_total = (
        exchange_capacity_mmol_per_100g * 1e-5
        * soil_bulk_density_g_m3 * soil_depth_m * conv_mol
    )
    cec_fractions = np.array(initial_exchange_fractions)
    initial_concentrations, exchange_constants = smew.f_CEC_to_conc(
        cec_fractions, initial_ph, soil_type, conv_mol, conv_al,
    )
    if wet_surface_model == "nonlinear":
        pore_diameters, pore_pdf = smew.soil_pore_pdf(
            soil_type, pore_model=pore_model,
        )
    else:
        pore_diameters, pore_pdf = None, None
    chemistry = smew.biogeochem_balance(
        n=porosity, s=moisture, L=leakage, T=transpiration, I=infiltration,
        v=vegetation, k_v=vegetation_capacity_g_m2, RAI=root_area_index,
        root_d=root_diameter_m, Zr=soil_depth_m, r_het=microbial_respiration,
        r_aut=root_respiration, D=diffusivity, temp_soil=temp_soil,
        pH_in=initial_ph, conc_in=initial_concentrations, f_CEC_in=cec_fractions,
        K_CEC=exchange_constants, CEC_tot=cec_total, Si_in=initial_silicon_umol_l,
        CaCO3_in=initial_caco3_umol_m2, MgCO3_in=initial_mgco3_umol_m2,
        M_rock_in=rock_mass_g_m2, t_app=application_day, mineral=minerals,
        rock_f_in=np.array(mineral_mass_fractions),
        d_in=np.array(particle_diameters_um) * 1e-6,
        psd_perc_in=np.array(particle_mass_fractions),
        SSA_in=np.nan,  # estimate surface area from particle size
        diss_f=dissolution_factor, dt=dt, conv_Al=conv_al, conv_mol=conv_mol,
        keyword_add=balance_background_inputs,
        keyword_ssa=wet_surface_model, pore_d_in=pore_diameters,
        pore_pdf_in=pore_pdf, rho_rock_in=rock_density_g_m3,
        mixalf_in=pore_particle_mixing,
    )
    return (chemistry,)


@app.cell
def _(
    chemistry,
    evaporation,
    infiltration,
    leakage,
    microbial_respiration,
    moisture,
    rain,
    root_respiration,
    runoff,
    soc,
    temp_soil,
    time_days,
    transpiration,
    vegetation,
):
    results = {
        "time_days": time_days,
        "rain": rain,
        "temp_soil": temp_soil,
        "s": moisture,
        "I": infiltration,
        "L": leakage,
        "T": transpiration,
        "E": evaporation,
        "Q": runoff,
        "v": vegetation,
        "SOC": soc,
        "r_het": microbial_respiration,
        "r_aut": root_respiration,
    }
    for _name in (
        "pH", "Alk", "Ca", "Mg", "K", "Na", "Si", "DIC", "CO2_air",
        "f_Ca", "f_Mg", "f_K", "f_Na", "f_Al", "f_H", "M_rock", "EW",
        "wet_f",
    ):
        results[_name] = chemistry[_name]
    return (results,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 6. Inspect the results
    Click below to generate the plots from the calculated time series.
    """)
    return


@app.cell
def _(mo):
    plot_button = mo.ui.run_button(label="Generate plots")
    plot_button
    return (plot_button,)


@app.cell
def _(mo, plot_button, results):
    mo.stop(not plot_button.value)
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(3, 2, figsize=(10, 9), sharex=True)
    for _axis, _name, _label in zip(
        axes.flat,
        ("s", "pH", "Alk", "SOC", "M_rock", "Mg"),
        ("Moisture (fraction)", "pH", "Alkalinity (µmol/L)",
         "Organic carbon (g C/m³)", "Rock mass (g/m²)", "Magnesium (µmol/L)"),
    ):
        _axis.plot(results["time_days"], results[_name])
        _axis.set_ylabel(_label)
        _axis.set_xlabel("Time (days)")
        _axis.grid(alpha=0.25)
    figure.tight_layout()
    mo.output.replace(figure)
    plt.close(figure)
    return


if __name__ == "__main__":
    app.run()
