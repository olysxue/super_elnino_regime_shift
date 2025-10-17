
# Main Figure Scripts (NCL v6.6.2)

This folder contains NCL scripts used to generate **Main Figures 1–4** in the manuscript:

**"Super El Niño events drive climate Regime Shifts with Enhanced Risks under Global Warming"**, Accepted in *Nature Communications*.

All scripts were developed and tested using **NCL version 6.6.2**.

---

## Script Overview and Data Sources

| Script | Figure | Description | Required Data |
|--------|--------|-------------|---------------|
| `Figure1_new_change_points_number_and_SE_delta_probability.ncl` | Fig. 1 |Calculates the number of detected regime‐shift change points and the probability difference between Super El Niño (SE) years and non-SE years | Regime‐shift detection results (change points, probabilities) and Super El Niño event years|
| `Figure2_plot_CESM_PD_2xCO2_4xCO2_zonal_ccurrent_TIW_intensity.ncl` | Fig. 2 |Projected changes of tropical mean background state modulating in CESM-UHR | CESM-UHR monthly mean outputs under PD, 2xCO2, and 4xCO2: wind stress, temp, zonal currents |
| `Figure3_CESM_seasonal_zonal_current_TIW_intensity_response.ncl` | Fig. 3 | Seasonal evolution of background conditions and TIW intensity under different CO2 simulations in CESM-UHR | four TIW metrics in CESM-UHR: EKE, SST, SSH, and VVEL fields|
| `Figure4_TIW_EKE_budget_analysis.ncl` | Fig. 4 | TIW-EKE budget analysis under different CO2 simulations in CESM-UHR | Barotrpic energy conversion rate and Baroclinic energy conversion rate in CESM-UHR |

---

## How to Run

Each script can be executed independently in a terminal with NCL:

ncl Figure*********.ncl
