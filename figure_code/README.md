
# Main Figure Scripts (NCL v6.6.2)

This folder contains NCL scripts used to generate **Main Figures 1–4** in the manuscript:

**"Super El Niño events drive climate Regime Shifts with Enhanced Risks under Global Warming"**, Accepted in *Nature Communications*.

All scripts were developed and tested using **NCL version 6.6.2**.

---

## Script Overview and Data Sources

| Script | Figure | Description | Required Data |
|--------|--------|-------------|---------------|
| `Figure1_new_change_points_number_and_SE_delta_probability.ncl` | Fig. 1 |Calculates the number of detected regime‐shift change points and the probability difference between Super El Niño (SE) years and non-SE years | Regime‐shift detection results (change points, probabilities) and Super El Niño event years|
| `Figure2_yearly_plot_clm_shift_probability_table.ncl` | Fig. 2 | Plots the yearly probability of climate regime shifts, highlighting differences between Super El Niño and normal years | Yearly or decadal regime‐shift probability data derived from observational time series |
| `Figure3_plot_map_pattern_regime_shift_evidence.ncl` | Fig. 3 | Super El Niño–induced global regime-shift evidence. The script composites regime-shift probabilities and normalized anomalies (SST, SAT, and soil moisture) between Super El Niño (SE) and non-SE years, highlighting regions where climate variables experienced statistically significant mean-state shifts following SE events | Multi-variable gridded datasets (SST, SAT, soil moisture) and regime-shift detection outputs; Super El Niño event years list (e.g., 1982/83, 1997/98, 2015/16).|
| `Figure4_CESM2_LE_historical_VS_SSP585_prob.ncl` | Fig. 4 | Compares regime‐shift probabilities between CESM2-LE historical and SSP5-8.5 simulations, focusing on future changes in Super El Niño–related regime-shift likelihood | CESM2-LE ensemble outputs (historical and SSP5-8.5) with SE and non-SE regime-shift probability diagnostics |

---

## How to Run

Each script can be executed independently in a terminal with NCL:

ncl Figure*********.ncl
