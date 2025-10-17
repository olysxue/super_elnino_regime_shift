
# Supplemrntary Figure Scripts (NCL v6.6.2)

This folder contains NCL scripts used to generate the Supplementary Figures (1–21) from the manuscript:

**"Super El Niño events drive climate Regime Shifts with Enhanced Risks under Global Warming"**, Accepted in *Nature Communications*.

All scripts were developed and tested with **NCL version 6.6.2**. Each script corresponds to one Extended Data Figure in the paper.

---


All scripts were developed and tested using **NCL version 6.6.2+**.

---

## Script Overview and Data Sources

| Script | Figure | Description | Required Data |
|:--|:--:|:--|:--|
| `S1_plot_mean_regular_ElNino_LaNina_prob.ncl` | Fig. S1 | Calculates global climatological and ENSO-related regime-shift probabilities for SST, SAT, and soil-moisture datasets. | Observational datasets (ERSSTv5, HadISST, ERA5, CPC) from 1948–2022. |
| `S2_plot_SST_different_datasets_probability1_1948_2022.ncl` | Fig. S2 | Maps increased CRS probabilities due to Super El Niño across multiple SST datasets (ERSSTv5, HadISST, COBE2, Kaplan, etc.). | SST datasets, detrended and regridded (2.5° – 1°). |
| `S3_plot_SAT_different_datasets_probability_1948_2022.ncl` | Fig. S3 | Maps increased CRS probabilities due to Super El Niño across multiple SAT datasets (NCEP, ERA5, CRUTEM4, GISS, etc.). | Global SAT datasets (land + reanalysis). |
| `S4_plot_Soil_different_datasets_probability_1948-2022.ncl` | Fig. S4 | Maps increased CRS probabilities due to Super El Niño across multiple soil-moisture datasets. | CPC, ERA5-Land, GLDAS, GLEAM, MERRA-2, etc. |
| `S5_plot_SST_different_datasets_regime_shift_probability_table_new.ncl` | Fig. S5 | Plots annual CRS probabilities (1948–2022) for SST datasets and compares SE vs non-SE years. | Annual CRS probabilities from multiple SST datasets. |
| `S6_plot_SAT_different_datasets_regime_shift_probility_table_1948-2022_new.ncl` | Fig. S6 | Plots annual CRS probabilities for SAT datasets and highlights enhanced probabilities during Super El Niño events. | SAT CRS probability datasets. |
| `S7_plot_soil_different_datasets_regime_shift_probability_table_new.ncl` | Fig. S7 | Plots annual CRS probabilities for soil-moisture datasets and multi-dataset mean differences between SE and non-SE years. | Soil-moisture CRS probabilities from 10 datasets. |
| `S8_plot_clm_shift_probability_table_L6.ncl` | Fig. S8 | Computes regime-shift probabilities using cut-off L = 6 and p = 0.05; compares SE vs non-SE years. | Annual CRS probability data (SST, SAT, soil moisture). |
| `S9_plot_clm_shift_probability_table_L8.ncl` | Fig. S9 | Same as S8 but L = 8. | Same as above. |
| `S10_plot_clm_shift_probability_table_L12.ncl` | Fig. S10 | Same as S8 but L = 12. | Same as above. |
| `S11_plot_clm_shift_probability_table_p01.ncl` | Fig. S11 | Same as S8 but p = 0.01. | Same as above. |
| `S12_plot_clm_shift_probability_table_p001.ncl` | Fig. S12 | Same as S8 but p = 0.001. | Same as above. |
| `S13_plot_GODAS_SSTA_MLD_depth_Super_El_Nino_central_pacific_reemergence.ncl` | Fig. S13 | Diagnoses the reemergence mechanism of SE-induced SST regime shifts in the Central-North Pacific using mixed-layer heat budget. | GODAS/ERA5 temperature, MLD, and wind data. |
| `S14_plot_NCEP_soil_moisture_precipitation_MAM.ncl` | Fig. S14 | Composites MAM soil-moisture and DJF precipitation anomalies during SE years. | NCEP precipitation + multi-source soil-moisture datasets. |
| `S15_plot_NCEP_soil_moisture_precipitation_JJA.ncl` | Fig. S15 | Composites JJA soil-moisture and MAM precipitation anomalies during SE decaying phases. | Same datasets as S14. |
| `S16_plot_soil_moisture_memory.ncl` | Fig. S16 | Calculates 1-year autocorrelation of soil-moisture anomalies to quantify memory and its relation to SE-induced CRS patterns. | ERA5-Land, GLDAS, GLEAM, MERRA-2 datasets. |
| `S17_plot_CESM2_LE_sat_sst_soil_hist_ssp_rs_prob.ncl` | Fig. S17 | Calculates CRS probabilities in CESM2-LE historical and SSP3-7.0 simulations; evaluates SE vs non-SE differences. | CESM2-LE monthly outputs (SST, SAT, soil moisture). |
| `S18_plot_CESM2_LE_regular_ENSO_hist.ncl` | Fig. S18 | Climatological and ENSO-composite changes in CRS probabilities under CESM2-LE historical runs. | CESM2-LE historical outputs (SST, SAT, soil moisture). |
| `S19_plot_CESM2_LE_regular_ENSO_SSP.ncl` | Fig. S19 | Same as S18 but for SSP3-7.0 simulations. | CESM2-LE SSP3-7.0 outputs. |
| `S20_plot_CMIP6_sat_sst_soil_hist_ssp_rs_prob.ncl` | Fig. S20 | Compares CRS probabilities in CMIP6 historical and SSP5-8.5 simulations and quantifies SE impact. | 48 CMIP6 models (listed in Table S1). |
| `S21_SE_arctic_Antarctic_SIC_SIE_regime_shift.ncl` | Fig. S21 | Analyzes SE-induced CRS probabilities in Arctic and Antarctic sea ice concentration and extent. | Satellite SIC/SIE datasets (1979–2022). |

---

## How to Run

Each script can be executed independently in a terminal with NCL:

ncl S1_plot_mean_regular_ElNino_LaNina_prob.ncl

---
