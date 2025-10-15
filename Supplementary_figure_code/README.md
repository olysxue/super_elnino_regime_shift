
# Extended Data Figure Scripts (NCL v6.6.2)

This folder contains NCL scripts used to generate the Extended Data Figures (1–10) from the manuscript:

**"Super El Niño events drive climate Regime Shifts with Enhanced Risks under Global Warming"**, published in *Nature Communications*.

All scripts were developed and tested with **NCL version 6.6.2**. Each script corresponds to one Extended Data Figure in the paper.

---

## Script Overview

| Script | Corresponding Figure | Description | Required Data |
|--------|----------------------|-------------|---------------|
| `Extended_Figure1_TIW_CMEMS_CESM_E3SMv1v2.ncl` | Extended Data Fig. 1 | Comparison of observed (CMEMS) and modeled (CESM-UHR, E3SMv1/v2) TIW metrics (EKE, VVEL, SST, SSH) | TIW metrics in CESM-UHR |
| `Extended_Figure2_E3SMv1_v2_TIW_intensity.ncl` | Extended Data Fig. 2 | Changes in TIW-EKE, SST, and SSH between PD and 4xCO2 in E3SMv1 and v2 | TIW metrics in E3SMv1/v2 |
| `Extended_Figure3_GLORYS12_CESM_E3SM_comparison.ncl` | Extended Data Fig. 3 | Comparison of climatological SST and zonal current structures in models and GLORYS12 reanalysis | monthly TEMP and UVEL for GLORYS12, CESM-UHR and E3SM|
| `Extended_Figure4_CMIP6_highresMIP_TIW_intensity_change.ncl` | Extended Data Fig. 4 | Projected changes in SST-based TIW intensity across CMIP6 HighResMIP models | SST-based TIW metric in HighresMIP |
| `Extended_Figure5_E3SMv1_E3SMv2_mean_state.ncl` | Extended Data Fig. 5 | Projected mean-state changes in SST, wind stress and zonal currents in E3SM | monthly mean temp, wind stress and UVEL in E3SM |
| `Extended_Figure6_CESM_mean_state_seasonal_change.ncl` | Extended Data Fig. 6 | Seasonal changes in SST, winds, and currents in CESM-UHR under PD and 4xCO2 forcing | monthly mean TEMP, wind and UVEL in CESM-UHR |
| `Extended_Figure7_E3SMv1_E3SMv2_TIW_intensity_seasonal_cycle.ncl` | Extended Data Fig. 7 | Seasonal evolution of TIW intensity metrics (EKE, SST, SSH) in E3SMv1 and v2 | TIW metrics in E3SM |
| `Extended_Figure8_EKE_budget_seasonal_evolution.ncl` | Extended Data Fig. 8 | Seasonal evolution of barotropic and baroclinic energy conversion terms | TIW-EKE budget terms in CESM-UHR |
| `Extended_Figure9_mean_heat_momentum_flux_convergence.ncl` | Extended Data Fig. 9 | long-term mean TIW-induced heat flux and momentum flux convergence in CESM-UHR| TIW-induced feedback in CESM-UHR |
| `Extended_Figure10_different_EKE_calculation_method.ncl` | Extended Data Fig. 10 | Comparison between direct and indirect TIW-EKE calculation methods | Total KE, Large scale KE and TIW-EKE in CESM-UHR and E3SM |

---

## Usage

Each script is self-contained and can be run independently:

ncl Extended_Figure*******.ncl
