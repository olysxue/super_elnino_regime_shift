
# **"Super El Niño events drive climate Regime Shifts with Enhanced Risks under Global Warming"**, Accepted in Nature Communications (2025).

## 📘 Overview

This repository contains the analysis and visualization codes used in the study of **Super El Niño–induced Climate Regime Shifts (CRSs)**, which quantifies how extreme ENSO events trigger and amplify abrupt transitions in global climate states.

The study integrates:
- **Observational datasets** (e.g., ERSST, HadISST, ERA5, GLDAS)  
- **Large-ensemble simulations** (e.g., CESM2-LE, CMIP6)  
- **Statistical regime-shift detection algorithms** (based on the Sequential t-test Analysis of Regime Shifts, *STARS*)  

All scripts are written primarily in **NCL** (NCAR Command Language), with supporting examples in **MATLAB**.

---

## 📂 Repository Structure

| Folder | Description |
|--------|--------------|
| **`figure_code/`** | Contains NCL scripts for generating the **main figures** in the paper (Figs. 1–4). Each script reproduces a specific figure using observational and model datasets. |
| **`Supplementary_figure_code/`** | Includes additional scripts for **supplementary figures (S1–S21)**. These codes compute diagnostics such as regime-shift probabilities in multiple reanalysis datasets and associated physical mechanisms. |
| **`regime_shift_detection/`** | Provides representative examples of the **regime-shift detection framework**, implemented in both NCL and MATLAB. These demonstrate the Sequential t-test Analysis of Regime Shifts (STARS) algorithm used to identify persistent climate mean-state changes. |

---
## ⚙️ Workflow

1. **Data Preparation**
   - Input  data for SST, Surface Air temperature and surface soil moisture for observational and reanalysis products (e.g., HadISST, ERSSTv5, ERA5, ERA5-Land, CPC, GLDAS) and model simulations (CESM2-LE, CMIP6).
   - All datasets are **detrended** and **regridded** to 1° × 1° resolution, and expressed as anomalies** prior to detection.

2. **Regime Shift Detection**
   - Conducted using the **STARS algorithm** (Rodionov, 2004) with configurable parameters:
     - Window length `L` (6–12 years)
     - Significance level `p` (0.05, 0.01, 0.001)
     - Gridwise detection for gridded datasets
      
3. **Regime Shift Occurrence Probability Calculation**
   - After detection, regime shift probabilities are calculated for:
   - Total CRS frequency at each grid  
   - Super El Niño–induced CRS probability
   - Latitude-weighted global CRS probability time series

4. **Visualization**
   - Scripts in `figure_code/` and `Supplementary_figure_code/` generate maps, time series, and probability plots for Super El Niño–induced regime-shift patterns.

---

## 🛠️ Dependencies

- [**NCL**](https://www.ncl.ucar.edu/) (tested with version **6.6.2**)  
- [**CDO**](https://code.mpimet.mpg.de/projects/cdo) (Climate Data Operators)
- [**NCO**](https://nco.sourceforge.net/) (NetCDF Operators; used for variable merging, dimension trimming, etc.)
- [**MATLAB R2021a**]

---
## 📧 Contact
For questions or assistance with the code and data, please contact:

Aoyun Xue
Postdoctoral Researcher
📩 aoyunxue@ucsb.edu
University of California, Santa Barbara
