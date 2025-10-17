
# **"Super El Niño events drive climate Regime Shifts with Enhanced Risks under Global Warming"**, Accepted in Nature Communications.


## 🗂️ Repository Structure

This repository contains the analysis and visualization codes used in the study of **Super El Niño–induced climate regime shifts (CRSs)**.  
The scripts are organized by purpose — **main figures**, **supplementary analyses**, and **regime-shift detection algorithms** — and written primarily in **NCL** with supporting examples in **MATLAB**.

---

## 📂 Repository Structure

| Folder | Description |
|--------|--------------|
| **`figure_code/`** | Contains NCL scripts for generating the **main figures** in the paper (Figs. 1–4). Each script reproduces a specific figure using observational and model datasets. |
| **`Supplementary_figure_code/`** | Includes additional scripts for **supplementary figures (S1–S21)**. These codes compute diagnostics such as regime-shift probabilities in multiple reanalysis datasets and associated physical mechanisms. |
| **`regime_shift_detection/`** | Provides representative examples of the **regime-shift detection framework**, implemented in both NCL and MATLAB. These demonstrate the Sequential t-test Analysis of Regime Shifts (STARS) algorithm used to identify persistent climate mean-state changes. |

---
## ⚙️ Workflow Overview

1. **Data Preparation**
   - Input data include observational and reanalysis products (e.g., HadISST, ERSSTv5, ERA5, ERA5-Land, CPC, GLDAS) and model simulations (CESM2-LE, CMIP6).
   - All datasets are **detrended** and **regridded** to 1° × 1° resolution before analysis.

2. **Regime Shift Detection**
   - Conducted using the **STARS algorithm** (Rodionov, 2004) with configurable parameters:
     - Window length `L` (6–12 years)
     - Significance level `p` (0.05, 0.01, 0.001)
   - Implemented in both NCL (for global fields) and MATLAB (for 1-D validation).

3. **Visualization**
   - Scripts in `figure_code/` and `Supplementary_figure_code/` generate maps, time series, and probability plots for Super El Niño–induced regime-shift patterns.

---

## 📦 Data Availability

▶ **Raw model output sources:**
---

## 🛠️ Dependencies

- [**NCL**](https://www.ncl.ucar.edu/) (tested with version **6.6.2**)  
- [**CDO**](https://code.mpimet.mpg.de/projects/cdo) (Climate Data Operators)
- [**NCO**](https://nco.sourceforge.net/) (NetCDF Operators; used for variable merging, dimension trimming, etc.)

---

## 📦 Citation

If you use this code, please cite:
---
## 📧 Contact
For questions or assistance with the code and data, please contact:

Aoyun Xue
Postdoctoral Researcher
📩 aoyunxue@ucsb.edu
University of California, Santa Barbara
