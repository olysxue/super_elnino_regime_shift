
# **"Super El Niño events drive climate Regime Shifts with Enhanced Risks under Global Warming"**, Accepted in Nature Communications.

This repository contains the code supporting the manuscript:

## 🗂️ Repository Structure

This repository contains the analysis and visualization codes used in the study of **Super El Niño–induced climate regime shifts (CRSs)**.  
The scripts are organized by purpose — **main figures**, **supplementary analyses**, and **regime-shift detection algorithms** — and written primarily in **NCL** with supporting examples in **MATLAB**.

All codes are designed for reproducibility of key results shown in the main and supplementary figures of the paper.

---

## 📂 Repository Structure

| Folder | Description |
|--------|--------------|
| **`figure_code/`** | Contains NCL scripts for generating the **main figures** in the paper (Figs. 1–4). Each script reproduces a specific figure using observational and model datasets. |
| **`Supplementary_figure_code/`** | Includes additional scripts for **supplementary figures (S1–S21)**. These codes compute diagnostics such as regime-shift probabilities in multiple reanalysis datasets and associated physical mechanisms. |
| **`regime_shift_detection/`** | Provides representative examples of the **regime-shift detection framework**, implemented in both NCL and MATLAB. These demonstrate the Sequential t-test Analysis of Regime Shifts (STARS) algorithm used to identify persistent climate mean-state changes. |
| **`README.md`** | This file — overview of the repository and workflow. |

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
