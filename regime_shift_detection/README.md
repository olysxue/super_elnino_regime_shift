
# 📁 Regime Shift Detection

## 🔍 Climate Regime Shift Detection Steps

The regime shift detection algorithm identifies statistically significant and persistent changes in the mean state of a climate variable (e.g., SST, SAT, soil moisture or sea ice).  
The method follows the **Sequential t-test Analysis of Regime Shifts (STARS)** framework (Rodionov 2004) with optional bootstrap confidence testing.

---

### **Step 1: Input Data**
- Load an anomaly time series (monthly or annual) of the target variable.  
- The series should be detrended.  
- Input can be either:
  - 1D (time series of regional mean), or  
  - 3D (time–lat–lon field for gridwise detection).

This method enables the identification of climate regime shifts (CRSs) in a purely data-driven manner, without requiring prior assumptions regarding the timing or number of shifts.

---

### **Step 2: Apply Pre-Filtering (Optional)**
- Smooth or band-pass filter the data to remove high-frequency variability (e.g., 2–8 day or interseasonal noise).
- Filtering is optional but recommended when focusing on decadal-scale CRSs.
---

### **Step 3: Define a Moving Window**
- Specify a window length `L` (usually 10 years for climate applications).  
- The algorithm compares the mean of the **current regime** with that of a **candidate new regime** `L` points ahead.
- 
---

### **Step 4: Perform Sequential t-Tests**
For each time t:
1. Compute the mean of two segments of length `L` (before and after t).  
2. Apply a two-sample Student’s t-test to test for a significant difference in means.  
3. If |t-statistic| > critical value (p < threshold), a potential shift is detected.

---

### **Step 5: Confirm Persistence**
- The detected shift must persist for at least `L` subsequent points with the same sign of anomaly.  
- Otherwise it is classified as a transient fluctuation and discarded.

#### 🧮 Regime Shift Index (RSI)
To verify persistence quantitatively, the **Regime Shift Index (RSI)** is calculated as the **cumulative sum of normalized anomalies** relative to the previous regime mean:

$$
\mathrm{RSI}(t)=\sum_{i=1}^{t}\frac{x_i-\bar{x}_{\mathrm{ref}}}{\sigma_{\mathrm{ref}}}
$$

**where:**
- **x<sub>i</sub>** — value of the variable at time step *i* (e.g., SST anomaly, SAT anomaly, soil moisture anomaly)  
- **x̄<sub>ref</sub>** — mean value of the previous (reference) regime  
- **σ<sub>ref</sub>** — standard deviation within the previous regime, used for normalization  
- **t** — time index (in years, months, or time steps, depending on input)  
- **RSI(t)** — cumulative index indicating the direction and strength of the regime shift  

A shift is confirmed when the RSI crosses zero and remains non-zero for at least cut-off length `L`, indicating a sustained departure from the prior mean state.  

**Interpretation**
- **Positive RSI →** persistent shift toward higher mean values
- **Negative RSI →** persistent shift toward lower mean values  
- The **magnitude** of RSI reflects both the strength and duration of the deviation.  

RSI complements the t-test by capturing the cumulative persistence of anomalies, making it especially useful for identifying decadal or multi-decadal transitions.

---

### **Step 6: Update Regime Mean**
- Once a significant and persistent shift is confirmed,  
  - Update the reference mean to the new regime.  
  - Record the shift year, mean difference, and significance level.

---

### **Typical Output**
| Variable | Description |
|-----------|--------------|
| `shift_years` | Years when significant regime changes occur |
| `mean_before`, `mean_after` | Means of two adjacent regimes |
| `t_value`, `p_value` | Test statistics |
| `persistence` | Boolean indicating sustained shift |

---

### **Example Interpretation**
A 10-year moving-window test applied to 1958–2024 SST anomalies may identify persistent regime shifts in:
- **1976–77:** Pacific climate transition,  
- **1998:** Post-Super El Niño cooling phase,  
- **2016:** Recent Super El Niño-induced shift toward a new warming regime.
---
# 📈 Regime Shift Occurrence Probability Calculation

To enhance the statistical robustness of the regime shift analysis, Climate Regime Shifts (CRSs) were detected in both **observational** and **climate model** datasets, for:
- **Annual mean**, and  
- **Four meteorological seasons:**  
  - **DJF** (December–February)  
  - **MAM** (March–May)  
  - **JJA** (June–August)  
  - **SON** (September–November)

---

### **1️⃣ Grid-level CRS Probability**

For each grid cell, the total probability of CRS occurrence is computed as:

$$
\mathrm{Prob}_t(x,y) = \frac{n(x,y)}{T}
$$

where:
- **n(x,y)** — number of years in which a regime shift was detected at grid location *(x, y)*  
- **T** — total number of years in the time series  

---

### **2️⃣ Super El Niño–Induced CRS Probability**

To quantify the specific influence of **Super El Niño (SE)** events, the CRS probability due to Super El Niño is defined as:

$$
\mathrm{Prob}_{SE}(x,y) = \frac{n_{SE}(x,y)}{T_{SE}}
$$

where:
- **n<sub>SE</sub>(x,y)** — number of regime shifts detected during Super El Niño development or following years  
- **T<sub>SE</sub>** — total number of Super El Niño years (e.g., 1982–83, 1997–98, 2015–16 in observations)

The enhanced CRS probability associated with Super El Niño is calculated as:

$$
\Delta \mathrm{Prob}(x,y) = \mathrm{Prob}_{SE}(x,y) - \mathrm{Prob}_t(x,y)
= \frac{n_{SE}(x,y)}{T_{SE}} - \frac{n(x,y)}{T}
$$

This provides a spatial and seasonal map of **Super El Niño–induced regime shift enhancement**.

---

### **3️⃣ Ensemble Aggregation (CESM2-LE)**

In CESM2-LE, regime-shift detection is performed for each of the 100 ensemble members.  
The occurrence probabilities are then **aggregated across all members** to improve robustness and reduce internal variability noise.

---

### **4️⃣ Global and Latitudinal Integration**

To obtain an integrated **global probability** view, the annual CRS probability is computed with **latitude weighting**:

$$
\mathrm{Prob}(i) = 
\frac{\int_{-\pi/2}^{\pi/2} \cos\phi \, n(i,\phi) \, d\phi}
{\int_{-\pi/2}^{\pi/2} \cos\phi \, N(\phi) \, d\phi}
$$

where:
- **n(i,φ)** — number of CRSs detected at latitude φ during year *i*  
- **N(φ)** — number of grid cells at latitude φ  
- **cosφ weighting** — ensures accurate area representation across latitudes

This latitude-weighted integration provides a globally consistent estimate of the temporal evolution of CRS occurrence probabilities.

---

### 🧮 Summary
- Detect regime shifts (annual and seasonal) for each grid cell  
- Calculate total and SE-specific CRS probabilities  
- Aggregate ensemble members for robustness  
- Compute latitude-weighted global CRS probability time series  

These steps together quantify **how and where Super El Niño events statistically increase the probability of climate regime shifts**.

---


## 📂 Directory Overview


This directory provides representative examples of the **regime shift detection framework** applied in the *Super El Niño–induced Climate Regime Shifts* study.  
It contains both **NCL** and **MATLAB** implementations of the **Sequential t-test Analysis of Regime Shifts (STARS)** algorithm (Rodionov 2004), illustrating how the detection was performed for different datasets and variables.

These are **demonstration files** — not the full analysis suite — intended to help users reproduce the key detection workflow and adapt it for their own applications.


| File | Language | Purpose |
|------|-----------|----------|
| `regime_shift_detection_sst_seasonally_ersst.ncl` | NCL | Example workflow for detecting SST regime shifts in the ERSST dataset. |
| `rsi_seasonally_sst_cesm2_lens.ncl` | NCL | Seasonal SST regime-shift detection in CESM2-LENS. |
| `rsi_seasonally_sat_cesm2_lens.ncl` | NCL | Regime-shift detection for surface air temperature (SAT) in CESM2-LENS. |
| `rsi_seasonally_soilw_cesm2_lens2.ncl` | NCL | Surface Soil moisture regime-shift detection example for CESM2-LENS. |
| `calculate_prob_SAT.ncl` | NCL | Calculates regime-shift probabilities for Super El Niño vs non-SE years after detection verification. |
| `stars.m` | MATLAB | Core STARS algorithm implementation for 1-D time series validation. |
| `stars_ERA5_land_soil_moisture.m` | MATLAB | MATLAB example applied to ERA5-Land soil-moisture anomalies. |

---
## ⚙️ Workflow Summary

1. **Detection** – Apply the STARS algorithm to identify significant and persistent mean-state shifts for each grid point or region.   
3. **Probability Calculation** –  
   - Count the **number of regime-shift occurrences** at each grid cell.   
   - Compute probabilities separately for:  
     - **Super El Niño (SE) years** (e.g., 1982/83, 1997/98, 2015/16)  
     - **Regular El Niño years**  
     - **Climatological mean (non-El Niño) years**  
   - Output includes spatial maps and temporal evolution of regime shift shift frequency and SE vs non-SE differences.
---

## 🧠 Notes

- Users can modify the window length (**L**) and significance level (**p**) within each script header.  
All scripts are **modular** and can be easily adapted to different climate variables or datasets.
---


### **References**

## 📚 References

**(1)** **Rodionov, S. N.** (2004). *A sequential algorithm for testing climate regime shifts.*  
**Geophysical Research Letters**, **31**, L09204. [https://doi.org/10.1029/2004GL019448](https://doi.org/10.1029/2004GL019448)

**(2)** **Rodionov, S. N.** (2005a). *A brief overview of the regime shift detection methods.*  
In: **Large-Scale Disturbances (Regime Shifts) and Recovery in Aquatic Ecosystems: Challenges for Management Toward Sustainability**,  
V. Velikova and N. Chipev (Eds.), UNESCO–ROSTE/BAS Workshop on Regime Shifts, 14–16 June 2005, Varna, Bulgaria, pp. 17–24.

**(3)** **Rodionov, S. N.** (2005b). *Detecting regime shifts in the mean and variance: Methods and specific examples.*  
In: **Large-Scale Disturbances (Regime Shifts) and Recovery in Aquatic Ecosystems**,  
V. Velikova and N. Chipev (Eds.), UNESCO–ROSTE/BAS Workshop on Regime Shifts, 14–16 June 2005, Varna, Bulgaria, pp. 68–72.

**(4)** **Rodionov, S. N.** (2006). *The use of prewhitening in climate regime shift detection.*  
**Geophysical Research Letters**, **31**, L12707. [https://doi.org/10.1029/2006GL025904](https://doi.org/10.1029/2006GL025904)

**(5)** **Rodionov, S. N.**, and **Overland, J. E.** (2005). *Application of a sequential regime shift detection method to the Bering Sea ecosystem.*  
**ICES Journal of Marine Science**, **62**, 328–332. [https://doi.org/10.1016/j.icesjms.2005.01.013](https://doi.org/10.1016/j.icesjms.2005.01.013)


---
### 🔗 Useful Websites
For detailed algorithm explanation and original source code, see:
- NOAA Bering Climate Regime Shift site: https://www.beringclimate.noaa.gov/regimes/
- STARS method documentation: https://www.afsc.noaa.gov/REFM/REEM/STARS/

---

