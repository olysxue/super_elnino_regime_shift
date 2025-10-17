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

### **References**

Rodionov, S.N., 2004: A sequential algorithm for testing climate regime shifts. Geophys. Res. Lett., 31, L09204, doi:10.1029/2004GL019448.

Rodionov, S.N., 2005a: A brief overview of the regime shift detection methods. In: Large-Scale Disturbances (Regime Shifts) and Recovery in Aquatic Ecosystems: Challenges for Management Toward Sustainability, V. Velikova and N. Chipev (Eds.), UNESCO-ROSTE/BAS Workshop on Regime Shifts, 14-16 June 2005, Varna, Bulgaria, 17-24.

Rodionov, S.N., 2005b: Detecting regime shifts in the mean and variance: Methods and specific examples. In: Large-Scale Disturbances (Regime Shifts) and Recovery in Aquatic Ecosystems: Challenges for Management Toward Sustainability, V. Velikova and N. Chipev (Eds.), UNESCO-ROSTE/BAS Workshop on Regime Shifts, 14-16 June 2005, Varna, Bulgaria, 68-72.

Rodionov, S.N., 2006: The use of prewhitening in climate regime shift detection, Geophys. Res. Lett., 31, L12707, doi:101029/2006GL025904.

Rodionov, S.N., and J.E. Overland, 2005: Application of a sequential regime shift detection method to the Bering Sea ecosystem. ICES J. Mar. Sci., 62: 328-332.

Rodionov, S.N., J.E. Overland, and N.A. Bond 2004: Detecting Change in the Bering Sea Ecosystem. A PowerPoint file (1.4 Mb) of the presentation at the Symposium on Marine Science in Alaska, Jan 12-14, 2004, Anchorage, Alaska

---
### 🔗 Useful Websites
For detailed algorithm explanation and original source code, see:
- NOAA Bering Climate Regime Shift site: https://www.beringclimate.noaa.gov/regimes/
- STARS method documentation: https://www.afsc.noaa.gov/REFM/REEM/STARS/

---

