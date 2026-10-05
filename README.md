# find_streaks

`find_streaks` is a Python module designed to detect readout and pileup streaks in Chandra X-ray Observatory ACIS observations.

Given a Level-2 FITS event file (`acisf..._evt2.fits`), the pipeline searches every active CCD chip for streaks, matches results with archival Chandra Data Archive (CDA) records, saves diagnostic plots, and outputs a summary CSV table.

The module contains several geometry and math utilities, with `find_streaks()` serving as the primary caller function that initiates the entire detection sequence.


## Method & Detection Logic
The find_streaks logic a two tier spatial test on each CCD to flag suspected streaks. The calculations are performed on the event file's sky coordinates ('[x,y]') after dither correction has been applied, concentrating streaks into consistent, continuous lines.

### 1. Orientation & Slicing:
First, the borders of the CCD are estimated and two non parallel vector, representing the two axes of the borders are chosen. Because streaks traverse the physical readout columns along the $\text{CHIPY}$ direction, spatial band slicing is performed along the axis coinciding with the $\text{CHIPX}$ direction.
Slicing is executed by computing the norm of the vector projection of each event onto the chosen axis vector (sharing the top-left corner as the origin). Events are then divided into discrete bands according to the scalar value of their projection norm.

### 2. Primary Check
The band containing the maximum photon count is identified and assessed against the background expectation. If the peak exceeds the background mean by more than $n_1 \cdot \sigma$ (where $\sigma$ is the standard deviation across non-peak bands), it is flagged as a streak candidate and forwarded to the secondary check.

### 3. Secondary Check
In order to rule out maximas caused by bright sources, the secondary check is performed oved two halves of the image along the readout direction and the same asessment is performed independantly along the two halves, according to a second threshold $n_2$. Both halves must independently show a statistically significant peak ($> n_2 \cdot \sigma$) at the identical projection index (sharing the same whole band).

### 4. Adaptive Second Iteration
If a narrow streak falls on the boundary between two adjacent bands, its counts may be split across bands, leading to a false negative. To mitigate this edge case, an adaptive second pass runs with a rescaled bandwidth if either check fails. The second iteration triggers at most once per chip, scaling the bandwidth by `second_iter_factor` (default: `0.8`).

Flow chart describing the function's decision tree:

```text
[Raw FITS Events for CCD]
            │
            ▼
1. Geometric Alignment (Estimate CCD corners & match CHIPX/CHIPY readout axes)
            │
            ▼
2. Spatial Projection (Divide photon coordinates into slices of width = bandwidth)
            │
            ▼
3. Primary Check: Is there an anomalous count peak?
   (Peak band count > background mean + n1 * sigma)
            │
            ├── NO  ──► [Adaptive 2nd Iteration: rerun at second_iter_factor * bandwidth] ──► Final: False
            │
            └── YES ──► 4. Bilateral Check: Does the peak cross the ENTIRE chip?
                        (Bisect chip along the orthogonal axis into Side 0 & Side 1)
                        Do BOTH halves share a coincident peak > n2 * sigma?
                              │
                              ├── YES ──► Final: TRUE
                              └── NO  ──► [Adaptive 2nd Iteration: rerun at second_iter_factor * bandwidth] ──► Final: FALSE
```

## Dependencies
This module is verified on Python 3.11 – 3.13 and relies on numpy, pandas, astropy and matplotlib as well as reference catalog file for streaked observations from the Chandra Data Archive - 'cda_flagged_list.txt' (included) **to be put in the same repository as the module's .py file**.


## Function Parameters
The module's primary function is:

find_streaks(

event_path,

plot_path, 

bandwidth=50, 

n1=10, 

n2=4, 

second_iter_factor=0.8, 

second=False

)

### Quickstart

Note that the only two positional arguments, not set to a default value are the 'event_path' argument that recieves the path for a .fits file and the 'plot_path' argument for the desired path for the output plots to be saved. Other argument for the find_streaks function are set to default value that have been tested as typical for a standard run.

### Parameters

**bandwidth** - Width of a single band along the chipx-parallel axis in sky pixels.

**n1** - Threshold distance of the maximum band's photon count from the mean count per band, measured in standard deviations ($\sigma$) for whole bands.

**n2** - Outlier threshold ($\sigma$) required independently on both chip halves during the second (split-band) verification.

**second_iter_factor** - Multiplier applied to bandwith if a second iteration is needed to rule out cases where the band's border runs exactly along the streak.

**second** - An internal recursion tracker, leave as False unless you wish to avoid a second iteration.


## Output
For a given observation, the module saves and returns a verdict_obsid_<OBSID>.csv table with the following columns:

**ccd (int)** - ACIS chip identifier

**is_streak (bool)** - Final verdict regarding suspected streak presence in the CCD.

**sigma_distance (float)** - Outlier distance of the maximum band in standard deviations, assessed against $n_1 \cdot \sigma$.

**sigma_distance_side1/2 (float)** -  Outlier distance in standard deviations on each half of the chip, assessed against $n_2 \cdot \sigma$.

**mean_photon_count (float)** - Mean background count across primary slices (excluding candidate peak).

**mean_photon_count_side1/2 (float)** - Mean background count on secondary half-bands for each side (excluding the candidate peak).

**labeled_as_streaked_in_cda (bool)** - Flag matching with cases allready flagged as streaked in CDA.

**The verdict DataFrame is saved to disk as well as directly returned by the function**
Aditionally, two figures are saved in the plot_path repository:

**`<OBSID>_first_check_ccd_<CCD>.png`** - Diagnostics from the primary check, containing three subplots from left to right:
1. Histogram of photon counts - the background mean is marked with a vertical dashed line, and the sigma_distance value presented in red.
2. Photon count per band - 1D slice profile of photon counts across the chip, the peak band is marked with a dashed red vertical line and the projection index specified in the bottom axis.
3. spatial map of the photons on chip, in [x,y] coordinates, including the bands' borders in purple and the candidate band marked in red.

example:

<img width="1800" height="600" alt="3956_first_check_ccd_7" src="https://github.com/user-attachments/assets/397ec792-8c81-46c8-8476-3c809c2c3fa2" />

**`<OBSID>_second_check_ccd_<CCD>.png`** - Split-band check diagnostics for each half of the detector separately (Side 1 on the left, Side 2 on the right). Displays band count profiles on top, count distribution histograms in the middle, and a spatial map with split-band boundaries on the bottom.

example:

<img width="1200" height="1800" alt="3956_second_check_7" src="https://github.com/user-attachments/assets/907ab7d2-b37a-40a1-bded-d548034442d6" />



The module creates a folder hierarchy under plot_path as such:
```text
plot_path/
└── obsid_<OBSID>/
    ├── verdict_obsid_<OBSID>.csv
    ├── ccd_0/
    │   ├── <OBSID>_first_check_ccd_<CCD>.png
    │   └── <OBSID>_second_check_ccd_<CCD>.png   (only if Check 1 passed)
    ├── ccd_1/
    └── ...
```



