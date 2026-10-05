# find_streaks
find_streaks is a python module designed to detect readout and pileup streak contamination in Chandra X-ray Observatory ACIS observations.
Given a Level-2 FITS event file (`acisf..._evt2.fits`), the pipeline searches every active CCD chip for streaks, matches results with archival Chandra Data Archive (CDA) records, saves diagnostic plots, and outputs a summary CSV table.
The module includes several utility functions but find_streaks is the main caller function that initiates the entire logical sequence.


# Method of conclusion
The find_streaks logic uses a spatial two tier test on each CCD to flag suspected streaks. The calculations are done on the event file's [x,y] coordinates after the dithering is corrected, which concentrates the streaks into a more consistant and clear line.

**Orientation & Coordinate Transformation:**
First, the borders of the CCD are estimated and two non parallel vector, representing the two axis of the borders are chosen. Streaks in the data show up along the columns of the pixels on the CCD - along the axis of the chipy and so the band slicing will be done along the border axis coinciding with the chipx axis.
The slicing itself is done by calculating the norm of the projection of each point on the cosen axis vector (the point - as a vector with a shared origin with the axis vector). The points are then devided to bands according to the value of their projection's norm.

**Primary Check**
The band containing a maximum photon count is chosen and is asessed againts the mean photon count per band - if it exeedes $n_1 \cdot \sigma$, sigma being the standart deviation, it is flagged as a canidate and passes on to the secondary check.

**Secondary Check**
In order to rule out maximas caused by bright sources, the secondary check is performed oved two halves of the image along the readout direction and the same asessment is performed independantly along the two halves, according to a second threshold n2. Both halves must independently show a statistically significant peak ($> n_2 \cdot \sigma$) at the identical projection index.

**Adaptive Second Iteration**
Streaks may fall on the border between two bands resulting in a false negative. In order to rule that out a second iteration with an augmented bandwidth has been implemented. The second iteration triggers if one of the two former checks fails, and will only be performed once. The bandwidth will be multipled by a factor - the function's second_iter_factor argument which is defaulted to 0.8.

Flow chart describing the function's decision tree:

[Raw FITS Events for CCD]
            │
            ▼
1. Geometric Alignment (Estimate CCD corners & match CHIPX/CHIPY readout axes)
            │
            ▼
2. Spatial Projection (Bin photon coordinates into slices of width = bandwidth)
            │
            ▼
3. Primary Check: Is there an anomalous count peak?
   (Peak band count > background mean + n1 * sigma)
            │
            ├── NO  ──► [Adaptive 2nd Pass: rerun at second_iter_factor * bandwidth] ──► Final: False
            │
            └── YES ──► 4. Bilateral Check: Does the peak cross the ENTIRE chip?
                        (Bisect chip along the orthogonal axis into Side 0 & Side 1)
                        Do BOTH halves share a coincident peak > n2 * sigma?
                              │
                              ├── YES ──► Final: TRUE
                              └── NO  ──► [Adaptive 2nd Pass: rerun at second_iter_factor * bandwidth] ──► Final: FALSE


# Dependencies
This module is verified on Python 3.11 – 3.13 and relies on numpy, pandas, astropy and matplotlib as well as reference catalog file for streaked observations from the Chandra Data Archive - "cda_flagged_list.txt" (included) to be put in the same repository as the module's .py file.


# Function Parameters
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

***Quickstart***

Note that the only two positional arguments, not set to a default value are the 'event_path' argument that recieves the path for a .fits file and the 'plot_path' argument for the desired path for the output plots to be saved. Other argument for the find_streaks function are set to default value that have been tested as typical for a standart run.

**bandwidth** - Slice width of the image along the chipy-parallel axis in sky pixels.

**n1** - The threshold distance of the maximal band's photon count from the image's mean photon count per band, measured in number of standart deviations. This outlier is measured for the whole bands in the first asessment.

**n2** - A similar outlier for the second assesment on both chip halves during the second, split band verification.

**second_iter_factor** - Multiplier applied to bandwith if a second iteration of the function is needed to rule out cases where the band's border runs exactly along the streak.

**second** - An internal recursion tracker, leave as False unless you wish to avoid a second iteration.


# Output
For a given observation, the module outputs a verdict_obsid_<OBSID>.csv table with the following columns:

**ccd (int)** - ACIS chip identifier

**is_streak (bool)** - Function's verdict regarding suspection of streak presence in the image

**sigma_distance (flaot)** - The distance in number of standart deviations of the maximum band that is asessed againts $n_1 \cdot \sigma$.

**sigma_distance_side1/2 (flaot)** - The distance in number of standart deviations of the maximum half band in each side that is asessed againts $n_2 \cdot \sigma$.

**mean_photon_count (flaot)** - Mean background count across primary slices (excluding candidate peak).

**mean_photon_count_side1/2 (flaot)** - Background mean photon count on secondary, half-bands on each side (excluding canidate peak).

**labeled_as_streaked_in_cda (bool)** - Flag matching with cases allready flagged in cda.

Aditionally, two figures are saved in the plot_path repository:

**<OBSID>_first_check_ccd_<CCD>.png** - Diagnostics from the primary check, containing three subplots from left to right:
1. Histogram of photon counts - the background mean is marked with a vertical dashed line, and the sigma_distance value presented in red.
2. Photon count per band - 1D slice profile across the chip, the peak band is marked with a dashed red vertical line.
3. Spacial map of the photons on chip, in [x,y] coordinates, including the bands' borders in purple and the canidate band marked in red.

**<OBSID>_second_check_ccd_<CCD>.png** - The same plots are presented for each half of the data seperately - side 1 on the left and side 2 on the right. The photon graph per band is on top, Hisograms below it, and a spacial map with the spilt bands' borders on the bottom:

The module creates a folder hierarchy under plot_path as such:
screening_output/
└── obsid_<OBSID>/
    ├── verdict_obsid_<OBSID>.csv
    ├── ccd 0/
    │   ├── <OBSID> first check.png
    │   └── <OBSID> second check.png   (only if Check 1 passed)
    ├── ccd 1/
    └── ...





