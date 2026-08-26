# streak-detection
This function takes an event file and returns True/False in accordance to wether or not it concluded that a pileup streak is present within the data.
Function header is:
find_streaks(event,bandwidth,n1,n2)


**Function Arguments**


**event (pd.DataFrame)** - this is the argument for the eventfile in the form of a pandas dataframe


**bandwitdth (int)** - the desired bandwidth by which the search is conducted, see search method


**n1 (int)** - the distance, in number of standart deviations, required to define a band that is suspected of containing a streak, counting number of photons in each band.


**n2 (int)** - similar to n2 but for the second check that is run on two sides of the same split bands, split along the orthogonal axis.


**Method of conclusion**


1. The borders of the CCD are estimated - and two non parallel vector representing two borders are chosen
2. Streaks in the data - wether caused by pileup or by readout error, show up along the columns of the pixels on the CCD - along the axis of the chipy. For that reason we choose the vector that is parallel to the chipx axis, on which we devide the data into bands parallel to the chipy axis.
3. The projection of each point upon each side is calculated, and the points are devided into bands according to the size of their projection on the side parallel to the chipx axis
4. The bands containing the absolute maximal number of photons is marked as suspected of streak
5. **test 1:**
The standart deviation of the photon counts of the bands, excluding the suspected one is calculated. The condition upon the band is weighed upon is:
Wether or not the number of photons in the suspected band is larger than the mean (excluding the suspected band) in more that n1 times the std. If not the function returns False.
6.**test 2:**
Is conducted if test 1 is passed, the bands are split into two sides along the orthogonal axis to that of the bands and the above calculation is done for each of the two sides seperately, over a distance of n2 times the std of each side. If the condition is met on both sides - for the same bands as the one that was suspected earlier, the function will return True.


