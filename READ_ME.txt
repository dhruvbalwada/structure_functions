----------------------------------------------------------------------
GLAD 15-minuted interpolated, low-pass filtered drifter trajectories
----------------------------------------------------------------------
Time period:  July 20 through October 22, 2012
----------------------------------------------------------------------
ASCII data file:  GLAD_15min_filtered.dat  (155 MB)
----------------------------------------------------------------------
Contact:  Ed Ryan, RSMAS, University of Miami (eryan@rsmas.miami.edu)
----------------------------------------------------------------------
This dataset was created by the Consortium for Advanced Research on
Transport of Hydrocarbon in the Environment (CARTHE). This research was
made possible by a grant from BP/The Gulf of Mexico Research Initiative.
----------------------------------------------------------------------
ASCII file format
----------------------------------------------------------------------

Column 1:  drifter ID string (like CARTHE_XXX, XXX=3 digit integer)
Column 2:  date (yyyy-mm-dd)
Column 3:  time (HH:MM:SS.SS)
Column 4:  latitude (decimal degrees)
Column 5:  longitude (decimal degrees)
Column 6:  estimated position error (meters) 
Column 7:  u (east-west) velocity (m/sec)
Column 8:  v (north-south) velocity (m/sec)
Column 9:  estimated velocity error (m/sec)

----------------------------------------------------------------------
Summary
----------------------------------------------------------------------

297 trajectories from near-surface CODE-type ocean drifters (drogued at a
depth of one meter) tracked in real-time using SPOT GPS units, launched in
the northern Gulf of Mexico near DeSoto Canyon in July 2012 as part of the
CARTHE Grand Lagrangian Deployment (GLAD) experiment.  Most of these
drifters were launched as triplets (separated by roughly 100 meters at
launch) in an attempt to measure multi-scale near surface dispersion.

Positions are low-pass filtered (one hour cutoff period) and interpolated to
uniform 15 minute intervals starting on whole hours over the period July 20
through October 22, 2012.  No temperature or salinity sensors were attached
to these drifters.

----------------------------------------------------------------------
How these trajectories were processed:
----------------------------------------------------------------------

GLAD CODE drifter positions were reported in real-time roughly every five
minutes via the Globalstar satellite network.  Positions reported by the
SPOT handheld GPS units inside each drifter were nominally accurate to
within seven meters. Occasionally, position errors much larger than this
nominal value are found.  These may be due to errors in GPS baselines (due
to poor satellite reception), data dropout, and, in some cases, drifters
that were picked up by small boats.

Each drifter record ends when the drifter was known to be picked up by a
boat, when the signal was lost for more than 24 hours, or when the drifter
traveled more than 80 km in a 12 hours period (implying a mean speed of 1.85
m/sec over 12 hours). 

For each record, positions that imply an instantaneous drifter
speed greater than 3 m/sec were deleted.  Also, positions that imply the
drifter track rotated through more than 360 compass degrees within three
hours were deleted.

Next, outliers were identified as positions that were more than 100 meters
away from estimated positions at the same times interpolated from a set of
both past and future positions.  These outliers were deleted.

All valid positions were then interpolated to uniform, five-minute time
intervals using spline interpolation.  These five-minute records were used
to compute finite difference estimates of u (east-west) and v (north-south)
velocities.

Finally, the five-minute position and velocity records were filtered
using a Butterworth fourth-order low-pass filter with a one hour period
cutoff.  These low-pass filtered records were then interpolated to uniform
15-minute intervals beginning on whole hours.  

Crude estimates of position error are also provided for each position in the
final 15-minute interval records. The error was set to 10 meters (a nominal
value) for positions recorded at times for which there were no gaps in the
five-minute raw data records over the previous hour. The error was increased
according to the ratio of the average sample time interval for the previous
12 raw data positions to the nominal sample interval of 5 minutes,
multiplied by the 10 meter nominal error value.  Velocity errors were
computed by dividing the sum of two consecutive position errors by the
corresponding time interval.  
