CAFE Observations Preparation
-----------------------------

Code to prepare for CAFE observations

CAFE-SNR.py estimates the required exposure times and total
observing time including overheads for obvservations with the
CAFE spectrogrpah at Calar Alto, Spain.

mkfinders.py makes 3' by 3' finding charts (N up, E left)
using pyds9 - this should be replaced with Paul's new gnerator

![CAFE SNR](img/cafe-snr.png)

CAFE Data Reduction
-------------------

My ristretto.py code has been retired before completion since the release of CERES.
All reductions for CAFE should be done using CERES. See the CERES repo README for
instructions on how to run it.

CAFE Housekeeping Scripts
-------------------------

parseCafeHeaders.py - a script to check what targets have been observed

parseCafeSkyApertures.py - a script to analyse manually placed sky apertures for automating the procedure


Contributors
------------

James McCormac
