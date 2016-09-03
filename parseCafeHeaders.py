"""
Rename this script when finished!
Goal is to look through all the new data and check
what was observed. Look for quads and count spectra etc.

First step: Get all spectra associated with a swasp_id
Second step: Then start populating the database table

CAFE headers are pretty sparse and in a funny format.
The real values seem to be in the comments after the
actual header value but can be accessed from the value
if converted correctly.

This script shows how to get useful info from the headers
when checking what has been done during a run. Below is
a key showing the most important header entries

Key:
    IMAGETYP: IMAGETYP
        Science
        Bias
        Flat
        Calibration
    UTSTART:  UT_START (format?)
    EXPTIME:  EXPTIME (sec)
    RA:       RA (arcsec)
    DEC:      DEC (arcsec)
    AIRMASS:  AIRMASS
    OBJECT:   OBJECT
        Name
        [arc] ThAr
        [flat] Hal
        [Bias]
    LAT:      37.2246
    LON:      -2.5463
    ELV:      2168

Time system:
    DATE - (EXPTIME/2)
"""
import glob as g
from astropy.io import fits

t=g.glob('*.fits')
for i in t:
    h=fits.open(i)[0].header
    if h['IMAGETYP'] == 'Science':
        line="{} {} {} {} {} {} {} {}".format(i,
                                              h['IMAGETYP'],
                                              h['OBJECT'],
                                              h['RA'],
                                              h['DEC'],
                                              h['UT_START'],
                                              h['EXPTIME'],
                                              h['AIRMASS'])
        print(line)
