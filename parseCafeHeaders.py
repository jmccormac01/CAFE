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
import os
import glob as g
from astropy.io import fits
import pymysql
from collections import defaultdict
from contextlib import contextmanager
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import (
    EarthLocation,
    SkyCoord
    )

# connect to database
db = pymysql.connect(host='localhost', db='eblm', password='mysqlpassword')

# observatory location
OLAT = 37.2236
OLON = -2.5463
ELEV = 2168.
OBSERVATORY = EarthLocation(lat=OLAT*u.deg, lon=OLON*u.deg, height=ELEV*u.m)

def getLightTravelTimes(ra, dec, time2corr):
    """
    Get the light travel times to the helio and
    bary centres

    TODO: Finish docstring
    """
    target = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    ltt_bary = time2corr.light_travel_time(target)
    ltt_helio = time2corr.light_travel_time(target, 'heliocentric')
    return ltt_bary, ltt_helio

def matchSwaspId(object_id, image_id):
    """
    Match the object ID to the swasp_id
    """
    qry="""
        SELECT swasp_id, epoch, period
        FROM eblm_parameters
        WHERE swasp_id LIKE '%{}%'
        """.format(object_id)
    swasp_id, epoch, period = [], [], []
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            swasp_id.append(row[0])
            epoch.append(float(row[1]))
            period.append(float(row[2]))
    if len(swasp_id) == 0:
        print('{} - WARNING: NO MATCH FOR {}'.format(image_id, object_id))
        return None, None, None
    elif len(swasp_id) > 1:
        print('{} - WARNING: MULTIPLE MATCHES FOR {}'.format(image_id, object_id))
        i = 0
        for s, e, p in zip(swasp_id, epoch, period):
            print(i, s, e, p)
            i += 1
        choice = 1E6
        while choice > len(swasp_id) or choice < 0:
            choice = int(input('Pick a row for the right object: '))
        return swasp_id[choice], epoch[choice], period[choice]
    else:
        print('{} - {} matches {} {} {}'.format(image_id,
                                                object_id,
                                                swasp_id[0],
                                                epoch[0],
                                                period[0]))
        return swasp_id[0], epoch[0], period[0]

def getPhaseObserved(dateobs, ra, dec, epoch, period):
    """
    Determine the phase of each observation
    """
    jd_mid = Time(dateobs, scale='utc', format='isot', location=OBSERVATORY) - (exptime*u.second/2.)
    ltt_bary, ltt_helio = getLightTravelTimes(ra, dec, jd_mid)
    time_bary = jd_mid.tdb + ltt_bary
    time_helio = jd_mid.utc + ltt_helio
    phase = ((time_helio.jd - (epoch+2450000))/period)%1
    return phase

@contextmanager
def cd(path):
    """
    Improved os.chdir()
    """
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)

if __name__ == '__main__':
    # defdict to store results
    observed_times = defaultdict(list)
    # loop over the directories and get all info at once
    data_directories = g.glob('1607*')
    for directory in data_directories:
        print('\nMoving to {}\n'.format(directory))
        with cd(directory):
            t = g.glob('*.fits')
            for i in t:
                h=fits.open(i)[0].header
                object_id = h['OBJECT']
                dateobs = h['DATE']
                exptime = float(h['EXPTIME'])
                imagetyp = h['IMAGETYP']
                ra = float(h['RA'])/3600.
                dec = float(h['DEC'])/3600.
                if imagetyp == 'Science' and not object_id.startswith('HD'):
                    # match swasp_ids
                    swasp_id, epoch, period = matchSwaspId(object_id, i)
                    # get the phases observed
                    phase = getPhaseObserved(dateobs, ra, dec, epoch, period)
                    print('{0:s} {1:s} Phase = {2:.2f}'.format(i, dateobs, phase))
                    observed_times[swasp_id].append([i,dateobs,round(phase,2)])
    # print a summary
    for obj in sorted(observed_times.keys()):
        print("\n{}:".format(obj))
        for spec in sorted(observed_times[obj]):
            print(spec)

