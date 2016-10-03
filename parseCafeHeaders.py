"""
Rename this script when finished!
Goal is to look through all the new data and check
what was observed. Look for quads and count spectra etc.

First step: Get all spectra associated with a swasp_id
Update the OBJECT keyword if necessary

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
import argparse as ap
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

def argParse():
    description = """
                Script to summarise CAFE data headers and
                make some small modifications for reduction
                """
    parser = ap.ArgumentParser(description=description)
    parser.add_argument('--update_object',
                        help='update the object keyword with swasp_id',
                        action='store_true')
    parser.add_argument('--reffile',
                        help='output a per night reffile.txt file',
                        action='store_true')
    return parser.parse_args()

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


def matchRvStandardId(object_id, image_id):
    """
    Get RV standard info from database
    """
    qry="""
        SELECT ra_hms,
               dec_dms,
               pm_ra_mas_y,
               pm_dec_mas_y,
               spec_type
        FROM rv_standards
        WHERE object_id = '{}'
        """.format(object_id)
    ra, dec, pra, pdec, spec_type = [], [], [], [], []
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            ra.append(row[0])
            dec.append(row[1])
            pra.append(float(row[2]))
            pdec.append(float(row[3]))
            spec_type.append(row[4])
    if len(ra) == 0:
        print('{} - WARNING: NO MATCH FOR {}'.format(image_id, object_id))
        return None, None, None, None, None
    else:
        return ra[0], dec[0], pra[0], pdec[0], spec_type[0]

def matchSwaspId(object_id, image_id):
    """
    Match the object ID to the swasp_id
    """
    qry="""
        SELECT swasp_id,
               epoch,
               period,
               paramfit_spec_type,
               ra_hms,
               dec_dms,
               pm_ra_mas_y,
               pm_dec_mas_y
        FROM eblm_parameters
        WHERE swasp_id LIKE '%{}%'
        """.format(object_id)
    swasp_id, epoch, period, spec_type  = [], [], [], []
    ra, dec, pra, pdec = [], [], [], []
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            swasp_id.append(row[0])
            epoch.append(float(row[1]))
            period.append(float(row[2]))
            spec_type.append(row[3])
            ra.append(row[4])
            dec.append(row[5])
            pra.append(float(row[6]))
            pdec.append(float(row[7]))
    if len(swasp_id) == 0:
        print('{} - WARNING: NO MATCH FOR {}'.format(image_id, object_id))
        return None, None, None, None, None, None, None, None
    elif len(swasp_id) > 1:
        print('{} - WARNING: MULTIPLE MATCHES FOR {}'.format(image_id, object_id))
        i = 0
        for s, e, p, st in zip(swasp_id, epoch, period, spec_type):
            print(i, s, round(e, 6), round(p, 6), st)
            i += 1
        choice = 1E6
        while choice > len(swasp_id) or choice < 0:
            choice = int(input('Pick a row for the right object: '))
        return swasp_id[choice], epoch[choice], period[choice], spec_type[choice], ra[choice], dec[choice], pra[choice], pdec[choice]
    else:
        print('{} - {} matches {} {} {} {}'.format(image_id,
                                                   object_id,
                                                   swasp_id[0],
                                                   round(epoch[0], 6),
                                                   round(period[0], 6),
                                                   spec_type[0]))
        return swasp_id[0], epoch[0], period[0], spec_type[0], ra[0], dec[0], pra[0], pdec[0]

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

def getCcfMask(st):
    """
    Determine best CCF mask for RVs
    based on stars spectral type
    """
    if st.lower().startswith('f') or st.lower().startswith('g'):
        return 'G2'
    elif st.lower().startswith('k'):
        return 'K5'
    elif st.lower().startswith('m'):
        return 'M2'
    else:
        return None

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
    args =argParse()
    # defdict to store results
    observed_times = defaultdict(list)
    spec_types, coords = {}, {}
    prop_motions, ccf_masks = {}, {}
    # loop over the directories and get all info at once
    data_directories = g.glob('160*')

    print('-------------------------')
    print('----PER NIGHT SUMMARY----')
    print('-------------------------')

    for directory in data_directories:
        print('\nMoving to {}\n'.format(directory))
        with cd(directory):
            # get a list of spectra in the night's directory
            templist = g.glob('*.fits')
            # keep a reffile per directory
            # this will be used by CERES to reduce each night's data
            reffile = []
            for image in templist:
                # grab a bunch of important header keywords
                h = fits.open(image)[0].header
                object_id = h['OBJECT']
                dateobs = h['DATE']
                exptime = float(h['EXPTIME'])
                imagetyp = h['IMAGETYP']
                ra = float(h['RA'])/3600.
                dec = float(h['DEC'])/3600.
                # only look for matching swasp objects for phasing
                if imagetyp == 'Science' and not object_id.startswith('HD'):
                    # match swasp_ids
                    swasp_id, epoch, period, spec_type, ra, dec, pra, pdec = matchSwaspId(object_id, image)
                    # get the phases observed
                    phase = getPhaseObserved(dateobs, ra, dec, epoch, period)
                    print('{0:s} {1:s} Phase = {2:.2f}'.format(image, dateobs, phase))
                    observed_times[swasp_id].append([image, dateobs, round(phase, 2)])
                    if swasp_id not in spec_types:
                        spec_types[swasp_id] = spec_type
                    if swasp_id not in coords:
                        coords[swasp_id] = '{},{}'.format(ra, dec)
                    if swasp_id not in prop_motions:
                        prop_motions[swasp_id] = '{},{}'.format(pra, pdec)
                    if swasp_id not in ccf_masks:
                        ccf_masks[swasp_id] = getCcfMask(spec_type)

                    # append the frame info to the reffile
                    reffile.append("{},{},{},1,{},4".format(swasp_id,
                                                            coords[swasp_id],
                                                            prop_motions[swasp_id],
                                                            ccf_masks[swasp_id]))

                    # update the OBJECT keyword with the swasp_id
                    if args.update_object:
                        data = fits.getdata(image)
                        h['OBJECT'] = swasp_id
                        print('{} {} --> {}'.format(image, object_id, swasp_id))
                        fits.writeto(image, data, header=h, clobber=True)

                # add the standards to the reffile also
                if imagetyp == 'Science' and object_id.startswith('HD'):
                    ra, dec, pra, pdec, spec_type = matchRvStandardId(object_id, image)
                    if object_id not in spec_types:
                        spec_types[object_id] = spec_type
                    if object_id not in coords:
                        coords[object_id] = '{},{}'.format(ra, dec)
                    if object_id not in prop_motions:
                        prop_motions[object_id] = '{},{}'.format(pra, pdec)
                    if object_id not in ccf_masks:
                        ccf_masks[object_id] = getCcfMask(spec_type)
                    reffile.append("{},{},{},1,{},4".format(object_id,
                                                           coords[object_id],
                                                           prop_motions[object_id],
                                                           ccf_masks[object_id]))

            # print the reffile
            print('\n')
            for line in set(reffile):
                print(line)
            # output the reffile if required
            if args.reffile:
                outfile = open('reffile.txt','w')
                for line in set(reffile):
                    outfile.write('{}\n'.format(line))
                outfile.close()

    # print a per object summary
    # the top line is the line needed for reffile
    # for this object
    print('\n--------------------------')
    print('----PER OBJECT SUMMARY----')
    print('--------------------------')
    for obj in sorted(observed_times.keys()):
        print("\n{},{},{},{},1,{},4:".format(obj, coords[obj], prop_motions[obj], spec_types[obj], ccf_masks[obj]))
        for spec in sorted(observed_times[obj]):
            print(spec)

