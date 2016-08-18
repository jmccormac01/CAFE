import os
import time
import pyds9
import pymysql 
import argparse as ap
from astropy.coordinates import SkyCoord
import astropy.units as u

# connect to the database
db = pymysql.connect(host='localhost', db='eblm', password='mysqlpassword')

def argParse():
    description = "Code to make find charts for CAFE obs"
    status_choices = ["OBSERVE",
                    "CONTINUE",
                    "CONTINUE-CAFE",
                    "CONTINUE-STABILIZED",
                    "PHASE_COVERAGE",
                    "ANALYSE_BLENDS"]
    parser = ap.ArgumentParser(description=description)
    parser.add_argument('status_flag',
                        choices = status_choices,
                        help='Status flag to ID objects in DB')
    return parser.parse_args()

def getTargetsRaDec(status):
    objects = {}
    qry = """
        SELECT
        swasp_id
        FROM
        eblm_parameters
        WHERE
        current_status = '{0:s}'
        """.format(status)
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            obj = row[0]
            ra = "{0:s}:{1:s}:{2:s}".format(obj[7:9],obj[9:11],obj[11:16])
            dec = "{0:s}:{1:s}:{2:s}".format(obj[16:19],obj[19:21],obj[21:25])
            objects[obj] = (ra, dec)
    return objects

if __name__ == "__main__":
    args = argParse()
    objects = getTargetsRaDec(args.status_flag)
    d = pyds9.DS9()
    time.sleep(5)

    outdir = '/Users/James/Dropbox/Observing/INTObs/EBLM/201608_10/finders/'

    for obj in sorted(objects):
        outfile = '{0:s}{1:s}.jpeg'.format(outdir, obj)
        while not os.path.exists(outfile):
            ra = objects[obj][0]
            dec = objects[obj][1]
            c= SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg), frame='icrs')
            print(obj, objects[obj])
            d.set('frame clear all')
            d.set('frame 1')
            d.set('dsseso size 3 3 arcmin')
            d.set('dsseso coord {0:s} {1:s}'.format(ra, dec))
            d.set('frame center')
            d.set('zoom 4')
            d.set('regions system wcs')
            d.set('regions format ds9')
            r = 'fk5; circle({0:.6f},{1:.6f},30")'.format(c.ra.deg, c.dec.deg)
            print(r)
            # BODGEEEEE!!!!! *facepalm*
            d.set('regions', 'fk5; circle({0:.6f},{1:.6f},30")'.format(c.ra.deg, c.dec.deg))
            try:
                d.set('saveimage {0:s} 75'.format(outfile))
            except ValueError:
                print('Failed to save {}, trying again...'.format(outfile))
                try:
                    d.set('saveimage {0:s} 75'.format(outfile))
                except ValueError:
                    print('Failed to save {} again, starting over...'.format(outfile))
                    continue
                # x=input()
        else:
            print('Found finder for {0:s}'.format(obj))
