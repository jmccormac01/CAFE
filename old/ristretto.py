"""
Script to reduce data from the CAFE spectrograph
"""
import sys
import os
import time
import glob as g
from contextlib import contextmanager
from ccdproc import (
    CCDData,
    ImageFileCollection,
    combine,
    subtract_bias,
    flat_correct
    )
from astropy.io import fits
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    SkyCoord
    )

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

# pylint: disable=invalid-name
# pylint: disable=redefined-builtin
# pylint: disable=no-member
# pylint: disable=redefined-outer-name
# pylint: disable=superfluous-parens

###############################
########## EDIT HERE ##########
###############################

# file locations
LINELIST = "/Users/James/Documents/LineLists/thar.dat"
LOGINCL = "/Users/James/"
# filenames
MASTER_FLAT = 'master_flat.fits'
MASTER_FLAT_FLAT = 'master_flat_flattened.fits'
MASTER_BIAS = 'master_bias.fits'
# ccd params
GAIN = 'CCDGAIN'
RDNOISE = 'CCDRON'
SATURATION = '65000'
DISPAXIS = 1
# image filename extension
# use wildcard to find them all
IMAGE_EXTENSION = '*.fits'
# header keyword map
BIAS_KEYWORD = 'Bias'
FLAT_KEYWORD = 'Flat'
SCIENCE_KEYWORD = 'Science'
ARC_KEYWORD = 'Calibration'
RA_KEYWORD = 'RA'
DEC_KEYWORD = 'DEC'
DATEOBS_KEYWORD = 'DATE'
EXPTIME_KEYWORD = 'EXPTIME'
OBJECT_KEYWORD = 'OBJECT'
# observatory location
OLAT = 37.2236
OLON = -2.5463
ELEV = 2168.
OBSERVATORY = EarthLocation(lat=OLAT*u.deg, lon=OLON*u.deg, height=ELEV*u.m)

##############################

def makeDirs():
    """
    Make directories for original and reduced data
    """
    if not os.path.isdir('original'):
        os.mkdir('original')
    if not os.path.isdir('reduced'):
        os.mkdir('reduced')

def copyFiles():
    """
    Copy the data
    """
    currentdir = os.getcwd()
    contents = os.listdir(currentdir)
    if not os.path.isdir('original'):
        os.mkdir('original')
    else:
        print('original directory exists, it shouldn\'t, exiting...')
        sys.exit(1)
    for i in contents:
        os.system('mv {} original/'.format(i))
        print("Moved {} to original folder".format(i))
    if not os.path.isdir('reduced'):
        os.mkdir('reduced')
    else:
        print('reduced directory exists, it shouldn\'t, exiting...')
        sys.exit(1)
    os.chdir('original/')
    for i in g.glob(IMAGE_EXTENSION):
        os.system("cp {} ../reduced/".format(i))
        print("Copied {} to reduction folder".format(i))
    os.chdir('../reduced/')

def getImageList():
    """
    Return ImageFileCollection object for the current
    working directory
    """
    return ImageFileCollection('.')

def checkHeader(filename, hdr):
    """
    Add DISPAXIS to the header object if it is missing
    """
    try:
        hdr['DISPAXIS']
    except KeyError:
        print('{}: Adding DISPAXIS = {}'.format(filename, DISPAXIS))
        hdr['DISPAXIS'] = DISPAXIS
    return hdr

def trimFiles():
    """
    CAFE data contains a useless region along the first
    120 columns, kill it.

    Update any missing header keywords at before writing file
    """
    t = g.glob(IMAGE_EXTENSION)
    for filename in t:
        print('Trimming image {0:s} to {0:s}[:,120:]'.format(filename))
        data, hdr = fits.getdata(filename, header=True)
        hdr = checkHeader(filename, hdr)
        fits.writeto(filename, data[:, 120:], header=hdr, clobber=True)

def renameFiles(images):
    """
    Rename the files so they are easier to scan by eye
    """
    for f in images.files_filtered(imagetyp=FLAT_KEYWORD):
        os.rename(f, "f_s_{}".format(f))
    for f in images.files_filtered(imagetyp=BIAS_KEYWORD):
        os.rename(f, "b_s_{}".format(f))
    for f in images.files_filtered(imagetyp=ARC_KEYWORD):
        os.rename(f, "a_s_{}".format(f))
    for f in images.files_filtered(imagetyp=SCIENCE_KEYWORD):
        os.rename(f, "i_s_{}".format(f))

def makeMasterBias(images):
    """
    Make a master bias using all biases found in
    images object

    TODO: Finish docstring
    """
    try:
        master_bias = CCDData.read(MASTER_BIAS, unit=u.adu)
        return master_bias
    except FileNotFoundError:
        bias_list = []
        for f in images.files_filtered(imagetyp=BIAS_KEYWORD):
            print(f)
            ccd = CCDData.read(f, unit=u.adu)
            bias_list.append(ccd)
        try:
            master_bias = combine(bias_list, method='median')
            master_bias.write(MASTER_BIAS, clobber=True)
            return master_bias
        except IndexError:
            return None

def subtractMasterBias(images, master_bias):
    """
    Subtract the master bias from all the science and
    calibration data

    TODO: Finish docstring
    """
    for filename in images.files_filtered(imagetyp=SCIENCE_KEYWORD):
        ccd = CCDData.read(filename, unit=u.adu)
        ccd = subtract_bias(ccd, master_bias)
    for filename in images.files_filtered(imagetyp=ARC_KEYWORD):
        ccd = CCDData.read(filename, unit=u.adu)
        ccd = subtract_bias(ccd, master_bias)
    for filename in images.files_filtered(imagetyp=FLAT_KEYWORD):
        ccd = CCDData.read(filename, unit=u.adu)
        ccd = subtract_bias(ccd, master_bias)

def makeMasterFlat(images):
    """
    Make a master flat from the already bias corrected flats

    TODO: Finish docstring
    """
    try:
        master_flat = CCDData.read(MASTER_FLAT, unit=u.adu)
        return master_flat
    except FileNotFoundError:
        # empty list for the flats
        flat_list = []
        # create the master flat field
        print('Combining flats')
        for f in images.files_filtered(imagetyp=FLAT_KEYWORD):
            print(f)
            ccd = CCDData.read(f, unit=u.adu)
            flat_list.append(ccd)
        try:
            master_flat = combine(flat_list, method='median')
            master_flat.write(MASTER_FLAT, clobber=True)
            return master_flat
        except IndexError:
            print('There are no flats, skipping...')
            master_flat = None

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

def importIraf():
    """
    Import IRAF and return it for others to use

    Echelle reductions use the noao>imred>echelle
    packages
    """
    with cd(LOGINCL):
        from pyraf import iraf
        iraf.noao(_doprint=0)
        iraf.imred(_doprint=0)
        iraf.echelle(_doprint=0)
        time.sleep(2)
    return iraf

def identifyOrders(iraf, images):
    """
    Function to identify the location of the echelle orders
    """
    # set up apall for echelle order identification
    iraf.apall.setParam('apertur', '')
    iraf.apall.setParam('format', 'echelle')
    iraf.apall.setParam('reference', '')
    iraf.apall.setParam('profile', '')
    # operation mode
    iraf.apall.setParam('interac', 'yes')
    iraf.apall.setParam('find', 'yes')
    iraf.apall.setParam('recen', 'no')
    iraf.apall.setParam('resize', 'yes')
    iraf.apall.setParam('edit', 'yes')
    iraf.apall.setParam('trace', 'yes')
    iraf.apall.setParam('fittrac', 'yes')
    iraf.apall.setParam('extract', 'no')
    iraf.apall.setParam('extras', 'no')
    iraf.apall.setParam('review', 'no')
    iraf.apall.setParam('line', '904')
    iraf.apall.setParam('nsum', '5')
    # default aperture parameters
    iraf.apall.setParam('lower', '-5')
    iraf.apall.setParam('upper', '5')
    # default background parameters
    iraf.apall.setParam('b_funct', 'chebyshev')
    iraf.apall.setParam('b_order', '2')
    iraf.apall.setParam('b_sampl', '-9:-7,5:7')
    iraf.apall.setParam('b_naver', '-3')
    iraf.apall.setParam('b_niter', '3')
    iraf.apall.setParam('b_low_r', '3')
    iraf.apall.setParam('b_high', '3')
    iraf.apall.setParam('b_grow', '0')
    # aperture centering parameters
    iraf.apall.setParam('width', '5')
    iraf.apall.setParam('radius', '10')
    iraf.apall.setParam('threshold', '800')
    # automatic finding and centering parameters
    iraf.apall.setParam('nfind', '78')
    iraf.apall.setParam('minsep', '5')
    iraf.apall.setParam('order', 'increasing')
    # resizing parameters
    iraf.apall.setParam('aprecen', '')
    iraf.apall.setParam('npeaks', 'INDEF')
    iraf.apall.setParam('shift', 'yes')
    # resizing parameters
    iraf.apall.setParam('llimit', 'INDEF')
    iraf.apall.setParam('ulimit', 'INDEF')
    iraf.apall.setParam('ylevel', '0.1')
    iraf.apall.setParam('peak', 'yes')
    iraf.apall.setParam('bkg', 'yes')
    iraf.apall.setParam('r_grow', '0')
    iraf.apall.setParam('avglimi', 'no')
    # tracing parameters
    iraf.apall.setParam('t_nsum', '3')
    iraf.apall.setParam('t_step', '10')
    iraf.apall.setParam('t_nlost', '4')
    iraf.apall.setParam('t_funct', 'spline3')
    iraf.apall.setParam('t_order', '3')
    iraf.apall.setParam('t_sampl', '*')
    iraf.apall.setParam('t_naver', '1')
    iraf.apall.setParam('t_niter', '0')
    iraf.apall.setParam('t_low_r', '3')
    iraf.apall.setParam('t_high_', '3')
    iraf.apall.setParam('t_grow', '0')
    # extraction parameters
    iraf.apall.setParam('backgro', 'none')
    iraf.apall.setParam('skybox', '1')
    iraf.apall.setParam('weights', 'none')
    iraf.apall.setParam('pfit', 'fit1d')
    iraf.apall.setParam('clean', 'no')
    iraf.apall.setParam('saturat', SATURATION)
    iraf.apall.setParam('readnoi', RDNOISE)
    iraf.apall.setParam('gain', GAIN)
    iraf.apall.setParam('lsigma', '4.')
    iraf.apall.setParam('usigma', '4.')
    iraf.apall.setParam('nsubaps', '1')
    iraf.apall.saveParList(filename="apall_orderdef.pars")

    print('Finding last flat to use as order location reference')
    flat_list = []
    for f in images.files_filtered(imagetyp=FLAT_KEYWORD):
        flat_list.append(f)
    if len(flat_list) == 0:
        print('No order reference frame, aborting...')
        sys.exit(1)
    order_ref_frame = sorted(flat_list)[-1]
    print('Found {} as the order reference'.format(order_ref_frame))
    print('Identify all the orders in this frame.')
    print('The first CAFE order starts around row 75')
    print('REMEMBER TO SAVE THE APERTURES TO THE DATABASE WHEN FINISHED')
    iraf.apall(input=order_ref_frame)
    return order_ref_frame

def setupFlatApertures(iraf, order_ref_frame):
    """
    Tailor the aperture model for the flats
    """
    iraf.apall.setParam('referen', order_ref_frame)
    # set up apall for echelle order identification
    iraf.apall.setParam('format', 'echelle')
    # operation mode
    iraf.apall.setParam('interac', 'yes')
    iraf.apall.setParam('find', 'no')
    iraf.apall.setParam('recen', 'no')
    iraf.apall.setParam('resize', 'yes')
    iraf.apall.setParam('edit', 'yes')
    iraf.apall.setParam('trace', 'no')
    iraf.apall.setParam('fittrac', 'no')
    iraf.apall.setParam('extract', 'no')
    iraf.apall.setParam('extras', 'no')
    iraf.apall.setParam('review', 'no')
    iraf.apall(input=MASTER_FLAT)

def flattenFlat():
    """
    Flatten the blaze shape in the dispersion direction
    of the flat fields
    """
    iraf.apflatten.setParam('output', MASTER_FLAT_FLAT.split('.')[0])
    iraf.apflatten.setParam('reference', '')
    iraf.apflatten.setParam('interac', 'yes')
    iraf.apflatten.setParam('find', 'no')
    iraf.apflatten.setParam('recente', 'no')
    iraf.apflatten.setParam('resize', 'no')
    iraf.apflatten.setParam('edit', 'no')
    iraf.apflatten.setParam('trace', 'no')
    iraf.apflatten.setParam('fittrac', 'no')
    iraf.apflatten.setParam('flatten', 'yes')
    iraf.apflatten.setParam('fitspec', 'yes')
    iraf.apflatten.setParam('line', 'INDEF')
    iraf.apflatten.setParam('nsum', '1')
    iraf.apflatten.setParam('thresho', '10')
    iraf.apflatten.setParam('pfit', 'fit2d')
    iraf.apflatten.setParam('clean', 'no')
    iraf.apflatten.setParam('saturat', SATURATION)
    iraf.apflatten.setParam('readnoi', RDNOISE)
    iraf.apflatten.setParam('gain', GAIN)
    iraf.apflatten.setParam('lsigma', '3.')
    iraf.apflatten.setParam('usigma', '3.')
    iraf.apflatten.setParam('functio', 'spline3')
    iraf.apflatten.setParam('order', '3')
    iraf.apflatten.setParam('sample', '*')
    iraf.apflatten.setParam('naverag', '1')
    iraf.apflatten.setParam('niterat', '5')
    iraf.apflatten.setParam('low_rej', '2.')
    iraf.apflatten.setParam('high_rej', '3.')
    iraf.apflatten.setParam('grow', '0.')
    iraf.apflatten.saveParList(filename="apflatten.pars")
    iraf.apflatten(input=MASTER_FLAT)

def getFlattenedFlat():
    """
    Read in the newly flattened flat for all to flat...
    """
    return CCDData.read(MASTER_FLAT_FLAT, unit=u.adu)

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

def correctData(filename, flattened_flat, filetype):
    """
    The data has been bias corrected already,
    now apply the flattened flat correction
    """
    with fits.open(filename) as fitsfile:
        hdr = fitsfile[0].header
        if filetype == 'science':
            half_exptime = float(hdr[EXPTIME_KEYWORD])/2.
            dateobs = hdr[DATEOBS_KEYWORD]
            ra = float(hdr[RA_KEYWORD])/3600.
            dec = float(hdr[DEC_KEYWORD])/3600.
            time_end = Time(dateobs,
                            scale='utc',
                            format='isot',
                            location=OBSERVATORY)
            # correct to mid exposure time
            jd_mid = time_end - half_exptime*u.second
            ltt_bary, ltt_helio = getLightTravelTimes(ra, dec, jd_mid)
            time_bary = jd_mid.tdb + ltt_bary
            time_helio = jd_mid.utc + ltt_helio
            hdr['BJD-MID'] = time_bary.jd
            hdr['HJD-MID'] = time_helio.jd
            hdr['JD-MID'] = jd_mid.jd
            hdr['UT-MID'] = jd_mid.isot
    ccd = CCDData.read(filename, unit=u.adu)
    ccd = flat_correct(ccd, flattened_flat)
    fits.writeto(filename, ccd, hdr, clobber=True)

# TODO: Decide on how to treat (or not) scattered light

def parseHeader(filename, keyword):
    """
    Parse a fits header to return a specific keyword
    """
    with fits.open(filename) as fitsfile:
        value = fitsfile[0].header[keyword]
    return value

def getArcTimes(images):
    """
    Get a dictionary of arc times to use for
    determining which arcs go with each science
    spectrum
    """
    arc_times = {}
    for f in images.files_filtered(imagetyp=ARC_KEYWORD):
        dateobs = parseHeader(f, DATEOBS_KEYWORD)
        arc_times[f] = Time(dateobs,
                            scale='utc',
                            format='isot',
                            location=OBSERVATORY)
    return arc_times

def findMatchingArcPair(arc_times, filename):
    """
    Figure out which two arcs sandwhich the science
    spectrum. Add those values to the header using
    the REFSPEC1 and REFSPEC2 keywords
    """
    dateobs = Time(parseHeader(filename, DATEOBS_KEYWORD),
                   scale='utc',
                   format='isot',
                   location=OBSERVATORY)
    # split the arcs into before and after
    diff_b, diff_a = [], []
    for key in sorted(arc_times):
        diff = (dateobs-arc_times[key]).sec
        if diff > 0:
            diff_b.append(key)
        else:
            diff_a.append(key)
    # keep track of arc sandwhich
    arclist = []
    # check for the BEFORE arc
    try:
        before = diff_b[-1]
        arclist.append(before)
    except IndexError:
        print('WARNING: NO BEFORE ARC FOUND FOR {}'.format(filename))
    # check for the AFTER arc
    try:
        after = diff_a[0]
        arclist.append(after)
    except IndexError:
        print('WARNING: NO AFTER ARC FOUND FOR {}'.format(filename))
    if len(arclist) == 0:
        print('FATAL: NO BEFORE/AFTER ARC FOUND FOR {}'.format(filename))
        sys.exit(1)
    else:
        refspec_factor = round(1./len(arclist), 1)
    # edit the filename headers for update the REFSPEC1/2 keywords
    # this is used for identify to wavelength calibration
    for i in range(0, len(arclist)):
        refspec = "{} {}".format(arclist[i].split(".fits")[0], refspec_factor)
        print('{}: REFSPEC{} {}'.format(filename, i+1, arclist[i]))
        iraf.hedit(images=filename,
                   fields="REFSPEC{}".format(i+1),
                   value=refspec,
                   add="yes",
                   verify="no",
                   show="yes")
    return arclist

def extractSpectra(filename, arclist, order_ref_frame, global_arclist):
    """
    Extract echelle spectra using IRAF interactively

    TODO:
        Interpolate across the two arcs either side to
        get the most accurate wavelength solution
        Finish docstring
        Add method of using super arc for inital
        identify
    """
    # set up apall for echelle order identification
    iraf.apall.setParam('apertur', '')
    iraf.apall.setParam('format', 'echelle')
    iraf.apall.setParam('reference', order_ref_frame)
    iraf.apall.setParam('profile', '')
    # operation mode
    iraf.apall.setParam('interac', 'yes')
    iraf.apall.setParam('find', 'no')
    iraf.apall.setParam('recen', 'no')
    iraf.apall.setParam('resize', 'no')
    iraf.apall.setParam('edit', 'no')
    iraf.apall.setParam('trace', 'no')
    iraf.apall.setParam('fittrac', 'no')
    iraf.apall.setParam('extract', 'yes')
    iraf.apall.setParam('extras', 'yes')
    iraf.apall.setParam('review', 'yes')
    iraf.apall.setParam('line', '904')
    iraf.apall.setParam('nsum', '5')
    # default aperture parameters
    iraf.apall.setParam('lower', '-5')
    iraf.apall.setParam('upper', '5')
    # default background parameters
    iraf.apall.setParam('b_funct', 'chebyshev')
    iraf.apall.setParam('b_order', '2')
    iraf.apall.setParam('b_sampl', '-9:-7,5:7')
    iraf.apall.setParam('b_naver', '-3')
    iraf.apall.setParam('b_niter', '3')
    iraf.apall.setParam('b_low_r', '3')
    iraf.apall.setParam('b_high', '3')
    iraf.apall.setParam('b_grow', '0')
    # aperture centering parameters
    iraf.apall.setParam('width', '5')
    iraf.apall.setParam('radius', '10')
    iraf.apall.setParam('threshold', '800')
    # automatic finding and centering parameters
    iraf.apall.setParam('nfind', '78')
    iraf.apall.setParam('minsep', '5')
    iraf.apall.setParam('order', 'increasing')
    # resizing parameters
    iraf.apall.setParam('aprecen', '')
    iraf.apall.setParam('npeaks', 'INDEF')
    iraf.apall.setParam('shift', 'yes')
    # resizing parameters
    iraf.apall.setParam('llimit', 'INDEF')
    iraf.apall.setParam('ulimit', 'INDEF')
    iraf.apall.setParam('ylevel', '0.1')
    iraf.apall.setParam('peak', 'yes')
    iraf.apall.setParam('bkg', 'yes')
    iraf.apall.setParam('r_grow', '0')
    iraf.apall.setParam('avglimi', 'no')
    # tracing parameters
    iraf.apall.setParam('t_nsum', '3')
    iraf.apall.setParam('t_step', '10')
    iraf.apall.setParam('t_nlost', '4')
    iraf.apall.setParam('t_funct', 'spline3')
    iraf.apall.setParam('t_order', '3')
    iraf.apall.setParam('t_sampl', '*')
    iraf.apall.setParam('t_naver', '1')
    iraf.apall.setParam('t_niter', '0')
    iraf.apall.setParam('t_low_r', '3')
    # extraction parameters
    iraf.apall.setParam('backgro', 'fit')
    iraf.apall.setParam('skybox', '6')
    iraf.apall.setParam('weights', 'variance')
    iraf.apall.setParam('pfit', 'fit2d')
    iraf.apall.setParam('clean', 'no')
    iraf.apall.setParam('saturat', SATURATION)
    iraf.apall.setParam('readnoi', RDNOISE)
    iraf.apall.setParam('gain', GAIN)
    iraf.apall.setParam('lsigma', '3.')
    iraf.apall.setParam('usigma', '3.')
    iraf.apall.setParam('nsubaps', '1')
    iraf.apall.saveParList(filename="apall_extraction.pars")

    # extract the science spectrum
    with fits.open(filename) as fitsfile:
        target_id = fitsfile[0].header[OBJECT_KEYWORD]
    print("Extracting spectrum of {} from {}".format(target_id, filename))
    print("Check aperture and background. Change if required")
    print("AP: m = mark aperture, d = delete aperture")
    print("SKY: s = mark sky, t = delete sky, f = refit")
    print("q = continue\n")
    iraf.apall(input=filename)
    print("Spectrum extracted!\n")

    # extract the arcs
    for arc in arclist:
        if arc not in global_arclist:
            iraf.apall(input=arc,
                       reference=order_ref_frame,
                       interac="no",
                       extract="yes",
                       extras="no",
                       review="no",
                       backgro="none",
                       weights="none",
                       clean="no")
            print("Arc {} extracted".format(arc))
            global_arclist.append(arc)

    # wavelength calibration

if __name__ == '__main__':
    print('\n\n---------------------------------------------------------')
    print('------------- Spectroscopy by RISTRETTO ----------------')
    print('---------------------------------------------------------\n')

    junk_yn = raw_input("Has the junk been removed from current folder? (y/n): ")
    if junk_yn != 'y':
        print('Remove all files that should not be analysed, then restart')
        sys.exit(1)

    # make a local backup copy of the data
    copyFiles()
    # get a list of images in current directory
    images = getImageList()
    # trim the CAFE images
    trimFiles()
    # rename the images to make them easier to read
    renameFiles(images)
    # get new filenames
    images = getImageList()
    # make a master bias
    master_bias = makeMasterBias(images)
    # bias correct the images
    subtractMasterBias(images, master_bias)
    # make a master flat
    master_flat = makeMasterFlat(images)
    # identify the aperture locations in a flat field
    iraf = importIraf()
    order_ref_frame = identifyOrders(iraf, images)
    # get aperture model for flat field
    setupFlatApertures(iraf, order_ref_frame)
    # flatten the flat
    flattenFlat()
    # read in the flattened flat
    flattened_flat = getFlattenedFlat()
    # correct the arcs
    for filename in images.files_filtered(imagetyp=ARC_KEYWORD):
        correctData(filename, flattened_flat, 'arc')
    # correct the science spectra
    for filename in images.files_filtered(imagetyp=SCIENCE_KEYWORD):
        correctData(filename, flattened_flat, 'science')
    # get the arc times for file association
    arc_times = getArcTimes(images)
    # make a global list to track the extracted arcs
    # some arcs are REFSPEC for several science images
    global_arclist = []
    for filename in images.files_filtered(imagetyp=SCIENCE_KEYWORD):
        arclist = findMatchingArcPair(arc_times, filename)
        extractSpectra(filename, arclist, order_ref_frame, global_arclist)

    #cleanCalibs()
    # extract wavelength calibrated and continuum normalised spectra
    #extractSpectra()
    # round up the reduced spectra
    #roundUpSpectra()
    # remove all the intermediate data products
    #eraseIntermediateProducts()
    # log the spectra to the database
    #logSpectraToDb()

