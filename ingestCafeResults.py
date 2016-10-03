"""
Script to ingest the results files from cafepipe.py

THIS IS NOT FINISHED AND HAD NOT BEEN TESTED

Results file has the following structure:
    0-  Object name
    1-  MBJD
    2-  RV
    3-  error in RV
    4-  Bisector span
    5-  error in bisector span
    6-  instrument
    7-  pipeline
    8-  resolving power
    9-  Efective Temperture
    10- log(g)
    11- [Fe/H]
    12- v*sin(i)
    13- value of the continuum normalized CCF at it lowest point
    14- standard deviation of the gaussian fitted to the CCF
    15- Exposure time
    16- Signal to noise ratio at ~5150A
    17- path to the CCF plot file
"""
import sys
import pymysql

if __name__ == '__main__':
    db = pymysql.connect(host='localhost',
                         db='eblm',
                         password='mysqlpassword')
    f = open('results.txt').readlines()
    for line in f:
        ls = line.rstrip().split()
        if len(ls) != 18:
            print('ERROR: Wrong number of columns in results.txt')
            sys.exit(1)
        obj = ls[0]
        if obj.startswith('1SWASP'):
            swasp_id = obj
        else:
            swasp_id = None
        bjd_mid = ls[1]
        mask_rv = ls[2]
        mask_rv_err = ls[3]
        bisector = ls[4]
        bisector_err = ls[5]
        mask_ccf_height = ls[13]
        mask_ccf_fwhm = ls[14]
        snr_5150 = ls[16]
        pdf_name = ls[17].split('/')[-1]
        image_id = '{}.fits'.format(pdf_name.split('.')[0])
        mask = pdf_name.split('.')[1].split('_')[2]
        qry = """
            INSERT INTO eblm_cafe (
                image_id,
                swasp_id,
                object_name,
                bjd_mid,
                mask,
                mask_rv,
                mask_rv_err,
                mask_ccf_height,
                mask_ccf_fwhm,
                bisector,
                bisector_err,
                snr5150
                )
            VALUES (
                '{}',
                '{}',
                '{}',
                {},
                '{}',
                {},
                {},
                {},
                {},
                {},
                {},
                {}
            """.format(image_id,
                       swasp_id,
                       obj,
                       bjd_mid,
                       mask,
                       mask_rv,
                       mask_rv_err,
                       mask_ccf_height,
                       mask_ccf_fwhm,
                       bisector,
                       bisector_err,
                       snr_5150)
        print(qry)
