# CAFE SNR calculator
from numpy.polynomial import Polynomial as P
from collections import defaultdict
import matplotlib.pyplot as pl
import argparse as ap
import numpy as np
import seaborn 
import math

def argParse():
    parser=ap.ArgumentParser()
    parser.add_argument('--plot',choices=['save','show'],help='Plot the results?')
    parser.add_argument('--predict',help='object flag to predict exptime',default='CONTINUE-STABILIZED')
    parser.add_argument('--snr',type=float,help='SNR ratio for calculations',default=30)
    return parser.parse_args()

args=argParse()
E=0.08 # 0.16 mags per airmass - therefore assuming airmass of 1.5 for most obs
V=np.arange(8,14.0,0.5)
t=np.array([30,120,300,900,2700])

# overheads
readout=90
acq_guide=300
arcs=2*(120+90)+60 # 2x120s arcs + 2x readout + 1 min buffer
n_target_spectra=15 # target number of spectra to estimate the total time requested

# exposure times obtained manually from 
# CAFE SN calculator
snr=defaultdict(list)
snr[8.]   = np.array([23.381,46.762,73.937,128.063,221.812]) 
snr[8.5]  = np.array([18.572,37.144,58.731,101.724,176.192])
snr[9.]   = np.array([15.306,30.612,48.402,83.835,145.206])
snr[9.5]  = np.array([12.158,24.316,38.447,66.592,115.341])
snr[10.]  = np.array([9.657,19.315,30.540,52.896,91.619])
snr[10.5] = np.array([7.671,15.342,24.258,42.017,72.776])
snr[11.]  = np.array([6.093,12.817,19.269,33.375,57.807])
snr[11.5] = np.array([4.840,9.680,15.306,26.511,45.918])
snr[12.]  = np.array([3.844,7.689,12.158,21.058,36.474])
snr[12.5] = np.array([3.05,6.108,9.657,16.727,28.972])
snr[13.]  = np.array([2.426,4.852,7.671,13.287,23.014])
snr[13.5] = np.array([1.923,3.854,6.093,10.554,18.280])

seaborn.axes_style("darkgrid")
pl.rc('legend',**{'fontsize':10})
coeffs_store=defaultdict(list)
fig,ax=pl.subplots(1,figsize=(10,10))
ax.set_title('CAFE')
for i in sorted(snr.keys()):
    coeffs=np.polyfit(t,snr[i],2)
    coeffs_store[i]=coeffs
    tn=np.arange(30,3060,60)
    besty=np.polyval(coeffs,tn)
    ax.plot(t,snr[i],'o',tn,besty,'k--')

if args.predict:
    import pymysql
    db=pymysql.connect(host='localhost',db='eblm')
    qry="SELECT swasp_id,Vmag FROM eblm_parameters WHERE current_status='%s'" % (args.predict)
    
    total_time=0
    with db.cursor() as cur:
        cur.execute(qry)
        for row in cur:
            swasp_id=row[0]
            vmag=float(row[1])
            diff=V-vmag
            n=np.where(abs(diff)==min(abs(diff)))[0]                
            p=P.fit(t,snr[V[n[0]]],2)
            t1,t2=(p-args.snr).roots()
            if t1<t2:
                predicted=t1
            else:
                predicted=t2    
            # check for texp>2700, scale to right number of
            # spectra to get required SNR
            if predicted > 2700:
                snr_max=snr[V[n[0]]]
                n_spectra=(args.snr/snr[V[n[0]]][-1])**2
                predicted=2700
            elif predicted < 0:
                predicted = 60
                n_spectra=1
            else:
                n_spectra=1
            print("%s %.2f %d x %.2fs" % (swasp_id,vmag,n_spectra,predicted))
            total_time=total_time+((n_spectra*(predicted+readout+arcs))+acq_guide)*n_target_spectra
    print("Total time needed w/ overheads: %ds [%.2f hours]" % (total_time,(total_time/3600.)))

ax.set_ylabel('SNR')
ax.set_xlabel('Exposure Time (s)')
ax.legend(('8.0','fit','8.5','fit','9.0','fit','9.5','fit','10.0','fit','10.5','fit','11.0','fit','11.5','fit','12.0','fit','12.5','fit','13.0','fit','13.5','fit'),loc='upper left')
if args.plot == 'save':
    fig.savefig('CAFE-SNR.png',dpi=200)
else:
    pl.show()







