"""
Script to deteremine the average best sky apertures
for CAFE data.

I took a manual sweep through all the apertures in an
example master_flat.fits image (see master_calibrations
folder in CAFE data directory). Parse this to look for 
the best average values
"""

master_database = '/Users/James/Dropbox/data/cafe/master_calibrations/database/apmaster_flat'

f = open(master_database).readlines()
l_u, l_l = [], []
r_u, r_l = [], []
for i in f:
    if 'sample' in i and ':' in i:
        blh, left, right = i.split()
        l_l.append(float(left.split(':')[0]))
        l_u.append(float(left.split(':')[1]))
        r_l.append(float(right.split(':')[0]))
        r_u.append(float(right.split(':')[1]))

print('Left:')
print('{}:{}'.format(round(max(l_l),0), round(max(l_u),0)))
print('Right:')
print('{}:{}'.format(round(min(r_l),0), round(min(r_u),0)))
