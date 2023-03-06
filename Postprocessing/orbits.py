#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: angel
"""

from netCDF4 import Dataset
import os
import shutil
from datetime import datetime

f = 's5p_l2_ch4_0017_00000.nc'
output = './orbits/'    
    
fh = Dataset(f, mode='r')
orbit = f[:-3].split('_')[4]
    
scan = fh.groups['instrument'].variables['scanline'][:]
grpi = fh.groups['instrument'].variables['ground_pixel'][:]
xch4 = fh.groups['target_product'].variables['xch4_corrected'][:]
time = fh.groups['instrument'].variables['time'][:]
Y = time[len(time)-1][0]
m = time[len(time)-1][1]
d = time[len(time)-1][2]
H = time[len(time)-1][3]
M = time[len(time)-1][4]
S = time[len(time)-1][5]
if (scan[len(scan)-1]+1)*(grpi[len(grpi)-2]+2) == len(xch4):
    l_grpi = grpi[len(grpi)-2]+2
    l_scan = scan[len(scan)-1]+1
elif (len(xch4)/215).is_integer() and Y*m*d > 0:
    l_grpi = 215
    l_scan = len(xch4)/215
else:
    print('corrupted file!')
    os.remove(f)
    exit()
if H > 9 and H < 15:
    date = datetime(Y,m,d,H,M,S)
    shutil.move(f, output)
    print(f, date)
else:
    os.remove(f)
