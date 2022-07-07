#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 13:48:53 2020

@author: angel
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator
import glob
from datetime import datetime

path = '/home/angel/Documents/iag-usp/tropomi/lorenteetal2021/1004_2304_2018/'
output = '/home/angel/Documents/iag-usp/tropomi/lorenteetal2021/figs_sron/'
nc_file = sorted(glob.glob(path+'*.nc'))

#lon_0 = 5.66
#lat_0 = 51.98

for i in nc_file:    
    
    print(i)
    fh   = Dataset(i, mode='r')
    
    scan = fh.groups['instrument'].variables['scanline'][:]
    grpi = fh.groups['instrument'].variables['ground_pixel'][:]
    xch4 = fh.groups['target_product'].variables['xch4'][:]
    qa   = fh.groups['diagnostics'].variables['qa_value'][:]
    lats = fh.groups['instrument'].variables['latitude_center'][:]
    lons = fh.groups['instrument'].variables['longitude_center'][:]
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
        continue

    date = datetime(Y,m,d,H,M,S)
    scan_2d = np.zeros((l_scan,l_grpi))
    grpi_2d = np.zeros((l_scan,l_grpi))
    xch4_2d = np.zeros((l_scan,l_grpi))
    qa_2d   = np.zeros((l_scan,l_grpi))
    lats_2d = np.zeros((l_scan,l_grpi))
    lons_2d = np.zeros((l_scan,l_grpi))
    
    for j in range(l_scan):
        for k in range(l_grpi):
            scan_2d[j][k] = scan[j*l_grpi+k]
            grpi_2d[j][k] = grpi[j*l_grpi+k]
 #           if qa_2d[j][k] > 0.5:
            xch4_2d[j][k] = xch4[j*l_grpi+k]
            lats_2d[j][k] = lats[j*l_grpi+k]
            lons_2d[j][k] = lons[j*l_grpi+k]
            
    lons_2d = np.clip(lons_2d, -180.0, 180.0)
    lats_2d = np.clip(lats_2d, -90.0, 90.0)

    levels = MaxNLocator(nbins=25).tick_values(1700,1950)
    cmap = colors.LinearSegmentedColormap.from_list(
        'white11red', ['#FFFFFF','paleturquoise','lightskyblue','aqua',
        'lightgreen','limegreen','yellow','goldenrod','darkorange',
        'red','darkred'])   
  #  cmap = cm.get_cmap("Spectral_r",lut=25)
  #  cmap.set_bad("w")
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    
    fig = plt.figure(figsize=(25, 15))
    ax = fig.add_subplot(111, frame_on=False)
    ax.set_xticks([])
    ax.set_yticks([])
  #  m = Basemap(width = 3600000, height = 3600000,               # 3600000
  #            projection = 'stere', lat_0 = lat_0, lon_0 = lon_0, resolution = 'l')
    m = Basemap(projection='cyl', resolution='l', llcrnrlat=-90.0, urcrnrlat=90.0,
                llcrnrlon=-180.0, urcrnrlon=180.0)
    m.drawcoastlines(linewidth=1.0)
    m.drawcountries(linewidth=1.0)
    m.drawmapboundary()  
    m.drawparallels(np.arange(-90., 120., 30), labels=[1, 0, 0, 0], fontsize=20)
    m.drawmeridians(np.arange(-180., 181., 30), labels=[0, 0, 0, 1], fontsize=20)
    cs = m.contourf(lons_2d, lats_2d, xch4_2d, levels=levels, cmap=cmap, norm=norm)
    
    plt.title(date, loc='left', fontsize=30)
    cbar = m.colorbar(cs, location='bottom', pad="5%", ticks=levels[::5], ax=ax)        
    cbar.set_label("xch4 [ppb]", fontsize=30)
    cbar.ax.tick_params(labelsize=30)
    fig.savefig(output+"xch4_"+str(date)+'.png', bbox_inches='tight')
    plt.close(fig)
