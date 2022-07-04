#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 11:17:56 2022

@author: angel
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator

df = xr.open_dataset("/home/angel/Documents/iag-usp/modis/methane_project/v6.0_CH4_ALL_2018.0.1x0.1.nc")

df["lon_adj"] = xr.where(df["lon"] > 180, df["lon"]-360, df["lon"])
df = df.sortby(df.lon_adj)

lon_0 = 5.66
lat_0 = 51.98
lat = df.lat.values
lon = df.lon_adj.values

x, y = np.meshgrid(lon,lat)

ch4 = df.emis_tot.values
ch4[ch4 == 0] = -999.
levels = MaxNLocator(nbins=25).tick_values(0,5)

cmap = colors.LinearSegmentedColormap.from_list(
'white11red', ['#FFFFFF','paleturquoise','lightskyblue','aqua',
               'lightgreen','limegreen','yellow','goldenrod','darkorange',
               'red','darkred'])    

fig = plt.figure(figsize=(12, 4))
ax = fig.add_subplot(111)
for axis in ['top','bottom','left','right']:
     ax.spines[axis].set_linewidth(3.0)
m = Basemap(width = 3600000, height = 3600000,
             projection = 'stere', lat_0 = lat_0, lon_0 = lon_0, resolution = 'l')
m.drawcoastlines(linewidth=1.0)
m.drawstates(linewidth=1.0)  
m.drawcountries(linewidth=1.0)
m.drawmapboundary()
x1, y1 = m(x, y)
cs = m.pcolormesh(x1, y1, ch4[4]*10**10, cmap='YlOrBr', vmin=0, vmax=5)

plt.title("All sources", loc='left', fontsize=15)
cbar = m.colorbar(cs, location='bottom', ax=ax, pad="2%", ticks=levels[::5])
cbar.set_label("$Kg$ $m^{-2}s^{-1}$ [$x10^{-10}$]", labelpad=2, fontsize=15)
cbar.ax.tick_params(labelsize=15) 

plt.savefig("edgar_ch4_all.png",bbox_inches='tight')
#plt.show()
