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
from matplotlib.ticker import MaxNLocator

df = xr.open_dataset("/Users/au710474/Documents/scripts/EDGAR/v6.0_CH4_TOTAL_NO_FIRES_2018.0.1x0.1.nc")

df["lon_adj"] = xr.where(df["lon"] > 180, df["lon"]-360, df["lon"])
df = df.sortby(df.lon_adj)

lon_0 = 5.66
lat_0 = 51.98
lat = df.lat.values
lon = df.lon_adj.values

x, y = np.meshgrid(lon,lat)

ch4 = df.emis_tot.values
ch4[ch4 == 0] = np.nan
levels = MaxNLocator(nbins=20).tick_values(1,10) 

fig = plt.figure(figsize=(15, 5))
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
cs = m.contourf(x1, y1, ch4[4]*10**10, cmap='YlOrBr', levels=levels)

plt.title("$Kg$ $m^{-2}s^{-1}$ [$x10^{-10}$]", loc='right', fontsize=12)
cbar = m.colorbar(cs, location='right', ax=ax, pad="2%", ticks=levels[::6])
#cbar.set_label("$Kg$ $m^{-2}s^{-1}$ [$x10^{-10}$]", labelpad=2, fontsize=15)
cbar.ax.tick_params(labelsize=15) 

plt.savefig("edgar_ch4_anthro.png",bbox_inches='tight')
plt.show()
