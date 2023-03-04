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

df_tot = xr.open_dataset("/home/angel/Documents/iag-usp/modis/methane_project/v6.0_CH4_TOTAL_NO_FIRES_2018.0.1x0.1.nc")
df_agr = xr.open_dataset("/home/angel/Documents/iag-usp/modis/methane_project/v6.0_CH4_AGRICULTURE_2018.0.1x0.1.nc")
df_was = xr.open_dataset("/home/angel/Documents/iag-usp/modis/methane_project/v6.0_CH4_WASTE_2018.0.1x0.1.nc")
df_ene = xr.open_dataset("/home/angel/Documents/iag-usp/modis/methane_project/v6.0_CH4_ENERGY_2018.0.1x0.1.nc")

df_tot["lon_adj"] = xr.where(df_tot["lon"] > 180, df_tot["lon"]-360, df_tot["lon"])
df_tot = df_tot.sortby(df_tot.lon_adj)
df_agr["lon_adj"] = xr.where(df_agr["lon"] > 180, df_agr["lon"]-360, df_agr["lon"])
df_agr = df_agr.sortby(df_agr.lon_adj)
df_was["lon_adj"] = xr.where(df_was["lon"] > 180, df_was["lon"]-360, df_was["lon"])
df_was = df_was.sortby(df_was.lon_adj)
df_ene["lon_adj"] = xr.where(df_ene["lon"] > 180, df_ene["lon"]-360, df_ene["lon"])
df_ene = df_ene.sortby(df_ene.lon_adj)

lon_0 = 5.66
lat_0 = 51.98
pinf = float('+inf')
lat = df_tot.lat.values
lon = df_tot.lon_adj.values

x, y = np.meshgrid(lon,lat)

ch4_tot = df_tot.emis_tot.values
ch4_tot[ch4_tot == 0] = np.nan
ch4_agr = df_agr.emis_tot.values
ch4_agr[ch4_agr == 0] = np.nan
ch4_was = df_was.emis_tot.values
ch4_was[ch4_was == 0] = np.nan
ch4_ene = df_ene.emis_tot.values
ch4_ene[ch4_ene == 0] = np.nan
ch4 = np.array([ch4_agr, ch4_was, ch4_ene, ch4_tot])
levels = MaxNLocator(nbins=20).tick_values(1,5) 

titulos = ['(a) Agriculture', '(b) Waste', '(c) Energy', '(d) All sources (except fires)']

fig, axes = plt.subplots(2, 2, figsize=(28,28))
c=0
for n_row in range(2):
    for n_col in range(2):    
        m = Basemap(width = 3600000, height = 3600000,
                    projection = 'stere', lat_0 = lat_0, lon_0 = lon_0, resolution = 'l', ax=axes[n_row,n_col])
        m.drawcoastlines(linewidth=1)
        m.drawstates(linewidth=1)  
        m.drawcountries(linewidth=1)
        m.drawmapboundary()
        x1, y1 = m(x, y)
        x1[x1 == pinf] = 0
        x1[x1 == pinf] = 0
        cs = m.contourf(x1, y1, ch4[c][4]*10**10, cmap='YlOrRd', levels=levels, extend='max')
        axes[n_row,n_col].set_title(titulos[c], loc='left', fontsize=13)
        c=c+1
plt.subplots_adjust(bottom=0.12, right=0.38, top=0.33)
clb=plt.colorbar(cs, ax=axes.ravel().tolist(),
                 shrink=0.7, pad=0.05, ticks=[1, 2, 3, 4, 5])
clb.set_label(label='$Kg$ $m^{-2}s^{-1}$ [$x10^{-10}$]', size=15, weight='bold')
clb.ax.tick_params(labelsize=15)

plt.show()
