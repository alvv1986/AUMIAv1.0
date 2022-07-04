#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:39:04 2020

@author: angel
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.geoaxes
from netCDF4 import Dataset
from shapely.geometry.polygon import LinearRing
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

wrf_sp_d01 = Dataset('/home/angel/Documents/iag-usp/modis/Domain/Methane/geo_em.d01.nc', mode='r')
wrf_sp_d02 = Dataset('/home/angel/Documents/iag-usp/modis/Domain/Methane/geo_em.d02.nc', mode='r')

# Loading coordinates
xlat1 = wrf_sp_d01['XLAT_M'][0]
xlon1 = wrf_sp_d01['XLONG_M'][0]
topo1 = wrf_sp_d01['HGT_M'][0]
xlat2 = wrf_sp_d02['XLAT_M'][0]
xlon2 = wrf_sp_d02['XLONG_M'][0]
topo2 = wrf_sp_d02['HGT_M'][0]

# Defining d02 domain
corner_lats_2 = wrf_sp_d02.getncattr('corner_lats')
corner_lons_2 = wrf_sp_d02.getncattr('corner_lons')
d02_lat = [corner_lats_2[0], corner_lats_2[1], corner_lats_2[2], corner_lats_2[3]]
d02_lon = [corner_lons_2[0], corner_lons_2[1], corner_lons_2[2], corner_lons_2[3]]
d02_sqr = LinearRing(list(zip(d02_lon, d02_lat)))

fig = plt.figure(figsize=(10, 9))
ax = plt.axes(projection = ccrs.Orthographic(central_latitude=49.1, central_longitude=15.1))

ax.add_feature(cfeature.BORDERS)
ax.coastlines('10m', zorder=15)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)

cbax = ax.pcolormesh(xlon1, xlat1, topo1,
                   transform=ccrs.PlateCarree(),
                   cmap="terrain")

ax.add_geometries([d02_sqr], ccrs.PlateCarree(), facecolor='none', 
                  edgecolor="red", linewidth=1.75, zorder=20)
ax.text(xlon2.min(), xlat2.max() + 0.30, "D02", fontsize=16,
        transform=ccrs.PlateCarree(), color="red")
axins = inset_axes(ax, width="32%", height="32%", loc="upper left", 
                   axes_class=cartopy.mpl.geoaxes.GeoAxes, 
                   axes_kwargs=dict(
                       map_projection=cartopy.crs.Orthographic(
                           central_longitude=19.0, central_latitude=50.0)
                       ))
axins.stock_img()
axins.coastlines()
axins.add_feature(cfeature.OCEAN)
axins.add_feature(cfeature.LAND)

def add_sub_region_box_1(ax):
    """ """
    plt.pcolormesh(xlon1, xlat1, topo1, transform=ccrs.PlateCarree(), cmap='terrain')
    return ax

add_sub_region_box_1(axins)

cb = plt.colorbar(cbax, ax=ax, label = 'HGT (m)', shrink=0.8)
cb.ax.tick_params(labelsize=15)
cb.set_label(label="HGT (m)", labelpad=-40, rotation=0, y=1.05, fontsize=15)
plt.tick_params(labelsize=20)
plt.savefig('/home/angel/Documents/iag-usp/modis/Domain/Methane/domains.png', dpi=300,  bbox_inches="tight")
#plt.show()