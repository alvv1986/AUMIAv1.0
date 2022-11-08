#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 12:39:04 2020

@author: angel
"""

import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.geoaxes
from netCDF4 import Dataset
from shapely.geometry.polygon import LinearRing
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.colors

class FixPointNormalize(matplotlib.colors.Normalize):
    """ 
    Inspired by https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    Subclassing Normalize to obtain a colormap with a fixpoint 
    somewhere in the middle of the colormap.
    This may be useful for a `terrain` map, to set the "sea level" 
    to a color in the blue/turquise range. 
    """
    def __init__(self, vmin=None, vmax=None, sealevel=0, col_val = 0.21875, clip=False):
        # sealevel is the fix point of the colormap (in data units)
        self.sealevel = sealevel
        # col_val is the color value in the range [0,1] that should represent the sealevel.
        self.col_val = col_val
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.sealevel, self.vmax], [0, self.col_val, 1]
        return np.ma.masked_array(np.interp(value, x, y))

# Combine the lower and upper range of the terrain colormap with a gap in the middle
# to let the coastline appear more prominently.
# inspired by https://stackoverflow.com/questions/31051488/combining-two-matplotlib-colormaps
#colors_undersea = plt.cm.terrain(np.linspace(0, 0.17, 56))
#colors_land = plt.cm.terrain(np.linspace(0.25, 1, 200))
colors_undersea = plt.cm.terrain(np.linspace(0, 0.06, 57))
colors_land = plt.cm.terrain(np.linspace(0.25, 1.0, 200))

# combine them and build a new colormap
colors = np.vstack((colors_undersea, colors_land))
cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors)

norm = FixPointNormalize(sealevel=0,vmax=2000,vmin=0)

wrf_sp_d01 = Dataset('/home/angel/Documents/iag-usp/tropomi/Domains/geo_em.d01.nc', mode='r')
wrf_sp_d02 = Dataset('/home/angel/Documents/iag-usp/tropomi/Domains/geo_em.d02.nc', mode='r')

# Loading coordinates
xlat1 = wrf_sp_d01['XLAT_M'][0]
xlon1 = wrf_sp_d01['XLONG_M'][0]
topo1 = wrf_sp_d01['HGT_M'][0]
landmask1 = wrf_sp_d01['LANDMASK'][0]
xlat2 = wrf_sp_d02['XLAT_M'][0]
xlon2 = wrf_sp_d02['XLONG_M'][0]
topo2 = wrf_sp_d02['HGT_M'][0]

# Defining d02 domain
corner_lats_2 = wrf_sp_d02.getncattr('corner_lats')
corner_lons_2 = wrf_sp_d02.getncattr('corner_lons')
d02_lat = [corner_lats_2[0], corner_lats_2[1], corner_lats_2[2], corner_lats_2[3]]
d02_lon = [corner_lons_2[0], corner_lons_2[1], corner_lons_2[2], corner_lons_2[3]]
d02_sqr = LinearRing(list(zip(d02_lon, d02_lat)))

fig = plt.figure(figsize=(11, 11))
ax = plt.axes(projection = ccrs.Orthographic(central_latitude=49.1, central_longitude=15.1))

ax.plot(14.1226, 52.1663, 'r^', markersize=8, transform=ccrs.Geodetic())   # Lindenberg
ax.text(14.4, 52.1663, '1', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(8.6360, 45.8147, 'ro', markersize=8, transform=ccrs.Geodetic())    # Ispra
ax.text(9.0, 45.8147, '2', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(2.1125, 47.9647, 'ro', markersize=8, transform=ccrs.Geodetic())    # Trainou
ax.text(2.6, 47.4, '3', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(5.5036, 48.5619, 'r^', markersize=8, transform=ccrs.Geodetic())    # Observatoire pérenne de l'environnement
ax.text(5.8, 48.5619, '4', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(15.0800, 49.572, 'r^', markersize=8, transform=ccrs.Geodetic())    # Křešín u Pacova
ax.text(15.3, 49.572, '5', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(11.0246, 47.8011, 'ro', markersize=8, transform=ccrs.Geodetic())   # Hohenpeissenberg
ax.text(11.3, 47.8011, '6', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(24.2947, 61.8474, 'ro', markersize=8, transform=ccrs.Geodetic())   # Hyytiälä
ax.text(24.6, 61.8474, '7', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(2.142, 48.7227, 'ro', markersize=8, transform=ccrs.Geodetic())     # Saclay
ax.text(2.5, 48.7227, '8', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(10.535, 51.8088, 'ro', markersize=8, transform=ccrs.Geodetic())    # Torfhaus
ax.text(10.9, 51.8088, '9', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(2.9658, 45.7719, 'ro', markersize=8, transform=ccrs.Geodetic())    # Puy de Dôme
ax.text(3.3, 45.7719, '10', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')
ax.plot(10.6999, 44.19369, 'ro', markersize=8, transform=ccrs.Geodetic())    # Puy de Dôme
ax.text(11.0, 44.19369, '11', transform=ccrs.Geodetic(), fontsize=13, weight = 'bold')

ax.add_feature(cfeature.BORDERS)
ax.coastlines('10m', zorder=15)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)

cbax = ax.pcolormesh(xlon1, xlat1, topo1,
                   transform=ccrs.PlateCarree(),
                   cmap=cut_terrain_map, norm=norm)

ax.add_geometries([d02_sqr], ccrs.PlateCarree(), facecolor='none', 
                  edgecolor="red", linewidth=1.5, zorder=20)
ax.text(xlon2.max() - 0.5, xlat2.min(), "D02", fontsize=20,
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
    plt.pcolormesh(xlon1, xlat1, topo1, transform=ccrs.PlateCarree(), cmap=cut_terrain_map, norm=norm)
    return ax

add_sub_region_box_1(axins)

cb = plt.colorbar(cbax, ax=ax, shrink=0.8)
cb.ax.tick_params(labelsize=15)
cb.set_label(label="HGT (m)", labelpad=-40, rotation=0, y=1.085, fontsize=20)
plt.tick_params(labelsize=22)
#plt.savefig('/Users/au710474/Documents/scripts/WRFSP5/Dom/domains.png', dpi=300,  bbox_inches="tight")
plt.show()
