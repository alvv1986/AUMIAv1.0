#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 16:42:32 2022

@author: angel
"""

import numpy as np
import os
import xesmf as xe
import xarray as xr
from netCDF4 import Dataset
from shapely.geometry import Polygon
from datetime import datetime
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import matplotlib.colors as colors
import shutil

levels = MaxNLocator(nbins=25).tick_values(1700,1950)
cmap = colors.LinearSegmentedColormap.from_list(
    'white11red', ['#FFFFFF','paleturquoise','lightskyblue','aqua',
    'lightgreen','limegreen','yellow','goldenrod','darkorange',
    'red','darkred'])   
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp)
    y = np.where(x<xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    y = np.where(x>xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]), y)
    return y

in_tropo = '/home/angel/tropomi/sron_s5p/'
in_wrf = '/home/angel/tropomi/wrfinput/'
out_tropo = '/home/angel/tropomi/regridded/'
out_xe = '/home/angel/tropomi/regridded/xe/'
figs = '/home/angel/tropomi/figs/'

f_tropo = sorted(os.listdir(in_tropo))

################################################# grid_out ###############################################
ds_disk = xr.open_dataset(in_wrf+"wrfinput_d01.nc")
dat = ds_disk["LANDMASK"].values[0]
lon_out = ds_disk.XLONG.values[0]
lat_out = ds_disk.XLAT.values[0]

############################################### corner points ############################################
lon_out_b = np.zeros((lat_out.shape[0]-1,lon_out.shape[1]-1))
lat_out_b = np.zeros((lat_out.shape[0]-1,lon_out.shape[1]-1))
for i in range(lat_out.shape[0]-1):
    for j in range(lon_out.shape[1]-1):
        lon1 = lon_out[i:i+2,j:j+2]
        lat1 = lat_out[i:i+2,j:j+2]
        p1 = Polygon([(lon1[0,0],lat1[0,0]),(lon1[0,1],lat1[0,1]),
                      (lon1[1,1],lat1[1,1]),(lon1[1,0],lat1[1,0])])
        lon_out_b[i,j] = p1.centroid.x
        lat_out_b[i,j] = p1.centroid.y
        
lon2 = np.zeros((lon_out_b.shape[0]+2,lon_out_b.shape[1]+2))
lat2 = np.zeros((lon_out_b.shape[0]+2,lon_out_b.shape[1]+2))

for i in range(lat_out_b.shape[0]):
    x = np.arange(1,lat2.shape[1]-1)
    p2 = lat_out_b[i,:]
    xtest = np.array((lat_out_b.shape[1],lat_out_b.shape[1]+1))
    xtest1 = np.array((0,1))
    lat2[i+1,1:-1] = p2
    lat2[i+1,lat2.shape[1]-1:lat2.shape[1]] = extrap(xtest,x,p2)[1]
    lat2[i+1,0:1] = extrap(xtest1,x,p2)[0]
    
for i in range(lat2.shape[1]):
    x = np.arange(1,lat2.shape[0]-1)
    p2 = lat2[1:-1,i]
    xtest = np.array((lat_out_b.shape[0],lat_out_b.shape[0]+1))
    xtest1 = np.array((0,1))
    lat2[lat2.shape[0]-1:lat2.shape[0],i] = extrap(xtest,x,p2)[1]
    lat2[0:1,i] = extrap(xtest1,x,p2)[0]

for i in range(lat_out_b.shape[0]):
    x = np.arange(1,lat2.shape[1]-1)
    p2 = lon_out_b[i,:]
    xtest = np.array((lat_out_b.shape[1],lat_out_b.shape[1]+1))
    xtest1 = np.array((0,1))
    lon2[i+1,1:-1] = p2
    lon2[i+1,lat2.shape[1]-1:lat2.shape[1]] = extrap(xtest,x,p2)[1]
    lon2[i+1,0:1] = extrap(xtest1,x,p2)[0]
    
for i in range(lat2.shape[1]):
    x = np.arange(1,lat2.shape[0]-1)
    p2 = lon2[1:-1,i]
    xtest = np.array((lat_out_b.shape[0],lat_out_b.shape[0]+1))
    xtest1 = np.array((0,1))
    lon2[lat2.shape[0]-1:lat2.shape[0],i] = extrap(xtest,x,p2)[1]
    lon2[0:1,i] = extrap(xtest1,x,p2)[0]

out_pol = Polygon([(lon_out[0,0],lat_out[0,0]),(lon_out[0,lon_out.shape[1]-1],lat_out[0,lat_out.shape[1]-1]),
                   (lon_out[lon_out.shape[0]-1,lon_out.shape[1]-1],lat_out[lat_out.shape[0]-1,lat_out.shape[1]-1]),
                   (lon_out[lon_out.shape[0]-1,0],lat_out[lat_out.shape[0]-1,0])])

################################################# grid_in ################################################
for file_name in f_tropo:   
    
    print(file_name)
    file_name = file_name.strip()
    fh = Dataset(in_tropo+file_name, mode='r')
    orbit = file_name[:-3].split('_')[4]
    
    scan = fh.groups['instrument'].variables['scanline'][:]
    grpi = fh.groups['instrument'].variables['ground_pixel'][:]
    xch4 = fh.groups['target_product'].variables['xch4_corrected'][:]
    aker = fh.groups['target_product'].variables['xch4_column_averaging_kernel'][:]
    apri = fh.groups['target_product'].variables['ch4_profile_apriori'][:]
    qa   = fh.groups['diagnostics'].variables['qa_value'][:]
    lats = fh.groups['instrument'].variables['latitude_center'][:]
    lons = fh.groups['instrument'].variables['longitude_center'][:]
    time = fh.groups['instrument'].variables['time'][:]
    sp = fh.groups['meteo'].variables['surface_pressure'][:]
    dp = fh.groups['meteo'].variables['dp'][:]
    alts = fh.groups['meteo'].variables['altitude_levels'][:]
    sa = fh.groups['meteo'].variables['surface_altitude'][:]
    dry_air = fh.groups['meteo'].variables['dry_air_subcolumns'][:]
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
    aker_2d = np.zeros((12,l_scan,l_grpi))
    apri_2d = np.zeros((12,l_scan,l_grpi))
    sp_2d = np.zeros((l_scan,l_grpi))
    dp_2d = np.zeros((l_scan,l_grpi))
    lp_2d = np.zeros((13,l_scan,l_grpi))
    alt_2d = np.zeros((13,l_scan,l_grpi))
    sa_2d = np.zeros((l_scan,l_grpi))
    dry_air_2d = np.zeros((12,l_scan,l_grpi))
    lats_2d = np.zeros((l_scan,l_grpi))
    lons_2d = np.zeros((l_scan,l_grpi))
    
    for j in range(l_scan):
        for k in range(l_grpi):
            scan_2d[j][k] = scan[j*l_grpi+k]
            grpi_2d[j][k] = grpi[j*l_grpi+k]
            qa_2d[j][k] = qa[j*l_grpi+k]
            if qa_2d[j][k] == 1:    
                xch4_2d[j][k] = xch4[j*l_grpi+k]
                sp_2d[j][k] = sp[j*l_grpi+k]
                dp_2d[j][k] = dp[j*l_grpi+k]
                sa_2d[j][k] = sa[j*l_grpi+k]
                alt_2d[0][j][k] = alts[j*l_grpi+k][12]
            else:
                xch4_2d[j][k] = np.nan
                sp_2d[j][k] = np.nan
                dp_2d[j][k] = np.nan
                sa_2d[j][k] = np.nan
                alt_2d[0][j][k] = np.nan
            lats_2d[j][k] = lats[j*l_grpi+k]
            lons_2d[j][k] = lons[j*l_grpi+k]
            
    lp_2d[0] = sp_2d
    for l in range(12):
        for j in range(l_scan):
            for k in range(l_grpi):
                scan_2d[j][k] = scan[j*l_grpi+k]
                grpi_2d[j][k] = grpi[j*l_grpi+k]
                qa_2d[j][k] = qa[j*l_grpi+k]
                if qa_2d[j][k] == 1:    
                    aker_2d[l][j][k] = aker[j*l_grpi+k][l]
                    apri_2d[l][j][k] = apri[j*l_grpi+k][l]
                    dry_air_2d[l][j][k] = dry_air[j*l_grpi+k][l]
                    lp_2d[l+1][j][k] = sp_2d[j][k]-(l+1)*dp_2d[j][k]
                    alt_2d[l+1][j][k] = alts[j*l_grpi+k][11-l]
                else:
                    aker_2d[l][j][k] = np.nan
                    apri_2d[l][j][k] = np.nan
                    dry_air_2d[l][j][k] = np.nan
                    lp_2d[l+1][j][k] = np.nan
                    alt_2d[l+1][j][k] = np.nan
            
    lons_2d = np.clip(lons_2d, -180.0, 180.0) ##corner
    lats_2d = np.clip(lats_2d, -90.0, 90.0) ##corner

############################################### corner points ############################################
    lon_b = np.zeros((lats_2d.shape[0]-1,lons_2d.shape[1]-1))
    lat_b = np.zeros((lats_2d.shape[0]-1,lons_2d.shape[1]-1))
    for i in range(lats_2d.shape[0]-1):
        for j in range(lons_2d.shape[1]-1):
            lon1 = lons_2d[i:i+2,j:j+2]
            lat1 = lats_2d[i:i+2,j:j+2]
            if str(np.count_nonzero(np.isnan(lon1))) != '0':
                lon_b[i,j] = np.nan
                lat_b[i,j] = np.nan
            else:               
                p1 = Polygon([(lon1[0,0],lat1[0,0]),(lon1[0,1],lat1[0,1]),
                              (lon1[1,1],lat1[1,1]),(lon1[1,0],lat1[1,0])])
                lon_b[i,j] = p1.centroid.x
                lat_b[i,j] = p1.centroid.y
                
################################################ regridding ##############################################
    grid_in = {'lon': lon_b, 'lat': lat_b,
               'lon_b': lons_2d, 'lat_b': lats_2d}
    
    grid_out = {'lon': lon_out, 'lat': lat_out,
               'lon_b': lon2, 'lat_b': lat2}
        
    file_xe = file_name[:-3]+'_xe_v17_qa.nc'
    if os.path.exists(file_xe):
        regridder = xe.Regridder(grid_in, grid_out, 'bilinear', reuse_weights=True, filename=file_xe)
    else:
        regridder = xe.Regridder(grid_in, grid_out, 'bilinear', filename=file_xe)
    shutil.move(file_xe, out_xe)
    
    data_xch4 = regridder(xch4_2d[0:-1,0:-1])
    data_lpre = regridder(lp_2d[:,0:-1,0:-1])
    data_alts = regridder(alt_2d[:,0:-1,0:-1])
    data_aker = regridder(aker_2d[:,0:-1,0:-1])
    data_apri = regridder(apri_2d[:,0:-1,0:-1])
    data_dry_air = regridder(dry_air_2d[:,0:-1,0:-1])
    data_xch4[data_xch4 <= 0.0] = np.nan
    data_lpre[data_lpre <= 0.0] = np.nan
    data_alts[data_alts <= 0.0] = np.nan
    data_aker[data_aker <= 0.0] = np.nan
    data_apri[data_apri <= 0.0] = np.nan
    data_dry_air[data_dry_air <= 0.0] = np.nan
    
    df = xr.Dataset(
        data_vars = dict(
            xch4 = (["y","x"], data_xch4),
            p_levs = (["z","y","x"], data_lpre),
            h_alts = (["z","y","x"], data_alts),
            a_prio = (["nz","y","x"], data_apri),
            a_kern = (["nz","y","x"], data_aker),
            dry_air = (["nz","y","x"], data_dry_air),),
        coords = dict(
            lon = (["y","x"],lon_out),
            lat = (["y","x"],lat_out),),
        attrs = dict(time = str(date)),)  

    df.to_netcdf(out_tropo+file_name[:-3]+"_bilinear_v17_qa.nc")

    fig = plt.figure(figsize=(25, 15))
    ax = fig.add_subplot(111, frame_on=False)
    ax.set_xticks([])
    ax.set_yticks([])
    m1 = Basemap(projection='cyl', resolution='l', llcrnrlat=10.0, urcrnrlat=90.0,
                llcrnrlon=-50.0, urcrnrlon=60.0)
    m1.drawcoastlines(linewidth=1.0)
    m1.drawcountries(linewidth=1.0)
    m1.drawmapboundary()  
    m1.drawparallels(np.arange(-90., 120., 30), labels=[1, 0, 0, 0], fontsize=20)
    m1.drawmeridians(np.arange(-180., 181., 30), labels=[0, 0, 0, 1], fontsize=20)
    
    cs = m1.contourf(lons_2d, lats_2d, xch4_2d, levels=levels, cmap=cmap, norm=norm) 
    cs1 = m1.contourf(lon_out, lat_out, dat, alpha=0.6)
    cs2 = m1.contourf(lon_out, lat_out, data_xch4, cmap=cmap, norm=norm)  
    plt.title(date, loc='left', fontsize=30)
    cbar = m1.colorbar(cs, location='bottom', pad="5%", ticks=levels[::5], ax=ax)        
    cbar.set_label('xch4 [ppb]', fontsize=30)
    cbar.ax.tick_params(labelsize=30)
    fig.savefig(figs+'xch4_'+str(date)+'_'+orbit+'_v17_qa.png', bbox_inches='tight')   
    plt.close(fig)