#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:46:11 2022

@author: angel
"""

import numpy as np
import os
from netCDF4 import Dataset
from wrf import interplevel
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
import matplotlib.colors as colors
from scipy import stats

path_s5p = '/home/angel/tropomi/regridded/'
path_mlev = '/home/angel/tropomi/mlevels/'
path_colu = '/home/angel/tropomi/mcolumn/'
figs = '/home/angel/tropomi/monthly_maps/'

lon_0 = 5.66
lat_0 = 51.98
pinf = float('+inf')

cmap = {'(a) SRON XCH$_{4}$ [ppb]': colors.LinearSegmentedColormap.from_list(
        'white11red', ['#FFFFFF','paleturquoise','lightskyblue','aqua',
               'lightgreen','limegreen','yellow','goldenrod','darkorange',
               'red','darkred']),
        '(b) WRF-GHG XCH$_{4}$ [ppb]': colors.LinearSegmentedColormap.from_list(
        'white11red', ['#FFFFFF','paleturquoise','lightskyblue','aqua',
               'lightgreen','limegreen','yellow','goldenrod','darkorange',
               'red','darkred']),
        '(c) \u0394XCH$_{4,\\frac{model-satellite}{satellite}}$ [%]': colors.LinearSegmentedColormap.from_list(
        'white11red', ['darkblue','blue','deepskyblue','skyblue','#FFFFFF',
               'lightsalmon','darksalmon','red','darkred']),
        '(d) SRON XCH$_{4}$ [ppb]': colors.LinearSegmentedColormap.from_list(
        'white11red', ['#FFFFFF','paleturquoise','lightskyblue','aqua',
               'lightgreen','limegreen','yellow','goldenrod','darkorange',
               'red','darkred']),
        '(e) WRF-GHG XCH$_{4}$ [ppb]': colors.LinearSegmentedColormap.from_list(
        'white11red', ['#FFFFFF','paleturquoise','lightskyblue','aqua',
               'lightgreen','limegreen','yellow','goldenrod','darkorange',
               'red','darkred']),
        '(e) \u0394XCH$_{4,\\frac{model-satellite}{satellite}}$ [%]': colors.LinearSegmentedColormap.from_list(
        'white11red', ['darkblue','blue','deepskyblue','skyblue','#FFFFFF',
               'lightsalmon','darksalmon','red','darkred'])}

lev = {'(a) SRON XCH$_{4}$ [ppb]': MaxNLocator(nbins=30).tick_values(1700,2000), 
       '(b) WRF-GHG XCH$_{4}$ [ppb]': MaxNLocator(nbins=30).tick_values(1700,2000), 
       '(c) \u0394XCH$_{4,\\frac{model-satellite}{satellite}}$ [%]': MaxNLocator(nbins=50).tick_values(-5,5),
       '(d) SRON XCH$_{4}$ [ppb]': MaxNLocator(nbins=30).tick_values(1700,2000), 
       '(e) WRF-GHG XCH$_{4}$ [ppb]': MaxNLocator(nbins=30).tick_values(1700,2000), 
       '(f) \u0394XCH$_{4,\\frac{model-satellite}{satellite}}$ [%]': MaxNLocator(nbins=50).tick_values(-5,5)}

f_s5p = sorted(os.listdir(path_s5p))
f_mlev = sorted(os.listdir(path_mlev))
f_colu = sorted(os.listdir(path_colu))

s5p_all_0 = np.array([])
mlev_all_0 = np.array([])
dif_mlev_0 = np.array([])
colu_all_0 = np.array([])
dif_colu_0 = np.array([])
s5p_stats = np.array([])
mlev_stats = np.array([])
dif_mlev_stats = np.array([])
colu_stats = np.array([])
dif_colu_stats = np.array([])

for s5p, mlev, colu in zip(f_s5p, f_mlev, f_colu):
    
    print(s5p, mlev)
    s5p = s5p.strip()
    fh_s5p = Dataset(path_s5p+s5p, mode='r')
    mlev = mlev.strip()
    fh_mlev = Dataset(path_mlev+mlev, mode='r')
    colu = colu.strip()
    fh_colu = Dataset(path_colu+colu, mode='r')
    
    date = mlev[:-3].split('_')[2]
    lon = fh_s5p['lon'][:]
    lat = fh_s5p['lat'][:]
    
    ch4_s5p = fh_s5p['xch4'][:]
    plevs_s5p = fh_s5p['p_levs'][:]
    plevs_s5p = plevs_s5p*100.
    alts_s5p = fh_s5p['h_alts'][:]
    akern_s5p = fh_s5p['a_kern'][:]
    aprio_s5p = fh_s5p['a_prio'][:]
    dry_air_s5p = fh_s5p['dry_air'][:]
    
    P = fh_mlev['P'][0][:]
    PB = fh_mlev['PB'][0][:]
    Pabs = P+PB
    PSFC = fh_mlev['PSFC'][0][:]
    PTOP = fh_mlev['P_TOP'][0]
    ch4_wrf_ant_mlev = fh_mlev['CH4_ANT'][0][:]
    ch4_wrf_bio_mlev = fh_mlev['CH4_BIO'][0][:]
    ch4_wrf_bbu_mlev = fh_mlev['CH4_BBU'][0][:]
    ch4_wrf_bck_mlev = fh_mlev['CH4_BCK'][0][:]
    ch4_wrf_all_mlev = (ch4_wrf_ant_mlev+ch4_wrf_bio_mlev+ch4_wrf_bbu_mlev+ch4_wrf_bck_mlev)*1000.
    ch4_wrf_all_colu = fh_colu['CH4_ALL'][0][:]
    
    ch4_interp = np.zeros((plevs_s5p.shape[0],plevs_s5p.shape[1],plevs_s5p.shape[2]))
    for k in range(ch4_interp.shape[0]):
        ch4_interp[k] = interplevel(ch4_wrf_all_mlev,Pabs,plevs_s5p[k])
 
    dh_s5p = np.zeros((alts_s5p.shape[0]-1,alts_s5p.shape[1],alts_s5p.shape[2]))
    for k in range(dh_s5p.shape[0]):
        dh_s5p[k] = alts_s5p[k+1]-alts_s5p[k]
    
    Pb = np.zeros((plevs_s5p.shape[0]-1,plevs_s5p.shape[1],plevs_s5p.shape[2]))
    Pa = np.zeros((plevs_s5p.shape[0]-1,plevs_s5p.shape[1],plevs_s5p.shape[2]))
    Pc = np.zeros((2,plevs_s5p.shape[1],plevs_s5p.shape[2]))
    Pb[0] = plevs_s5p[0]
    Pa[0] = (plevs_s5p[0]+plevs_s5p[1])*0.5
    Pc[0] = (Pa[0]+Pb[0])*0.5
    for k in range(plevs_s5p.shape[0]-3):
        Pb[k+1] = Pa[k]
        Pa[k+1] = 2*plevs_s5p[k+1]-Pb[k+1]  
    Pb[11] = Pa[10]
    Pa[11] = plevs_s5p[12]
    Pc[1] = (Pa[11]+Pb[11])*0.5
     
    ch4_smooth = np.zeros((ch4_interp.shape[0]-1,ch4_interp.shape[1],ch4_interp.shape[2]))
    ch4_limits = np.zeros((2,ch4_interp.shape[1],ch4_interp.shape[2]))
    ch4_limits[0] = interplevel(ch4_wrf_all_mlev,Pabs,Pc[0])
    ch4_limits[1] = interplevel(ch4_wrf_all_mlev,Pabs,Pc[1])
    A1 = np.zeros((ch4_smooth.shape[0],ch4_smooth.shape[1],ch4_smooth.shape[2]))
    A2 = np.zeros((ch4_smooth.shape[0],ch4_smooth.shape[1],ch4_smooth.shape[2]))
    for k in range(ch4_smooth.shape[0]):
        if k==0:
            A1[k] = akern_s5p[k]*ch4_limits[0]
            A2[k] = (np.ones((119,119))-akern_s5p[k])*aprio_s5p[k]*(10**9)/(dry_air_s5p[k]*dh_s5p[k])
            ch4_smooth[k] = A1[k]+A2[k]             
        elif k>0 and k<11:
            A1[k] = akern_s5p[k]*ch4_interp[k]
            A2[k] = (np.ones((119,119))-akern_s5p[k])*aprio_s5p[k]*(10**9)/(dry_air_s5p[k]*dh_s5p[k])
            ch4_smooth[k] = A1[k]+A2[k]
        else:
            A1[k] = akern_s5p[k]*ch4_limits[1]
            A2[k] = (np.ones((119,119))-akern_s5p[k])*aprio_s5p[k]*(10**9)/(dry_air_s5p[k]*dh_s5p[k])
            ch4_smooth[k] = A1[k]+A2[k]
    
    xch4 = np.zeros((ch4_interp.shape[0]-1,ch4_interp.shape[1],ch4_interp.shape[2]))
    for k in range(xch4.shape[0]):
        xch4[k] = ch4_smooth[k]*(Pb[k]-Pa[k])/(plevs_s5p[0]-plevs_s5p[12])

    xch4_col = np.nansum(xch4, axis=0)
    
    xch4_col = np.where(np.isnan(ch4_s5p), np.nan, xch4_col)
    ch4_wrf_all_colu = np.where(np.isnan(ch4_s5p), np.nan, ch4_wrf_all_colu)

    dif_mlev_s5p = xch4_col-ch4_s5p
    dif_mlev_s5p = (dif_mlev_s5p/ch4_s5p)*100
    dif_colu_s5p = ch4_wrf_all_colu-ch4_s5p
    dif_colu_s5p = (dif_colu_s5p/ch4_s5p)*100
    
    obs_1d = np.ravel(ch4_s5p)
    mlev_1d = np.ravel(xch4_col)
    colu_1d = np.ravel(ch4_wrf_all_colu)
    dif_mlev_1d = np.ravel(dif_mlev_s5p)
    dif_colu_1d = np.ravel(dif_colu_s5p)
    
    obs_1d_nonan = obs_1d[np.logical_not(np.isnan(obs_1d))]
    mlev_1d_nonan = mlev_1d[np.logical_not(np.isnan(mlev_1d))]
    dif_mlev_nonan = dif_mlev_1d[np.logical_not(np.isnan(dif_mlev_1d))]
    colu_1d_nonan = colu_1d[np.logical_not(np.isnan(colu_1d))]
    dif_colu_nonan = dif_colu_1d[np.logical_not(np.isnan(dif_colu_1d))]

    s5p_all_0 = np.concatenate([s5p_all_0, obs_1d])
    mlev_all_0 = np.concatenate([mlev_all_0, mlev_1d])
    dif_mlev_0 = np.concatenate([dif_mlev_0, dif_mlev_1d])
    colu_all_0 = np.concatenate([colu_all_0, colu_1d])
    dif_colu_0 = np.concatenate([dif_colu_0, dif_colu_1d])    
    
    s5p_stats = np.concatenate([s5p_stats, obs_1d_nonan])
    mlev_stats = np.concatenate([mlev_stats, mlev_1d_nonan])
    dif_mlev_stats = np.concatenate([dif_mlev_stats, dif_mlev_nonan])
    colu_stats = np.concatenate([colu_stats, colu_1d_nonan])
    dif_colu_stats = np.concatenate([dif_colu_stats, dif_colu_nonan])  
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(obs_1d_nonan, mlev_1d_nonan)
    x1=np.linspace(1750, 1950, 1000)
    y1=slope*x1+intercept        
    plt.plot(obs_1d_nonan, mlev_1d_nonan, 'r*', linewidth=1)
    plt.plot(x1,y1,'-',color='r',linewidth=2)
    plt.plot([1750, 1950], [1750, 1950], 'r--', linewidth=1)
    plt.xlim(1750, 1950)
    plt.ylim(1750, 1950)
    plt.show()
    
s5p_all_1 = np.reshape(s5p_all_0, (len(f_s5p), np.shape(ch4_s5p)[0], np.shape(ch4_s5p)[1]))
mlev_all_1 = np.reshape(mlev_all_0, (len(f_mlev), np.shape(xch4_col)[0], np.shape(xch4_col)[1]))
dif_mlev_1 = np.reshape(dif_mlev_0, (len(f_mlev), np.shape(dif_mlev_s5p)[0], np.shape(dif_mlev_s5p)[1]))
colu_all_1 = np.reshape(colu_all_0, (len(f_colu), np.shape(ch4_wrf_all_colu)[0], np.shape(ch4_wrf_all_colu)[1]))
dif_colu_1 = np.reshape(dif_colu_0, (len(f_colu), np.shape(dif_colu_s5p)[0], np.shape(dif_colu_s5p)[1]))

s5p_all_2 = np.nanmean(s5p_all_1, axis=0)
mlev_all_2 = np.nanmean(mlev_all_1, axis=0)
dif_mlev_2 = np.nanmean(dif_mlev_1, axis=0)
colu_all_2 = np.nanmean(colu_all_1, axis=0)
dif_colu_2 = np.nanmean(dif_colu_1, axis=0)

ch4_1 = {'(a) SRON XCH$_{4}$ [ppb]': s5p_all_2, 
         '(b) Sth. WRF-GHG XCH$_{4}$ [ppb]': mlev_all_2, 
         '(c) \u0394XCH$_{4,\\frac{model-satellite}{satellite}}$ [%]': dif_mlev_2,
         '(e) SRON XCH$_{4}$ [ppb]': s5p_all_2, 
         '(f) WRF-GHG XCH$_{4}$ [ppb]': colu_all_2, 
         '(g) \u0394XCH$_{4,\\frac{model-satellite}{satellite}}$ [%]': dif_colu_2}

white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)

fig = plt.figure(figsize=(40,30))
axs = fig.subplots(2,4,gridspec_kw={"width_ratios":[1,1,1,1],"height_ratios":[0.9,0.9]})

y = np.zeros((2,s5p_stats.shape[0]))
y[0] = mlev_stats 
y[1] = colu_stats
lbs = ['(d)', '(h)']

c = 0
for row_n in range(axs.shape[0]): 
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(s5p_stats, y[row_n])
    x1=np.linspace(1750, 1950, 1000)
    y1=slope*x1+intercept    
    ab = np.vstack([s5p_stats, y[row_n]])
    z = stats.gaussian_kde(ab)(ab)
    textstr = '\n'.join((
        r'$n=%.0f$' % (len(s5p_stats),),
        r'$R=%.2f$' % (r_value, ),
        ))
    cb = axs[row_n,3].scatter(s5p_stats, y[row_n], marker = ".", c=z, s=90, cmap='jet', vmin=0.00020, vmax=0.00160)
    cbar = plt.colorbar(cb, ax=axs[row_n,3], ticks = [0.0003, 0.0005, 0.0007, 0.0009, 0.0011, 0.0013, 0.0015])
    cbar.ax.set_yticklabels(['3', '5', '7', '9', '11', '13', '15'])
    cbar.ax.tick_params(labelsize=18)
    axs[row_n,3].plot(x1,x1,'--',color='black',linewidth=1)
    axs[row_n,3].plot(x1,y1,'-',color='r',linewidth=2)
    axs[row_n,3].axis([1750, 1950, 1750, 1950])    
    
    for axis in ['top','bottom','left','right']:
        axs[row_n,3].spines[axis].set_linewidth(2.0)
    if row_n == 0:
        axs[row_n,3].set_xticks([int(k) for k in np.linspace(1750,1950,3)])
        axs[row_n,3].set_yticks([int(k) for k in np.linspace(1750,1950,3)])
    else:
        axs[row_n,3].set_xticks([round(k,1) for k in np.linspace(1750,1950,3)])
        axs[row_n,3].set_yticks([round(k,1) for k in np.linspace(1750,1950,3)])

    axs[row_n,3].tick_params(axis='both', which='major', labelsize=18)
    axs[row_n,3].text(0.14, 1.002, 'Density [x10$^{-4}$]', transform=axs[row_n,3].transAxes, 
                      ha='left', va='bottom', fontsize=25)
    axs[row_n,3].text(0.0005, 1.015, lbs[row_n], transform=axs[row_n,3].transAxes, 
                      ha='left', va='bottom', fontsize=25)
    axs[row_n,3].text(0.04, 0.98, r'$n=%.0f$' % (len(s5p_stats),), 
    transform=axs[row_n,3].transAxes, ha='left', va='top', fontsize=25)
    axs[row_n,3].text(0.04, 0.91, r'$r$ = {:.2f}'.format(np.corrcoef(s5p_stats, y[row_n])[0, 1]), 
    transform=axs[row_n,3].transAxes, ha='left', va='top', fontsize=25)
    
    for col_n in range(axs.shape[1]-1):
        for axis in ['top','bottom','left','right']:
            axs[row_n,col_n].spines[axis].set_linewidth(4.0)                
        norm = BoundaryNorm(list(lev.values())[c], ncolors=list(cmap.values())[c].N, clip=True) 
        m = Basemap(width = 3600000, height = 3600000,
                    projection = 'stere', lat_0 = lat_0, lon_0 = lon_0, resolution = 'l', ax=axs[row_n,col_n])
        m.drawcoastlines(linewidth=1.5)
        m.drawstates(linewidth=1.5)  
        m.drawcountries(linewidth=1.5)
        m.drawmapboundary() 
        x1, y1 = m(lon, lat)
        x1[x1 == pinf] = 0
        y1[y1 == pinf] = 0
        cs = m.contourf(x1, y1, list(ch4_1.values())[c], levels=list(lev.values())[c], cmap=list(cmap.values())[c], norm=norm)
        if col_n == 2:
            cbar = m.colorbar(cs, location='bottom', ax=axs[row_n,col_n], pad="2%", ticks=list(lev.values())[col_n][::5])
        else:
            cbar = m.colorbar(cs, location='bottom', ax=axs[row_n,col_n], pad="2%", ticks=list(lev.values())[col_n][::6])
        cbar.ax.tick_params(labelsize=18)
        axs[row_n,col_n].set_title(list(ch4_1.keys())[c], loc='left', fontsize=25)
        c = c+1
    
plt.subplots_adjust(top=0.5, bottom=0.1, left=0, right=0.6, hspace=0.2, wspace=0.2)
fig.savefig(figs+'ch4_mean_MaytoAug.png', bbox_inches='tight')
plt.show() 

from math import sqrt

MBE_mlev = np.mean(mlev_stats-s5p_stats)
MSE_mlev = np.square(np.subtract(s5p_stats, mlev_stats)).mean()
RMSE_mlev = sqrt(MSE_mlev)

MBE_colu = np.mean(colu_stats-s5p_stats)
MSE_colu = np.square(np.subtract(s5p_stats, colu_stats)).mean()
RMSE_colu = sqrt(MSE_colu) 
    
