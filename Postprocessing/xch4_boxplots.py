#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 14:46:11 2022

@author: angel
"""

import numpy as np
import pandas as pd
import os
from netCDF4 import Dataset
from wrf import interplevel
import matplotlib.pyplot as plt
import seaborn as sns

mon = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
mon_lb = ['Jan19', 'Feb19', 'Mar19', 'Apr18', 'May18', 'Jun18', 'Jul18', 'Aug18', 'Sep18',
          'Oct18', 'Nov18', 'Dec18']

path_s5p = '/home/angel/tropomi/regridded/'
path_mlev = '/home/angel/tropomi/mlevels/'
path_colu = '/home/angel/tropomi/column/'
figs = '/home/angel/tropomi/boxplots/'

pinf = float('+inf')

s5p_all_0 = np.array([])
mlev_all_0 = np.array([])
colu_all_0 = np.array([])
s5p_stats = np.array([])
mlev_stats = np.array([])
colu_stats = np.array([])
xps_stats = np.array([])
hue_s5p_stats = np.array([])
hue_colu_stats = np.array([])
hue_mlev_stats = np.array([])

c = 0
for m in range(len(mon)):

    c = c+1
    f_s5p = sorted(os.listdir(path_s5p+mon[m]))
    f_mlev = sorted(os.listdir(path_mlev+mon[m]))
    f_colu = sorted(os.listdir(path_colu+mon[m]))

    for s5p, mlev, colu in zip(f_s5p, f_mlev, f_colu):
    
        print(s5p, mlev)
        s5p = s5p.strip()
        fh_s5p = Dataset(path_s5p+mon[m]+'/'+s5p, mode='r')
        mlev = mlev.strip()
        fh_mlev = Dataset(path_mlev+mon[m]+'/'+mlev, mode='r')
        colu = colu.strip()
        fh_colu = Dataset(path_colu+mon[m]+'/'+colu, mode='r')
        
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
        colu_1d_nonan = colu_1d[np.logical_not(np.isnan(colu_1d))]
        xps_1d = np.ones(obs_1d_nonan.shape[0])
        xps_1d = xps_1d*c
        hue_1d = np.ones(obs_1d_nonan.shape[0])
        hue_colu = hue_1d*5        # wrf code
        hue_s5p = hue_1d*6         # s5p code
        hue_mlev = hue_1d*7        # wrf sth code
        
        dif_mlev_nonan = dif_mlev_1d[np.logical_not(np.isnan(dif_mlev_1d))]
        dif_colu_nonan = dif_colu_1d[np.logical_not(np.isnan(dif_colu_1d))]

        s5p_all_0 = np.concatenate([s5p_all_0, obs_1d])
        mlev_all_0 = np.concatenate([mlev_all_0, mlev_1d])
        colu_all_0 = np.concatenate([colu_all_0, colu_1d])   
    
        s5p_stats = np.concatenate([s5p_stats, obs_1d_nonan])
        mlev_stats = np.concatenate([mlev_stats, mlev_1d_nonan])
        colu_stats = np.concatenate([colu_stats, colu_1d_nonan])
        xps_stats = np.concatenate([xps_stats, xps_1d])
        hue_s5p_stats = np.concatenate([hue_s5p_stats, hue_s5p])
        hue_colu_stats = np.concatenate([hue_colu_stats, hue_colu])
        hue_mlev_stats = np.concatenate([hue_mlev_stats, hue_mlev])        
        
        s5p_wrf_stats = np.concatenate([colu_stats, s5p_stats])
        s5p_wrf_stats = np.concatenate([s5p_wrf_stats, mlev_stats])         # 1
        xps_mon_stats = np.concatenate([xps_stats, xps_stats])
        xps_mon_stats = np.concatenate([xps_mon_stats, xps_stats])          # 2
        hue_mon_stats = np.concatenate([hue_colu_stats, hue_s5p_stats])
        hue_mon_stats = np.concatenate([hue_mon_stats, hue_mlev_stats])     # 3
        
    A = np.zeros((3,s5p_wrf_stats.shape[0]))
    A[0] = s5p_wrf_stats 
    A[1] = xps_mon_stats
    A[2] = hue_mon_stats    
    A = np.transpose(A)
    df_A = pd.DataFrame(A, columns=['XCH4','Xp','SRC'])
    fig, ax = plt.subplots(figsize=(25,10), dpi=80)
    bp = sns.boxplot(x='Xp', y='XCH4', data=df_A, hue='SRC')
    ax.set_xticklabels(mon_lb, fontsize=26)  
    ax.set_xlabel('')
    ax.set_xlim([-0.60, 11.60])
    ax.set_ylim([1700, 2000])
    ax.set_yticklabels(['1700', '1750', '1800', '1850', '1900', '1950', '2000'],fontsize=24)  
    ax.set_ylabel(r' XCH$_{4}$ [ppb]',fontsize=24)
    ax.vlines(3.5, 1700, 2000, linestyles='dashed', colors='black', alpha=0.75)
    ax.vlines(7.5, 1700, 2000, linestyles='dashed', colors='black', alpha=0.75)
    lb = [r' WRF-GHG', r' SRON', r' Sth. WRF-GHG']
    bp.legend(title='', fontsize=24, loc='upper left', borderaxespad=0.2)
    n = 0
    for i in lb:
        bp.legend_.texts[n].set_text(i)
        n += 1
    plt.show()
    plt.close(fig) 
