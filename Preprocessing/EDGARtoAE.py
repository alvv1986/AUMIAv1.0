#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xarray as xr
import numpy as np
from pathlib import Path

def add_date_datesec(edgar_file, year=2018):
    '''
    Add the 'date' and 'datesec' variables into
    edgar emission file

    Parameters
    ----------
    edgar_file : str
        edgar monthly emission file.
    year : int, optional
        year of the emission file. The default is 2015

    Returns
    -------
    edgar_date: xarray Dataset
        emission file with 'date' and 'datesec' variable
    '''
    month = int(edgar_file.split("2018_")[1].split("_")[0])
    time = np.array([month - 1])
    date = np.array([year * 10000 + month * 100 + 1])

    edgar_date = xr.open_dataset(edgar_file)
    var_name = list(edgar_date.data_vars)[0]
    edgar_date = edgar_date.rename({var_name:"emis_tot"})

    edgar_date["time"] = time
    edgar_date["date"] = xr.DataArray(date, dims=["time"])
    edgar_date["datesec"] = xr.DataArray(np.array([0]), dims=["time"])
    return edgar_date

def concat_sector_by_month(sector_path):
    '''
    Load and concat monthly emission into one netcdf

    Parameters
    ----------
    sector_path : str
        path of pollutant emission sector folder

    Returns
    -------
    year_emiss: xarray Dataset
        montly emission into one dataset

    '''
    sector = Path(sector_path)
    sector_files = [str(file) for file in list(sector.glob("*.nc"))]
    ds = {int(f.split("2018_")[1].split("_")[0]): add_date_datesec(f) 
            for f in sector_files}
    ds_sorted = dict(sorted(ds.items()))
    year_emiss = xr.concat(list(ds_sorted.values()), dim="time")
    return year_emiss

def join_pol_by_sector(pol_path, total=False):
    '''
    Read sector emissions and save them in one dict. If total = True
    it returns the total emission

    Parameters
    ----------
    pol_path : str
        Polutant emissions folder where sectors are
    total : Bool, optional
        sum all the sector emission. The default is False
    
    Returns
    -------
    pol_emi : xarray Dataset
       Total emission (sum of all sources)
    pol_by_sector: dictionary
        Dictionary with sectors monthtly emission in one dataset
    '''
    pol_folder = Path(pol_path)
    pol_sectors = [str(folder) for folder in pol_folder.glob("*/")]
    pol_by_sector = {sector.split("/")[-1]: concat_sector_by_month(sector) 
                     for sector in pol_sectors}
    if total:
        pol_emi = sum(pol_by_sector.values())
        pol_emi["date"] = pol_by_sector[list(pol_by_sector.keys())[1]].date
        return pol_emi
    else:
        return pol_by_sector

def group_sectors(pol_by_sector, sectors):
    '''
    Group emissions sector in a single group (i.e. IPCC 96)

    Parameters
    ----------
    pol_by_sector : dict
        Dictionary with sectors monthly emission in onde dataset
    sectors : list
        List of sectors to merge
    '''
    group = {sector: pol_by_sector[sector] for sector in sectors 
             if sector in pol_by_sector.keys()}
    group_total = sum(group.values())
    group_total["date"] = pol_by_sector[list(pol_by_sector.keys())[1]].date
    return group_total

def write_netcdf_toAE(emiss, pollutant, sec, prefix="v6.0_",
                      suffix=".0.1x0.1.nc"):
    '''
    Export emission dataset into netcdf

    Parameters
    ----------
    emiss : xarray Dataset
        Emission dataset to export to netcdf
    pollutatnt : str
        Name of emitted pollutant
    sec : str
        Sector       
    prefix : str, optional
        prefix used in anthro emiss. Default is "v6.0_"
    suffix : str, optional
        suffix to add to the file name. Default is ".0.1x0.1.nc"
    '''
    file_name = (prefix + pollutant + "_" + sec + "_2018" +  suffix)
    print(file_name)
    emiss['emis_tot'].attrs['units'] = 'kg m-2 s-1'
#    emiss['emis_tot'].attrs['_FillValue'] = 0.
    emiss.to_netcdf(file_name, unlimited_dims={"time": True})

# Group sectors 
energy      = ["ENE", "REF_TRF", "IND", "RCO", "PRO_COAL", "PRO", "PRO_OIL", "PRO_GAS"]
transport   = ["TRO_noRES", "TNR_Other"]
aviation    = ["TNR_Aviation_CDS", "TNR_Aviation_CRS", "TNR_Aviation_LTO"]
shipping    = ["TNR_Ship"]
industrial  = ["CHE", "IRO"]
agriculture = ["ENF", "MNM", "AWB", "AGS"]
waste       = ["SWD_LDF", "SWD_INC", "WWT"]
fires       = ["FFF"]
tot_no_fire = ["ENE", "REF_TRF", "IND", "RCO", "PRO_COAL", "PRO", "PRO_OIL", "PRO_GAS",        
               "TRO_noRES", "TNR_Other",
               "TNR_Aviation_CDS", "TNR_Aviation_CRS", "TNR_Aviation_LTO",
               "TNR_Ship",
               "CHE", "IRO",
               "ENF", "MNM", "AWB", "AGS",
               "SWD_LDF", "SWD_INC", "WWT"]

pol_path = "/home/angel/tropomi/CH4"

ch4_all = join_pol_by_sector(pol_path, total=True)
write_netcdf_toAE(ch4_all, "CH4", "ALL")

ch4_sec = join_pol_by_sector(pol_path)
ch4_energy      = group_sectors(ch4_sec, energy)
ch4_transport   = group_sectors(ch4_sec, transport)
ch4_aviation    = group_sectors(ch4_sec, aviation)
ch4_shipping    = group_sectors(ch4_sec, shipping)
ch4_industrial  = group_sectors(ch4_sec, industrial)
ch4_agriculture = group_sectors(ch4_sec, agriculture)
ch4_waste       = group_sectors(ch4_sec, waste)
ch4_fires       = group_sectors(ch4_sec, fires)
ch4_tot_no_fire = group_sectors(ch4_sec, tot_no_fire)

write_netcdf_toAE(ch4_energy, "CH4", "ENERGY")
write_netcdf_toAE(ch4_transport, "CH4", "TRANSPORT")
write_netcdf_toAE(ch4_aviation, "CH4", "AVIATION")
write_netcdf_toAE(ch4_shipping, "CH4", "SHIPPING")
write_netcdf_toAE(ch4_industrial, "CH4", "INDUSTRIAL")
write_netcdf_toAE(ch4_agriculture, "CH4", "AGRICULTURE")
write_netcdf_toAE(ch4_waste, "CH4", "WASTE")
write_netcdf_toAE(ch4_fires, "CH4", "FIRES")
write_netcdf_toAE(ch4_tot_no_fire, "CH4", "TOTAL_NO_FIRES")
