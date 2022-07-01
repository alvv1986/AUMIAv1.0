# MethaneAU
This reposity proposes a modeling tool for methane inversion over Europe. 

1. Creating a CH4 a-priori emission file in the proper WRF netcdf file format

- Download EDGAR CH4 data 
    
    For this example we consider CH4 monthly gridmaps from https://edgar.jrc.ec.europa.eu/dataset_ghg60#p2. Click to expand on 
the option Annual sector-specific gridmaps (1970-2018) and montlhy sector-specific gridmaps (2000-2018), and then on CH4 for each sector, if 
available.

- 24 different sectors are available for 2018

    "ENE", "REF_TRF", "IND", "RCO", "PRO_COAL", "PRO", "PRO_OIL", "PRO_GAS", "TRO_noRES", "TNR_Other", "TNR_Aviation_CDS", "TNR_Aviation_CRS", 
"TNR_Aviation_LTO", "TNR_Ship", "CHE", "IRO", "ENF", "MNM", "AWB", "AGS", "SWD_LDF", "SWD_INC", "WWT" and "FFF"

    Each sector being assigned a folder with the same name, and containing 12 nc files (monthly). Set ``pol_path`` in the script ``EDGARtoAE.py`` 
(``CH4`` in this example contains all the 24 directories describe in step 2), e.g. for ENE (Power Industry) we should have:

    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_1_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_2_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_3_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_4_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_5_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_6_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_7_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_8_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_9_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_10_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_11_ENE.0.1x0.1.nc
    
    /home/angel/Documents/iag-usp/modis/methane_project/CH4/ENE/v6.0_CH4_2018_12_ENE.0.1x0.1.nc

    and something similar for the other sectors

- Run the script EDGARtoAE.py

- Run the WRF-GHG WPS for the period of interest

- Now you should be ready to run the anthro_emis

2. Interpolating CAMS CH4 fields to the WRF-GHG initial and boundary conditions 

- Download CAMS data 
 
    CAMS methane (chemistry) and surface pressure fields can be obtained from https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-reanalysis-eac4?tab=form.

- 
