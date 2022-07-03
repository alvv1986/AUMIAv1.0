# MethaneAU
This reposity proposes a WRF-based modeling tool for CH4 inversion over Europe. Prior to proceed with the tasks below, select a study period with availability of TROPOMI CH4 fields. The forward and backward models are the WRF-GHG and WRF-STILT, while the a-priori emissions are taken from the EDGAR model version 6. IC/CC for the forward modeling are based on ECMWF CAMS EAC4 fields.

1. Run the WRF WPS for a given study period using ECMWF ERA5 fields

2. Creating a CH4 a-priori emission file in the proper WRF netcdf file format

For this example we consider CH4 monthly gridmaps (https://edgar.jrc.ec.europa.eu/dataset_ghg60#p2). Click to expand on the option Annual sector-specific gridmaps (1970-2018) and montlhy  sector-specific gridmaps (2000-2018), and then on CH4 for each sector, if available. 24 different sectors are available for 2018: 

    ENE: Power industry
    
    REF_TRF: Oil refineries and transformation industry
    
    IND: Combustion for manufacturing
    
    RCO: Energy for Buildings
    
    PRO_COAL: Fuel exploitation coal
    
    PRO: Fuel exploitation
    
    PRO_OIL: Fuel exploitation oil
    
    PRO_GAS: Fuel exploitation gas
    
    TRO_noRES: Road transportation no resuspension
    
    TNR_Other: Railways, pipelines, off-road transport
    
    TNR_Aviation_CDS: Aviation climbing & descent
    
    TNR_Aviation_CRS: Aviation cruise
    
    TNR_Aviation_LTO: Aviation landing & take-off 
    
    TNR_Ship: Shipping
    
    CHE: Chemical processes
    
    IRO: Iron and steel production
    
    ENF: Enteric fermentation
    
    MNM: Manure management
    
    AWB: Agricultural waste burning
    
    AGS: Agricultural soils burning
    
    SWD_LDF: Solid waste landfills
    
    SWD_INC: Solid waste incineration
    
    WWT: Waste water handling
    
    FFF: Fosil fuel fires

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

with something similar for the other sectors. Run the scripts ``EDGARtoAE.py`` to create a EDGAR like , and ``CH4.py`` to make a quick visualization. Now you should be ready to run the anthro_emis, ``./anthro_emis < anthro_ghg.inp``

3. Interpolating background CH4 global concentrations to the WRF-GHG initial and boundary conditions 

Background methane concentrations for Europe are obtained from CAMS global reanalysis EAC4 fields (https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-reanalysis-eac4?tab=form). In addition to methane (chemistry), CAMS EAC4 surface pressure fields are also required for interpolation.

- Create a txt file (``levels.txt``) containing the model levels according to the selected CAMS product, L60 in this example. For CAMS EAC4 products, the model levels (1 to 60) can be obtanied from https://confluence.ecmwf.int/display/UDOC/L60+model+level+definitions.

- Set up the parameters indir, outdir and ab in the ncl script MACC_BC2MOZART_CH4.ncl, and then run it by typing ``ncl fileid="03062022" MACC_BC2MOZART_CH4.ncl``. Name fileid whatever you want. This script is a modified version of the original one at https://confluence.ecmwf.int/pages/viewpage.action?pageId=174865233.

- Set your mozbc namelist file according to the example in GHG.inp and then run mozbc, ``./mozbc < mozbc_ghg.inp``.

4. Run the WRF-GHG model


