# AUMIAv1.0
This methodology proposes a satellite-based tool for the quantification of CH<sub>4</sub> emissions over Europe. The WRF-based models WRF-GHG and WRF-STILT were selected for the forward and backward modeling; the WRF-STILT implementation is currently under development. The CH<sub>4</sub> space observations are based on TROPOMI measurements, in this case the SRON RemoTeC-S5P XCH4 scientific product version 17 which represents the most recent improvements to the TROPOMI operational product. The a-priori emissions are taken from the EDGAR model version 6. IC/BC for the forward modeling are based on ERA5 fields for meteorology and on CAM-chem fields for CH<sub>4</sub> concentration. Prior to proceed with the tasks below, make sure you have the NCAR utilities anthro_emis, fire_emis and mozbc properly installed.

1. Run the WRF WPS for a given study period using ECMWF ERA5 fields

Download, via the Climate Data Store Application Program Interface (cdsapi), met fields for both surface and vertical levels by running ``GetERA5-sl.py`` and ``GetERA5-pl.py``. Then, run ``geogrid.exe``, ``ungrib.exe`` and ``metgrid.exe`` as usually. For information on the grid configuration look up the ``namelist.wps`` file.

2. Create a CH<sub>4</sub> a-priori emission file in the proper WRF netcdf file format

For this example we consider CH<sub>4</sub> monthly gridmaps (https://edgar.jrc.ec.europa.eu/dataset_ghg60#p2). Click to expand on the option Annual sector-specific gridmaps (1970-2018) and montlhy sector-specific gridmaps (2000-2018), and then on CH<sub>4</sub> for each sector, if available. 24 different sectors are available for 2018: 

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

Each sector is assigned a folder with the same name, and contains 12 nc files (monthly). Set ``pol_path`` in the script ``EDGARtoAE.py`` 
(``CH4`` in this example contains all the 24 directories describe above), e.g. for Power Industry (ENE) we should have:

    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_1_ENE.0.1x0.1.nc

    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_2_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_3_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_4_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_5_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_6_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_7_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_8_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_9_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_10_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_11_ENE.0.1x0.1.nc
    
    /home/angel/tropomi/CH4/ENE/v6.0_CH4_2018_12_ENE.0.1x0.1.nc

with something similar for the other sectors. Run the scripts ``EDGARtoAE.py`` to write the data in the proper WRF data file format, and ``edgarv6_ch4.py`` to make a quick visualization of the emissions data. Now you should be ready to run the ``anthro_emis`` by typing ``./anthro_emis < anthro_ghg.inp``

3. Create the fire emission files

Fire emission data are taken from FINN (https://rda.ucar.edu/datasets/ds312.9/). The emission files are created using the ``fire_emis`` utility. Set the file ``finn_ghg.inp`` accordingly and then run the fire_emis by typing ``./anthro_emis < anthro_ghg.inp`` 

4. Interpolate background CH<sub>4</sub> global concentrations to the WRF-GHG initial and boundary conditions 

Background methane concentrations for Europe are taken from CAM-chem (https://www.acom.ucar.edu/cam-chem/cam-chem.shtml). Set your mozbc namelist file according to the example in ``mozbc_ghg.inp`` and then run mozbc by typing ``./mozbc < mozbc_ghg.inp``.

5. Run the WRF-GHG model 

If the emission files for anthro (from sectors other than biomass burning) and biomass burning sources are all ready to use, then the script bash ``run_wrf.sh`` can be used for automation, for example by typing ``sbatch run_wrf.sh`` on Lumi. 

6. Postprocessing routines

- Set the first and last Sentinel-5 Precursor (S5P) orbits in ``orbit_filter.csh`` and then run it for downloading TROPOMI XCH<sub>4</sub> fields from https://ftp.sron.nl/open-access-data-2/TROPOMI/tropomi/ch4/18_17/
- Set the paths and dates in the file ``1.extractvars.csh`` to extract the model parameters of interest
- Set the paths and dates in the file ``4.pick.column.v1.csh`` to calculate simulated XCH<sub>4</sub> concentrations without smoothing. This step is optional and needs NCL to be run.
- Set and run the file ``regridding.py`` to... 
- Run the file ``xch4_maps.py`` for a quick view of the data...
