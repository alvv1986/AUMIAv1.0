#!/bin/bash -l
#SBATCH --job-name=test_wrf
#SBATCH --account=project_465000083
#SBATCH --time=02:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=8
#SBATCH --partition=standard
ulimit -s unlimited
############################# env vars and paths ##############################
export HDF5_DISABLE_VERSION_CHECK=1
DIR_WRF=$EBROOTWRF/test/em_real
DIR_MET=$EBROOTWPS/era5/0304_2304_2018
DIR_MOZ=/users/anliduvi/mozbc_gfortran
DIR_OUT=/scratch/project_465000083/wrfout/
###############################################################################
zero=0
# Defining starting year and month for loop 
year=2018
month=4
cd $DIR_WRF
rm rsl.* wrfrst* wrfinput_d01 wrfbdy_d01
# Defining the end of the month loop:
while [ $month -le 4 ]
do
########### Set the start and end days for each month to loop over ############
if [ $month -eq 1 ];then cmonth=JAN; idx=01; fdx=31; fi
if [ $month -eq 2 ];then cmonth=FEB; idx=01; fdx=28; fi  # Check if leap-year
if [ $month -eq 3 ];then cmonth=MAR; idx=01; fdx=31; fi
if [ $month -eq 4 ];then cmonth=APR; idx=10; fdx=19; fi
if [ $month -eq 5 ];then cmonth=MAY; idx=01; fdx=31; fi
if [ $month -eq 6 ];then cmonth=JUN; idx=01; fdx=30; fi
if [ $month -eq 7 ];then cmonth=JUL; idx=01; fdx=31; fi
if [ $month -eq 8 ];then cmonth=AUG; idx=01; fdx=31; fi
if [ $month -eq 9 ];then cmonth=SEP; idx=01; fdx=30; fi
if [ $month -eq 10 ];then cmonth=OCT; idx=01; fdx=31; fi
if [ $month -eq 11 ];then cmonth=NOV; idx=01; fdx=30; fi
if [ $month -eq 12 ];then cmonth=DEC; idx=01; fdx=31; fi
###################### Run a 3-day WRF-GHG simulation #########################
day=`expr $idx \* 1`
while [ $day -le $fdx ]
do
cmonth=$month
if [ $month -lt 10 ]; then cmonth=$zero$month; fi
cday=$day
if [ $day -lt 10 ]; then cday=$zero$day; fi
echo $year $cmonth $cday
cday1=`date -u +%d -d "${year}-${cmonth}-${cday} 4 day" `
cmonth1=`date -u +%m -d "${year}-${cmonth}-${cday} 4 day" `
year1=`date -u +%Y -d "${year}-${cmonth}-${cday} 4 day" `
cday2=`date -u +%d -d "${year}-${cmonth}-${cday} -1 day" `
cmonth2=`date -u +%m -d "${year}-${cmonth}-${cday} -1 day" `
year2=`date -u +%Y -d "${year}-${cmonth}-${cday} -1 day" `
mkdir $DIR_OUT/${year}${cmonth}${cday}
echo ${year}${cmonth}${cday} ${year1}${cmonth1}${cday1}
cp namelist.input.opera.prepchem namelist.input
cat << Eof > namelist.input
&time_control
run_days                            = 0,
run_hours                           = 0,
run_minutes                         = 0,
run_seconds                         = 0,
start_year                          = $year, 2018,
start_month                         = $cmonth,   07,
start_day                           = $cday,   31,
start_hour                          = 00,   00,
end_year                            = $year1, 2018,
end_month                           = $cmonth1,   08,
end_day                             = $cday1,   07,
end_hour                            = 00,   00,
interval_seconds                    = 21600
input_from_file                     = .true.,.true.,
history_interval                    = 60,  60,
history_outname                     = "${DIR_OUT}/${year}${cmonth}${cday}/wrfout_d<domain>_<date>",
frames_per_outfile                  = 1, 1,
restart                             = .false.,
restart_interval                    = 7200,
io_form_history                     = 2
io_form_restart                     = 2
io_form_input                       = 2
io_form_boundary                    = 2
io_form_auxinput5                   = 0
io_form_auxinput12                  = 2
auxinput5_inname                    = 'wrfchemi_d<domain>_<date>'
auxinput12_inname                   = 'wrf_chem_input'
auxinput5_interval_m                = 60, 60, 60
frames_per_auxinput5                = 12
frames_per_auxinput12               = 1
force_use_old_data                  = .True.,
/

&domains
time_step                           = 180,
time_step_fract_num                 = 0,
time_step_fract_den                 = 1,
max_dom                             = 1,
e_we                                = 120,    67,
e_sn                                = 120,    61,
e_vert                              = 40,     35,
eta_levels                          = 1.000, 0.998, 0.996, 0.993, 0.990,
                                      0.986, 0.982, 0.977, 0.972, 0.966,
                                      0.960, 0.953, 0.945, 0.930, 0.905,
                                      0.875, 0.840, 0.800, 0.750, 0.700,
                                      0.650, 0.600, 0.550, 0.500, 0.450,
                                      0.400, 0.350, 0.300, 0.250, 0.200,
                                      0.150, 0.127, 0.106, 0.088, 0.070,
                                      0.055, 0.040, 0.026, 0.013, 0.000
!dzstretch_s                         = 1.1
p_top_requested                     = 5000,
num_metgrid_levels                  = 38,
num_metgrid_soil_levels             = 4,
dx                                  = 30000,
dy                                  = 30000,
grid_id                             = 1,     2,
parent_id                           = 0,     1,
i_parent_start                      = 1,     61,
j_parent_start                      = 1,     65,
parent_grid_ratio                   = 1,     3,
parent_time_step_ratio              = 1,     3,
feedback                            = 1,
smooth_option                       = 0
/

&physics
!physics_suite                       = 'CONUS'
mp_physics                          = 4,    -1,
cu_physics                          = 93,    -1,
ra_lw_physics                       = 1,    -1,
ra_sw_physics                       = 1,    -1,
bl_pbl_physics                      = 1,    -1,
sf_sfclay_physics                   = 1,    -1,
sf_surface_physics                  = 2,    -1,
radt                                = 15,    15,
bldt                                = 0,     0,
cudt                                = 0,     0,
cu_diag                             = 1,
icloud                              = 0,
num_land_cat                        = 21,
sf_urban_physics                    = 0,     0,
fractional_seaice                   = 1,
/

&fdda
/

&dynamics
hybrid_opt                          = 2,
w_damping                           = 0,
diff_opt                            = 2,      2,
km_opt                              = 4,      4,
diff_6th_opt                        = 0,      0,
diff_6th_factor                     = 0.12,   0.12,
base_temp                           = 290.
damp_opt                            = 3,
zdamp                               = 5000.,  5000.,
dampcoef                            = 0.2,    0.2,
khdif                               = 0,      0,
kvdif                               = 0,      0,
non_hydrostatic                     = .true., .true.,
moist_adv_opt                       = 1,      1,
scalar_adv_opt                      = 1,      1,
chem_adv_opt                        = 1,      1,
tracer_adv_opt                      = 1,      1,
gwd_opt                             = 1,      0,
/

&bdy_control
spec_bdy_width                      = 5,
specified                           = .true.
/

&grib2
/

&chem
kemit                               = 1,
chem_opt                            = 17,        2,
bioemdt                             = 15,       30,
photdt                              = 15,       30,
chemdt                              = 0.,       2.,
!frames_per_emissfile                = 36,
io_style_emissions                  = 1,
emiss_inpt_opt                      = 16,        1,
emiss_opt                           = 17,        3,
chem_in_opt                         = 1,        0,
phot_opt                            = 1,        1,
gas_drydep_opt                      = 0,        0,
aer_drydep_opt                      = 0,        0,
bio_emiss_opt                       = 17,        0,
dust_opt                            = 0,
dmsemis_opt                         = 0,
seas_opt                            = 0,
gas_bc_opt                          = 1,        1,
gas_ic_opt                          = 1,        1,
aer_bc_opt                          = 1,        1,
aer_ic_opt                          = 1,        1,
gaschem_onoff                       = 1,        1,
aerchem_onoff                       = 0,        1,
wetscav_onoff                       = 0,        0,
cldchem_onoff                       = 0,        0,
vertmix_onoff                       = 1,        1,
chem_conv_tr                        = 1,        1,
biomass_burn_opt                    = 0,        0,
plumerisefire_frq                   = 30,       30,
aer_ra_feedback                     = 0,        0,
have_bcs_chem                       = .false., .false.,
vprm_opt                            = "VPRM_table_EUROPE",
wpeat                               = 0.05,
wflood                              = 0.19,
term_opt                            = "CH4_termite_OW",
/

&namelist_quilt
nio_tasks_per_group = 0,
nio_groups = 1,
/
Eof
ln -s $DIR_OUT/${year2}${cmonth2}${cday2}/wrfout_d01_${year}-${cmonth}-${cday}_00:00:00 wrf_chem_input_d01
echo 'Running real.exe...'
$DIR_WRF/real.exe
rm -f rsl.*
cd $DIR_MOZ
rm met_em.d01* wrfinput_d01 wrfbdy_d01
ln -sf $DIR_MET/met_em.d01* .
ln -s $DIR_WRF/wrfinput_d01 .
ln -s $DIR_WRF/wrfbdy_d01 .
./mozbc < GHG.inp
cd $DIR_WRF
cp wrfinput_d01 wrfinput_d01_${year}${cmonth}${cday}
cp wrfbdy_d01 wrfbdy_d01_${year}${cmonth}${cday}
cat << Eof > namelist.input
&time_control
run_days                            = 0,
run_hours                           = 0,
run_minutes                         = 0,
run_seconds                         = 0,
start_year                          = $year, 2018,
start_month                         = $cmonth,   07,
start_day                           = $cday,   31,
start_hour                          = 00,   00,
end_year                            = $year1, 2018,
end_month                           = $cmonth1,   08,
end_day                             = $cday1,   07,
end_hour                            = 00,   00,
interval_seconds                    = 21600
input_from_file                     = .true.,.true.,
history_interval                    = 60,  60,
history_outname                     = "${DIR_OUT}/${year}${cmonth}${cday}/wrfout_d<domain>_<date>",
frames_per_outfile                  = 1, 1,
restart                             = .false.,
restart_interval                    = 7200,
io_form_history                     = 2
io_form_restart                     = 2
io_form_input                       = 2
io_form_boundary                    = 2
io_form_auxinput5                   = 2
io_form_auxinput12                  = 2
auxinput5_inname                    = 'wrfchemi_d<domain>_<date>'
auxinput12_inname                   = 'wrf_chem_input'
auxinput5_interval_m                = 60, 60, 60
frames_per_auxinput5                = 12
frames_per_auxinput12               = 1
force_use_old_data                  = .True.,
/

&domains
time_step                           = 180,
time_step_fract_num                 = 0,
time_step_fract_den                 = 1,
max_dom                             = 1,
e_we                                = 120,    67,
e_sn                                = 120,    61,
e_vert                              = 40,     35,
eta_levels                          = 1.000, 0.998, 0.996, 0.993, 0.990,
                                      0.986, 0.982, 0.977, 0.972, 0.966,
                                      0.960, 0.953, 0.945, 0.930, 0.905,
                                      0.875, 0.840, 0.800, 0.750, 0.700,
                                      0.650, 0.600, 0.550, 0.500, 0.450,
                                      0.400, 0.350, 0.300, 0.250, 0.200,
                                      0.150, 0.127, 0.106, 0.088, 0.070,
                                      0.055, 0.040, 0.026, 0.013, 0.000
!dzstretch_s                         = 1.1
p_top_requested                     = 5000,
num_metgrid_levels                  = 38,
num_metgrid_soil_levels             = 4,
dx                                  = 30000,
dy                                  = 30000,
grid_id                             = 1,     2,
parent_id                           = 0,     1,
i_parent_start                      = 1,     61,
j_parent_start                      = 1,     65,
parent_grid_ratio                   = 1,     3,
parent_time_step_ratio              = 1,     3,
feedback                            = 1,
smooth_option                       = 0
/

&physics
!physics_suite                       = 'CONUS'
mp_physics                          = 4,    -1,
cu_physics                          = 93,    -1,
ra_lw_physics                       = 1,    -1,
ra_sw_physics                       = 1,    -1,
bl_pbl_physics                      = 1,    -1,
sf_sfclay_physics                   = 1,    -1,
sf_surface_physics                  = 2,    -1,
radt                                = 15,    15,
bldt                                = 0,     0,
cudt                                = 0,     0,
cu_diag                             = 1,
icloud                              = 0,
num_land_cat                        = 21,
sf_urban_physics                    = 0,     0,
fractional_seaice                   = 1,
/

&fdda
/

&dynamics
hybrid_opt                          = 2,
w_damping                           = 0,
diff_opt                            = 2,      2,
km_opt                              = 4,      4,
diff_6th_opt                        = 0,      0,
diff_6th_factor                     = 0.12,   0.12,
base_temp                           = 290.
damp_opt                            = 3,
zdamp                               = 5000.,  5000.,
dampcoef                            = 0.2,    0.2,
khdif                               = 0,      0,
kvdif                               = 0,      0,
non_hydrostatic                     = .true., .true.,
moist_adv_opt                       = 1,      1,
scalar_adv_opt                      = 1,      1,
chem_adv_opt                        = 1,      1,
tracer_adv_opt                      = 1,      1,
gwd_opt                             = 1,      0,
/

&bdy_control
spec_bdy_width                      = 5,
specified                           = .true.
/

&grib2
/

&chem
kemit                               = 1,
chem_opt                            = 17,        2,
bioemdt                             = 15,       30,
photdt                              = 15,       30,
chemdt                              = 0.,       2.,
!frames_per_emissfile                = 36,
io_style_emissions                  = 1,
emiss_inpt_opt                      = 16,        1,
emiss_opt                           = 17,        3,
chem_in_opt                         = 1,        0,
phot_opt                            = 1,        1,
gas_drydep_opt                      = 0,        0,
aer_drydep_opt                      = 0,        0,
bio_emiss_opt                       = 17,        0,
dust_opt                            = 0,
dmsemis_opt                         = 0,
seas_opt                            = 0,
gas_bc_opt                          = 1,        1,
gas_ic_opt                          = 1,        1,
aer_bc_opt                          = 1,        1,
aer_ic_opt                          = 1,        1,
gaschem_onoff                       = 1,        1,
aerchem_onoff                       = 0,        1,
wetscav_onoff                       = 0,        0,
cldchem_onoff                       = 0,        0,
vertmix_onoff                       = 1,        1,
chem_conv_tr                        = 1,        1,
biomass_burn_opt                    = 0,        0,
plumerisefire_frq                   = 30,       30,
aer_ra_feedback                     = 0,        0,
have_bcs_chem                       = .true., .false.,
vprm_opt                            = "VPRM_table_EUROPE",
wpeat                               = 0.05,
wflood                              = 0.19,
term_opt                            = "CH4_termite_OW",
/

&namelist_quilt
nio_tasks_per_group = 0,
nio_groups = 1,
/
Eof
echo 'Running wrf.exe...'
srun $DIR_WRF/wrf.exe
sync
day=`expr $day + 1`
done
month=`expr $month + 1`
done
echo "Successful completion of WRF-GHG"
exit