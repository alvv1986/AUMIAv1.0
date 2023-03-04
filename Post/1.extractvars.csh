#!/bin/csh
setenv HDF5_DISABLE_VERSION_CHECK 1
set id = 2018
#set inputdirectory = /scratch/project_465000083/wrfout/0304_2704_2018/wrfout_tropo
set inputdirectory = /scratch/project_465000083/wrfout/test/cam-chem_4dayrun_nofire_ysu
#set outputdirectory = /scratch/project_465000083/regrid/mlevels/2018
set outputdirectory = /scratch/project_465000083/regrid/test/cam-chem_4dayrun_nofire_ysu
foreach mm (04)
foreach day (01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)
   foreach hr (00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23)

      foreach dom (d01)
       set fin = {$inputdirectory}/wrfout_{$dom}_{$id}-{$mm}-{$day}_{$hr}:00:00
       set fon2 = {$outputdirectory}/mlevels_{$dom}_{$id}-{$mm}-{$day}-{$hr}Z.nc

processing:
        if (-e $fin) then
	  echo "processing ... $fin"
          ncks -O -v Times,XLONG,XLAT,T,PSFC,P_TOP,P,PB,CH4_ANT,CH4_BIO,CH4_BCK,CH4_BBU,CH4_TST $fin $fon2
	else if (-e $fin.gz) then
		echo "gunzipping ... $fin"
		gunzip -c -f $fin.gz > $fin
                goto processing
        else
                echo "No $fin"
	endif
	        if (-e $fin.gz) rm -f $fin
      end
   end
end
end
echo "Successfully Done!"
