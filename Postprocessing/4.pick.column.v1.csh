#!/bin/csh
setenv HDF5_DISABLE_VERSION_CHECK 1
set echo

foreach year (2018)
 foreach mon (04)
  foreach day (01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)
   foreach hour (00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23)

        setenv FIN /home/angel/tropomi/mlevels/mlevels_d01_${year}-${mon}-${day}-${hour}Z.nc
        if (-e $FIN) then
        setenv FON /home/angel/tropomi/mcolumn/column_d01_${year}-${mon}-${day}-${hour}Z.nc  #

ncl << EOF

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

	fn      = getenv("FIN")
        fon     = getenv("FON")
        fi      = addfile(fn,"r")
        system("rm -f "+fon)
        fo      = addfile(fon,"c")
        filedimdef(fo,"Time",-1,True)

;================================= begin calculation ==================================

        CH4_ANT           = fi->CH4_ANT
        CH4_BIO           = fi->CH4_BIO
        CH4_BCK           = fi->CH4_BCK
        CH4_BBU           = fi->CH4_BBU
        P                 = fi->P
        PB                = fi->PB
        Pabs              = P+PB
        Pabs!0            = "Time"
        Pabs!1            = "bottom_top"
        Pabs!2            = "south_north"
        Pabs!3            = "west_east"
        PSFC              = fi->PSFC
        PTOP              = fi->P_TOP
        Times             = fi->Times
        XLAT              = fi->XLAT
        XLONG             = fi->XLONG
 
        CH4_ANT           = CH4_ANT*1000.
        CH4_BIO           = CH4_BIO*1000.
        CH4_BCK           = CH4_BCK*1000.
        CH4_BBU           = CH4_BBU*1000.
        CH4_ALL           = CH4_ANT+CH4_BIO+CH4_BCK+CH4_BBU
        dsizes            = dimsizes(CH4_ALL)
        CH4_LY_ANT        = new((/dsizes(0),dsizes(1),dsizes(3),dsizes(2)/),float)
        CH4_LY_BIO        = new((/dsizes(0),dsizes(1),dsizes(3),dsizes(2)/),float)
        CH4_LY_BCK        = new((/dsizes(0),dsizes(1),dsizes(3),dsizes(2)/),float)
        CH4_LY_BBU        = new((/dsizes(0),dsizes(1),dsizes(3),dsizes(2)/),float)
        CH4_LY_ALL        = new((/dsizes(0),dsizes(1),dsizes(3),dsizes(2)/),float)
        CH4_LY_ANT!0      = "Time"
        CH4_LY_ANT!1      = "bottom_top"
        CH4_LY_ANT!2      = "south_north"
        CH4_LY_ANT!3      = "west_east"
        CH4_LY_BIO!0      = "Time"
        CH4_LY_BIO!1      = "bottom_top"
        CH4_LY_BIO!2      = "south_north"
        CH4_LY_BIO!3      = "west_east"
        CH4_LY_BCK!0      = "Time"
        CH4_LY_BCK!1      = "bottom_top"
        CH4_LY_BCK!2      = "south_north"
        CH4_LY_BCK!3      = "west_east"
        CH4_LY_BBU!0      = "Time"
        CH4_LY_BBU!1      = "bottom_top"
        CH4_LY_BBU!2      = "south_north"
        CH4_LY_BBU!3      = "west_east"
        CH4_LY_ALL!0      = "Time"
        CH4_LY_ALL!1      = "bottom_top"
        CH4_LY_ALL!2      = "south_north"
        CH4_LY_ALL!3      = "west_east"
        Pb                = new((/dsizes(1),dsizes(2),dsizes(3)/),float)
        Pa                = new((/dsizes(1),dsizes(2),dsizes(3)/),float)
        do j = 0, dsizes(2)-1
         do i = 0, dsizes(3)-1
          Pb(0,j,i) = PSFC(0,j,i)
          Pa(0,j,i) = 2*Pabs(0,0,j,i)-Pb(0,j,i)
         end do
        end do
        do k = 1, dsizes(1)-1
         do j = 0, dsizes(2)-1
          do i = 0, dsizes(3)-1
           Pb(k,j,i) = (/Pa(k-1,j,i)/)
           Pa(k,j,i) = 2*Pabs(0,k,j,i)-Pb(k,j,i)           
          end do
         end do
        end do
        do k = 0, dsizes(1)-1
         do j = 0, dsizes(2)-1
          do i = 0, dsizes(3)-1
           CH4_LY_ANT(0,k,j,i) = CH4_ANT(0,k,j,i)*(Pb(k,j,i)-Pa(k,j,i))/(PSFC(0,j,i)-PTOP)
           CH4_LY_BIO(0,k,j,i) = CH4_BIO(0,k,j,i)*(Pb(k,j,i)-Pa(k,j,i))/(PSFC(0,j,i)-PTOP)
           CH4_LY_BCK(0,k,j,i) = CH4_BCK(0,k,j,i)*(Pb(k,j,i)-Pa(k,j,i))/(PSFC(0,j,i)-PTOP)
           CH4_LY_BBU(0,k,j,i) = CH4_BBU(0,k,j,i)*(Pb(k,j,i)-Pa(k,j,i))/(PSFC(0,j,i)-PTOP)
           CH4_LY_ALL(0,k,j,i) = CH4_ALL(0,k,j,i)*(Pb(k,j,i)-Pa(k,j,i))/(PSFC(0,j,i)-PTOP)
          end do
         end do
        end do 

        fo->CH4_ANT       = dim_sum_Wrap(CH4_LY_ANT(Time|:,south_north|:,west_east|:,bottom_top|:))
        fo->CH4_BIO       = dim_sum_Wrap(CH4_LY_BIO(Time|:,south_north|:,west_east|:,bottom_top|:))
        fo->CH4_BCK       = dim_sum_Wrap(CH4_LY_BCK(Time|:,south_north|:,west_east|:,bottom_top|:))
        fo->CH4_BBU       = dim_sum_Wrap(CH4_LY_BBU(Time|:,south_north|:,west_east|:,bottom_top|:))
        fo->CH4_ALL       = dim_sum_Wrap(CH4_LY_ALL(Time|:,south_north|:,west_east|:,bottom_top|:))
        fo->Times         = Times
        fo->XLAT          = XLAT
        fo->XLONG         = XLONG
        fo@history        = "Angel V., May 2022"

end
EOF

	echo "$year - $mon - $day - $hour Z is ok"
        else  
        echo "No $FIN"
        endif
	end
    end
end

    echo "Successfully Done!"
end
