load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
setfileoption("nc","Format","LargeFile") 

if (.not. isvar("fileid")) then      ; is fileid on command line?
    print("parameter fileid missing")
    exit
end if

indir    = "/scratch/project_465000083/cams"
outdir   = "/users/anliduvi/public_box/ECMWF2WRFChem"
f_ch4    = indir+"/levtype_ml.nc"
f_sp     = indir+"/levtype_sfc.nc"
grg_ch4  = addfile(f_ch4,"r")
grg_sp   = addfile(f_sp,"r")

; output file
system ("rm "+outdir+"/test"+fileid+".nc")      ; remove pre-existing file
out  = addfile(outdir+"/test"+fileid+".nc","c")

date              = toint(cd_calendar(grg_ch4->time, -2))
date@long_name    = "current date as 6 digit integer (YYYYMMDD)"
date!0            = "time"
date&time         = grg_ch4->time&time

datesec1          = cd_calendar(grg_ch4->time, 0)
datesec           = toint(datesec1(:,3) * 3600 + datesec1(:,4) * 60 + datesec1(:,5))
datesec@units     = "s"
datesec@long_name = "seconds to complete current date"
datesec!0         = "time"
datesec&time      = grg_ch4->time&time

time              = cd_convert(grg_ch4->time, "days since 0000-01-01 00:00:00")
time@long_name    = "simulation time"
time@calendar     = "gregorian"
time!0            = "time"
time&time         = grg_ch4->time&time

lev1              = grg_ch4->level
lev               = ispan(0,dimsizes(lev1)-1,1)
lat               = grg_ch4->latitude
lon               = grg_ch4->longitude
setfileoption(out,"DefineMode",True)
fAtt              = True
fAtt@TITLE        = "CAMS global reanalysis (EAC4): https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-reanalysis-eac4?tab=form"
fileattdef(out,fAtt)
dimNames          = (/"time","lev","lat","lon"/)
dimSizes          = (/dimsizes(time),dimsizes(lev),dimsizes(lat),dimsizes(lon)/)
dimUnlim          = (/True,False,False,False/)
filedimdef(out,dimNames,dimSizes,dimUnlim)

out->time         = time
out->date         = date
out->datesec      = datesec

ab                = readAsciiTable("/scratch/project_465000083/cams/levels.txt",9,"float",1)
hyam              = ab(1:dimsizes(lev1),1)
hybm              = ab(1:dimsizes(lev1),2)
hyam!0            = "lev"
hybm!0            = "lev"
hyam              = hyam / 100000.0

out->hyam         = hyam
out->hybm         = hybm

; surface pressure
pres              = short2flt(grg_sp->sp)
;pres              = lonFlip(pres)
pres!0            = "time"
pres!1            = "lat"
pres!2            = "lon"
pres@units        = "PA"
out->PS           = pres

p0                = todouble(100000)
p0@long_name      = "reference pressure"
p0@units          = "Pa"
out->P0           = p0

; methane
ch4               = short2flt(grg_ch4->ch4_c)
ch4               = ch4*(28.97/16.04)
;ch4               = lonFlip(ch4)
ch4!0             = "time"
ch4!1             = "lev"
ch4!2             = "lat"
ch4!3             = "lon"
ch4@units         = "mol/mol"
out->CH4          = ch4
