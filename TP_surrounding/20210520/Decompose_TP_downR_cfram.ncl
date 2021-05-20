
begin

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)

;=============================================================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili1 = "surface_thermal_downwards_mon_1x1_1979_2019.nc"
fili2 = "surface_thermal_downwards_clear_sky_mon_1x1_1979_2019.nc"
; fili3 = "surface_net_thermal_mon_1x1_1979_2019.nc"
; fili4 = "surface_net_thermal_clear_sky_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili1, "r")
Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

f2 = addfile(diri+fili2, "r")
; f3 = addfile(diri+fili3, "r")
; f4 = addfile(diri+fili4, "r")

STRD = short2flt(f1->strd(iYYYY,:,:))
STRD = STRD/86400.
STRDC = short2flt(f2->strdc(iYYYY,:,:))
STRDC = STRDC/86400.
; STR  = short2flt(f3->str(iYYYY,:,:))
; STR  = STR/86400.
; STRC  = short2flt(f4->strc(iYYYY,:,:))
; STRC  = STRC/86400.

STRD_CRE = STRD - STRDC
copy_VarCoords(STRD, STRD_CRE)
printVarSummary(STRD_CRE)

strd_cre = month_to_season(STRD_CRE, season)
strd_cre&time = year
printVarSummary(strd_cre)

cloud = runave_n_Wrap(strd_cre, 10, 0, 0)
printVarSummary(cloud)
cloud_region = wgt_areaave_Wrap(cloud(:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
printVarSummary(cloud_region)
; print(cloud_region)

;================================================================================
nn = 31
nlat = dimsizes(cloud&latitude)
nlon = dimsizes(cloud&longitude)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210513_clearsky/"
nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)


UNDEF = -999.
x = new((/nn,6,nvars,nlat,nlon/),float,UNDEF)
x!0 = "time"
x!1 = "process"
x!2 = "var"
x!3 = "lat"
x!4 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x@long_name = "forcing"
x@units     = "W/m2"


do nnn = 0, nn-1
    fili1 = dirir2+"base_radsfc_ranc_"+nnn+".grd"
    fili2 = dirir2+"co2_radsfc_ranc_"+nnn+".grd"
    fili3 = dirir2+"o3_radsfc_ranc_"+nnn+".grd"
    fili4 = dirir2+"wv_radsfc_ranc_"+nnn+".grd"
    fili5 = dirir2+"ta_radsfc_ranc_"+nnn+".grd"
    fili6 = dirir2+"warm_radsfc_ranc_"+nnn+".grd"

    do nv = 0, nvars-1
        x(nnn,0,nv,:,:) = fbindirread(fili1,nv,(/nlat,nlon/),"float")
        x(nnn,1,nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
        x(nnn,2,nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
        x(nnn,3,nv,:,:) = fbindirread(fili4,nv,(/nlat,nlon/),"float")
        x(nnn,4,nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
        x(nnn,5,nv,:,:) = fbindirread(fili6,nv,(/nlat,nlon/),"float")
    end do










end do




end