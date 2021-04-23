
begin

yrStrt = 1979
yrLast = 2019

diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili0 = "surface_thermal_downwards_mon_1x1_1979_2019.nc"
fili1 = "surface_thermal_downwards_clear_sky_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili0, "r")
f2 = addfile(diri+fili1, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

strd = f1->strd(iYYYY,:,:)
strdc = f2->strdc(iYYYY,:,:)

lat = strd&latitude
lon = strd&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)
;===================================================================================================
nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210423_clearsky/"
fili2 = diri2+"baseline_radsfc_ranc_31.grd"
fili3 = diri2+"wv_radsfc_ranc_31.grd"
fili4 = diri2+"t_radsfc_ranc_31.grd"
fili5 = diri2+"albedo_radsfc_ranc_31.grd"
fili6 = diri2+"co2_radsfc_ranc_31.grd"
fili7 = diri2+"solar_radsfc_ranc_31.grd"
fili8 = diri2+"o3_radsfc_ranc_31.grd"
fili9 = diri2+"warm_radsfc_ranc_31.grd"

UNDEF = -999.
x = new((/nvars,nlat,nlon/),float,UNDEF)
x!0 = "var"
x!1 = "lat"
x!2 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x@long_name = "forcing"
x@units     = "W/m2"

wv = x
t = x
al = x
co2 = x
solar = x
o3 = x
y = x

do nv = 0,nvars-1
    x(nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
    wv(nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
    t(nv,:,:) = fbindirread(fili4,nv,(/nlat,nlon/),"float")
    al(nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
    co2(nv,:,:) = fbindirread(fili6,nv,(/nlat,nlon/),"float")
    solar(nv,:,:) = fbindirread(fili7,nv,(/nlat,nlon/),"float")
    o3(nv,:,:) = fbindirread(fili8,nv,(/nlat,nlon/),"float")
    y(nv,:,:) = fbindirread(fili9,nv,(/nlat,nlon/),"float")
end do
printVarSummary(x)
printMinMax(x(0,:,:),1)


nx = 2
lines = 

