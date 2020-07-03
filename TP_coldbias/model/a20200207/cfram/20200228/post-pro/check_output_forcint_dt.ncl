begin

;================read metadata=============================================
diri1 = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_Ctrl_FAMIP/a20191206/"
fili1 = "Ctrl_FAMIP.cam.h0.T.1979-2005.nc"
f1 = addfile(diri1+fili1,"r")
lon = f1->lon
lat = f1->lat
lev = f1->level

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)
delete(f1)
;===========================================================================
nvname = (/"solarin","swdn_sfc","swup_sfc","tsurf","huss","ps","o3","ta",\
           "hus","cld_amt","liq_wat","ice_wat"/)
nvars = dimsizes(nvname)

x = new((/6,nlev,nlat,nlon/),float,-999.)
x!0 = "var"
x!1 = "lev"
x!2 = "lat"
x!3 = "lon"
x&var = nvname(6:11)
x&lat = lat
x&lon = lon
x&lev = lev
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x&lev@units = "hPa"

y = new((/6,nlat,nlon/),float,-999.)
y!0 = "var"
y!1 = "lat"
y!2 = "lon"
y&var = nvname(0:5)
y&lat = lat
y&lon = lon
y&lat@units = "degrees_north"
y&lon@units = "degrees_east"

diri = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200207/cfram/20200228/"
fili = "baseline_input_1.dat"

do i = 0, 5
    y(i,:,:) = fbindirread(diri+fili, i, (/nlat,nlon/), "float")
end do
do nv = 0,5
   do nl = 0,nlev-1
      x(nv,nl,:,:) = fbindirread(diri+fili,nv*nlev+nl+6,(/nlat,nlon/),"float")
   end do
end do
printVarSummary(x)
printVarSummary(y)

ytest = 50
xtest = 120
print(x(:,14,{ytest},{xtest}))
print(y(:,{ytest},{xtest}))
