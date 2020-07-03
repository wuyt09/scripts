
begin

plev = (/10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000/)

nvname = (/"albedo","wv","cloud","cloud_sw","cloud_lw","co2","o3","solar",\
           "dyn","atm_dyn","ocean_dyn","shflx","lhflx"/)  ; ,"error"

f = addfile("./Ctrl_FAMIP.cam.h0.2005-12.nc", "r")

lat = tofloat(f->lat)
lon = tofloat(f->lon)
level = tofloat(f->lev)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = 27
nvars = dimsizes(nvname)

hyam = f->hyam
hybm = f->hybm
PS   = f->PS
p0 = 1000.

ps = new((/nvars,nlat,nlon/), float)
do i = 0, nvars-1
    ps(i,:,:) = PS
end do
printVarSummary(ps)


lev = new(nlev,float,"No_FillValue")
lev(0:nlev-2) = level
lev(nlev-1) = 1001

setfileoption("bin","ReadByteOrder","Native")

UNDEF = -999.
x = new((/nvars,nlev,nlat,nlon/),float,UNDEF)
x!0 = "var"
x!1 = "lev"
x!2 = "lat"
x!3 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lev = lev
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x&lev@units = "hPa"
x@long_name = "Delta T"
x@units     = "K"

y = x
y@long_name = "forcing"
y@units     = "W/m2"

fili2 = "../partial_T_1.grd"
fili3 = "../forcing_1.grd"

do nv = 0,nvars-1
   do nl = 0,nlev-1
      x(nv,nl,:,:) = fbindirread(fili2,nv*nlev+nl,(/nlat,nlon/),"float")
      y(nv,nl,:,:) = fbindirread(fili3,nv*nlev+nl,(/nlat,nlon/),"float")
   end do
end do

if (any(isnan_ieee(x))) then   ; if x got a element is Not A Number, then the Fill_value
  value = -999.
  replace_ieeenan (x,value,0)
  x@_FillValue = value
  replace_ieeenan (y,value,0)
  y@_FillValue = value
end if
printVarSummary(x)
printVarSummary(y)


x26 = x(:,0:25,:,:)
y26 = y(:,0:25,:,:)

; type of interpolation: 1 = linear, 2 = log, 3 = loglog
    interp = 2
; is extrapolation desired if data is outside the range of PS
    extrap = False

xp = vinth2p(x26, hyam, hybm, plev, ps, interp, p0, 1, extrap)
yp = vinth2p(y26, hyam, hybm, plev, ps, interp, p0, 1, extrap)
printVarSummary(yp)


dt27 = 


end