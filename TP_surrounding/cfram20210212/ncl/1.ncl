
begin

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEVs = 1000
LEVe = 300

latstr_east = 35
latend_east = 55
lonstr_east = 80
lonend_east = 120

latstr_west = 40
latend_west = 60
lonstr_west = 10
lonend_west = 50


;=============================================================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili = "t_mon_1x1_1979_2019.nc"

f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t = short2flt(f->t(iYYYY,:,:,:))
printVarSummary(t)


filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,:,:))
; printVarSummary(ps)

plev = tofloat(t&level)
plev = plev*100
plev!0 = "lev"
plev&lev = t&level
plev@units = "Pa"
; print(plev)

dp = dpres_plevel_Wrap(plev, ps, min(plev), 0)
; printVarSummary(dp)

tdp = t*dp
copy_VarCoords(t, tdp)
; printVarSummary(tdp)

t_tropo = dim_sum_n_Wrap(tdp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
copy_VarCoords(t(:,0,:,:), t_tropo)

t_tropo_JJA = month_to_season(t_tropo, season)
t_series = wgt_areaave_Wrap(t_tropo_JJA(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
; printVarSummary(t_series)
print(t_series)

t_base = avg(t_series(0:9))
t_new = t_series - t_base
print(t_new)
t_runmean = runave(t_new, 10, 0)
t_runmean!0 = "year"
t_runmean&year = year
print(t_runmean)
;=============================================================================================================

nvname = (/"albedo","wv","cloud","cloud_sw","cloud_lw","co2","o3","solar",\
           "dyn","atm_dyn","ocean_dyn","shflx","lhflx"/)  ; ,"error"
nn = 31

lat  = f->latitude ;(::-1)
lon  = f->longitude
level  = f->level
lev = new(38, integer, "No_FillValue")
lev(0:nlev-2) = level
lev(nlev-1) = 1001
print(lev)

delete(f)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)

setfileoption("bin","ReadByteOrder","Native")

UNDEF = -999.
x = new((/nn,nvars,nlev,nlat,nlon/),float,UNDEF)
x!0 = "series"
x!1 = "var"
x!2 = "lev"
x!3 = "lat"
x!4 = "lon"
x&series = ispan(1983, 2013, 1)
x&var = nvname
x&lat = lat
x&lon = lon
x&lev = lev
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x&lev@units = "hPa"
x@long_name = "forcing"
x@units     = "W/m2"


end