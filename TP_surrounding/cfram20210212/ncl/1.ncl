
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


;================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili = "t_mon_1x1_1979_2019.nc"

f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t = short2flt(f->t(iYYYY,:,:,:))
printVarSummary(t)

lat  = f->latitude ;(::-1)
lon  = f->longitude
lev  = f->level
print(lev)

delete(f)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)



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
print(t_runmean)
end