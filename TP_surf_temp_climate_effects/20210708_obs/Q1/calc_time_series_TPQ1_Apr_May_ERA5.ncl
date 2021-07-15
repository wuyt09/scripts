begin

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEVs = 600
LEVe = 300

latstr = 25
latend = 40
lonstr = 75
lonend = 105

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)


;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/thermal_equation/"
fili = "q1.monmean.1979.2019.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t = short2flt(f->q1(iYYYY,::-1,::-1,:))
time = f->time(iYYYY)
lat = t&lat
lon = t&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))
printVarSummary(ps)

;===================================================================
; diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
; fili = "q1.mon.mean.1979.2016.ERAI.nc"
; f = addfile(diri+fili, "r")

; Time = f->time
; YYYY = cd_calendar(Time, -1)/100
; iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

; t = short2flt(f->q1(iYYYY,::-1,::-1,:))
; time = f->time(iYYYY)
; lat = t&lat
; lon = t&lon
; nlat = dimsizes(lat)
; nlon = dimsizes(lon)

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))
printVarSummary(ps)


plev = tofloat(t&lev)
plev = plev*100
plev!0 = "lev"
plev&lev = t&lev
plev@units = "Pa"
dp = dpres_plevel_Wrap(plev, ps, min(plev), 0)
printVarSummary(dp)
tdp = t*dp
copy_VarCoords(t, tdp)
printVarSummary(tdp)

ts_obs = dim_sum_n_Wrap(tdp, 1)
copy_VarCoords(t(:,0,:,:), ts_obs)


; ;==============================================
; diri = "/home/yangsong3/wangzq/DATA/jra/"
; fili = "Q1I.jra55.mon.mean.nc"
; f = addfile(diri+fili, "r")

; TIME = f->time
; YYYY = cd_calendar(TIME, -1)/100
; iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
; ts_obs = f->Q1(iYYYY,:,:)
; ; printVarSummary(q1)
; lon = ts_obs&lon
; lat = ts_obs&lat
; nlon = dimsizes(lon)
; nlat = dimsizes(lat)
;==============================================
; ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, lon, lat, 0)
TS_obs = mask(ts_obs, TOPO.ge.1500, True)
copy_VarCoords(ts_obs, TS_obs)
TS_obs_sea = new((/3,nyear,nlat,nlon/), float)
TS_obs_sea(0,:,:,:) = TS_obs(2::12,:,:)
TS_obs_sea(1,:,:,:) = TS_obs(3::12,:,:)
TS_obs_sea(2,:,:,:) = TS_obs(4::12,:,:)
copy_VarCoords(ts_obs(0,:,:), TS_obs_sea(0,0,:,:))
TS_obs_sea!1 = "year"
TS_obs_sea&year = year
; TS_obs_sea = dtrend_msg_n(year, TS_obs_sea, False, False, 1)
TP_obs = wgt_areaave_Wrap(TS_obs_sea(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs)

TP_obs_std = dim_standardize_n_Wrap(TP_obs, 1, 1)
printVarSummary(TP_obs_std)

write_table("time_series_TP_Q1_Mar_Apr_May_ERA5.txt", "w", [/"YEAR","Mar","Apr","May"/], "%4s %4s %4s %4s")
write_table("time_series_TP_Q1_Mar_Apr_May_ERA5.txt", "a", [/year,TP_obs_std(0,:),TP_obs_std(1,:),TP_obs_std(2,:)/], "%0.4i %4.2f %4.2f %4.2f %4.2f")


end