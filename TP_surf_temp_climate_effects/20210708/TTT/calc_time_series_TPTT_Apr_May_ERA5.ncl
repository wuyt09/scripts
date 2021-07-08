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
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili = "t_mon_1x1_1979_2020.nc" ;"temp.monthly.197901-201812.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t = short2flt(f->t(iYYYY,::-1,::-1,:))
time = f->time(iYYYY)
lat = t&latitude
lon = t&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))
printVarSummary(ps)


plev = tofloat(t&level)
plev = plev*100
plev!0 = "lev"
plev&lev = t&level
plev@units = "Pa"
dp = dpres_plevel_Wrap(plev, ps, min(plev), 0)
printVarSummary(dp)
tdp = t*dp
copy_VarCoords(t, tdp)
printVarSummary(tdp)

ts_obs = dim_sum_n_Wrap(tdp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
copy_VarCoords(t(:,0,:,:), ts_obs)

ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
TS_obs = mask(ts_obs, topo.ge.1500, True)
copy_VarCoords(ts_obs, TS_obs)
TS_obs_sea = new((/3,nyear,nlat,nlon/), float)
TS_obs_sea(0,:,:,:) = TS_obs(2::12,:,:)
TS_obs_sea(1,:,:,:) = TS_obs(3::12,:,:)
TS_obs_sea(2,:,:,:) = TS_obs(4::12,:,:)
copy_VarCoords(ts_obs(0,:,:), TS_obs_sea(0,0,:,:))
TS_obs_sea!1 = "year"
TS_obs_sea&year = year
TP_obs = wgt_areaave_Wrap(TS_obs_sea(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs)

TP_obs_std = dim_standardize_n_Wrap(TP_obs, 1, 1)
printVarSummary(TP_obs_std)

write_table("time_series_TP_TT_Mar_Apr_May_ERA5.txt", "w", [/"YEAR","Mar","Apr","May"/], "%4s %4s %4s %4s")
write_table("time_series_TP_TT_Mar_Apr_May_ERA5.txt", "a", [/year,TP_obs_std(0,:),TP_obs_std(1,:),TP_obs_std(2,:)/], "%0.4i %4.2f %4.2f %4.2f %4.2f")


end