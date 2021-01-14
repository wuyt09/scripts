
begin

lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

season = (/"DJF","MAM","JJA","SON"/)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "skt_mon_1x1_1979_2019.nc"
f1 = addfile(diri+fili, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
ts_obs = short2flt(f1->skt(iYYYY,::-1,:))
time = f1->time(iYYYY)
ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
printVarSummary(ts_obs)

TS_obs = mask(ts_obs, topo.ge.1500, True)
copy_VarCoords(ts_obs, TS_obs)
TS_obs_sea = month_to_seasonN(TS_obs, season)
TP_obs = wgt_areaave_Wrap(TS_obs_sea(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs)

std_value = dim_stddev_n_Wrap(TP_obs, 1)
printVarSummary(std_value)
print(season+"  "+sprintf("%4.2f", std_value))

; TP_obs_std = dim_standardize_n_Wrap(TP_obs, 1, 1)
; printVarSummary(TP_obs_std)

; write_table("time_series_TP_skt_"+season(0)+"_"+season(1)+"_ERA5.txt", "w", [/"YEAR","DJF","MAM"/], "%4s %4s %4s")
; write_table("time_series_TP_skt_"+season(0)+"_"+season(1)+"_ERA5.txt", "a", [/year,TP_obs_std(0,:),TP_obs_std(1,:)/], "%0.4i %4.2f %4.2f %4.2f")

; corr = escorc(TP_obs_std(0,:), TP_obs_std(1,:))


end