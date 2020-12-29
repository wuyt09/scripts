
begin

lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

season = (/"DJF","MAM"/)

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

TP_obs_std = dim_standardize_n_Wrap(TP_obs, 1, 1)
printVarSummary(TP_obs_std)

write_table("time_series_TP_skt_"+season(0)+"_"+season(1)+"_ERA5.txt", "w", [/"YEAR","DJF","MAM"/], "%4s %4s %4s")
write_table("time_series_TP_skt_"+season(0)+"_"+season(1)+"_ERA5.txt", "a", [/year,TP_obs_std(0,:),TP_obs_std(1,:)/], "%0.4i %4.2f %4.2f %4.2f")

corr = escorc(TP_obs_std(0,:), TP_obs_std(1,:))
;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "calc_time_series_TP_skt_"+season(0)+"_"+season(1)+"_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

Colors = (/"black","red","blue","black","red","blue"/)
Patterns = (/0,0,0,2,2,2/)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@trXMinF = yrStrt
res@trXMaxF = yrLast

res@gsnYRefLine = (/0,1,-1/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

res@gsnRightString = "Corr: "+sprintf("%4.2f", corr)
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, year, TP_obs_std, res)
draw(plot)

end