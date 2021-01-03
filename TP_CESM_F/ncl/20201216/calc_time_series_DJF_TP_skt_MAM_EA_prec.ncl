
begin


lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40


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
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "skt_mon_1x1_1979_2019.nc"
f1 = addfile(diri+fili, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
time = f1->time(iYYYY)
ts_obs = short2flt(f1->skt(iYYYY,::-1,:))
printVarSummary(ts_obs)
nlon = dimsizes(ts_obs&longitude)
nlat = dimsizes(ts_obs&latitude)

ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
TP_obs = mask(ts_obs, topo.ge.1500, True)
copy_VarCoords(ts_obs, TP_obs)
; TP_obs_sea = ts_obs(4::12,{latstr:latend},{lonstr:lonend})
TS_obs_sea = month_to_season(TP_obs, "DJF")
TS_obs_sea = dim_standardize_n_Wrap(TS_obs_sea, 1, 0)
TP_obs_sea = wgt_areaave_Wrap(TS_obs_sea(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
;===================================================================

diri2 = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
filip = "tp_mon_1x1_1979_2019.nc"

fp = addfile(diri2+filip, "r")
Time = fp->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
prec = short2flt(fp->tp(iYYYY,::-1,:))
prec = prec*1000
prec@units = "mm/day"

prec = dtrend_msg_n(time, prec, False, False, 0)
printVarSummary(prec)

prec_sea = month_to_season(prec, "MAM")
prec_sea = dim_standardize_n_Wrap(prec_sea, 1, 0)
prec_sea_EA = wgt_areaave_Wrap(prec_sea(:,{21:31},{108:120}), 1, 1, 0)
printVarSummary(prec_sea_EA)
;===================================================================

lines = new((/2,nyear/), float)
lines(0,:) = TP_obs_sea
lines(1,:) = prec_sea_EA
corr = escorc(lines(0,:), lines(1,:))
;======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "calc_time_series_DJF_TP_skt_MAM_EA_prec"
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
res@trYMinF = -2
res@trYMaxF = 2

res@gsnYRefLine = (/0,1,-1/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

res@gsnRightString = "Corr: "+sprintf("%4.2f", corr)
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, year, lines, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"DJF TP skt","MAM EA prec"/)

;============================================================
    gx = (/1980,1982/)
    text_pos = 1983
    gy1 = (/1.8,1.8/)
    interval = 0.3
    gy2 = gy1 - interval
    gy3 = gy2 - interval

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = 0
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = 2
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)

    draw(plot(0))
end