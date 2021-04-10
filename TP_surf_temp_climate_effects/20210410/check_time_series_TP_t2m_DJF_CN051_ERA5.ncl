
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
topo1 = ftopo->topo({14.75:55.25},{69.75:140.25})
printVarSummary(topo1)
topo2 = ftopo->topo(::-1,:)
;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

x = new((/4,nyear/), float)
;===================================================================
; CN051
;===================================================================
diri = "/home/yangsong3/data-observation/CN05.1/"
fili = "CN05.1_Tm_1961_2018_month_025x025.nc"
f1 = addfile(diri+fili, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
time = f1->time(iYYYY)
ts_obs = short2flt(f1->tm(iYYYY,:,:))
printVarSummary(ts_obs)

nlon = dimsizes(ts_obs&lon)
nlat = dimsizes(ts_obs&lat)

ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
TOPO = linint2_Wrap(topo1&lon, topo1&lat, topo1, False, ts_obs&lon, ts_obs&lat, 0)
TS_obs = mask(ts_obs, TOPO.ge.1500, True)
copy_VarCoords(ts_obs, TS_obs)
TS_obs_sea = month_to_seasonN(ts_obs, season)
printVarSummary(TS_obs_sea)

TP_obs = wgt_areaave_Wrap(TS_obs_sea(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs)

TP_obs_std = dim_standardize_n_Wrap(TP_obs, 1, 1)
printVarSummary(TP_obs_std)

corr = escorc(TP_obs_std(0,:), TP_obs_std(1,:))

x(0:1,:) = TP_obs_std
; write_table("time_series_TP_t2m_DJF_MAM_CN051.txt", "w", [/"YEAR","DJF","MAM"/], "%4s %4s %4s")
; write_table("time_series_TP_t2m_DJF_MAM_CN051.txt", "a", [/year,TP_obs_std(0,:),TP_obs_std(1,:)/], "%0.4i %4.2f %4.2f %4.2f")

ts_mask = linint2_Wrap(ts_obs&lon, ts_obs&lat, ts_obs, False, topo2&lon, topo2&lat, 0)
delete([/diri,fili,f1,Time,YYYY,iYYYY,time,ts_obs,nlon,nlat,ts_obs,TOPO,TS_obs,TS_obs_sea,TP_obs,TP_obs_std,corr/])
;===================================================================
; ERA5
;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "t2m_mon_1x1_1979_2019.nc"
f1 = addfile(diri+fili, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
time = f1->time(iYYYY)
ts_obs = short2flt(f1->t2m(iYYYY,::-1,:))
printVarSummary(ts_obs)
nlon = dimsizes(ts_obs&longitude)
nlat = dimsizes(ts_obs&latitude)

ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
ts_obs = mask(ts_obs, .not.ismissing(ts_mask), True)
TS_obs = mask(ts_obs, topo2.ge.1500, True)
copy_VarCoords(ts_obs, TS_obs)
TS_obs_sea = month_to_seasonN(ts_obs, season)
printVarSummary(TS_obs_sea)

TP_obs = wgt_areaave_Wrap(TS_obs_sea(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs)

TP_obs_std = dim_standardize_n_Wrap(TP_obs, 1, 1)
printVarSummary(TP_obs_std)

corr = escorc(TP_obs_std(0,:), TP_obs_std(1,:))

x(2:3,:) = TP_obs_std
; write_table("time_series_TP_t2m_DJF_MAM_ERA5.txt", "w", [/"YEAR","DJF","MAM"/], "%4s %4s %4s")
; write_table("time_series_TP_t2m_DJF_MAM_ERA5.txt", "a", [/year,TP_obs_std(0,:),TP_obs_std(1,:)/], "%0.4i %4.2f %4.2f %4.2f")


;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210410/"
figname = "check_time_series_TP_t2m_DJF_MAM_CN051_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

Colors = (/"black","red","black","red"/)
Patterns = (/0,0,2,2/)

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
res@trYMaxF = 3
res@trYMinF = -3

res@gsnYRefLine = (/0,1,-1/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

; res@gsnRightString = "Corr: "+sprintf("%4.2f", corr)
; res@gsnLeftString = "CN05.1"
res@gsnStringFontHeightF = 0.03
res@tiYAxisString = ""

plot(0) = gsn_csm_xy(wks, year, x, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"CN05.1 DJF","CN05.1 MAM","ERA5 DJF","ERA5 MAM"/)

;============================================================
    gx = (/1980,1982/)
    text_pos = 1982.5
    gx2 = (/1990,1992/)
    text_pos2 = 1992.5

    gy1 = (/2.8,2.8/)
    interval = 0.5
    gy2 = gy1 - interval
    gy3 = gy1 ;- interval
    gy4 = gy2 ;- interval

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx2, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos2, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx2, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos2, gy4, res_text)

draw(plot)

end