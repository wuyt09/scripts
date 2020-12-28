
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
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

ts_obs = new((/nyear*nm,72,128/), float,-99.)

diriobs = "/home/yangsong3/data-observation/stationChina2000/Grid/t2m/"

do i = 0, nyear-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        filiobs = "SURF_CLI_CHN_TEM_MON_GRID_0.5-MEAN-"+yr+mn+".txt"
        x = asciiread(diriobs+filiobs, -1, "float")
        ts_obs(i*12+j,:,:) = reshape(x(6:), (/72,128/))
    end do
end do
printVarSummary(ts_obs)

ts_obs = ts_obs(:,::-1,:)
time = ispan(1, nyear*12, 1)
ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)


lat = fspan(18.25, 53.75, 72)
lon = fspan(72.25, 135.75, 128)
nlon = dimsizes(lon)
nlat = dimsizes(lat)
; print(lat)
; print(lon)

ts_obs!0 = "time"
ts_obs!1 = "lat"
ts_obs!2 = "lon"
ts_obs&time = ispan(1, nyear*nm, 1)
ts_obs&lat = lat
ts_obs&lon = lon
ts_obs&lat@units = "degrees_north"
ts_obs&lon@units = "degrees_east"
printVarSummary(ts_obs)

topox1 = ftopo->topo({18:54},{72:136})
lonx = topox1&lon
latx = topox1&lat
printVarSummary(topox1)
topox2 = linint2_Wrap(lonx, latx, topox1, False, lon, lat, 0)
print(topox2({35},{90}))

TS_obs = mask(ts_obs, topox2.ge.1500, True)
copy_VarCoords(ts_obs, TS_obs)
printVarSummary(TS_obs)
; TS_obs_sea = month_to_seasonN(TS_obs, season)
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

write_table("time_series_TP_t2m_Mar_Apr_May_OBS.txt", "w", [/"YEAR","Mar","Apr","May"/], "%4s %4s %4s %4s")
write_table("time_series_TP_t2m_Mar_Apr_May_OBS.txt", "a", [/year,TP_obs_std(0,:),TP_obs_std(1,:),TP_obs_std(2,:)/], "%0.4i %4.2f %4.2f %4.2f %4.2f")


;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "calc_time_series_TP_t2m_Mar_Apr_May"
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

res@gsnRightString = ""
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, year, TP_obs_std, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"March","April","May"/)

;============================================================
    gx = (/1980,1982/)
    text_pos = 1983
    gy1 = (/2.8,2.8/)
    interval = 0.5
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
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = 2
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

draw(plot)

end