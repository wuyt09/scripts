
begin

lonstr = 0
lonend = 130
latstr = 10
latend = 60

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

latstr_east = 35
latend_east = 45
lonstr_east = 105
lonend_east = 140

latstr_west = 30
latend_west = 55
lonstr_west = -10
lonend_west = 60

latstr_south = 5 ;10
latend_south = 25
lonstr_south = 70 ;50
lonend_south = 90 ;100

latstr_north = 45
latend_north = 60
lonstr_north = 85
lonend_north = 120

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

region = (/"West","East","South","North","TP"/)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

fa      = addfile ("/home/yangsong3/wuyt/data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc","r")
; lsdata = doubletofloat(a->lsm(0,:,:))
lsdata = fa->lsm(0,::-1,:)
printVarSummary(lsdata)
;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "t2m_mon_1x1_1979_2019.nc"

; diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
; fili = "2_metre_temperature.mon.mean.nc"

; diri = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
; fili  = "T2m.mon.mean.7916.r1.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t2m = short2flt(f->t2m(iYYYY,::-1,:))
printVarSummary(t2m)

t2m_clm = month_to_season(t2m, season)
printVarSummary(t2m_clm)

lines_west = wgt_areaave_Wrap(t2m_clm(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
lines_east = wgt_areaave_Wrap(t2m_clm(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
lines_north = wgt_areaave_Wrap(t2m_clm(:,{latstr_north:latend_north},{lonstr_north:lonend_north}), 1., 1., 0)

lsm = linint2_Wrap(lsdata&longitude, lsdata&latitude, lsdata, True, t2m&longitude, t2m&latitude, 0)
printVarSummary(lsm)
delete(lsdata)
lsm1 = conform_dims(dimsizes(t2m_clm), lsm, (/1,2/))
x_clm = where(lsm1.ge.0.8, t2m_clm, t2m@_FillValue)
copy_VarCoords(t2m_clm, x_clm)
lines_south = wgt_areaave_Wrap(x_clm(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)

x_mask = mask(t2m_clm, topo.ge.1500, True)
copy_VarCoords(t2m_clm, x_mask)
lines_tp = wgt_areaave_Wrap(x_mask(:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
; lines_tp = wgt_areaave_Wrap(t2m_clm(:,:,{20:latstr_tp},{lonstr_tp:lonstr_tp}), 1., 1., 0)


lines = new((/5,nyear/), float)
lines(0,:) = lines_west - lines_west(0)
lines(1,:) = lines_east - lines_east(0)
lines(2,:) = lines_south - lines_south(0)
lines(3,:) = lines_north - lines_north(0)
lines(4,:) = lines_tp - lines_tp(0)
copy_VarCoords(t2m_clm(:,0,0), lines(0,:))
lines!0 = "region"
lines&region = region
printVarSummary(lines)

lines_runave = runave_n_Wrap(lines, 10, 0, 1)
printVarSummary(lines_runave)

; rc = regCoef_n(year, lines, 0, 2)
; copy_VarCoords(lines(:,:,0), rc)
; printVarSummary(rc)
; exit
lines_trend = lines
rc = lines(:,0)
do j = 0, 4
    x_med = regline(year, lines(j,:))
    rc(j) = x_med
    lines_trend(j,:) = year*x_med+x_med@yintercept
end do
printVarSummary(lines_trend)

; tval = onedtond(rc@tval, dimsizes(rc))
; df   = onedtond(rc@nptxy, dimsizes(rc))
; b    = tval
; b    = 0.5
; prob = betainc(df/(df+tval^2), df/2.0, b)
; copy_VarCoords(rc, prob)

; rc@long_name = "regression coefficient"
; prob@long_name = "probability"


trend = rc*10
copy_VarCoords(rc, trend)
trend@units = "K/decade"
printVarSummary(trend)

; pp = new((/3,5,nyear/), float)
; pp(0,:,:) = lines
; pp(1,:,:) = lines_runave
; pp(2,:,:) = lines_trend
; copy_VarCoords(lines, pp(0,:,:))
pp = new((/2*5,nyear/), float)
pp(0:4,:) = lines
; pp(5:9,:) = lines_runave

;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/20210415/"
figname = "time_series_"+season+"_surfT_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(1, graphic)


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.8
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = year
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@tmXTOn = False
; res@tiXAxisString = "Month"

res@gsnYRefLine = 0
res@gsnYRefLineColor = "gray25"
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineThicknessF = 8


; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-20, 20, 5)
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYLPrecision = 4
res@tmYROn = False
res@tiYAxisString = ""
res@gsnStringFontHeightF = 0.025

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"forestgreen","red","orange2","blue","black","forestgreen","red","orange2","blue","black"/)
; Patterns = (/0,16,16,0,0/)
Patterns = (/0,0,0,0,0,16,16,16,16,16/)

res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

; res@gsnRightString = "trend: "+sprintf("%4.2f", trend(j))+" K/decade"
res@gsnLeftString = season+" T2m trend ("+yrStrt+"-"+yrLast+")"
plot(0) = gsn_csm_xy(wks, year, pp, res)


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = region +" "+sprintf("%4.2f", trend)+" K/decade"

    res_lines@gsLineDashPattern = 0

;============================================================
    gx = (/1980,1982/)
    text_pos = 1982.5
    gy1 = (/2.8,2.8/)
    interval = 0.4
    gy2 = gy1 - interval
    gy3 = gy2 - interval
    gy4 = gy3 - interval
    gy5 = gy4 - interval
    ; gy6 = (/8,8/)

    res_lines@gsLineThicknessF = 10
    res_lines@gsLineColor = Colors(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ; res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ; res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    ; res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
    res_lines@gsLineColor = Colors(4)
    ; res_lines@gsLineDashPattern = Patterns(4)
    dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)

draw(plot(0))
end