
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

fili1 = "./TTT/time_series_TP_TT_Mar_Apr_May_ERA5.txt"
fili2 = "./Q1/time_series_TP_Q1_Mar_Apr_May_ERA5_2015.txt"

X1 = asciiread(fili1, (/nyear,4/), "float")
X2 = asciiread(fili2, (/nyear,4/), "float")


lines1 = new((/2,nyear/), float)
lines2 = lines1
lines1(0,:) = X1(:,2)
lines1(1,:) = X2(:,2)
lines2(0,:) = X1(:,3)
lines2(1,:) = X2(:,3)
;============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210708_obs/"
figname = "draw_TT_Q1_time_series"
Colors = (/"green4","goldenrod2"/)

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 2
plot = new(nplot, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPattern = 0

res@trXMinF = min(year)
res@trXMaxF = max(year)

res@gsnRightString = ""
res@gsnStringFontHeightF = 0.025

res@gsnYRefLine = (/1.,0.7,-0.7,-1./)
res@gsnYRefLineColor = "gray54"
res@gsnYRefLineThicknesses = 3
res@gsnYRefLineDashPattern = 2
; res@gsnLeftString = "TTT ("+levs+"-"+levt+" hPa)"

res@gsnLeftString = "Apr"
plot(0) = gsn_csm_xy(wks, year, lines1, res)
res@gsnLeftString = "May"
plot(1) = gsn_csm_xy(wks, year, lines2, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"TPTT","TPQ1"/)

; ;============================================================
    interval = 0.3
    gx = (/1985,1988/)
    text_pos = 1989
    gy1 = (/-2.0,-2.0/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    gy4 = gy3-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    ;res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/1,2/), pres)

end