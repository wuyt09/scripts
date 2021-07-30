begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

yrStrt = 1951
yrLast = 2007

year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"


diri = "/home/yangsong3/wuyt/sysu/scripts/model_precip/APHRO/"
fili = "APHRO_MA_025deg_V1003R1.monmean.1951-2007.nc"
f = addfile(diri+fili, "r")

precip = f->precip
printVarSummary(precip)

prec = month_to_season(precip, season)
printVarSummary(prec)


pr_TP = wgt_areaave_Wrap(prec(:,{28:40},{85:105}), 1., 1., 0)
pr_IS = wgt_areaave_Wrap(prec(:,{18:28},{72:85}), 1., 1., 0)

pp = new((/4,nyear/), float)
pp(0,:) = pr_TP - pr_TP(0)
pp(1,:) = pr_IS - pr_IS(0)
pp(2,:) = runave_n(pp(0,:), 10, 0, 0)
pp(3,:) = runave_n(pp(1,:), 10, 0, 0)

corr_tp = regCoef_n(year, pr_TP, 0, 0)
corr_is = regCoef_n(year, pr_IS, 0, 0)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey40"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "trend_pr_CETP_NorIn_time_series_JJA_APHRO_MA"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = year
; res@tmXBLabels = year
res@tmXTOn = False
; res@tiXAxisString = "Month"

res@gsnYRefLine = 0
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineThicknessF = 8
; res@trYMinF = -3
; res@trYMaxF = 3
; res@tmYROn = False
; res@tiYAxisString = "Precip std(mm/day)"
; res@gsnLeftString = "Corr: "+sprintf("%4.2f", corr6)
; res@gsnRightString = "Sig: "+sprintf("%4.0f", 100*(1-p6))+"%"
res@gsnStringFontHeightF = 0.015

Colors = (/"black","red","black","red"/)
Patterns = (/0,0,16,16/)
res@xyLineThicknessF = 8
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, year, pp, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"CE-TP: "+sprintf("%4.2f", corr_tp),"North-India: "+sprintf("%4.2f", corr_is),"CE-TP: 10-yr runnung ave","North-India: 10-yr runnung ave"/)

;============================================================
    gx = (/1955,1960/)
    text_pos = 1961
    interval = 0.4
    gy1 = (/3.6,3.6/)
    gy2 = gy1 - interval

    gx2 = (/1975,1980/)
    text_pos2 = 1981
    gy3 = gy1
    gy4 = gy2

    res_lines@gsLineThicknessF = 15
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
    dum31 = gsn_add_polyline(wks, plot(0), gx2, gy1, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos2, gy1, res_text)
    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx2, gy2, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos2, gy2, res_text)

draw(plot)


end