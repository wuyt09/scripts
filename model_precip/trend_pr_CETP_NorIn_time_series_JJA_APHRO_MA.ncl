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

; prec = month_to_season(precip, season)
prec = precip(6::12,:,:) + precip(7::12,:,:)
copy_VarCoords(precip(6::12,:,:), prec)
printVarSummary(prec)


pr_tp = wgt_areaave_Wrap(prec(:,{26:33},{85:100}), 1., 1., 0)
; pr_IS = wgt_areaave_Wrap(prec(:,{18:28},{72:85}), 1., 1., 0)
pr_is = wgt_areaave_Wrap(prec(:,{15:30},{70:80}), 1., 1., 0)


rc_tp = regline(year, pr_tp)
rc_is = regline(year, pr_is)

trend_tp = year*rc_tp+rc_tp@yintercept
trend_is = year*rc_is+rc_is@yintercept

line_tp = runave_n_Wrap(pr_tp, 10, 0, 0)
line_is = runave_n_Wrap(pr_is, 10, 0, 0)

lines = new((/4,nyear/), float)
lines(0,:) = pr_tp
lines(1,:) = pr_is
lines(2,:) = line_tp
lines(3,:) = line_is


corr = escorc(pr_tp, pr_is)
corr_10 = escorc(line_tp, line_is)
pr_tp_dtrend = dtrend_msg(year, pr_tp, False, False)
pr_is_dtrend = dtrend_msg(year, pr_is, False, False)
pr_tp_dtrend_10 = runave_n_Wrap(pr_tp_dtrend, 10, 0, 0)
pr_is_dtrend_10 = runave_n_Wrap(pr_is_dtrend, 10, 0, 0)
corr_dtrend_10 = escorc(pr_tp_dtrend_10, pr_is_dtrend_10)
print(sprintf("%4.2f", corr_dtrend_10))
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
res@gsnStringFontHeightF = 0.025

Colors = (/"black","red","blue","purple"/)
Patterns = (/0,0,0,0/)
res@xyLineThicknessF = 8
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@gsnLeftString = "Corr: "+sprintf("%4.2f", corr)
res@gsnRightString = "10-yr runing ave corr: "+sprintf("%4.2f", corr_10)
plot(0) = gsn_csm_xy(wks, year, lines, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"CE-TP: "+sprintf("%4.2f", rc_tp*10)+"mm/day/decade","North-India: "+sprintf("%4.2f", rc_is*10)+"mm/day/decade","CE-TP: 10-yr runnung ave","North-India: 10-yr runnung ave"/)

;============================================================
    gx = (/1955,1960/)
    text_pos = 1961
    interval = 0.8
    gy1 = (/9,9/)
    gy2 = gy1 - interval

    gx2 = (/1975,1980/)
    text_pos2 = 1981
    gy3 = gy2 - interval
    gy4 = gy3 - interval

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
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)

draw(plot)


end