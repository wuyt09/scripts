
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"

begin

yrStrt = 1871
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)
year@units = "hours since 1-1-1 00:00:0.0"
; print(year)
nyear = dimsizes(year)
; print(nyear)


diri = "/home/yangsong3/wuyt/sysu/scripts/model_precip/IndianRainfall/"
fili1 = "Northeast_India_rainfall.txt" ; blue
fili2 = "Northwest_India_rainfall.txt" ; red
fili3 = "Westcentral_India_rainfall.txt" ; purple

x1_med = asciiread(diri+fili1, -1, "float")
x1 = onedtond(x1_med(13:), (/nyear,18/))
print(x1(nyear-1,:))
x1 = x1/10

x2_med = asciiread(diri+fili2, -1, "float")
x2 = onedtond(x2_med(13:), (/nyear,18/))
print(x2(nyear-1,:))
x2 = x2/10

x3_med = asciiread(diri+fili3, -1, "float")
x3 = onedtond(x3_med(13:), (/nyear,18/))
print(x3(nyear-1,:))
x3 = x3/10


pr_TP = (x1(:,7) + x1(:,8))/2/31
pr_is1 = (x2(:,7) + x2(:,8))/2/31
pr_is2 = (x3(:,7) + x3(:,8))/2/31
pr_IS = (pr_is1+pr_is2)/2
pr_TP!0 = "year"
pr_TP&year = year
copy_VarCoords(pr_TP, pr_IS)

pr_tp = pr_TP ;({yrs:yre})
pr_is = pr_IS ;({yrs:yre})
printVarSummary(pr_tp)

rc_tp = regline(year, pr_tp)
rc_is = regline(year, pr_is)

trend_tp = year*rc_tp+rc_tp@yintercept
trend_is = year*rc_is+rc_is@yintercept

line_tp = runave_n_Wrap(pr_tp, 10, 0, 0)
line_is = runave_n_Wrap(pr_is, 10, 0, 0)


lines = new((/4,nyr/), float)
lines(0,:) = pr_tp
lines(1,:) = pr_is
lines(2,:) = line_tp
lines(3,:) = line_is
; lines!0 = "region"
; lines!1 = "month"
; lines!2 = "year"
; lines&year = 

print(rc_tp)
print(rc_is)

corr = escorc(pr_tp, pr_is)
corr_10 = escorc(line_tp, line_is)

pr_tp_dtrend = dtrend_msg(year, pr_tp, False, False)
pr_is_dtrend = dtrend_msg(year, pr_is, False, False)
pr_tp_dtrend_10 = runave_n_Wrap(pr_tp_dtrend, 10, 0, 0)
pr_is_dtrend_10 = runave_n_Wrap(pr_is_dtrend, 10, 0, 0)
corr_dtrend_10 = escorc(pr_tp_dtrend_10, pr_is_dtrend_10)
print(sprintf("%4.2f", corr_dtrend_10))
;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "trend_pr_IMD"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)
; dum11 = plot
; ; dum12 = plot
; dum21 = plot
; ; dum22 = plot
; dum31 = plot
; ; dum32 = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = yrs
res@trXMaxF = yre
; res@tmXBMode = "Explicit"
; res@tmXBValues = year
; res@tmXBLabels = year
res@tmXTOn = False
; res@tiXAxisString = "Month"

res@gsnYRefLine = (/-0.368,0,0.368/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineDashPatterns = (/16,0,16/)
res@gsnYRefLineThicknessF = 5
; res@trYMinF = -0.8
; res@trYMaxF = 0.8
; res@tmYROn = False
; res@tiYAxisString = "Precip std(mm/day)"
; res@gsnLeftString = "Corr: "+sprintf("%4.2f", corr6)
; res@gsnRightString = "Sig: "+sprintf("%4.0f", 100*(1-p6))+"%"
res@gsnStringFontHeightF = 0.02


Colors = (/"blue","red","black","purple"/)
Patterns = (/0,0,0,0/)
res@xyLineThicknessF = 8
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@gsnLeftString = "Corr: "+sprintf("%4.2f", corr)
res@gsnRightString = "10-yr runing ave corr: "+sprintf("%4.2f", corr_10)
plot(0) = gsn_csm_xy(wks, yr, lines, res)
;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"Northeast: "+sprintf("%4.2f", rc_tp*10)+" mm/day/decade","West: "+sprintf("%4.2f", rc_is*10)+" mm/day/decade","Northeast: 10-yr runing ave","West: 10-yr runing ave"/)

;============================================================
    gx = (/1880,1890/)
    text_pos = 1892
    interval = 1.5
    gy1 = (/17.5,17.5/)
    gy2 = gy1 - interval

    gx2 = (/1950,1960/)
    text_pos2 = 1962
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