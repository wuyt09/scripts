
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

x2_med = asciiread(diri+fili2, -1, "float")
x2 = onedtond(x2_med(13:), (/nyear,18/))
print(x2(nyear-1,:))

x3_med = asciiread(diri+fili3, -1, "float")
x3 = onedtond(x3_med(13:), (/nyear,18/))
print(x3(nyear-1,:))

pr_tp = new((/nyear,4/), float)
pr_is1 = pr_tp
pr_is2 = pr_tp
pr_is3 = pr_tp
pr_tp(:,0:2) = x1(:,6:8)
pr_tp(:,3) = (x1(:,7) + x1(:,8))/2
pr_is1(:,0:2) = x2(:,6:8)
pr_is1(:,3) = (x2(:,7) + x2(:,8))/2
pr_is2(:,0:2) = x3(:,6:8)
pr_is2(:,3) = (x3(:,7) + x3(:,8))/2
pr_is3 = (pr_is1+pr_is2)/2
printVarSummary(pr_tp)

rc_tp = regCoef_n(year, pr_tp, 0, 0)
rc_is1 = regCoef_n(year, pr_is1, 0, 0)
rc_is2 = regCoef_n(year, pr_is2, 0, 0)
printVarSummary(pr_tp)
; print(sprintf("%4.2f", rc_tp))
; print(sprintf("%4.2f", rc_is1))
; print(sprintf("%4.2f", rc_is2))

cor1 = new((/4,126/), float)
cor2 = cor1
cor3 = cor1
do i = 0, 3
    cor1(i,:) = run_cor(pr_tp(:,i), pr_is1(:,i), year, 21)
    cor2(i,:) = run_cor(pr_tp(:,i), pr_is2(:,i), year, 21)
    cor3(i,:) = run_cor(pr_tp(:,i), pr_is3(:,i), year, 21)
end do
printVarSummary(cor1)

lines = new((/3,4,126/), float)
lines(0,:,:) = cor1
lines(1,:,:) = cor2
lines(2,:,:) = cor3
; lines!0 = "region"
; lines!1 = "month"
; lines!2 = "year"
; lines&year = 

yr = ispan(yrStrt+10, yrLast-10, 1)
title = (/"June","July","August","July+August"/)
;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "runcor_pr_IMD"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(4, graphic)
dum11 = plot
; dum12 = plot
dum21 = plot
; dum22 = plot
dum31 = plot
; dum32 = plot

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

res@gsnYRefLine = (/-0.368,0,0.368/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineDashPatterns = (/16,0,16/)
res@gsnYRefLineThicknessF = 5
res@trYMinF = -0.8
res@trYMaxF = 0.8
; res@tmYROn = False
; res@tiYAxisString = "Precip std(mm/day)"
; res@gsnLeftString = "Corr: "+sprintf("%4.2f", corr6)
; res@gsnRightString = "Sig: "+sprintf("%4.0f", 100*(1-p6))+"%"
res@gsnStringFontHeightF = 0.025


Colors = (/"green4","purple","red"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 8
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    ; Texts  = (/"Northeast & Northwest","Northeast & Centralwest","Northeast & 2 regions"/)
    Texts  = (/"Northeast & Northwest","Northeast & Centralwest","Northeast & 2 regions"/)

;============================================================
    gy1 = (/0.65,0.65/)

    yr0 = 1872
    yrinterval = 10
    gx = (/yr0,yr0+yrinterval/)
    text_pos = yr0+yrinterval+1

    yr0 = yr0+55
    gx2 = (/yr0,yr0+yrinterval/)
    text_pos2 = yr0+yrinterval+1

    yr0 = yr0+55
    gx3 = (/yr0,yr0+yrinterval/)
    text_pos3 = yr0+yrinterval+1
    
    res_lines@gsLineThicknessF = 8


do i = 0, 3
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_xy(wks, yr, lines(:,i,:), res)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11(i) = gsn_add_polyline(wks, plot(i), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(i), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21(i) = gsn_add_polyline(wks, plot(i), gx2, gy1, res_lines)
    dum22 = gsn_add_text(wks, plot(i), Texts(1), text_pos2, gy1, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31(i) = gsn_add_polyline(wks, plot(i), gx3, gy1, res_lines)
    dum32 = gsn_add_text(wks, plot(i), Texts(2), text_pos3, gy1, res_text)
end do



pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,2/), pres)

end