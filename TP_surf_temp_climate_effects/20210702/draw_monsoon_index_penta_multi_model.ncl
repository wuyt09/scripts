begin

case=(/"CTRL","NOQ1","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(case)

x = new((/ncase,73,5/), float)
do i = 0, ncase-1
    fili = "monsoon_index_pentad_"+case(i)+".txt"
    x_med = asciiread(fili, (/73,5/), "float")
    x(i,:,:) = x_med
end do
lines = x(:,:,1:4)

title = (/"WY","SEA","SA","EA"/)
;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210702/"
figname = "draw_monsoon_index_penta_multi_model"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(4, graphic)

Colors = (/"black","red","green4","purple"/)
Patterns = (/0,0,0,0/)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@trXMinF = 1
res@trXMaxF = 73
res@tmXBPrecision = 2
res@tiXAxisString = "Pentad"
; res@trYMaxF = 3
; res@trYMinF = -3

res@gsnYRefLine = (/0/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

; res@gsnRightString = "Corr: "+sprintf("%4.2f", corr)
; res@gsnLeftString = "CN05.1"
res@gsnStringFontHeightF = 0.03
res@tiYAxisString = ""

do i = 0, 3
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_xy(wks, ispan(1, 73, 1), lines(:,:,i), res)
end do

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = case

;============================================================
    gx = (/2,8/)
    text_pos = 8.5

    gy1 = (/28,28/)
    interval = 5
    gy2 = gy1 - interval
    gy3 = gy2 - interval
    gy4 = gy3 - interval

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

pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/2,2/), pres)
end