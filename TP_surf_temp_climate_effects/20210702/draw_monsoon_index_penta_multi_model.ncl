begin

case=(/"CTRL","NOQ1","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(case)

x = new((/ncase,4,73/), float)
do i = 0, ncase-1
    fili = "monsoon_index_pentad_"+case(i)+".txt"
    x(i,:,:) = asciiread(fili, (/4,73/), "float")
end do

title = (/"WY","WNP","Gosiwami","EA"/)
;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210702/"
figname = "draw_monsoon_index_penta_multi_model"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(4, graphic)

Colors = (/"black","red","green4","purple"/)
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

; res@trXMinF = yrStrt
; res@trXMaxF = yrLast
; res@trYMaxF = 3
; res@trYMinF = -3

; res@gsnYRefLine = (/0,1,-1/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

; res@gsnRightString = "Corr: "+sprintf("%4.2f", corr)
; res@gsnLeftString = "CN05.1"
res@gsnStringFontHeightF = 0.03
res@tiYAxisString = ""

do i = 0, 3
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_xy(wks, ispan(1, 73, 1), x(:,i,:), res)
end do

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = case

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

pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/2,2/), pres)
end