
begin
yrStrt = 1979
yrLast = 2017
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)


inx_name = (/"TP_t2m","TTGI","WYI","WFI"/)

fili = "./ISM_onset_pentad.txt"
x = asciiread(fili, (/nyear,4/), "float")
print(x)
x(38,2) = x@Fill_Value
x(18,2) = 27

corr = new(3, float)
do i = 1, 3
    corr(i-1) = escorc(x(:,0), x(:,i))
end do

print(corr)

x!0 = "var"
x!1 = "year"
y = x(year|:,var|:)
;=============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201116/"
figname = "calc_corr_time_series_ISM_onset_Peatad_TP_surf_T"



wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = Month
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; res@tmXTOn = False
; res@tiXAxisString = "Month"

res@trYMinF = 22
res@trYMaxF = 36
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(22, 36, 1)
res@tmYLLabels = "P"+res@tmYLValues
; res@tmYROn = False
; res@tiYAxisString = "Temperature (~S~o~N~C)"

;Colors = (/"black","deeppink2", "orange2", "green4", "black" /);
Colors = (/"red","blue","forestgreen","purple"/)
Patterns = (/0,0,0,0/)
res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

; plot(0) = gsn_csm_xy(wks, year, ptime, res)

res2 = True
res2@gsnDraw  = False
res2@gsnFrame = False
res2@gsnMaximize = True

res2@vpWidthF = 0.7
res2@vpHeightF = 0.4

res2@trXMinF = yrStrt
res2@trXMaxF = yrLast

res2@xyLineThicknessF = 10
res2@xyLineColor = "black"
res2@xyDashPattern = 0


plot(0) = gsn_csm_xy2(wks, year, y(1:3,:), y(0,:), res, res2)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = inx_name(1:dimsizes(inx_name)-1)

;============================================================
    gx = (/1983,1985/)
    text_pos = 1985.5

    interval = 1
    gy1 = (/35,35/)
    gy2 = gy1 - interval
    gy3 = gy2 - interval
    ; gy4 = gy3 - interval

    gx2 = (/1991,1993/)
    text_pos2 = 1993.5
    gy5 = gy1
    ; gy6 = (/8,8/)

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
    ; res_lines@gsLineColor = Colors(3)
    ; res_lines@gsLineDashPattern = Patterns(3)
    ; dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    ; dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
    res_lines@gsLineColor = "black"
    res_lines@gsLineDashPattern = 0
    dum51 = gsn_add_polyline(wks, plot(0), gx2, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), "TP_t2m", text_pos2, gy5, res_text)

    draw(plot(0))

end