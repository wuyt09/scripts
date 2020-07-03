
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
fili = "air.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->t(iYYYY,:,:,:))
printVarSummary(T)

t = clmMonTLLL(T)
printVarSummary(t)

t_12mn = wgt_areaave_Wrap(t(:,:,{25:40},{75:105}), 1, 1, 0)
ttt_12mn = t_12mn(level|:,month|:)
printVarSummary(ttt_12mn)
; dt = ttt_12mn
; do i = 1, 11
;     dt(i) = ttt_12mn(i) - ttt_12mn(i-1)
; end do
; dt(0) = ttt_12mn(0) - ttt_12mn(11)
; ;dt = abs(dt)
; print(dt)

; rc1 = regline(Month(2:5), dt(2:5))
; rc2 = regline(Month(7:10), dt(7:10))
; print("rc1 "+rc1+"  rc2 "+rc2) ;rc1 1.48554  rc2 -1.20139

; exit

;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/"
figname = "draw_monthly_level_tenmp_ERAI"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False

; res@mpOutlineOn = True
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend
; res@mpMinLatF = latstr
; res@mpMaxLatF = latend

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

; res@lbLabelBarOn = False
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 2.
; res@cnMaxLevelValF = 250.
; res@cnMinLevelValF = 220.
res@cnFillPalette = "CBR_coldhot"

; do i = 0, 11
;     res@gsnLeftString = "Month = "+(i+1)
;     plot(i) = gsn_csm_contour_map(wks, ttt(i,:,:), res)
; end do
plot = gsn_csm_pres_hgt(wks, ttt_12mn({1000:50},:), res)
draw(plot)
frame(wks)
; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; ;pres@gsnPanelMainString = "tropo (400hPa-200hPa) temp"

; gsn_panel(wks, plot, (/4,3/), pres)

; ;====================================================================
; delete(res)

; res = True
; res@gsnDraw  = False
; res@gsnFrame = False
; res@gsnMaximize = True

; res@vpWidthF = 0.7
; res@vpHeightF = 0.4

; res@trXMinF = 1
; res@trXMaxF = 12
; res@tmXBMode = "Explicit"
; res@tmXBValues = Month
; res@tmXBLabels = Month
; res@tmXTOn = False
; res@tiXAxisString = "Month"

; ; res@trYMinF = -16
; ; res@trYMaxF = 20
; ; res@tmYLMode = "Explicit"
; ; res@tmYLValues = ispan(-15, 20, 5)
; ; res@tmYLLabels = ispan(-15, 20, 5)
; res@tmYROn = False
; res@tiYAxisString = "Temperature (~S~o~N~C)"

; ;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
; Colors = (/"black","blue","forestgreen"/)
; Patterns = (/0,0,0/)
; res@xyLineThicknessF = 5
; res@xyLineColors = Colors
; res@xyDashPatterns = Patterns


; plot(0) = gsn_csm_xy(wks, Month, ttt_12mn, res)
; draw(plot(0))
; frame(wks)



end