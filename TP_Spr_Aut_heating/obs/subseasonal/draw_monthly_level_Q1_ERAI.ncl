
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)

diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
fili = "q1.mon.mean.1979.2016.ERAI.nc" ; ERAI 

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

q1 = f->q1(iYYYY,:,:,:)    ; units K/s
q1 = q1*86400
q1@units = "K/day"
printVarSummary(q1)

; Cp = 1004.07 ; J/(Kg*K)
; q1 = q1 * Cp ; units W/Kg

filisp = "/home/yangsong3/data-observation/ERA-interim/monthly/Surface_pressure.mon.mean.nc"
fsp = addfile(filisp, "r")

sp0 = short2flt(fsp->sp(iYYYY,:,:))
printVarSummary(sp0)
sp = linint2_Wrap(sp0&longitude, sp0&latitude, sp0, True, q1&lon, q1&lat, 0)

; q1 = mask(q1, sp.le., mvalue)
; ;----------------integration-----------------------------
; plev = q1&lev*100.
; psfc = sp
; ptop = min(plev)
; dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
; printVarSummary(dP)

; ;dp = conform_dims(dimsizes(q1), dP, 1)
; q1dp = q1 * dP
; copy_VarCoords(q1, q1dp)
; printVarSummary(q1dp)
; Q1 = dim_sum_n_Wrap(q1dp(:,{1000:1},:,:), 1)
; copy_VarCoords(q1(:,0,:,:), Q1)
; printVarSummary(Q1)
; Q1 = Q1/9.8  ; units W/m2
; print(q1(30:50,35,{33},{95}))
; exit

t = clmMonTLLL(q1)
printVarSummary(t)

t_12mn = wgt_areaave_Wrap(t(:,:,{25:40},{75:105}), 1, 1, 1)
;t_12mn = t(:,:,{33},{95})
ttt_12mn = t_12mn(lev|:,month|:)
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
figname = "draw_monthly_level_Q1_ERAI"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

; res@trXMinF = 1
; res@trXMaxF = 12
; res@tmXBMode = "Explicit"
; res@tmXBValues = Month
; res@tmXBLabels = Month
; res@tmXTOn = False
; res@tiXAxisString = "Month"

; res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.2
res@cnMaxLevelValF = 3.
res@cnMinLevelValF = -3.
res@cnFillPalette = "CBR_coldhot"
res@gsnLeftString = "ERAI_heating_rate"

plot = gsn_csm_pres_hgt(wks, ttt_12mn({1000:100},:), res)
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