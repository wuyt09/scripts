
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = ispan(1, 12, 1)

diri = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
fili = "air.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->air(iYYYY,:,:,:))
printVarSummary(T)

; filisp = "Surface_pressure.mon.mean.nc"
; fsp = addfile(diri+filisp, "r")
; sp = short2flt(fsp->sp(iYYYY,:,:))
; printVarSummary(sp)
;===========integration========================================
lev = T&level
plev = lev*100.
psfc = max(plev)
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!0 = "lev"
dP&lev = lev
printVarSummary(dP)

dp = conform_dims(dimsizes(T), dP, 1)
Tdp = T*dp
copy_VarCoords(T, Tdp)
TTT = dim_sum_n_Wrap(Tdp(:,{400:200},:,:), 1)/sum(dP({400:200}))
;dim_sum_n_Wrap(dP(:,{400:200},:,:),1)
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)


ttt = clmMonTLL(TTT)
printVarSummary(ttt)

ttt_12mn = wgt_areaave_Wrap(ttt(:,{25:40},{75:105}), 1, 1, 0)

dt = ttt_12mn
do i = 1, 11
    dt(i) = ttt_12mn(i) - ttt_12mn(i-1)
end do
dt(0) = ttt_12mn(0) - ttt_12mn(11)
;dt = abs(dt)
print(dt)

rc1 = regline(Month(2:5), dt(2:5))
rc2 = regline(Month(7:10), dt(7:10))
print("rc1 "+rc1+"  rc2 "+rc2) ; rc1 1.38291  rc2 -1.10443

;================================================================
; plot
;================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/"
figname = "draw_monthly_tropo_temp_NCEP2"

lonstr = 60 ;20
lonend = 110 ;220
latstr = 20
latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(12, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False

res@mpOutlineOn = True
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpMinLatF = latstr
res@mpMaxLatF = latend

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 2.
res@cnMaxLevelValF = 250.
res@cnMinLevelValF = 220.
res@cnFillPalette = "CBR_coldhot"

do i = 0, 11
    res@gsnLeftString = "Month = "+(i+1)
    plot(i) = gsn_csm_contour_map(wks, ttt(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "tropo (400hPa-200hPa) temp"

gsn_panel(wks, plot, (/4,3/), pres)

;====================================================================
delete(res)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = 1
res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = Month
res@tmXBLabels = Month
res@tmXTOn = False
res@tiXAxisString = "Month"

; res@trYMinF = -16
; res@trYMaxF = 20
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-15, 20, 5)
; res@tmYLLabels = ispan(-15, 20, 5)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","blue","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns
res@gsnLeftString = "NCEP2 tropo temp"


plot(0) = gsn_csm_xy(wks, Month, ttt_12mn, res)
draw(plot(0))
frame(wks)



end