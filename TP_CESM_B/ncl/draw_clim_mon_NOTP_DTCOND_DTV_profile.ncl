
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40

yrStat = 260
yrLast = 304 

diri = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"

fili1 = "TP_NOTOPO.plevel.DTCOND.mon.260-304.nc"
fili2 = "TP_NOTOPO.plevel.DTV.mon.260-304.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")
x1 = f1->DTCOND
x2 = f2->DTV
printVarSummary(x1)
x1 = x1*86400
x2 = x2*86400

lev = x1&lev
lev!0 = "level"
lev&level = lev

total = x1
total = x1+x2
printVarSummary(total)

;x11 = clmMonTLLL(x1)
;x22 = dim_avg_n_Wrap(x2, 0)
Total = clmMonTLLL(total)
printVarSummary(Total)

totalx = wgt_areaave_Wrap(Total(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)

;=========================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_mon_NOTP_DTCOND_DTV_profile"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 12
plot = new(nplot, graphic)

 res = True
 res@gsnDraw  = False
 res@gsnFrame = False

 res@gsnXRefLine = 0
 res@gsnXRefLineDashPattern = 2
 res@gsnXRefLineColor = "grey"
 res@gsnXRefLineThicknessF = 3
 res@trXMinF = -6.
 res@trXMaxF = 6
 res@trYMinF   = 100
 res@trYMaxF   = 1000
 res@trYLog = True
 
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)
 res@xyLineColors = "red"
 res@xyLineThicknessF = 3
 res@xyDashPattern = 0

 res@trYReverse = True
 
 res@gsnRightString = "K"
 res@gsnStringFontHeightF = 0.035

 do j = 0,11
    res@gsnLeftString = "Month = "+(j+1)
    plot(j) = gsn_csm_xy(wks, totalx(j,{1000:50}), lev({1000:50}), res)  ;
 end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "NO_TP (Sum of heating entire TP)"

gsn_panel(wks, plot, (/3,4/), pres)

end