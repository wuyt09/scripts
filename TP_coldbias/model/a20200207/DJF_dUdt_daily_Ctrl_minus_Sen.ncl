
begin
lonstr = 80
lonend = 105
latstr = 30
latend = 40
yrStrt = 1981
yrLast = 2005
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"
LEV = 200
dudt = new((/ncase,nyear,12,96,144/), float)

diri = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200207/dudt/"

do j = 0, ncase-1
    do k = 0, nyear-1
        nyr = year(k)
        filiu = casename(j)+".dudt.monmean."+nyr+".nc"

        fu = addfile(diri+filiu, "r")
        dudt(j,k,:,:,:) = fu->dudt(:,{LEV},:,:)
    end do
end do
printVarSummary(dudt)

dUdt = reshape(dudt, (/ncase,nyear*12,96,144/))
dUdt!0 = "case"
dUdt&case = casename
dUdt!1 = "month"
copy_VarCoords(dudt(0,0,0,:,:), dUdt(0,0,:,:))
printVarSummary(dUdt)

dUdTX = month_to_season(dUdt(month|:,case|:,lat|:,lon|:), season)
printVarSummary(dUdTX)
dUdT = dim_avg_n_Wrap(dUdTX, 0)
dUdT = dUdT*1.e+4

diff_dudt = dUdT(0,:,:)
diff_dudt = dUdT(1,:,:) - dUdT(0,:,:)
printMinMax(dUdT(:,{latstr:latend},{lonstr:lonend}), True)
printMinMax(diff_dudt({latstr:latend},{lonstr:lonend}), True)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_dUdt_daily_Ctrl_minus_Sen"

lonstr = 60
lonend = 120
latstr = 10
latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
res@cnFillPalette         = "BlueWhiteOrangeRed" ;"MPL_RdBu"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 0.5
; res@cnMinLevelValF = -2
; res@cnMaxLevelValF = 2
; res@cnLevelSpacingF = 0.2
res@gsnStringFontHeightF = 0.035
;res@gsnRightString       = "K"


res@gsnLeftString = "Ctrl run"
plot(0) = gsn_csm_contour_map(wks, dUdT(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "Sen run"
plot(1) = gsn_csm_contour_map(wks, dUdT(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.1
res@gsnLeftString = "diff"
plot(2) = gsn_csm_contour_map(wks, diff_dudt, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean: dU/dt"
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/1,3/), pres)

end