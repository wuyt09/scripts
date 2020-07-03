
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"

vt_sea = new((/ncase,25,181,360/), float)

do i = 0, ncase-1
    ;diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    ; fili1 = casename(i)+".cam.h0.T.1979-2005.nc"
    ; fili2 = casename(i)+".cam.h0.PS.1979-2005.nc"
    diri = "./"
    fili1 = casename(i)+".cam.h0.T.1979-2005.1x1.nc"
    fili2 = casename(i)+".cam.h0.PS.1979-2005.1x1.nc"
    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    
    TIME = f2->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    VT = f1->T(iYYYY,:,:,:)
    ; VT = VTX(:,1:7,:,:) ; 1000hPa-300hPa

    PS = f2->PS(iYYYY,:,:)
    lev = VT&level
    plev = lev
    plev = plev *100
    plev!0 = "lev"
    plev&lev = lev
    plev@units = "Pa"

    ptop = min(plev)
    psfc = PS
    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    ; dP(:,{600},{30},{90}) = dP@_FillValue
    printVarSummary(dP)
    dP = where(ismissing(VT), dP@_FillValue, dP)
    ; exit

    VTdP = VT*dP
    copy_VarCoords(VT, VTdP)
    levs = 1000
    levt = 450
    vt = dim_sum_n_Wrap(VTdP(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dP(:,{levs:levt},:,:), 1)
    copy_VarCoords(VTdP(:,0,:,:), vt)
    vt_sea(i,:,:,:) = month_to_season(vt, season)
end do
printVarSummary(vt_sea)

vt_ave = dim_avg_n_Wrap(vt_sea, 1)
vt_var = dim_variance_n_Wrap(vt_sea, 1)
prob = ttest(vt_ave(0,:,:), vt_var(0,:,:), 25, vt_ave(1,:,:), vt_var(1,:,:), 25, False, False)
copy_VarCoords(vt_ave(0,:,:), prob)
printVarSummary(prob)

vt_diff = vt_ave(0,:,:)
vt_diff = vt_ave(1,:,:) - vt_ave(0,:,:)
printVarSummary(vt_diff)

print(vt_ave(:,{38},{98}))
print(vt_diff({38},{98}))
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = season+"_AtmosT"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(4, graphic)
plot_topo = plot
plot_prob = plot

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
res@cnFillPalette         = "BlueDarkRed18"
res@cnFillMode            = "CellFill"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035

res@gsnRightString = "K"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 240
res@cnMaxLevelValF = 280
res@cnLevelSpacingF = 2.
res@gsnLeftString = "Ctrl"
plot(0) = gsn_csm_contour_map(wks, vt(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
; draw(plot(0))
; frame(wks)
res@gsnLeftString = "Sen"
plot(1) = gsn_csm_contour_map(wks, vt(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
; draw(plot(1))
; frame(wks)

;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "red"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -2
res@cnMaxLevelValF = 2
res@cnLevelSpacingF = 0.2
res@gsnLeftString = "Diff"
plot(2) = gsn_csm_contour_map(wks, vt_diff, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot_prob(2) = gsn_csm_contour(wks, prob, res3)
plot_prob(2) = ShadeLtContour(plot_prob(2),0.055,17)
overlay(plot(2), plot_prob(2))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean AtmosT" ; (600hPa-300hPa)(500hPa to 200hPa)
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/1,3/), pres)

end