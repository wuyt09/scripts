
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)
nyear = 25

season = "DJF"
cmed = new((/ncase,25,96,144/), float)
chgh = cmed

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.CLDMED.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.CLDHGH.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    CLDMED = f1->CLDMED(iYYYY,:,:)
    CLDHGH = f2->CLDHGH(iYYYY,:,:)

    cmed(i,:,:,:) = month_to_season(CLDMED, season)
    chgh(i,:,:,:) = month_to_season(CLDHGH, season)
end do
printVarSummary(cmed)

;=======================================================================
; ttest
;=======================================================================

cmed_ave = dim_avg_n_Wrap(cmed, 1)
cmed_var = dim_variance_n_Wrap(cmed, 1)
prob_cmed = ttest(cmed_ave(0,:,:), cmed_var(0,:,:), nyear, cmed_ave(1,:,:), cmed_var(1,:,:), nyear, False, False)
copy_VarCoords(cmed_ave(0,:,:), prob_cmed)
diff_cmed = cmed_ave(0,:,:)
diff_cmed = cmed_ave(1,:,:) - cmed_ave(0,:,:)

chgh_ave = dim_avg_n_Wrap(chgh, 1)
chgh_var = dim_variance_n_Wrap(chgh, 1)
prob_chgh = ttest(chgh_ave(0,:,:), chgh_var(0,:,:), nyear, chgh_ave(1,:,:), chgh_var(1,:,:), nyear, False, False)
copy_VarCoords(chgh_ave(0,:,:), prob_chgh)
diff_chgh = chgh_ave(0,:,:)
diff_chgh = chgh_ave(1,:,:) - chgh_ave(0,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_MedCloud_HighCloud"

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
res@cnFillPalette         = "WhiteBlue"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@pmLabelBarOrthogonalPosF = 0.15

res@gsnStringFontHeightF = 0.035


res@gsnRightString = "Ctrl"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.1
res@cnMaxLevelValF = 0.6
res@cnLevelSpacingF = 0.05
res@gsnLeftString = "Medium Cloud"
plot(0) = gsn_csm_contour_map(wks, cmed_ave(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.1
res@cnMaxLevelValF = 0.6
res@cnLevelSpacingF = 0.05
res@gsnLeftString = "High Cloud"
plot(1) = gsn_csm_contour_map(wks, chgh_ave(0,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

;==========================================================
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
 res3@cnFillColor  =  "green4"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

res@cnFillPalette = "MPL_RdYlBu"
res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnLeftString = "Medium Cloud"
plot(2) = gsn_csm_contour_map(wks, diff_cmed, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot_prob(2) = gsn_csm_contour(wks,prob_cmed,res3)
plot_prob(2) = ShadeLtContour(plot_prob(2),0.055,17)
overlay(plot(2), plot_prob(2))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnLeftString = "High Cloud"
plot(3) = gsn_csm_contour_map(wks, diff_chgh, res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))
plot_prob(3) = gsn_csm_contour(wks,prob_chgh,res3)
plot_prob(3) = ShadeLtContour(plot_prob(3),0.055,17)
overlay(plot(3), plot_prob(3))


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean"
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,2/), pres)

end