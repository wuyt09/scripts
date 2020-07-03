
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
cliq = new((/ncase,25,96,144/), float)
cice = cliq
qq = cliq
cc = cliq

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.CLOUD.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.CLDICE.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.CLDLIQ.1979-2005.nc"
    fili4 = casename(i)+".cam.h0.PS.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    CC = f1->CLOUD(iYYYY,:,:,:)
    CI = f2->CLDICE(iYYYY,:,:,:)
    CL = f3->CLDLIQ(iYYYY,:,:,:)
    PS = f4->PS(iYYYY,:,:)


    QU = CI
    QV = CL

    lev = CI&level
    plev = lev
    plev = plev *100
    plev@units = "Pa"

    ptop = min(plev)
    psfc = PS
    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)

    QUdP = QU
    QUdP = QU*dP
    QVdP = QV
    QVdP = QV*dP
    CCdP = CC
    CCdP = CC*dP

    ; QUX = dim_sum_n_Wrap(QUdP(:,{1000:200},:,:), 1)
    ; QVX = dim_sum_n_Wrap(QVdP(:,{1000:200},:,:), 1)
    ; CCX = dim_sum_n_Wrap(CCdP(:,{1000:200},:,:), 1)/dim_sum_n_Wrap(dP(:,{1000:200},:,:), 1)
    ; copy_VarCoords(CCdP(:,0,:,:), CCX)
    QUX = dim_sum_n_Wrap(QUdP, 1)
    QVX = dim_sum_n_Wrap(QVdP, 1)
    CCX = dim_sum_n_Wrap(CCdP, 1)  ;/dim_sum_n_Wrap(dP, 1)
    copy_VarCoords(CC(:,0,:,:), CCX)
    QUX = QUX/9.8
    QVX = QVX/9.8
    CCX = CCX/9.8

    cice(i,:,:,:) = month_to_season(QUX, season)
    cliq(i,:,:,:) = month_to_season(QVX, season)
    cc(i,:,:,:) = month_to_season(CCX, season)

    qq = cice
    qq = cliq+cice
    qq@long_name = "cloud liq+ice"

end do
printVarSummary(qq)

;=======================================================================
; ttest
;=======================================================================

qq_ave = dim_avg_n_Wrap(qq, 1)
qq_var = dim_variance_n_Wrap(qq, 1)
prob_qq = ttest(qq_ave(0,:,:), qq_var(0,:,:), nyear, qq_ave(1,:,:), qq_var(1,:,:), nyear, False, False)
copy_VarCoords(qq_ave(0,:,:), prob_qq)
diff_qq = qq_ave(0,:,:)
diff_qq = qq_ave(1,:,:) - qq_ave(0,:,:)

cc_ave = dim_avg_n_Wrap(cc, 1)
cc_var = dim_variance_n_Wrap(cc, 1)
prob_cc = ttest(cc_ave(0,:,:), cc_var(0,:,:), nyear, cc_ave(1,:,:), cc_var(1,:,:), nyear, False, False)
copy_VarCoords(cc_ave(0,:,:), prob_cc)
diff_cc = cc_ave(0,:,:)
diff_cc = cc_ave(1,:,:) - cc_ave(0,:,:)

cice_ave = dim_avg_n_Wrap(cice, 1)
cice_var = dim_variance_n_Wrap(cice, 1)
prob_cice = ttest(cice_ave(0,:,:), cice_var(0,:,:), nyear, cice_ave(1,:,:), cice_var(1,:,:), nyear, False, False)
copy_VarCoords(cice_ave(0,:,:), prob_cice)
diff_cice = cice_ave(0,:,:)
diff_cice = cice_ave(1,:,:) - cice_ave(0,:,:)

cliq_ave = dim_avg_n_Wrap(cliq, 1)
cliq_var = dim_variance_n_Wrap(cliq, 1)
prob_cliq = ttest(cliq_ave(0,:,:), cliq_var(0,:,:), nyear, cliq_ave(1,:,:), cliq_var(1,:,:), nyear, False, False)
copy_VarCoords(cliq_ave(0,:,:), prob_cliq)
diff_cliq = cliq_ave(0,:,:)
diff_cliq = cliq_ave(1,:,:) - cliq_ave(0,:,:)

; printMinMax(cc_ave(0,:,:), False)
; printMinMax(diff_cc, False)
; exit
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_cloudcover_Cice_Cliq"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)
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
res@cnMinLevelValF = 0.
res@cnMaxLevelValF = 1000.
res@cnLevelSpacingF = 100.
res@gsnLeftString = "Cloud Cover (whole column)"
plot(0) = gsn_csm_contour_map(wks, cc_ave(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.02
res@cnMaxLevelValF = 0.2
res@cnLevelSpacingF = 0.02
res@gsnLeftString = "Cloud Liq+Ice (whole column)"
plot(1) = gsn_csm_contour_map(wks, qq_ave(0,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.01
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnLeftString = "Cloud Ice (whole column)"
plot(2) = gsn_csm_contour_map(wks, cice_ave(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.01
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnLeftString = "Cloud Liq (whole column)"
plot(3) = gsn_csm_contour_map(wks, cliq_ave(0,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

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
res@cnMinLevelValF = -300.
res@cnMaxLevelValF = 300.
res@cnLevelSpacingF = 30.
res@gsnLeftString = "Cloud Cover (whole column)"
plot(4) = gsn_csm_contour_map(wks, diff_cc, res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))
plot_prob(4) = gsn_csm_contour(wks,prob_cc,res3)
plot_prob(4) = ShadeLtContour(plot_prob(4),0.055,17)
overlay(plot(4), plot_prob(4))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.02
res@cnMaxLevelValF = 0.02
res@cnLevelSpacingF = 0.002
res@gsnLeftString = "Cloud Liq+Ice (whole column)"
plot(5) = gsn_csm_contour_map(wks, diff_qq, res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))
plot_prob(5) = gsn_csm_contour(wks,prob_qq,res3)
plot_prob(5) = ShadeLtContour(plot_prob(5),0.055,17)
overlay(plot(5), plot_prob(5))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.01
res@cnMaxLevelValF = 0.01
res@cnLevelSpacingF = 0.001
res@gsnLeftString = "Cloud Ice"
plot(6) = gsn_csm_contour_map(wks, diff_cice, res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))
plot_prob(6) = gsn_csm_contour(wks,prob_cice,res3)
plot_prob(6) = ShadeLtContour(plot_prob(6),0.055,17)
overlay(plot(6), plot_prob(6))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.02
res@cnMaxLevelValF = 0.02
res@cnLevelSpacingF = 0.002
res@gsnLeftString = "Cloud Liq"
plot(7) = gsn_csm_contour_map(wks, diff_cliq, res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7), plot_topo(7))
plot_prob(7) = gsn_csm_contour(wks,prob_cliq,res3)
plot_prob(7) = ShadeLtContour(plot_prob(7),0.055,17)
overlay(plot(7), plot_prob(7))



pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean"
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,4/), pres)

end