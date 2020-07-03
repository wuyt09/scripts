
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
qu = new((/ncase,96,144/), float)
qv = qu
qq = qu
cc = qu

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

    QUX = dim_sum_n_Wrap(QUdP(:,{1000:200},:,:), 1)
    QVX = dim_sum_n_Wrap(QVdP(:,{1000:200},:,:), 1)
    CCX = dim_sum_n_Wrap(CCdP(:,{1000:200},:,:), 1)/dim_sum_n_Wrap(dP(:,{1000:200},:,:), 1)
    copy_VarCoords(CCdP(:,0,:,:), CCX)
    ; QUX = dim_sum_n_Wrap(QUdP, 1)
    ; QVX = dim_sum_n_Wrap(QVdP, 1)
    QUX = QUX/9.8
    QVX = QVX/9.8
    CCX = CCX/9.8

    QUX0 = month_to_season(QUX, season)
    QVX0 = month_to_season(QVX, season)
    CCX0 = month_to_season(CCX, season)

    qu(i,:,:) = dim_avg_n_Wrap(QUX0, 0)
    qv(i,:,:) = dim_avg_n_Wrap(QVX0, 0)
    cc(i,:,:) = dim_avg_n_Wrap(CCX0, 0)

    qq = qu
    qq = qu+qv
    qq@long_name = "cloud liq+ice"

end do
printVarSummary(qq)

; qq = qq*1000.
; qq@units = "g/s/m"
; copy_VarCoords(qu, qq)
; printVarSummary(qq)

; printMinMax(qu, 0)
; printMinMax(qv, 0)
; printMinMax(qq, 0)
; exit

qu_diff = qu(0,:,:)
qu_diff = qu(1,:,:) - qu(0,:,:)
qv_diff = qv(0,:,:)
qv_diff = qv(1,:,:) - qv(0,:,:)
qq_diff = qq(0,:,:)
qq_diff = qq(1,:,:) - qq(0,:,:)
cc_diff = cc(0,:,:)
cc_diff = cc(1,:,:) - cc(0,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_cloudcover"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)
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
res@cnFillPalette         = "WhiteBlue"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035


res@gsnRightString = "Ctrl"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.1
res@cnMaxLevelValF = 1.
res@cnLevelSpacingF = 0.1
res@gsnLeftString = "Cloud Cover (whole column)"
plot(0) = gsn_csm_contour_map(wks, cc(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.02
res@cnMaxLevelValF = 0.2
res@cnLevelSpacingF = 0.02
res@gsnLeftString = "Cloud Liq+Ice (whole column)"
plot(1) = gsn_csm_contour_map(wks, qq(0,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.01
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnLeftString = "Cloud Ice"
plot(2) = gsn_csm_contour_map(wks, qu(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.01
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnLeftString = "Cloud Liq"
plot(3) = gsn_csm_contour_map(wks, qv(0,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

;==========================================================

res@cnFillPalette = "MPL_RdYlBu"
res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnLeftString = "Cloud Cover (whole column)"
plot(4) = gsn_csm_contour_map(wks, cc_diff, res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.02
res@cnMaxLevelValF = 0.02
res@cnLevelSpacingF = 0.002
res@gsnLeftString = "Cloud Liq+Ice (whole column)"
plot(5) = gsn_csm_contour_map(wks, qq_diff, res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.01
res@cnMaxLevelValF = 0.01
res@cnLevelSpacingF = 0.001
res@gsnLeftString = "Cloud Ice"
plot(6) = gsn_csm_contour_map(wks, qu_diff, res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.01
res@cnMaxLevelValF = 0.01
res@cnLevelSpacingF = 0.001
res@gsnLeftString = "Cloud Liq"
plot(7) = gsn_csm_contour_map(wks, qv_diff, res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7), plot_topo(7))



pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean"
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,4/), pres)

end