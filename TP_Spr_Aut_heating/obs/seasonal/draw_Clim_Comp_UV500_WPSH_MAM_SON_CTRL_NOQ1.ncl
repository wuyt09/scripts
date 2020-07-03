
begin
;=================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;=================================================================

;==================================================================
; Composite
;==================================================================
 lonstr = 50
 lonend = 160
 latstr = -10
 latend = 50
 season = (/"MAM","SON"/)
 year = ispan(6, 25, 1)
 Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
 LEV = 500
 casename = (/"CTRL","NOQ1"/)
 U_clm = new((/dimsizes(casename),12,96,144/), float)
 V_clm = U_clm
 Z_clm = U_clm


do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    U = fu->U(:,{LEV},:,:) ;{latstr:latend},{lonstr:lonend}
    V = fv->V(:,{LEV},:,:) ;{latstr:latend},{lonstr:lonend}
    U_clm(i,:,:,:) = clmMonTLL(U)
    V_clm(i,:,:,:) = clmMonTLL(V)
    printVarSummary(U_clm)

    filiz = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
    fz = addfile(diri+filiz, "r")
    Z = fz->Z3(:,{LEV},:,:) ;{latstr:latend},{lonstr:lonend}
    Z_clm(i,:,:,:) = clmMonTLL(Z)

end do
;==================================================================
;==================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/seasonal/"
figname = "draw_Clim_Comp_UV500_WPSH_MAM_SON_CTRL_NOQ1"

wks = gsn_open_wks("pdf", figpath+figname)
plotc = new(8, graphic)
plotz = plotc
plotv = plotc


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnAddCyclic = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillOn = True
res@cnFillPalette = "sunshine_9lev"
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 3.
res@cnMaxLevelValF = 5900.
res@cnMinLevelValF = 5870.
res@lbLabelBarOn = False
;------------------------------------------
resZ = True
resZ@gsnDraw  = False
resZ@gsnFrame = False

resZ@cnFillOn = False
resZ@cnLinesOn = True
resZ@cnLineLabelsOn = True
resZ@cnLineLabelFontColor = "red"
resZ@cnLineColor = "green4"
resZ@cnLineThicknessF = 4
resZ@cnInfoLabelOn = True
resZ@cnInfoLabelOrthogonalPosF = -0.12
resZ@cnInfoLabelFontHeightF = 0.02
;res@cnLineDrawOrder = "PostDraw"
resZ@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
resZ@cnLineLabelInterval = 1
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0    ; -0.35          ; move ref vector up

  vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.045            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------
vcres@vcRefMagnitudeF         = 10.0  ; define vector ref mag
resZ@lbLabelBarOn = False
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 10.
; res@cnMaxLevelValF = 12550.
; res@cnMinLevelValF = 12400.
resZ@cnLevelSelectionMode = "ExplicitLevels" ; "ManualLevels"      ;use explicit levels
;resZ@cnLevels             = (/12500,12520,12540,12560,12575/) ;(/1500/) ;(/16750/) ;
resZ@cnLevels =(/12460,12480,12500,12520,12540/)

;---------------------画到同一张图------------------------------------------------
    resP                     = True         ; modify the panel plot
    resP@gsnDraw             = True
    resP@gsnFrame            = True
    resP@gsnMaximize         = True         ; large format
    resP@lbLabelAutoStride   = True         ; auto stride on labels
    resP@lbLabelFontHeightF = 0.012
    resP@pmLabelBarOrthogonalPosF     = -0.005  ; 调色标和图的间距
    resP@gsnPanelLabelBar = True
;--------------------------------------------------------------------------------  

do i = 0, 3
  res@gsnLeftString = Month(i+2)
  res@gsnRightString = "CTRL"
  plotc(i) = gsn_csm_contour_map_ce(wks, Z_clm(0,i+2,:,:), res)
  ;plotz(i) = gsn_csm_contour(wks, Z_clm(0,i+2,:,:), resZ)
  plotv(i) = gsn_csm_vector(wks, U_clm(0,i+2,:,:), V_clm(0,i+2,:,:), vcres)
  ;overlay(plotc(i), plotz(i))
  overlay(plotc(i), plotv(i))
  res@gsnLeftString = Month(i+2)
  res@gsnRightString = "NOQ1"
  plotc(i+4) = gsn_csm_contour_map_ce(wks, Z_clm(1,i+2,:,:), res)
  ;plotz(i+4) = gsn_csm_contour(wks, Z_clm(1,i+2,:,:), resZ) ;
  plotv(i+4) = gsn_csm_vector(wks, U_clm(1,i+2,:,:), V_clm(1,i+2,:,:), vcres)
  ;overlay(plotc(i+4), plotz(i+4))
  overlay(plotc(i+4), plotv(i+4))
end do

resP@gsnPanelMainString = "Spring"
;gsn_panel(wks, plotc(4:7), (/1,4/), resP)
gsn_panel(wks, plotc, (/2,4/), resP)


do i = 0, 3
  res@gsnLeftString = Month(i+7)
  res@gsnRightString = "CTRL"
  plotc(i) = gsn_csm_contour_map_ce(wks, Z_clm(0,i+7,:,:), res)
  ;plotz(i) = gsn_csm_contour(wks, Z_clm(0,i+7,:,:), resZ)
  plotv(i) = gsn_csm_vector(wks, U_clm(0,i+7,:,:), V_clm(0,i+7,:,:), vcres)
  ;overlay(plotc(i), plotz(i))
  overlay(plotc(i), plotv(i))
  res@gsnLeftString = Month(i+7)
  res@gsnRightString = "NOQ1"
  plotc(i+4) = gsn_csm_contour_map_ce(wks, Z_clm(1,i+7,:,:), res)
  ;plotz(i+4) = gsn_csm_contour(wks, Z_clm(1,i+7,:,:), resZ) ;
  plotv(i+4) = gsn_csm_vector(wks, U_clm(1,i+7,:,:), V_clm(1,i+7,:,:), vcres)
  ;overlay(plotc(i+4), plotz(i+4))
  overlay(plotc(i+4), plotv(i+4))
end do
resP@gsnPanelMainString = "Autumn"
gsn_panel(wks, plotc, (/2,4/), resP)
end