
begin

LEV = 500
nlat = 96
nlon = 144

lonstr = 50
lonend = 150
latstr = -10
latend = 50

casename = (/"CTRL","NOQ1"/)
ncase = dimsizes(casename)
LEV = 500
t = new((/ncase,12,96,144/), float)
u = t
v = t

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  tname = "F2000_"+casename(i)+".plevel.T.mon.0601-3012.nc"
  uname = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
  vname = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
  pname = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"
  ft = addfile(diri+tname, "r")
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  fp = addfile(diri+pname, "r")
  U = fu->U(:,{LEV},:,:)
  V = fv->V(:,{LEV},:,:)
  u(i,:,:,:) = clmMonTLL(U)
  v(i,:,:,:) = clmMonTLL(V)

  T = ft->T
  PS = fp->PS
  plev = T&lev*100
  plev!0 = "lev"
  plev&lev = T&lev
  psfc = PS
  ptop = min(plev)
  dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
  TdP = T*dP
  copy_VarCoords(T, TdP)
  tx = dim_sum_n_Wrap(TdP(:,{500:200},:,:), 1)/dim_sum_n_Wrap(dP(:,{500:200},:,:), 1)
  copy_VarCoords(T(:,0,:,:), tx)
  t(i,:,:,:) = clmMonTLL(tx)
end do
printVarSummary(t)

diff_t = t(0,:,:,:)
diff_t = t(0,:,:,:) - t(1,:,:,:)
diff_u = u(0,:,:,:)
diff_u = u(0,:,:,:) - u(1,:,:,:)
diff_v = v(0,:,:,:)
diff_v = v(0,:,:,:) - v(1,:,:,:)
printVarSummary(diff_t)

;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "UV500_TropoT_Apr_Sep"

wks = gsn_open_wks("ps", figpath+figname)

nplot = 18
plot = new(nplot, graphic)
plotv = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0  ;  -0.35           ; move ref vector up

  vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.025            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------
vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag
;res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 250.
res@cnMinLevelValF = 240.
;res@cnFillPalette = "BlueWhiteOrangeRed"
res@gsnRightString = "K"


do i = 0,5
  res@gsnRightString = Month(i+3)
  
  vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 1.
  res@cnMaxLevelValF = 250.
  res@cnMinLevelValF = 240.
  res@cnFillPalette = "ncl_default"

  res@gsnLeftString = "Ctrl"
  plot(i) = gsn_csm_contour_map(wks, t(0,i+3,:,:), res)
  plotv(i) = gsn_csm_vector(wks, u(0,i+3,::2,::2), v(0,i+3,::2,::2), vcres)
  overlay(plot(i), plotv(i))

  res@gsnLeftString = "Sen"
  plot(i+6) = gsn_csm_contour_map(wks, t(1,i+3,:,:), res)
  plotv(i+6) = gsn_csm_vector(wks, u(1,i+3,::2,::2), v(1,i+3,::2,::2), vcres)
  overlay(plot(i+6), plotv(i+6))

  vcres@vcRefMagnitudeF         = 2.0  ; define vector ref mag
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 2.
  res@cnMaxLevelValF = 10.
  res@cnMinLevelValF = -10.
  res@cnFillPalette = "CBR_coldhot"
  res@gsnLeftString = "diff"
  plot(i+12) = gsn_csm_contour_map(wks, diff_t(i+3,:,:), res)
  plotv(i+12) = gsn_csm_vector(wks, diff_u(i+3,::2,::2), diff_v(i+3,::2,::2), vcres)
  overlay(plot(i+12), plotv(i+12))
end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "UV500_TropoT (500-200hPa)"
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,6/), pres)


end