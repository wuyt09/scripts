
begin

LEV = 200
nlat = 96
nlon = 144

; lonstr = 50
; lonend = 180
; latstr = -10
; latend = 40
lonstr = 50
lonend = 150
latstr = 0 ;-10
latend = 50

casename = (/"CTRL","NOQ1"/)
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)

LEV = 500
z = new((/ncase,nyear,nmonth,96,144/), float)
u = z
v = z

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
  uname = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
  vname = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  fu = addfile(diri+uname, "r")
  fv = addfile(diri+vname, "r")
  Z = fz->Z3(:,{LEV},:,:)
  U = fu->U(:,{LEV},:,:)
  V = fv->V(:,{LEV},:,:)
  z(i,:,:,:,:) = reshape(Z, (/nyear,nmonth,96,144/))
  copy_VarCoords(Z(0,:,:), z(i,0,0,:,:))
  u(i,:,:,:,:) = reshape(U, (/nyear,nmonth,96,144/))
  v(i,:,:,:,:) = reshape(V, (/nyear,nmonth,96,144/))
end do
z!0 = "case"
z&case = casename
z!1 = "year"
z&year = year
z!2 = "month"
z&month = month
printVarSummary(z)
copy_VarCoords(z, u)
copy_VarCoords(z, v)

;------------ttest--------------------------------------
u_var = dim_variance_n_Wrap(u, 1)
v_var = dim_variance_n_Wrap(v, 1)
z_var = dim_variance_n_Wrap(z, 1)

u_ave = dim_avg_n_Wrap(u, 1)
v_ave = dim_avg_n_Wrap(v, 1)
z_ave = dim_avg_n_Wrap(z, 1)


probu = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(1,:,:,:), u_var(1,:,:,:), nyear, False, False)
probv = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(1,:,:,:), v_var(1,:,:,:), nyear, False, False)
probz = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(1,:,:,:), z_var(1,:,:,:), nyear, False, False)

diff_z = z_ave(0,:,:,:)
diff_z = z_ave(0,:,:,:) - z_ave(1,:,:,:)
diff_u = u_ave(0,:,:,:)
diff_u = u_ave(0,:,:,:) - u_ave(1,:,:,:)
diff_v = v_ave(0,:,:,:)
diff_v = v_ave(0,:,:,:) - v_ave(1,:,:,:)

sig = 0.05
diff_u = mask(diff_u, (probu.lt.sig).or.(probv.lt.sig), True)
diff_v = mask(diff_v, (probu.lt.sig).or.(probv.lt.sig), True)
diff_z = mask(diff_z, probz.lt.sig, True)
printVarSummary(diff_z)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "UV500_H500_Apr_Jul_"+casename(0)+"_"+casename(1)

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 6
plot = new(nplot, graphic)
plotv = plot
plotc = plot
plot_topo = plot
;----------------------------------------------------
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 16
 res_topo@cnLineColor            = "blue"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@lbLabelBarOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;=====================================================
    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "green4"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = True
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1

;=====================================================

;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0  ;   -0.35          ; move ref vector up

  vcres@vcGlyphStyle            = "CurlyVector" ; "LineArrow"    ; turn on curly vectors
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
;res@lbLabelBarOn = False

;res@gsnRightString = "mm/day"

  vcres@vcRefMagnitudeF         = 2.0  ; define vector ref mag
  res@cnFillPalette = "WhiteYellowOrangeRed"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 5.
  res@cnMaxLevelValF = 5900.
  res@cnMinLevelValF = 5880.
  ; cnres@cnLevelSelectionMode = "ExplicitLevels"
  ; cnres@cnLevels = (/12460,12470,12480,12490,12500,15210,12520,12540/);12520,12560
  cnres@cnLevelSelectionMode = "ManualLevels"
  cnres@cnLevelSpacingF = 10.
  cnres@cnMaxLevelValF = 5900.
  cnres@cnMinLevelValF = 5880.

  do i = 0,5
  res@gsnLeftString = Month(i+3)
  
  plot(i) = gsn_csm_contour_map(wks, z_ave(1,i+3,:,:), res)
  plotc(i) = gsn_csm_contour(wks, z_ave(0,i+3,:,:), cnres)
  overlay(plot(i), plotc(i))
  plotv(i) = gsn_csm_vector(wks, diff_u(i+3,:,:), diff_v(i+3,:,:), vcres)
  overlay(plot(i), plotv(i))
  plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
  overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = "UV200_H200"
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/2,3/), pres)


end