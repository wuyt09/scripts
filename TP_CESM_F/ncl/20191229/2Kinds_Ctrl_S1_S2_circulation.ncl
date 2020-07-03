
begin

lonstr = 0
lonend = 360
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)
season = (/"MAM","JJA","SON","DJF"/)

LEV = (/850,200/)
u = new((/ncase,12,2,96,144/), float)
v = u
precp = new((/ncase,12,96,144/), float)
z = precp
uvar = u
vvar = v
pvar = precp
zvar = z

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    filiz = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
    filic = "F2000_"+casename(i)+".PRECC.mon.0601-3012.nc"
    filil = "F2000_"+casename(i)+".PRECL.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    u1 = fu->U(:,{LEV},:,:)
    u_zonal = dim_avg_n_Wrap(u1, 2)
    u1 = u1 - conform_dims(dimsizes(u1), u_zonal, (/0,1,3/))
    v1 = fv->V(:,{LEV},:,:)
    v_zonal = dim_avg_n_Wrap(v1, 2)
    v1 = v1 - conform_dims(dimsizes(v1), v_zonal, (/0,1,3/))
    printVarSummary(u1)
    
    fz = addfile(diri+filiz, "r")
    z1 = fz->Z3(:,{300},:,:)
    z_zonal = dim_avg_n_Wrap(z1, 1)
    z1 = z1 - conform_dims(dimsizes(z1), z_zonal, (/0,2/))

    fc = addfile(diri+filic, "r")
    fl = addfile(diri+filil, "r")
    c1 = fc->PRECC
    l1 = fl->PRECL
    PRECP = c1
    PRECP = (c1+l1)*1000*86400
    
    u(i,:,:,:,:) = clmMonTLLL(u1)
    v(i,:,:,:,:) = clmMonTLLL(v1)
    precp(i,:,:,:) = clmMonTLL(PRECP)
    z(i,:,:,:) = clmMonTLL(z1)
end do
printVarSummary(u)

ux = u(0:1,:,:,:,:)
vx = v(0:1,:,:,:,:)
px = precp(0:1,:,:,:)
zx = z(0:1,:,:,:)
probu = ux
probv = vx
probp = px
probz = zx
s1 = 25
s2 = 25

ux(0,:,:,:,:) = u(0,:,:,:,:) - u(1,:,:,:,:)
ux(1,:,:,:,:) = u(1,:,:,:,:) - u(2,:,:,:,:)
; probu(0,:,:,:,:) = ttest(u(0,:,:,:,:), uvar(0,:,:,:,:), s1, u(1,:,:,:,:), uvar(1,:,:,:,:), s2, True, False)
; probu(1,:,:,:,:) = ttest(u(1,:,:,:,:), uvar(1,:,:,:,:), s1, u(2,:,:,:,:), uvar(2,:,:,:,:), s2, True, False)

vx(0,:,:,:,:) = v(0,:,:,:,:) - v(1,:,:,:,:)
vx(1,:,:,:,:) = v(1,:,:,:,:) - v(2,:,:,:,:)
; probv(0,:,:,:,:) = ttest(v(0,:,:,:,:), vvar(0,:,:,:,:), s1, v(1,:,:,:,:), vvar(1,:,:,:,:), s2, True, False)
; probv(1,:,:,:,:) = ttest(v(1,:,:,:,:), vvar(1,:,:,:,:), s1, v(2,:,:,:,:), vvar(2,:,:,:,:), s2, True, False)

px(0,:,:,:) = precp(0,:,:,:) - precp(1,:,:,:)
px(1,:,:,:) = precp(1,:,:,:) - precp(2,:,:,:)
; probp(0,:,:,:) = ttest(precp(0,:,:,:), pvar(0,:,:,:), s1, precp(1,:,:,:), pvar(1,:,:,:), s2, True, False)
; probp(1,:,:,:) = ttest(precp(1,:,:,:), pvar(1,:,:,:), s1, precp(2,:,:,:), pvar(2,:,:,:), s2, True, False)

zx(0,:,:,:) = z(0,:,:,:) - z(1,:,:,:)
zx(1,:,:,:) = z(1,:,:,:) - z(2,:,:,:)
; probz(0,:,:,:) = ttest(z(0,:,:,:), zvar(0,:,:,:), s1, z(1,:,:,:), zvar(1,:,:,:), s2, True, False)
; probz(1,:,:,:) = ttest(z(1,:,:,:), zvar(1,:,:,:), s1, z(2,:,:,:), zvar(2,:,:,:), s2, True, False)

printVarSummary(ux)

siglvl = 0.05
; ux = mask(ux, probu.le.siglvl .and. probv.le.siglvl, True)
; vx = mask(vx, probu.le.siglvl .and. probv.le.siglvl, True)
; px = mask(px, probp.le.siglvl, True)
; zx = mask(zx, probz.le.siglvl, True)
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191229/"

figname = "2Kind_Ctrl_S1_S2_circulation"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12*2
plot = new(nplot, graphic)
plot_topo = plot
plotv = plot
plot_prob = plot


   res_topo            = True
   res_topo@cnFillOn   = False
   res_topo@cnLinesOn  = True
   res_topo@gsnFrame   = False
   res_topo@gsnDraw    = False
   res_topo@cnLineLabelsOn         = False
   res_topo@cnLineThicknessF       = 3
   res_topo@cnLineColor            = "blue"
   res_topo@cnInfoLabelOn          = False
   res_topo@cnLevelSelectionMode   ="ExplicitLevels"
   res_topo@cnLevels               = (/0.9/)
   res_topo@gsnLeftString = ""
   res_topo@gsnRightString = ""


  res = True
  res@gsnDraw  = False
  res@gsnFrame = False

  res@mpOutlineOn = True
  res@mpFillOn = True
  res@mpGeophysicalLineThicknessF = 2
  res@mpGeophysicalLineColor = "grey45" ;"grey25"
  res@mpFillDrawOrder = "PreDraw"

  res@mpMinLatF = latstr
  res@mpMaxLatF = latend
  res@mpMinLonF = lonstr
  res@mpMaxLonF = lonend
  res@mpCenterLonF = (lonstr+lonend)/2

  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False

  res@gsnStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025
  res@pmLabelBarOrthogonalPosF = 0.15
  ;------------------------------------------

    vcres                         = True             ; plot mods desired
    vcres@gsnFrame                = False            ; so we can draw time stamp
    vcres@gsnDraw                 = False
    vcres@vcRefAnnoOrthogonalPosF =  -0.7             ; move ref vector up

    ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
    ;vcres@vcMinDistanceF          = 0.017
    vcres@vcGlyphStyle            = "LineArrow"
    vcres@vcMinDistanceF          = 0.025
    vcres@vcRefLengthF            = 0.04            ; define length of vec ref
    
    vcres@vcLineArrowThicknessF   = 3.0
    ;vecres@vcLineArrowHeadMinSizeF = 0.008
    ;vecres@vcLineArrowHeadMaxSizeF = 0.3

    vcres@vcRefAnnoString2On = False
    vcres@vcRefAnnoFontHeightF = 0.02
  ;------------------------------------------
  
    ;------------------------
    ; draw prob
    ;------------------------
    res3  = True
    res3@gsnDraw   = False
    res3@gsnFrame  = False

    res3@cnFillOn  = False 
    res3@cnLinesOn = False
    res3@cnLineLabelsOn = False
    res3@cnInfoLabelOn  = False
    res3@cnFillColor = "blue"

    res3@gsnRightString = ""

    res3@cnLevelSelectionMode = "ManualLevels"
    res3@cnMinLevelValF = 0.00
    res3@cnMaxLevelValF = 0.10
    res3@cnLevelSpacingF = 0.01

    opt  = True
    opt@gsnShadeFillType = "pattern"
    opt@gsnShadeLow = 17
    opt@gsnShadeFillScaleF = 0.8
    opt@gsnShadeFillDotSizeF = 0.003
  ;------------------------------------------

    res@lbLabelBarOn = False
    res@cnLevelSelectionMode = "ManualLevels"
    
  vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
  res@cnLevelSpacingF = 0.5
  res@cnMaxLevelValF = 5.
  res@cnMinLevelValF = -5.
  res@cnFillPalette = "CBR_drywet"

  res@gsnLeftString = "Thermal"
  do i = 0, dimsizes(month)-1
    res@gsnRightString = month(i)
    plot(2*i) = gsn_csm_contour_map(wks, px(0,i,:,:), res)
    plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(2*i) = gsn_csm_vector(wks, ux(0,i,{LEV(0)},:,:), vx(0,i,{LEV(0)},:,:), vcres)
    overlay(plot(2*i), plot_topo(2*i))
    overlay(plot(2*i), plotv(2*i))

    plot_prob(2*i) = gsn_csm_contour(wks, probp(0,i,:,:), res3)
    plot_prob(2*i) = gsn_contour_shade(plot_prob(2*i), siglvl, -999., opt)
    ;overlay(plot(2*i), plot_prob(2*i))
  end do
  
  res@gsnLeftString = "Dynamical"
  do i = 0, dimsizes(month)-1
    res@gsnRightString = month(i)
    plot(2*i+1) = gsn_csm_contour_map(wks, px(1,i,:,:), res)
    plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(2*i+1) = gsn_csm_vector(wks, ux(1,i,{LEV(0)},:,:), vx(1,i,{LEV(0)},:,:), vcres)
    overlay(plot(2*i+1), plot_topo(2*i+1))
    overlay(plot(2*i+1), plotv(2*i+1))

    plot_prob(2*i+1) = gsn_csm_contour(wks, probp(1,i,:,:), res3)
    plot_prob(2*i+1) = gsn_contour_shade(plot_prob(2*i+1), siglvl, -999., opt)
    ;overlay(plot(2*i+1), plot_prob(2*i+1))
  end do
  
  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelMainString = "UV"+LEV(0)+"_Precip"
  pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01
  pres@lbLabelFontHeightF = 0.01
  gsn_panel(wks, plot, (/12,2/), pres)

;-------------------------------------------------------------------

  res@gsnLeftString = "Thermal"
  vcres@vcRefMagnitudeF         = 10.0  ; define vector ref mag
  res@cnLevelSpacingF = 10.
  res@cnMaxLevelValF = 100.
  res@cnMinLevelValF = -100.
  res@cnFillPalette = "BlueWhiteOrangeRed"
  do i = 0, dimsizes(month)-1
    res@gsnRightString = month(i)
    plot(2*i) = gsn_csm_contour_map(wks, zx(0,i,:,:), res)
    plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(2*i) = gsn_csm_vector(wks, ux(0,i,{LEV(1)},:,:), vx(0,i,{LEV(1)},:,:), vcres)
    overlay(plot(2*i), plot_topo(2*i))
    overlay(plot(2*i), plotv(2*i))

    plot_prob(2*i) = gsn_csm_contour(wks, probz(0,i,:,:), res3)
    plot_prob(2*i) = gsn_contour_shade(plot_prob(2*i), siglvl, -999., opt)
    ;overlay(plot(2*i), plot_prob(2*i))
  end do

  res@gsnLeftString = "Dynamical"
  do i = 0, dimsizes(month)-1
    res@gsnRightString = month(i)
    plot(2*i+1) = gsn_csm_contour_map(wks, zx(1,i,:,:), res)
    plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(2*i+1) = gsn_csm_vector(wks, ux(1,i,{LEV(1)},:,:), vx(1,i,{LEV(1)},:,:), vcres)
    overlay(plot(2*i+1), plot_topo(2*i+1))
    overlay(plot(2*i+1), plotv(2*i+1))

    plot_prob(2*i+1) = gsn_csm_contour(wks, probz(1,i,:,:), res3)
    plot_prob(2*i+1) = gsn_contour_shade(plot_prob(2*i+1), siglvl, -999., opt)
    ;overlay(plot(2*i+1), plot_prob(2*i+1))
  end do
  
  pres@gsnPanelMainString = "UV"+LEV(1)+"_H500"
  pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01
  pres@lbLabelFontHeightF = 0.01
  gsn_panel(wks, plot, (/12,2/), pres)

;==================================================================
;==================================================================
ppx = px
uux = ux
vvx = vx
zzx = zx
do i = 0, 8
    ppx(:,i,:,:) = px(:,i+3,:,:)
    uux(:,i,:,:,:) = ux(:,i+3,:,:,:)
    vvx(:,i,:,:,:) = vx(:,i+3,:,:,:)
    zzx(:,i,:,:) = zx(:,i+3,:,:)
end do
do i = 9, 11
    ppx(:,i,:,:) = px(:,i-9,:,:)
    uux(:,i,:,:,:) = ux(:,i-9,:,:,:)
    vvx(:,i,:,:,:) = vx(:,i-9,:,:,:)
    zzx(:,i,:,:) = zx(:,i-9,:,:)
end do


p_warm = ppx(:,0,:,:)
p_warm = dim_avg_n_Wrap(ppx(:,0:5,:,:), 1)
p_cold = ppx(:,0,:,:)
p_cold = dim_avg_n_Wrap(ppx(:,6:11,:,:), 1)

z_warm = zzx(:,0,:,:)
z_warm = dim_avg_n_Wrap(zzx(:,0:5,:,:), 1)
z_cold = zzx(:,0,:,:)
z_cold = dim_avg_n_Wrap(zzx(:,6:11,:,:), 1)

u_warm = uux(:,0,:,:,:)
u_warm = dim_avg_n_Wrap(uux(:,0:5,:,:,:), 1)
u_cold = uux(:,0,:,:,:)
u_cold = dim_avg_n_Wrap(uux(:,6:11,:,:,:), 1)

v_warm = vvx(:,0,:,:,:)
v_warm = dim_avg_n_Wrap(vvx(:,0:5,:,:,:), 1)
v_cold = vvx(:,0,:,:,:)
v_cold = dim_avg_n_Wrap(vvx(:,6:11,:,:,:), 1)


res@gsnLeftString = "Thermal"
res@gsnRightString = "Warm season"
res@cnLevelSelectionMode = "ManualLevels"
vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.
res@cnFillPalette = "CBR_drywet"
plot(0) = gsn_csm_contour_map(wks, p_warm(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(0) = gsn_csm_vector(wks, u_warm(0,{LEV(0)},:,:), v_warm(0,{LEV(0)},:,:), vcres)
overlay(plot(0), plot_topo(0))
overlay(plot(0), plotv(0))

res@gsnLeftString = "Dynamical"
plot(1) = gsn_csm_contour_map(wks, p_warm(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(1) = gsn_csm_vector(wks, u_warm(1,{LEV(0)},:,:), v_warm(1,{LEV(0)},:,:), vcres)
overlay(plot(1), plot_topo(1))
overlay(plot(1), plotv(1))

res@gsnLeftString = "Thermal"
res@gsnRightString = "Cold season"
plot(2) = gsn_csm_contour_map(wks, p_cold(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2) = gsn_csm_vector(wks, u_cold(0,{LEV(0)},:,:), v_cold(0,{LEV(0)},:,:), vcres)
overlay(plot(2), plot_topo(2))
overlay(plot(2), plotv(2))

res@gsnLeftString = "Dynamical"
plot(3) = gsn_csm_contour_map(wks, p_cold(1,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(3) = gsn_csm_vector(wks, u_cold(1,{LEV(0)},:,:), v_cold(1,{LEV(0)},:,:), vcres)
overlay(plot(3), plot_topo(3))
overlay(plot(3), plotv(3))

pres@gsnPanelMainString = "UV"+LEV(0)+"_Precip"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.01
gsn_panel(wks, plot, (/2,2/), pres)






res@gsnLeftString = "Thermal"
res@gsnRightString = "Warm season"
res@cnLevelSelectionMode = "ManualLevels"
vcres@vcRefMagnitudeF         = 10.0  ; define vector ref mag
res@cnLevelSpacingF = 10.
res@cnMaxLevelValF = 100.
res@cnMinLevelValF = -100.
res@cnFillPalette = "BlueWhiteOrangeRed"
plot(0) = gsn_csm_contour_map(wks, z_warm(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(0) = gsn_csm_vector(wks, u_warm(0,{LEV(1)},:,:), v_warm(0,{LEV(1)},:,:), vcres)
overlay(plot(0), plot_topo(0))
overlay(plot(0), plotv(0))

res@gsnLeftString = "Dynamical"
plot(1) = gsn_csm_contour_map(wks, z_warm(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(1) = gsn_csm_vector(wks, u_warm(1,{LEV(1)},:,:), v_warm(1,{LEV(1)},:,:), vcres)
overlay(plot(1), plot_topo(1))
overlay(plot(1), plotv(1))

res@gsnLeftString = "Thermal"
res@gsnRightString = "Cold season"
plot(2) = gsn_csm_contour_map(wks, z_cold(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2) = gsn_csm_vector(wks, u_cold(0,{LEV(1)},:,:), v_cold(0,{LEV(1)},:,:), vcres)
overlay(plot(2), plot_topo(2))
overlay(plot(2), plotv(2))

res@gsnLeftString = "Dynamical"
plot(3) = gsn_csm_contour_map(wks, z_cold(1,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(3) = gsn_csm_vector(wks, u_cold(1,{LEV(1)},:,:), v_cold(1,{LEV(1)},:,:), vcres)
overlay(plot(3), plot_topo(3))
overlay(plot(3), plotv(3))

pres@gsnPanelMainString = "UV"+LEV(1)+"_H300"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.01
gsn_panel(wks, plot, (/2,2/), pres)

end