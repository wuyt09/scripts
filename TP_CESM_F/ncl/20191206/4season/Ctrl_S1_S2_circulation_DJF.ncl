
begin

lonstr = 40
lonend = 240
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)
season = "DJF"

LEV = (/850,200/)
u = new((/ncase,2,96,144/), float)
v = u
precp = new((/ncase,96,144/), float)
z = precp

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
    z1 = fz->Z3(:,{500},:,:)
    z_zonal = dim_avg_n_Wrap(z1, 1)
    z1 = z1 - conform_dims(dimsizes(z1), z_zonal, (/0,2/))

    fc = addfile(diri+filic, "r")
    fl = addfile(diri+filil, "r")
    c1 = fc->PRECC
    l1 = fl->PRECL
    PRECP = c1
    PRECP = (c1+l1)*1000*86400
    
    U1 = month_to_season(u1, season)
    u(i,:,:,:) = dim_avg_n_Wrap(U1, 0)
    V1 = month_to_season(v1, season)
    v(i,:,:,:) = dim_avg_n_Wrap(V1, 0)
    PRECP1 = month_to_season(PRECP, season)
    precp(i,:,:) = dim_avg_n_Wrap(PRECP1, 0)
    Z1 = month_to_season(z1, season)
    z(i,:,:) = dim_avg_n_Wrap(Z1, 0)
end do

ux = u(0:1,:,:,:)
vx = v(0:1,:,:,:)
px = precp(0:1,:,:)
zx = z(0:1,:,:)

ux(0,:,:,:) = u(0,:,:,:) - u(1,:,:,:)
ux(1,:,:,:) = u(1,:,:,:) - u(2,:,:,:)
vx(0,:,:,:) = v(0,:,:,:) - v(1,:,:,:)
vx(1,:,:,:) = v(1,:,:,:) - v(2,:,:,:)
px(0,:,:) = precp(0,:,:) - precp(1,:,:)
px(1,:,:) = precp(1,:,:) - precp(2,:,:)
zx(0,:,:) = z(0,:,:) - z(1,:,:)
zx(1,:,:) = z(1,:,:) - z(2,:,:)
printVarSummary(ux)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191206/4season/"


  figname = "Ctrl_S1_S2_circulation_"+season

  wks = gsn_open_wks("pdf", figpath+figname)
  nplot = 6
  plot = new(nplot, graphic)
  plot_topo = plot
  plotv = plot


   res_topo            = True
   res_topo@cnFillOn   = False
   res_topo@cnLinesOn  = True
   res_topo@gsnFrame   = False
   res_topo@gsnDraw    = False
   res_topo@cnLineLabelsOn         = False
   res_topo@cnLineThicknessF       = 5
   res_topo@cnLineColor            = "grey55"
   res_topo@cnLineDashPattern      = 16
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
  res@pmLabelBarOrthogonalPosF = 0.25
  ;------------------------------------------

    vcres                         = True             ; plot mods desired
    vcres@gsnFrame                = False            ; so we can draw time stamp
    vcres@gsnDraw                 = False
    vcres@vcRefAnnoOrthogonalPosF =  -1.0             ; move ref vector up

    ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
    ;vcres@vcMinDistanceF          = 0.017
    vcres@vcGlyphStyle            = "LineArrow"
        
    vcres@vcLineArrowThicknessF   = 3.0
    ;vecres@vcLineArrowHeadMinSizeF = 0.008
    ;vecres@vcLineArrowHeadMaxSizeF = 0.3

    vcres@vcRefAnnoString2On = False
    vcres@vcRefAnnoFontHeightF = 0.02
  ;------------------------------------------
  
    ;res@lbLabelBarOn = False
    res@cnLevelSelectionMode = "ManualLevels"
    res@gsnLeftString = "UV"+LEV(0)+"_Precip"
 
  res@gsnRightString = "Ctrl"
  vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag
  vcres@vcMinDistanceF          = 0.04
  vcres@vcRefLengthF            = 0.05            ; define length of vec ref

  res@cnLevelSpacingF = 2.
  res@cnMaxLevelValF = 20.
  res@cnMinLevelValF = 0.
  res@cnFillPalette = "CBR_wet"
  plot(0) = gsn_csm_contour_map(wks, precp(0,:,:), res)
  plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
  plotv(0) = gsn_csm_vector(wks, u(0,{LEV(0)},:,:), v(0,{LEV(0)},:,:), vcres)
  overlay(plot(0), plot_topo(0))
  overlay(plot(0), plotv(0))

  res@gsnRightString = "Thermal"
  vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
  res@cnLevelSpacingF = 0.5
  res@cnMaxLevelValF = 5.
  res@cnMinLevelValF = -5.
  res@cnFillPalette = "CBR_drywet"
  plot(1) = gsn_csm_contour_map(wks, px(0,:,:), res)
  plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
  plotv(1) = gsn_csm_vector(wks, ux(0,{LEV(0)},:,:), vx(0,{LEV(0)},:,:), vcres)
  overlay(plot(1), plot_topo(1))
  overlay(plot(1), plotv(1))

  res@gsnRightString = "Dynamical"
  vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
  res@cnLevelSpacingF = 0.5
  res@cnMaxLevelValF = 5.
  res@cnMinLevelValF = -5.
  res@cnFillPalette = "CBR_drywet"
  plot(2) = gsn_csm_contour_map(wks, px(1,:,:), res)
  plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
  plotv(2) = gsn_csm_vector(wks, ux(1,{LEV(0)},:,:), vx(1,{LEV(0)},:,:), vcres)
  overlay(plot(2), plot_topo(2))
  overlay(plot(2), plotv(2))
;-------------------------------------------------------------------
  
  res@gsnLeftString = "UV"+LEV(1)+"_H500"

  res@gsnRightString = "Ctrl"
  vcres@vcRefMagnitudeF         = 20.0  ; define vector ref mag
  res@cnLevelSpacingF = 5.
  res@cnMaxLevelValF = 420.
  res@cnMinLevelValF = 380.
  res@cnFillPalette = "WhiteYellowOrangeRed"
  plot(3) = gsn_csm_contour_map(wks, z(0,:,:), res)
  plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
  plotv(3) = gsn_csm_vector(wks, u(0,{LEV(1)},:,:), v(0,{LEV(1)},:,:), vcres)
  overlay(plot(3), plot_topo(3))
  overlay(plot(3), plotv(3))

  res@gsnRightString = "Thermal"
  vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag
  res@cnLevelSpacingF = 5.
  res@cnMaxLevelValF = 50.
  res@cnMinLevelValF = -50.
  res@cnFillPalette = "BlueWhiteOrangeRed"
  plot(4) = gsn_csm_contour_map(wks, zx(0,:,:), res)
  plot_topo(4) = gsn_csm_contour(wks, tpmask, res_topo)
  plotv(4) = gsn_csm_vector(wks, ux(0,{LEV(1)},:,:), vx(0,{LEV(1)},:,:), vcres)
  overlay(plot(4), plot_topo(4))
  overlay(plot(4), plotv(4))

  res@gsnRightString = "Dynamical"
  vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag
  res@cnLevelSpacingF = 5.
  res@cnMaxLevelValF = 50.
  res@cnMinLevelValF = -50.
  res@cnFillPalette = "BlueWhiteOrangeRed"
  plot(5) = gsn_csm_contour_map(wks, zx(1,:,:), res)
  plot_topo(5) = gsn_csm_contour(wks, tpmask, res_topo)
  plotv(5) = gsn_csm_vector(wks, ux(1,{LEV(1)},:,:), vx(1,{LEV(1)},:,:), vcres)
  overlay(plot(5), plot_topo(5))
  overlay(plot(5), plotv(5))

  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelMainString = season
  ;pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01

  gsn_panel(wks, plot, (/2,3/), pres)



end