
begin

lonstr = 50
lonend = 130
latstr = 0
latend = 50

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

LEV = (/925/)
u = new((/ncase,12,96,144/), float)
v = u
precp = new((/ncase,12,96,144/), float)

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    filic = "F2000_"+casename(i)+".PRECC.mon.0601-3012.nc"
    filil = "F2000_"+casename(i)+".PRECL.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    u1 = fu->U(:,{LEV},:,:)
    u_zonal = dim_avg_n_Wrap(u1, 1)
    u1 = u1 - conform_dims(dimsizes(u1), u_zonal, (/0,2/))
    v1 = fv->V(:,{LEV},:,:)
    v_zonal = dim_avg_n_Wrap(v1, 1)
    v1 = v1 - conform_dims(dimsizes(v1), v_zonal, (/0,2/))
    printVarSummary(u1)

    fc = addfile(diri+filic, "r")
    fl = addfile(diri+filil, "r")
    c1 = fc->PRECC
    l1 = fl->PRECL
    PRECP = c1
    PRECP = (c1+l1)*1000*86400
    
    u(i,:,:,:) = clmMonTLL(u1)
    v(i,:,:,:) = clmMonTLL(v1)
    precp(i,:,:,:) = clmMonTLL(PRECP)
end do

;P_band = dim_avg_n_Wrap(precp(:,:,:,{110:120}), 2)

ux = u(0:1,:,:,:)
vx = v(0:1,:,:,:)
px = precp(0:1,:,:,:)

ux(0,:,:,:) = u(0,:,:,:) - u(1,:,:,:)
ux(1,:,:,:) = u(1,:,:,:) - u(2,:,:,:)
vx(0,:,:,:) = v(0,:,:,:) - v(1,:,:,:)
vx(1,:,:,:) = v(1,:,:,:) - v(2,:,:,:)
px(0,:,:,:) = precp(0,:,:,:) - precp(1,:,:,:)
px(1,:,:,:) = precp(1,:,:,:) - precp(2,:,:,:)
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


  figname = "Ctrl_S1_S2_Rainband"

  wks = gsn_open_wks("pdf", figpath+figname)
  nplot = 6*ncase
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
  res@pmLabelBarOrthogonalPosF = 0.2
  ;------------------------------------------

    vcres                         = True             ; plot mods desired
    vcres@gsnFrame                = False            ; so we can draw time stamp
    vcres@gsnDraw                 = False
    vcres@vcRefAnnoOrthogonalPosF =  -0.9             ; move ref vector up

    ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
    ;vcres@vcMinDistanceF          = 0.017
    vcres@vcGlyphStyle            = "LineArrow"
        
    vcres@vcLineArrowThicknessF   = 3.0
    ;vecres@vcLineArrowHeadMinSizeF = 0.008
    ;vecres@vcLineArrowHeadMaxSizeF = 0.3

    vcres@vcRefAnnoString2On = False
    vcres@vcRefAnnoFontHeightF = 0.02
  ;------------------------------------------
      
  res@gsnRightString = "mm/day"
  vcres@vcRefMagnitudeF         = 8.0  ; define vector ref mag
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.05            ; define length of vec ref

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 2.
  res@cnMaxLevelValF = 20.
  res@cnMinLevelValF = 0.
  res@cnFillPalette = "CBR_wet"

  do i = 0, 5
    res@gsnLeftString = "Ctrl"
    res@gsnRightString = month(i+3)
      plot(i) = gsn_csm_contour_map(wks, precp(0,i+3,:,:), res)
      plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
      plotv(i) = gsn_csm_vector(wks, u(0,i+3,:,:), v(0,i+3,:,:), vcres)
      overlay(plot(i), plot_topo(i))
      overlay(plot(i), plotv(i))
  end do

  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelMainString = "UV_"+LEV+"_Precip"
  ;pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01

  gsn_panel(wks, plot, (/6,1/), pres)



end