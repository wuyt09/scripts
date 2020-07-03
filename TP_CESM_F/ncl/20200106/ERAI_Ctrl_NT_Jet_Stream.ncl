
begin
lonstr = 40
lonend = 180
latstr = 20
latend = 60

casename = (/"CTRL","NOTP"/) ;
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

LEV = (/200/)
u = new((/ncase,12,96,144/), float)
v = u

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    u1 = fu->U(:,{LEV},:,:)
    ;u_zonal = dim_avg_n_Wrap(u1, 1)
    ;u1 = u1 - conform_dims(dimsizes(u1), u_zonal, (/0,2/))
    v1 = fv->V(:,{LEV},:,:)
    ;v_zonal = dim_avg_n_Wrap(v1, 1)
    ;v1 = v1 - conform_dims(dimsizes(v1), v_zonal, (/0,2/))
    printVarSummary(u1)
    
    u(i,:,:,:) = clmMonTLL(u1)
    v(i,:,:,:) = clmMonTLL(v1)
end do

;speed = wind_speed(u, v)
speed = u
printVarSummary(speed)

;==================================================================

yrStrt = 1979
yrLast = 2016

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"

fu = addfile(diri+filiu, "r")
TIMEu = fu->time
YYYYu = cd_calendar(TIMEu, -1)/100
iYYYYu = ind(YYYYu.ge.yrStrt .and. YYYYu.le.yrLast)
uwnd = short2flt(fu->u(iYYYYu,{LEV},:,:)) ;{latstr:latend},{lonstr:lonend}
printVarSummary(uwnd)

fv = addfile(diri+filiv, "r")
TIMEv = fv->time
YYYYv = cd_calendar(TIMEv, -1)/100
iYYYYv = ind(YYYYv.ge.yrStrt .and. YYYYv.le.yrLast)
vwnd = short2flt(fv->v(iYYYYv,{LEV},:,:)) ;{latstr:latend},{lonstr:lonend}
printVarSummary(vwnd)

;speed_era = wind_speed(uwnd, vwnd)
speed_era = uwnd
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200106/"


  figname = "ERAI_Ctrl_NT_Jet_Stream"

  wks = gsn_open_wks("pdf", figpath+figname)
  nplot = 6
  plot = new(nplot, graphic)
  plot_topo = plot


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
  res@vpHeightF = 0.4
  res@vpWidthF = 0.7

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
  res@lbLabelBarOn = False

  res@gsnStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025
  res@pmLabelBarOrthogonalPosF = 0.2
  ;------------------------------------------
      
  res@gsnRightString = "m/s"

  ; res@cnLevelSelectionMode = "ExplicitLevels"
  ; res@cnLevels = (/0.5,1,2,4,6,8,10,13/)
  ; res@cnFillColors = (/2,3,4,6,7,9,10,11,12/)
  res@cnFillPalette = "WhiteYellowOrangeRed"
  
  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01
  pres@lbLabelFontHeightF = 0.015

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 5
  res@cnMaxLevelValF = 60.
  res@cnMinLevelValF = 30.
  do i = 0, 5
    res@gsnLeftString = month(i+3)
    plot(i) = gsn_csm_contour_map(wks, speed_era(i+3,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i), plot_topo(i))
  end do

  pres@gsnPanelMainString = "ERAI Jet Stream"
  gsn_panel(wks, plot, (/6,1/), pres)
  
  ; res@cnLevelSelectionMode = "ManualLevels"
  ; res@cnLevelSpacingF = 5
  ; res@cnMaxLevelValF = 30.
  ; res@cnMinLevelValF = 10.
  do i = 0, 5
    res@gsnLeftString = month(i+3)
    plot(i) = gsn_csm_contour_map(wks, speed(0,i+3,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i), plot_topo(i))
  end do

  pres@gsnPanelMainString = "Ctrl Jet Stream"
  gsn_panel(wks, plot, (/6,1/), pres)


  ; res@cnLevelSelectionMode = "ManualLevels"
  ; res@cnLevelSpacingF = 10
  ; res@cnMaxLevelValF = 40.
  ; res@cnMinLevelValF = 10.
  do i = 0, 5
    res@gsnLeftString = month(i+3)
    plot(i) = gsn_csm_contour_map(wks, speed(1,i+3,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
    overlay(plot(i), plot_topo(i))
  end do

  pres@gsnPanelMainString = "NT Jet Stream"
  gsn_panel(wks, plot, (/6,1/), pres)


end