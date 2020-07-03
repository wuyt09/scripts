
begin


lonstr = 40
lonend = 220
latstr = 0
latend = 60

casename = (/"CTRL","NOTP"/) ;"NOQ1","NOTP_NOQ1",
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)
season = (/"MAM","JJA","SON","DJF"/)
ns = dimsizes(season)

LEV = 200
u = new((/ncase,4,96,144/), float)
v = u
z = u

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    filit = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
    filips = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    u1 = fu->U(:,{LEV},:,:)
    v1 = fv->V(:,{LEV},:,:)
    U = month_to_seasonN(u1, season)
    u(i,:,:,:) = dim_avg_n_Wrap(U, 1)
    V = month_to_seasonN(v1, season)
    v(i,:,:,:) = dim_avg_n_Wrap(V, 1)

    ft = addfile(diri+filit, "r")
    z1 = ft->Z3(:,{500},:,:)
    Z = month_to_seasonN(z1, season)
    z(i,:,:,:) = dim_avg_n_Wrap(Z, 1)
end do

ux = u(0,:,:,:)
ux = u(1,:,:,:) - u(0,:,:,:)
vx = v(0,:,:,:)
vx = v(1,:,:,:) - v(0,:,:,:)
zx = z(0,:,:,:)
zx = z(1,:,:,:) - z(0,:,:,:)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200106/"
figname = "UV"+LEV+"_H500_Ctrl_TPforcing"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4*2
plot = new(nplot, graphic)
plot_topo = plot
plotv = plot


 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 16
 res_topo@cnLineColor            = "gray"
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
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up

  ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.035
  vcres@vcRefLengthF            = 0.04            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
;res@cnLevelSelectionMode = "ExplicitLevels"

do i = 0, ns-1

    res@gsnRightString = season(i)
    res@gsnLeftString = "Ctrl"
    res@cnLevelSpacingF = 10.
    res@cnMaxLevelValF = 5920.
    res@cnMinLevelValF = 5860.
    ; delete(res@cnLevels)
    ; delete(res@cnFillColors)
    ; res@cnLevels = (/5500,5600,5700,5800,5850,5870,5890,5900/)
    ; res@cnFillColors = (/18,50,82,114,130,136,142,145,131/)
    res@cnFillPalette = "WhiteYellowOrangeRed"
    vcres@vcRefMagnitudeF         = 20.0
    plot(i) = gsn_csm_contour_map(wks, z(0,i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(i) = gsn_csm_vector(wks, u(0,i,:,:), v(0,i,:,:), vcres)
    overlay(plot(i), plot_topo(i))
    overlay(plot(i), plotv(i))

    res@gsnLeftString = "NT - Ctrl"
    res@cnLevelSpacingF = 10
    res@cnMaxLevelValF = 100
    res@cnMinLevelValF = -100
    ; delete(res@cnLevels)
    ; delete(res@cnFillColors)
    ; res@cnLevels = (/-100,-80,-60,-40,-20,-10,10,20,50/)
    ; res@cnLineColors = (/2,17,34,51,68,85,102,128,144,192,208/)
    res@cnFillPalette = "BlueWhiteOrangeRed"
    vcres@vcRefMagnitudeF         = 5.0
    plot(i+4) = gsn_csm_contour_map(wks, zx(i,:,:), res)
    plot_topo(i+4) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(i+4) = gsn_csm_vector(wks, ux(i,:,:), vx(i,:,:), vcres)
    overlay(plot(i+4), plot_topo(i+4))
    overlay(plot(i+4), plotv(i+4))

end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelMainString = "UV"+LEV+"_Precip"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@lbLabelFontHeightF = 0.01

gsn_panel(wks, plot(0:3), (/1,4/), pres)
gsn_panel(wks, plot(4:7), (/1,4/), pres)






end