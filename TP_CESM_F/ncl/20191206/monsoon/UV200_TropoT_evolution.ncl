
begin

lonstr = 40
lonend = 150
latstr = 0
latend = 60

casename = (/"CTRL","NOQ1","NOTP_NOQ1","NOTP"/)
ncase = dimsizes(casename)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

LEV = 200
u = new((/ncase,12,96,144/), float)
v = u
t = u

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    filit = "F2000_"+casename(i)+".plevel.T.mon.0601-3012.nc"
    filips = "F2000_"+casename(i)+".PS.mon.0601-3012.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    u1 = fu->U(:,{LEV},:,:)
    v1 = fv->V(:,{LEV},:,:)
    ;printVarSummary(u1)
    u(i,:,:,:) = clmMonTLL(u1)
    v(i,:,:,:) = clmMonTLL(v1)
    ;printVarSummary(u)

    ft = addfile(diri+filit, "r")
    T = ft->T
    fps = addfile(diri+filips, "r")
    PS = fps->PS

    plev = T&lev*100
    plev!0 = "lev"
    plev&lev = T&lev
    plev@units = "Pa"
    psfc = PS
    ptop = min(plev)
    dp = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    Tdp = T
    Tdp = T*dp
    t1 = dim_sum_n_Wrap(Tdp(:,{400:200},:,:), 1)/dim_sum_n_Wrap(dp(:,{400:200},:,:), 1)
    copy_VarCoords(T(:,0,:,:), t1)
    t(i,:,:,:) = clmMonTLL(t1)
end do

ux = u(0:1,:,:,:)
vx = v(0:1,:,:,:)
tx = t(0:1,:,:,:)

ux(0,:,:,:) = u(0,:,:,:) - u(1,:,:,:)
ux(1,:,:,:) = u(1,:,:,:) - u(2,:,:,:)
vx(0,:,:,:) = v(0,:,:,:) - v(1,:,:,:)
vx(1,:,:,:) = v(1,:,:,:) - v(2,:,:,:)
tx(0,:,:,:) = t(0,:,:,:) - t(1,:,:,:)
tx(1,:,:,:) = t(1,:,:,:) - t(2,:,:,:)
printVarSummary(ux)


;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;==================================================================
;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191206/monsoon/"
figname = "UV"+LEV+"_TropoT_evolution"

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
vcres@vcRefMagnitudeF         = 6.0  ; define vector ref mag
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1
res@cnMaxLevelValF = 5
res@cnMinLevelValF = -5
res@cnFillPalette = "BlueWhiteOrangeRed"

do i = 0, 3

    res@gsnRightString = month(i+5)
    res@gsnLeftString = "Thermal"
    plot(i) = gsn_csm_contour_map(wks, tx(0,i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(i) = gsn_csm_vector(wks, ux(0,i,:,:), vx(0,i,:,:), vcres)
    overlay(plot(i), plot_topo(i))
    overlay(plot(i), plotv(i))

    res@gsnLeftString = "Dynamical"
    plot(i+4) = gsn_csm_contour_map(wks, tx(1,i,:,:), res)
    plot_topo(i+4) = gsn_csm_contour(wks, tpmask, res_topo)
    plotv(i+4) = gsn_csm_vector(wks, ux(1,i,:,:), vx(1,i,:,:), vcres)
    overlay(plot(i+4), plot_topo(i+4))
    overlay(plot(i+4), plotv(i+4))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV"+LEV+"_TropoT_difference"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/2,4/), pres)


end