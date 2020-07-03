
begin
LEV = 300

casename = (/"CTRL","NOTP","NOQ1","NOTP_NOQ1"/)
ncase = dimsizes(casename)
u = new((/ncase,12,96,144/), float)
v = u
z = u
season = (/"MAM","JJA","SON","DJF"/)
do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    filiz = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fz = addfile(diri+filiz, "r")
    U = fu->U(:,{LEV},:,:) ;240:299
    V = fv->V(:,{LEV},:,:) ;240:299
    Z = fz->Z3(:,{LEV},:,:) ;240:299
    printVarSummary(U)
    u(i,:,:,:) = clmMonTLL(U)
    v(i,:,:,:) = clmMonTLL(V)
    z(i,:,:,:) = clmMonTLL(Z)
delete([/diri,filiu,filiv,filiz,U,V,Z/])
end do
printVarSummary(u)

u_diff = u(0:3,:,:,:)
v_diff = u_diff
z_diff = u_diff

u_diff(0,:,:,:) = u(0,:,:,:) - u(2,:,:,:)
u_diff(1,:,:,:) = u(0,:,:,:) - u(1,:,:,:)
u_diff(2,:,:,:) = u(2,:,:,:) - u(3,:,:,:)
u_diff(3,:,:,:) = u_diff(0,:,:,:) + u_diff(2,:,:,:)
v_diff(0,:,:,:) = v(0,:,:,:) - v(2,:,:,:)
v_diff(1,:,:,:) = v(0,:,:,:) - v(1,:,:,:)
v_diff(2,:,:,:) = v(2,:,:,:) - v(3,:,:,:)
v_diff(3,:,:,:) = v_diff(0,:,:,:) + v_diff(2,:,:,:)
z_diff(0,:,:,:) = z(0,:,:,:) - z(2,:,:,:)
z_diff(1,:,:,:) = z(0,:,:,:) - z(1,:,:,:)
z_diff(2,:,:,:) = z(2,:,:,:) - z(3,:,:,:)
z_diff(3,:,:,:) = z_diff(0,:,:,:) + z_diff(2,:,:,:)

uwarm_diff = dim_avg_n_Wrap(u_diff(:,3:8,:,:), 1)
vwarm_diff = dim_avg_n_Wrap(v_diff(:,3:8,:,:), 1)
zwarm_diff = dim_avg_n_Wrap(z_diff(:,3:8,:,:), 1)

ucool_diff = u_diff(:,0,:,:)
ucool_diff = dim_avg_n_Wrap(u_diff(:,9:11,:,:), 1)+dim_avg_n_Wrap(u_diff(:,0:2,:,:), 1)
vcool_diff = v_diff(:,0,:,:)
vcool_diff = dim_avg_n_Wrap(v_diff(:,9:11,:,:), 1)+dim_avg_n_Wrap(v_diff(:,0:2,:,:), 1)
zcool_diff = z_diff(:,0,:,:)
zcool_diff = dim_avg_n_Wrap(z_diff(:,9:11,:,:), 1)+dim_avg_n_Wrap(z_diff(:,0:2,:,:), 1)

;==================================================================

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191130/"
figname = "UV"+LEV+"_H"+LEV+"_12month"

lonstr = 0 ;20
lonend = 360 ;220
latstr = 10
latend = 80

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 24
plot = new(nplot, graphic)
plotv = plot
plot_topo = plot

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
;res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF =  -0.4 ;-1.0             ; move ref vector up
  ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefMagnitudeF         = 10.0             ; define vector ref mag
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.015

;================================================================================

res@lbLabelBarOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 10.
res@cnMaxLevelValF = 70.
res@cnMinLevelValF = -70.
res@cnFillPalette = "BlueWhiteOrangeRed"

vcres@vcRefMagnitudeF         = 4.0             ; define vector ref mag
vcres@vcRefLengthF            = 0.035            ; define length of vec ref

; do i = 0, dimsizes(season)-1

; res@gsnRightString = season(i)
; res@gsnLeftString = "Thermal (CTRL - NOQ1)"
; plot(3*i) = gsn_csm_contour_map_ce(wks, z_diff(0,i,:,:), res)
; plot_topo(3*i) = gsn_csm_contour(wks, tpmask, res_topo)
; plotv(3*i) = gsn_csm_vector(wks, u_diff(0,i,:,:), v_diff(0,i,:,:), vcres)
; overlay(plot(3*i), plot_topo(3*i))
; overlay(plot(3*i), plotv(3*i))

; res@gsnLeftString = "Dynamical (CTRL - NOTP)"
; plot(3*i+1) = gsn_csm_contour_map_ce(wks, z_diff(1,i,:,:), res)
; plot_topo(3*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
; plotv(3*i+1) = gsn_csm_vector(wks, u_diff(1,i,:,:), v_diff(1,i,:,:), vcres)
; overlay(plot(3*i+1), plot_topo(3*i+1))
; overlay(plot(3*i+1), plotv(3*i+1))

; res@gsnLeftString = "Dynamical (NOTP - NOTP_NOQ1)"
; plot(3*i+2) = gsn_csm_contour_map_ce(wks, z_diff(2,i,:,:), res)
; plot_topo(3*i+2) = gsn_csm_contour(wks, tpmask, res_topo)
; plotv(3*i+2) = gsn_csm_vector(wks, u_diff(2,i,:,:), v_diff(2,i,:,:), vcres)
; overlay(plot(3*i+2), plot_topo(3*i+2))
; overlay(plot(3*i+2), plotv(3*i+2))
; end do

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelMainString = "UV"+LEV+"_H"+LEV
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = -0.01

; gsn_panel(wks, plot, (/4,3/), pres)


do i = 0, 11

res@gsnRightString = "Month="+(i+1)
res@gsnLeftString = "Thermal"
plot(2*i) = gsn_csm_contour_map_ce(wks, z_diff(1,i,:,:), res)
plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i) = gsn_csm_vector(wks, u_diff(1,i,:,:), v_diff(1,i,:,:), vcres)
overlay(plot(2*i), plot_topo(2*i))
overlay(plot(2*i), plotv(2*i))

res@gsnLeftString = "Dynamical"
plot(2*i+1) = gsn_csm_contour_map_ce(wks, z_diff(2,i,:,:), res)
plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i+1) = gsn_csm_vector(wks, u_diff(2,i,:,:), v_diff(2,i,:,:), vcres)
overlay(plot(2*i+1), plot_topo(2*i+1))
overlay(plot(2*i+1), plotv(2*i+1))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV"+LEV+"_H"+LEV
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/12,2/), pres)


leftsti = (/"CTRL-NOQ1","CTRL-NOTP","NOQ1-NOTP_NOQ1","CTRL-NOTP_NOQ1"/)
do i = 0, 3

res@gsnRightString = "Warm season (Apr-Sep)"
res@gsnLeftString = leftsti(i)
plot(2*i) = gsn_csm_contour_map_ce(wks, zwarm_diff(i,:,:), res)
plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i) = gsn_csm_vector(wks, uwarm_diff(i,:,:), vwarm_diff(i,:,:), vcres)
overlay(plot(2*i), plot_topo(2*i))
overlay(plot(2*i), plotv(2*i))

res@gsnRightString = "Cool season (Oct-Mar)"
res@gsnLeftString = leftsti(i)
plot(2*i+1) = gsn_csm_contour_map_ce(wks, zcool_diff(i,:,:), res)
plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i+1) = gsn_csm_vector(wks, ucool_diff(i,:,:), vcool_diff(i,:,:), vcres)
overlay(plot(2*i+1), plot_topo(2*i+1))
overlay(plot(2*i+1), plotv(2*i+1))

end do

gsn_panel(wks, plot, (/4,2/), pres)

end