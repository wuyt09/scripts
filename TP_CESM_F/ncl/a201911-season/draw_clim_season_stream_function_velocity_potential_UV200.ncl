
begin

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

filiu1 = "F2000_CTRL.plevel.U.mon.0601-3012.nc"
filiu2 = "F2000_NOQ1.plevel.U.mon.0601-3012.nc"
filiu3 = "F2000_NOTP_NOQ1.plevel.U.mon.0601-3012.nc"

filiv1 = "F2000_CTRL.plevel.V.mon.0601-3012.nc"
filiv2 = "F2000_NOQ1.plevel.V.mon.0601-3012.nc"
filiv3 = "F2000_NOTP_NOQ1.plevel.V.mon.0601-3012.nc"

filic1 = "F2000_CTRL.plevel.Z3.mon.0601-3012.nc"
filic2 = "F2000_NOQ1.plevel.Z3.mon.0601-3012.nc"
filic3 = "F2000_NOTP_NOQ1.plevel.Z3.mon.0601-3012.nc"


fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")

u1 = fu1->U(:,{200},:,:)
u2 = fu2->U(:,{200},:,:)
u3 = fu3->U(:,{200},:,:)
v1 = fv1->V(:,{200},:,:)
v2 = fv2->V(:,{200},:,:)
v3 = fv3->V(:,{200},:,:)
printVarSummary(u1)

dv11 = uv2sfvpF(u1, v1)
dv22 = uv2sfvpF(u2, v2)
dv33 = uv2sfvpF(u3, v3)
copy_VarCoords(u1, dv11(0,:,:,:))
copy_VarCoords(u2, dv22(0,:,:,:))
copy_VarCoords(u3, dv33(0,:,:,:))
printVarSummary(dv11)

dv1 = dv11(1,:,:,:)
dv2 = dv22(1,:,:,:)
dv3 = dv33(1,:,:,:)
; exit
; ud1 = dv1 ;new(dimsizes(dv1), typeof(dv1))
; vd1 = dv1 ;new(dimsizes(dv1), typeof(dv1))
; ud2 = dv2 ;new(dimsizes(dv2), typeof(dv2))
; vd2 = dv2 ;new(dimsizes(dv2), typeof(dv2))
; ud3 = dv3 ;new(dimsizes(dv3), typeof(dv3))
; vd3 = dv3 ;new(dimsizes(dv3), typeof(dv3))

; dv2uvf(dv1, ud1, vd1)
; dv2uvf(dv2, ud2, vd2)
; dv2uvf(dv3, ud3, vd3)
; printVarSummary(ud1)

;=========================================

Season = (/"MAM","JJA","SON","DJF"/)

dv1 = dv1*1.0e-6
dv2 = dv2*1.0e-6
dv3 = dv3*1.0e-6
x1 = month_to_seasonN(dv1, Season)
x11 = dim_avg_n_Wrap(x1, 1)
x2 = month_to_seasonN(dv2, Season)
x22 = dim_avg_n_Wrap(x2, 1)
x3 = month_to_seasonN(dv3, Season)
x33 = dim_avg_n_Wrap(x3, 1)

x_h = x11 - x22
x_m = x22 - x33
copy_VarCoords(x11, x_h)
copy_VarCoords(x22, x_m)
printVarSummary(x_h)

; ud1_Sea = month_to_seasonN(ud1, Season)
; vd1_Sea = month_to_seasonN(vd1, Season)
; ud2_Sea = month_to_seasonN(ud2, Season)
; vd2_Sea = month_to_seasonN(vd2, Season)
; ud3_Sea = month_to_seasonN(ud3, Season)
; vd3_Sea = month_to_seasonN(vd3, Season)
; ud1_sea = dim_avg_n_Wrap(ud1_Sea, 1)
; vd1_sea = dim_avg_n_Wrap(vd1_Sea, 1)
; ud2_sea = dim_avg_n_Wrap(ud2_Sea, 1)
; vd2_sea = dim_avg_n_Wrap(vd2_Sea, 1)
; ud3_sea = dim_avg_n_Wrap(ud3_Sea, 1)
; vd3_sea = dim_avg_n_Wrap(vd3_Sea, 1)

; uh_sea = ud1_sea
; uh_sea = ud1_sea - ud2_sea
; um_sea = ud2_sea
; um_sea = ud2_sea - ud3_sea

; vh_sea = vd1_sea
; vh_sea = vd1_sea - vd2_sea
; vm_sea = vd2_sea
; vm_sea = vd2_sea - vd3_sea

; printVarSummary(uh_sea)

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_season_stream_function_velocity_potential_UV200"

lonstr = 40 ;20
lonend = 180 ;220
latstr = -10
latend = 80

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 8
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
 res_topo@cnLineColor            = "green4"
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

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@lbLabelBarOn = False
res@cnFillPalette = "BlueWhiteOrangeRed"

; ;------------------------------------------

;   vcres                         = True             ; plot mods desired
;   vcres@gsnFrame                = False            ; so we can draw time stamp
;   vcres@gsnDraw                 = False
;   vcres@vcRefAnnoOrthogonalPosF = -1.0 ;   -0.35           ; move ref vector up
;   ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
;   ;vcres@vcMinDistanceF          = 0.017
;   vcres@vcGlyphStyle            = "LineArrow"
;   vcres@vcMinDistanceF          = 0.04
 
;   vcres@vcLineArrowThicknessF   = 3.0
;   ;vecres@vcLineArrowHeadMinSizeF = 0.008
;   ;vecres@vcLineArrowHeadMaxSizeF = 0.3

;   vcres@vcRefAnnoString2On = False
;   vcres@vcRefAnnoFontHeightF = 0.02
;   vcres@gsnLeftString = ""
;   vcres@gsnRightString = ""

; ;------------------------------------------

;   vcres@vcRefMagnitudeF         = 1.0             ; define vector ref mag
;   vcres@vcRefLengthF            = 0.025            ; define length of vec ref

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.2
res@cnMaxLevelValF = 3.
res@cnMinLevelValF = -3.

do i = 0, dimsizes(Season)-1

res@gsnRightString = Season(i)
res@gsnLeftString = "Thermal effects"
plot(i) = gsn_csm_contour_map(wks, x_h(i,:,:), res)
plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
;plotv(i) = gsn_csm_vector(wks, uh_sea(i,:,:), vh_sea(i,:,:), vcres)
overlay(plot(i), plot_topo(i))
;overlay(plot(i), plotv(i))

res@gsnLeftString = "Dynamical effects"
plot(i+4) = gsn_csm_contour_map(wks, x_m(i,:,:), res)
plot_topo(i+4) = gsn_csm_contour(wks, tpmask, res_topo)
;plotv(i+4) = gsn_csm_vector(wks, um_sea(i,:,:), vm_sea(i,:,:), vcres)
overlay(plot(i+4), plot_topo(i+4))
;overlay(plot(i+4), plotv(i+4))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV200_velocity_potential"
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/2,4/), pres)
end