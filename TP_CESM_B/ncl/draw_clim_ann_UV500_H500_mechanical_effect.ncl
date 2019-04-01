
begin

yrStat = 260
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_B2000_F19G16_CAM4/pro/"

filiu1 = "B2000_f19g16_CP_CTRL.plevel.U.mon.260-304.nc"
filiu2 = "TP_NOTOPO.plevel.U.mon.260-304.nc"
filiu3 = "TP_HEAT.plevel.U.mon.260-304.nc"

filiv1 = "B2000_f19g16_CP_CTRL.plevel.V.mon.260-304.nc"
filiv2 = "TP_NOTOPO.plevel.V.mon.260-304.nc"
filiv3 = "TP_HEAT.plevel.V.mon.260-304.nc"

filic1 = "B2000_f19g16_CP_CTRL.plevel.Z3.mon.260-304.nc"
filic2 = "TP_NOTOPO.plevel.Z3.mon.260-304.nc"
filic3 = "TP_HEAT.plevel.Z3.mon.260-304.nc"


fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")

u1 = fu1->U
u2 = fu2->U
u3 = fu3->U
v1 = fv1->V
v2 = fv2->V
v3 = fv3->V
printVarSummary(u1)

;================================================================
fc1 = addfile(diri1+filic1, "r")
fc2 = addfile(diri2+filic2, "r")
fc3 = addfile(diri3+filic3, "r")

p1 = fc1->Z3
p2 = fc2->Z3
p3 = fc3->Z3
printVarSummary(p1)
;==================================================================

pp1 = dim_avg_n_Wrap(p1, 0)
pp2 = dim_avg_n_Wrap(p2, 0)
pp3 = dim_avg_n_Wrap(p3, 0)

uu1 = dim_avg_n_Wrap(u1, 0)
uu2 = dim_avg_n_Wrap(u2, 0)
uu3 = dim_avg_n_Wrap(u3, 0)

vv1 = dim_avg_n_Wrap(v1, 0)
vv2 = dim_avg_n_Wrap(v2, 0)
vv3 = dim_avg_n_Wrap(v3, 0)
printVarSummary(uu1)

p_m = pp1 - pp3
p_t = pp1 - pp2
p_h = p_t - p_m
copy_VarCoords(pp1, p_m)
copy_VarCoords(pp1, p_t)
copy_VarCoords(pp1, p_h)

u_m = uu1 - uu3
u_t = uu1 - uu2
u_h = u_t - u_m
copy_VarCoords(uu1, u_m)
copy_VarCoords(uu1, u_t)
copy_VarCoords(uu1, u_h)

v_m = vv1 - vv3
v_t = vv1 - vv2
v_h = v_t - v_m
copy_VarCoords(vv1, v_m)
copy_VarCoords(vv1, v_t)
copy_VarCoords(vv1, v_h)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_ann_UV500_H500_mechanical_effect"

lonstr = 40
lonend = 180
latstr = 0
latend = 70

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 6
plot = new(nplot, graphic)
plotv = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF =  -0.5 ;-1.0             ; move ref vector up
  vcres@vcRefMagnitudeF         = 5.0             ; define vector ref mag
  vcres@vcRefLengthF            = 0.045            ; define length of vec ref
  vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  vcres@vcMinDistanceF          = 0.017
  ;vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False

;------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 200.
res@cnMaxLevelValF = 6000.
res@cnMinLevelValF = 5000.

res@cnFillPalette = "WhiteYellowOrangeRed"

res@gsnLeftString = "CTRL"
plot(0) = gsn_csm_contour_map(wks, pp1({500},:,:), res)
;plot(0) = gsn_csm_vector_scalar_map(wks, uu1({850},:,:), vv1({850},:,:), pp1, res)
plotv(0) = gsn_csm_vector(wks, uu1({500},:,:), vv1({500},:,:), vcres)
overlay(plot(0), plotv(0))

res@gsnLeftString = "NO_TP"
plot(2) = gsn_csm_contour_map(wks, pp2({500},:,:), res)
;plot(2) = gsn_csm_vector_scalar_map(wks, uu2({850},:,:), vv2({850},:,:), pp2, res)
plotv(2) = gsn_csm_vector(wks, uu2({500},:,:), vv2({500},:,:), vcres)
overlay(plot(2), plotv(2))

res@gsnLeftString = "NOTP_HEAT"
plot(4) = gsn_csm_contour_map(wks, pp3({500},:,:), res)
;plot(4) = gsn_csm_vector_scalar_map(wks, uu3({850},:,:), vv3({850},:,:), pp3, res)
plotv(4) = gsn_csm_vector(wks, uu3({500},:,:), vv3({500},:,:), vcres)
overlay(plot(4), plotv(4))

;res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 5.
res@cnMaxLevelValF = 20.
res@cnMinLevelValF = -20.
res@cnFillPalette = "BlueWhiteOrangeRed"

vcres@vcRefMagnitudeF         = 1.0             ; define vector ref mag

res@gsnLeftString = "Total (CTRL - NO_TP)"
plot(1) = gsn_csm_contour_map(wks, p_t({500},:,:), res)
;plot(1) = gsn_csm_vector_scalar_map(wks, u_t({850},:,:), v_t({850},:,:), p_t, res)
plotv(1) = gsn_csm_vector(wks, u_t({500},:,:), v_t({500},:,:), vcres)
overlay(plot(1), plotv(1))

res@gsnLeftString = "Mecha (CTRL - NOTP_HEAT)"
plot(3) = gsn_csm_contour_map(wks, p_m({500},:,:), res)
;plot(3) = gsn_csm_vector_scalar_map(wks, u_m({850},:,:), v_m({850},:,:), p_m, res)
plotv(3) = gsn_csm_vector(wks, u_m({500},:,:), v_m({500},:,:), vcres)
overlay(plot(3), plotv(3))

res@gsnLeftString = "Heat (Total - Mecha)"
plot(5) = gsn_csm_contour_map(wks, p_h({500},:,:), res)
;plot(5) = gsn_csm_vector_scalar_map(wks, u_h({850},:,:), v_h({850},:,:), p_h, res)
plotv(5) = gsn_csm_vector(wks, u_h({500},:,:), v_h({500},:,:), vcres)
overlay(plot(5), plotv(5))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV500_H500"

gsn_panel(wks, plot, (/3,2/), pres)
end