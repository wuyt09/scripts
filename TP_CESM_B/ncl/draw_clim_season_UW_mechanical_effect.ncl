
begin

yrStat = 260
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_B2000_F19G16_CAM4/pro/"

filiu1 = "B2000_f19g16_CP_CTRL.plevel.U.mon.260-304.nc"
filiu2 = "TP_NOTOPO.plevel.U.mon.260-304.nc"
filiu3 = "TP_HEAT.plevel.U.mon.260-304.nc"

filiv1 = "B2000_f19g16_CP_CTRL.plevel.OMEGA.mon.260-304.nc"
filiv2 = "TP_NOTOPO.plevel.OMEGA.mon.260-304.nc"
filiv3 = "TP_HEAT.plevel.OMEGA.mon.260-304.nc"

filic1 = "B2000_f19g16_CP_CTRL.plevel.T.mon.260-304.nc"
filic2 = "TP_NOTOPO.plevel.T.mon.260-304.nc"
filic3 = "TP_HEAT.plevel.T.mon.260-304.nc"


fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")

u1 = fu1->U
u2 = fu2->U
u3 = fu3->U
v1 = fv1->OMEGA
v2 = fv2->OMEGA
v3 = fv3->OMEGA
printVarSummary(u1)

;================================================================
p1 = v1
p2 = v2
p3 = v3
p1 = v1*(-1.)
p2 = v2*(-1.)
p3 = v3*(-1.) ; upward is positive value

v1 = v1*(-50.)
v2 = v2*(-50.)
v3 = v3*(-50.)
;==================================================================
u1x = dim_avg_n_Wrap(u1(:,:,{25:40},:), 2)
u2x = dim_avg_n_Wrap(u2(:,:,{25:40},:), 2)
u3x = dim_avg_n_Wrap(u3(:,:,{25:40},:), 2)
v1x = dim_avg_n_Wrap(v1(:,:,{25:40},:), 2)
v2x = dim_avg_n_Wrap(v2(:,:,{25:40},:), 2)
v3x = dim_avg_n_Wrap(v3(:,:,{25:40},:), 2)
p1x = dim_avg_n_Wrap(p1(:,:,{25:40},:), 2)
p2x = dim_avg_n_Wrap(p2(:,:,{25:40},:), 2)
p3x = dim_avg_n_Wrap(p3(:,:,{25:40},:), 2)

Season = (/"DJF","MAM","JJA","SON"/)
p11 = month_to_seasonN(p1x, Season)
p22 = month_to_seasonN(p2x, Season)
p33 = month_to_seasonN(p3x, Season)
u11 = month_to_seasonN(u1x, Season)
u22 = month_to_seasonN(u2x, Season)
u33 = month_to_seasonN(u3x, Season)
v11 = month_to_seasonN(v1x, Season)
v22 = month_to_seasonN(v2x, Season)
v33 = month_to_seasonN(v3x, Season)

uu1 = dim_avg_n_Wrap(u11, 1)
uu2 = dim_avg_n_Wrap(u22, 1)
uu3 = dim_avg_n_Wrap(u33, 1)

vv1 = dim_avg_n_Wrap(v11, 1)
vv2 = dim_avg_n_Wrap(v22, 1)
vv3 = dim_avg_n_Wrap(v33, 1)
printVarSummary(uu1)

pp1 = dim_avg_n_Wrap(p11, 1)
pp2 = dim_avg_n_Wrap(p22, 1)
pp3 = dim_avg_n_Wrap(p33, 1)

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
figname = "draw_clim_season_UW_mechanical_effect"

lonstr = 40
lonend = 180

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
plot = new(nplot, graphic)
plotv = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@vpHeightF = 0.3        ; Changes the aspect ratio
res@vpWidthF  = 0.55

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
res@cnFillDrawOrder = "PreDraw"

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;------------------------------------------
 res@trYMinF   = 100
 res@trYMaxF   = 1000
  
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)

  res@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up
  ;res@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;res@vcMinDistanceF          = 0.017
  res@vcGlyphStyle            = "LineArrow"
  res@vcMinDistanceF          = 0.03
  res@vcRefLengthF            = 0.025            ; define length of vec ref
  
  res@vcLineArrowThicknessF   = 3.0
  ;res@vcLineArrowHeadMinSizeF = 0.008
  ;res@vcLineArrowHeadMaxSizeF = 0.3
  res@vcVectorDrawOrder = "PostDraw" 
 
  res@vcRefAnnoString2On = False
  res@vcRefAnnoFontHeightF = 0.02

;------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.005
res@cnMaxLevelValF = 0.05
res@cnMinLevelValF = -0.05
res@cnFillPalette = "BlueWhiteOrangeRed"
res@vcRefMagnitudeF         = 2.             ; define vector ref mag
res@lbLabelBarOn = False

do i = 0, dimsizes(Season)-1

res@gsnRightString = Season(i)
res@gsnLeftString = "Total (CTRL - NO_TP)"
plot(i) = gsn_csm_pres_hgt_vector(wks, p_t(i,:,{lonstr:lonend}), u_t(i,:,{lonstr:lonend}), v_t(i,:,{lonstr:lonend}), res)
;gsn_csm_contour_map(wks, p_t(i,:,:), res)
;plot(1) = gsn_csm_vector_scalar_map(wks, u_t({850},:,:), v_t({850},:,:), p_t, res)
;plotv(i) = gsn_csm_vector(wks, u_t(i,{200},:,:), v_t(i,{200},:,:), vcres)
;overlay(plot(i), plotv(i))

res@gsnLeftString = "Mecha (CTRL - NOTP_HEAT)"
plot(i+4) = gsn_csm_pres_hgt_vector(wks, p_m(i,:,{lonstr:lonend}), u_m(i,:,{lonstr:lonend}), v_m(i,:,{lonstr:lonend}), res)
;plot(3) = gsn_csm_vector_scalar_map(wks, u_m({850},:,:), v_m({850},:,:), p_m, res)
;plotv(i+4) = gsn_csm_vector(wks, u_m(i,{200},:,:), v_m(i,{200},:,:), vcres)
;overlay(plot(i+4), plotv(i+4))


res@gsnLeftString = "Heat (Total - Mecha)"
plot(i+8) = gsn_csm_pres_hgt_vector(wks, p_h(i,:,{lonstr:lonend}), u_h(i,:,{lonstr:lonend}), v_h(i,:,{lonstr:lonend}), res)
;plot(5) = gsn_csm_vector_scalar_map(wks, u_h({850},:,:), v_h({850},:,:), p_h, res)
;plotv(i+8) = gsn_csm_vector(wks, u_h(i,{200},:,:), v_h(i,{200},:,:), vcres)
;overlay(plot(i+8), plotv(i+8))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UW (25-40N)"
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)
;================================================================================

res@lbLabelBarOn = False

;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.01
res@cnMaxLevelValF = 0.2
res@cnMinLevelValF = -0.2
res@cnFillPalette = "BlueWhiteOrangeRed"

res@vcRefMagnitudeF         = 10.0             ; define vector ref mag

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)
res@gsnLeftString = "CTRL"
plot(i) = gsn_csm_pres_hgt_vector(wks, pp1(i,:,{lonstr:lonend}), uu1(i,:,{lonstr:lonend}), vv1(i,:,{lonstr:lonend}), res)
;plotv(i) = gsn_csm_vector(wks, uu1(i,{200},:,:), vv1(i,{200},:,:), vcres)
;overlay(plot(i), plotv(i))

res@gsnLeftString = "NO_TP"
plot(i+4) = gsn_csm_pres_hgt_vector(wks, pp2(i,:,{lonstr:lonend}), uu2(i,:,{lonstr:lonend}), uu2(i,:,{lonstr:lonend}), res)
;plotv(i+4) = gsn_csm_vector(wks, uu2(i,{200},:,:), vv2(i,{200},:,:), vcres)
;overlay(plot(i+4), plotv(i+4))

res@gsnLeftString = "NOTP_HEAT"
plot(i+8) = gsn_csm_pres_hgt_vector(wks, pp3(i,:,{lonstr:lonend}), uu3(i,:,{lonstr:lonend}), vv3(i,:,{lonstr:lonend}), res)
;plotv(i+8) = gsn_csm_vector(wks, uu3(i,{200},:,:), vv3(i,{200},:,:), vcres)
;overlay(plot(i+8), plotv(i+8))
end do
pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UW (25-40N)"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)
end