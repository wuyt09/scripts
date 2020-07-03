
begin

yrStat = 260
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOSH_B2000_F19G16_CAM4/pro/"

filiu1 = "B2000_f19g16_CP_CTRL.plevel.U.mon.260-304.nc"
filiu2 = "TP_NOTOPO.plevel.U.mon.260-304.nc"
filiu3 = "TP_NOSH.plevel.U.mon.260-304.nc"

filiv1 = "B2000_f19g16_CP_CTRL.plevel.V.mon.260-304.nc"
filiv2 = "TP_NOTOPO.plevel.V.mon.260-304.nc"
filiv3 = "TP_NOSH.plevel.V.mon.260-304.nc"

filic1 = "B2000_f19g16_CP_CTRL.plevel.Z3.mon.260-304.nc"
filic2 = "TP_NOTOPO.plevel.Z3.mon.260-304.nc"
filic3 = "TP_NOSH.plevel.Z3.mon.260-304.nc"


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

Season = (/"DJF","MAM","JJA","SON"/)
p11 = month_to_seasonN(p1, Season)
p22 = month_to_seasonN(p2, Season)
p33 = month_to_seasonN(p3, Season)
u11 = month_to_seasonN(u1, Season)
u22 = month_to_seasonN(u2, Season)
u33 = month_to_seasonN(u3, Season)
v11 = month_to_seasonN(v1, Season)
v22 = month_to_seasonN(v2, Season)
v33 = month_to_seasonN(v3, Season)

pp1 = dim_avg_n_Wrap(p11, 0)
pp2 = dim_avg_n_Wrap(p22, 0)
pp3 = dim_avg_n_Wrap(p33, 0)

uu1 = dim_avg_n_Wrap(u11, 0)
uu2 = dim_avg_n_Wrap(u22, 0)
uu3 = dim_avg_n_Wrap(u33, 0)

vv1 = dim_avg_n_Wrap(v11, 0)
vv2 = dim_avg_n_Wrap(v22, 0)
vv3 = dim_avg_n_Wrap(v33, 0)
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
figname = "draw_clim_season_UV500_H500"

lonstr = 40
lonend = 180
latstr = -20
latend = 60

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
plot = new(nplot, graphic)
plotv = plot

plot2 = plot
plotv2 = plotv

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "red" ;"grey25"
res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

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
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefMagnitudeF         = 10.0             ; define vector ref mag
  vcres@vcRefLengthF            = 0.025            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02

;================================================================================

res@lbLabelBarOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 50.
res@cnMaxLevelValF = 6000.
res@cnMinLevelValF = 5500.
res@cnFillPalette = "WhiteBlueGreenYellowRed"

vcres@vcRefMagnitudeF         = 10.0             ; define vector ref mag

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)
res@gsnLeftString = "CTRL"
plot(i) = gsn_csm_contour_map(wks, pp1(i,{500},:,:), res)
plotv(i) = gsn_csm_vector(wks, uu1(i,{500},:,:), vv1(i,{500},:,:), vcres)
overlay(plot(i), plotv(i))

res@gsnLeftString = "NO_TP"
plot(i+4) = gsn_csm_contour_map(wks, pp2(i,{500},:,:), res)
plotv(i+4) = gsn_csm_vector(wks, uu2(i,{500},:,:), vv2(i,{500},:,:), vcres)
overlay(plot(i+4), plotv(i+4))

res@gsnLeftString = "NOSH"
plot(i+8) = gsn_csm_contour_map(wks, pp3(i,{500},:,:), res)
plotv(i+8) = gsn_csm_vector(wks, uu3(i,{500},:,:), vv3(i,{500},:,:), vcres)
overlay(plot(i+8), plotv(i+8))
end do
pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV500_H500"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)

;------------------------------------------
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 10.
res@cnMaxLevelValF = 80.
res@cnMinLevelValF = -80.
res@cnFillPalette = "BlueWhiteOrangeRed"

vcres@vcRefMagnitudeF         = 1.0             ; define vector ref mag

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)
res@gsnLeftString = "CTRL - NOTP"
plot(i) = gsn_csm_contour_map(wks, p_t(i,{500},:,:), res)
plotv(i) = gsn_csm_vector(wks, u_t(i,{500},:,:), v_t(i,{500},:,:), vcres)
overlay(plot(i), plotv(i))

res@gsnLeftString = "CTRL - NOSH"
plot(i+4) = gsn_csm_contour_map(wks, p_m(i,{500},:,:), res)
plotv(i+4) = gsn_csm_vector(wks, u_m(i,{500},:,:), v_m(i,{500},:,:), vcres)
overlay(plot(i+4), plotv(i+4))

res@gsnLeftString = "NOSH - NOTP"
plot(i+8) = gsn_csm_contour_map(wks, p_h(i,{500},:,:), res)
plotv(i+8) = gsn_csm_vector(wks, u_h(i,{500},:,:), v_h(i,{500},:,:), vcres)
overlay(plot(i+8), plotv(i+8))
end do
pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV500_H500"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)
end