
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

;================================================================
;fc1 = addfile("/home/yangsong3/data-model/CESM_CTRL/F2000_F19_CAM4_CTRL/pro/CTRL.cam.h0.Z3.1-110.nc", "r")
fc1 = addfile(diri1+filic1, "r")
fc2 = addfile(diri2+filic2, "r")
fc3 = addfile(diri3+filic3, "r")

p1 = fc1->Z3(:,{500},:,:)
p2 = fc2->Z3(:,{500},:,:)
p3 = fc3->Z3(:,{500},:,:)
printVarSummary(p1)
;==================================================================

Season = (/"MAM","JJA","SON","DJF"/)
p11 = month_to_seasonN(p1, Season)
p22 = month_to_seasonN(p2, Season)
p33 = month_to_seasonN(p3, Season)
u11 = month_to_seasonN(u1, Season)
u22 = month_to_seasonN(u2, Season)
u33 = month_to_seasonN(u3, Season)
v11 = month_to_seasonN(v1, Season)
v22 = month_to_seasonN(v2, Season)
v33 = month_to_seasonN(v3, Season)

pp1 = dim_avg_n_Wrap(p11, 1)
pp2 = dim_avg_n_Wrap(p22, 1)
pp3 = dim_avg_n_Wrap(p33, 1)

uu1 = dim_avg_n_Wrap(u11, 1)
uu2 = dim_avg_n_Wrap(u22, 1)
uu3 = dim_avg_n_Wrap(u33, 1)

vv1 = dim_avg_n_Wrap(v11, 1)
vv2 = dim_avg_n_Wrap(v22, 1)
vv3 = dim_avg_n_Wrap(v33, 1)
printVarSummary(uu1)

p_m = pp2 - pp3
p_h = pp1 - pp2
copy_VarCoords(pp1, p_m)
copy_VarCoords(pp1, p_h)

u_m = uu2 - uu3
u_h = uu1 - uu2
copy_VarCoords(uu1, u_m)
copy_VarCoords(uu1, u_h)

v_m = vv2 - vv3
v_h = vv1 - vv2
copy_VarCoords(vv1, v_m)
copy_VarCoords(vv1, v_h)

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "fig4_CLIM_SEA_THERMAL_DYNAMIC_MED_LEVEL_horizon"

lonstr = 40 ;20
lonend = 180 ;220
latstr = -10
latend = 70

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 8
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
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
  ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02

;================================================================================

res@lbLabelBarOn = False

;res@cnLevelSelectionMode = "ManualLevels"
;res@cnLevelSpacingF = 10.
;res@cnMaxLevelValF = 50.
;res@cnMinLevelValF = -50.
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/-40,-30,-20,-10,-5,5,10,20,30,40/)
res@cnFillColors = (/8,24,40,56,72,120,168,184,200,216,232/)

vcres@vcRefMagnitudeF         = 4.0             ; define vector ref mag
vcres@vcRefLengthF            = 0.035            ; define length of vec ref

do i = 0, dimsizes(Season)-1

res@gsnRightString = Season(i)
res@gsnLeftString = "Thermal effects" ;"CTRL - TP&NOQ1"
plot(i) = gsn_csm_contour_map_ce(wks, p_h(i,:,:), res)
plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(i) = gsn_csm_vector(wks, u_h(i,:,:), v_h(i,:,:), vcres)
overlay(plot(i), plot_topo(i))
overlay(plot(i), plotv(i))

res@gsnLeftString = "Dynamic effects"  ;"TP&NOQ1 - NOTP&NOQ1"
plot(i+4) = gsn_csm_contour_map_ce(wks, p_m(i,:,:), res)
plot_topo(i+4) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(i+4) = gsn_csm_vector(wks, u_m(i,:,:), v_m(i,:,:), vcres)
overlay(plot(i+4), plot_topo(i+4))
overlay(plot(i+4), plotv(i+4))

end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelMainString = "UV200_Hgt500"
pres@gsnPanelLabelBar = True
pres@pmLabelBarWidthF = 0.55
pres@lbLabelFontHeightF = 0.018
pres@pmLabelBarOrthogonalPosF = -0.01

draw(plot(7))
frame(wks)

gsn_panel(wks, plot, (/2,4/), pres)


end