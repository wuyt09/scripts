
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

u1 = fu1->U(:,{500},:,:)
u2 = fu2->U(:,{500},:,:)
u3 = fu3->U(:,{500},:,:)
v1 = fv1->V(:,{500},:,:)
v2 = fv2->V(:,{500},:,:)
v3 = fv3->V(:,{500},:,:)
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

;Season = (/"DJF","MAM","JJA","SON"/)
Season = "DJF"
pp11 = month_to_season(p1, Season)
pp1 = dim_avg_n_Wrap(pp11, 0)
pp22 = month_to_season(p2, Season)
pp2 = dim_avg_n_Wrap(pp22, 0)
pp33 = month_to_season(p3, Season)
pp3 = dim_avg_n_Wrap(pp33, 0)

uu11 = month_to_season(u1, Season)
uu1 = dim_avg_n_Wrap(uu11, 0)
uu22 = month_to_season(u2, Season)
uu2 = dim_avg_n_Wrap(uu22, 0)
uu33 = month_to_season(u3, Season)
uu3 = dim_avg_n_Wrap(uu33, 0)

vv11 = month_to_season(v1, Season)
vv1 = dim_avg_n_Wrap(vv11, 0)
vv22 = month_to_season(v2, Season)
vv2 = dim_avg_n_Wrap(vv22, 0)
vv33 = month_to_season(v3, Season)
vv3 = dim_avg_n_Wrap(vv33, 0)
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

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/JJA/"
figname = "draw_clim_JJA_CTRL_S1_S2_UV500_Hgt500"

lonstr = 20 ;20
lonend = 160 ;220
latstr = -20
latend = 70

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 5
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
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.025            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------
vcres@vcRefMagnitudeF         = 2.0  ; define vector ref mag
res@lbLabelBarOn = False
res@gsnRightString = "gpm"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 10.
res@cnMaxLevelValF = 50.
res@cnMinLevelValF = -50.
res@cnFillPalette = "BlueWhiteOrangeRed"

res@gsnLeftString = "thermal"
plot(0) = gsn_csm_contour_map(wks, p_h, res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(0) = gsn_csm_vector(wks, u_h, v_h, vcres)
overlay(plot(0), plot_topo(0))
overlay(plot(0), plotv(0))

res@gsnLeftString = "mechanical"
plot(1) = gsn_csm_contour_map(wks, p_m, res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(1) = gsn_csm_vector(wks, u_m, v_m, vcres)
overlay(plot(1), plot_topo(1))
overlay(plot(1), plotv(1))


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = Season+"_UV500_Hgt500_difference"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot(0:1), (/1,2/), pres)

;------------------------------------------
vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
;vcres@vcRefLengthF            = 0.035            ; define length of vec ref
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF  = 5900
res@cnMinLevelValF  = 5870
res@cnLevelSpacingF = 5
res@cnFillPalette = "WhiteYellowOrangeRed"


res@gsnLeftString = "CTRL"
plot(2) = gsn_csm_contour_map(wks, pp1, res)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2) = gsn_csm_vector(wks, uu1, vv1, vcres)
overlay(plot(2), plot_topo(2))
overlay(plot(2), plotv(2))

res@gsnLeftString = "NOQ1"
plot(3) = gsn_csm_contour_map(wks, pp2, res)
plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(3) = gsn_csm_vector(wks, uu2, vv2, vcres)
overlay(plot(3), plot_topo(3))
overlay(plot(3), plotv(3))

res@gsnLeftString = "NOTP_NOQ1"
plot(4) = gsn_csm_contour_map(wks, pp3, res)
plot_topo(4) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(4) = gsn_csm_vector(wks, uu3, vv3, vcres)
overlay(plot(4), plot_topo(4))
overlay(plot(4), plotv(4))


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = Season+"_UV500_Hgt500_climatology"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot(2:4), (/1,3/), pres)

end