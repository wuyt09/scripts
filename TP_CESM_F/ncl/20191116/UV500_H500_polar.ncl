
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

;==================================================================
; t test
;==================================================================
p1_var = dim_variance_n_Wrap(p11, 1)
u1_var = dim_variance_n_Wrap(u11, 1)
v1_var = dim_variance_n_Wrap(v11, 1)
p2_var = dim_variance_n_Wrap(p22, 1)
u2_var = dim_variance_n_Wrap(u22, 1)
v2_var = dim_variance_n_Wrap(v22, 1)
p3_var = dim_variance_n_Wrap(p33, 1)
u3_var = dim_variance_n_Wrap(u33, 1)
v3_var = dim_variance_n_Wrap(v33, 1)
s = 25 ; n years
prob_ph12 = ttest(pp1, p1_var, s, pp2, p2_var, s, False, False)
prob_pm23 = ttest(pp2, p2_var, s, pp3, p3_var, s, False, False)
prob_uh12 = ttest(uu1, u1_var, s, uu2, u2_var, s, False, False)
prob_um23 = ttest(uu2, u2_var, s, uu3, u3_var, s, False, False)
prob_vh12 = ttest(vv1, v1_var, s, vv2, v2_var, s, False, False)
prob_vm23 = ttest(vv2, v2_var, s, vv3, v3_var, s, False, False)

siglev = 0.05
p_m = pp2 - pp3
p_h = pp1 - pp2
; p_m = mask(p_m, prob_pm23.le.siglev, True)
; p_h = mask(p_h, prob_ph12.le.siglev, True)
copy_VarCoords(pp1, p_m)
copy_VarCoords(pp1, p_h)

u_m = uu2 - uu3
u_h = uu1 - uu2
; u_m = mask(u_m, prob_um23.le.siglev .and. prob_vm23.le.siglev, True)
; u_h = mask(u_m, prob_uh12.le.siglev .and. prob_vh12.le.siglev, True)
copy_VarCoords(uu1, u_m)
copy_VarCoords(uu1, u_h)

v_m = vv2 - vv3
v_h = vv1 - vv2
; v_m = mask(v_m, prob_um23.le.siglev .and. prob_vm23.le.siglev, True)
; v_h = mask(v_m, prob_uh12.le.siglev .and. prob_vh12.le.siglev, True)
copy_VarCoords(vv1, v_m)
copy_VarCoords(vv1, v_h)

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191116/"
figname = "UV500_H500_polar"

lonstr = 0 ;20
lonend = 360 ;220
latstr = 0
latend = 90

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
 res_topo@cnLineColor            = "blue"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnPolar = "NH"

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey45" ;"grey25"
;res@mpFillDrawOrder = "PreDraw"

; res@mpMinLatF = latstr
; res@mpMaxLatF = latend
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend
; res@mpCenterLonF = (lonstr+lonend)/2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnPolarLabelFontHeightF = 0.025
res@gsnStringFontHeightF = 0.035
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  ;vcres@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up
  ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefMagnitudeF         = 10.0             ; define vector ref mag
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.025

;================================================================================
;------------------------
; draw prob
;------------------------
res3  = True
res3@gsnDraw   = False
res3@gsnFrame  = False

;res3@cnFillOn  = False 
res3@cnLinesOn = False
res3@cnLineLabelsOn = False
res3@cnInfoLabelOn  = False
;res3@lbLabelBarOn   = False

res3@gsnRightString = ""

res3@cnLevelSelectionMode = "ManualLevels"
res3@cnMinLevelValF = 0.00
res3@cnMaxLevelValF = 0.10
res3@cnLevelSpacingF = 0.01

;res3@cnFillPattern = 17 ;(/17/)
res3@cnFillColor  =  "green4" ;""
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1.2
;================================================================================
res@lbLabelBarOn = False

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 10.
; res@cnMaxLevelValF = 70.
; res@cnMinLevelValF = -70.
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/-70,-50,-30,-20,-10,-5,-2,2,5,10,20,30,50,70/)
res@cnFillPalette = "BlueWhiteOrangeRed"

vcres@vcRefMagnitudeF         = 4.0             ; define vector ref mag
vcres@vcRefLengthF            = 0.035            ; define length of vec ref

do i = 0, dimsizes(Season)-1

res@gsnRightString = Season(i)
res@gsnLeftString = "Thermal"
plot(2*i) = gsn_csm_contour_map_polar(wks, p_h(i,:,:), res)
plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i) = gsn_csm_vector(wks, u_h(i,:,:), v_h(i,:,:), vcres)
overlay(plot(2*i), plot_topo(2*i))
overlay(plot(2*i), plotv(2*i))

res@gsnLeftString = "Dynamical"
plot(2*i+1) = gsn_csm_contour_map_polar(wks, p_m(i,:,:), res)
plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i+1) = gsn_csm_vector(wks, u_m(i,:,:), v_m(i,:,:), vcres)
overlay(plot(2*i+1), plot_topo(2*i+1))
overlay(plot(2*i+1), plotv(2*i+1))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV500_H500"
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01
pres@pmLabelBarWidthF = 0.5
gsn_panel(wks, plot, (/4,2/), pres)

end