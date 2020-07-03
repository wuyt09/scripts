
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

p_m = p22 - p33
p_h = p11 - p22
copy_VarCoords(p11, p_m)
copy_VarCoords(p11, p_h)

u_m = u22 - u33
u_h = u11 - u22
copy_VarCoords(u11, u_m)
copy_VarCoords(u11, u_h)

v_m = v22 - v33
v_h = v11 - v22
copy_VarCoords(v11, v_m)
copy_VarCoords(v11, v_h)

p_m1 = dim_avg_n_Wrap(p_m(:,0:9,:,:), 1)
p_m2 = dim_avg_n_Wrap(p_m(:,15:24,:,:), 1)
p_h1 = dim_avg_n_Wrap(p_h(:,0:9,:,:), 1)
p_h2 = dim_avg_n_Wrap(p_h(:,15:24,:,:), 1)

u_m1 = dim_avg_n_Wrap(u_m(:,0:9,:,:), 1)
u_m2 = dim_avg_n_Wrap(u_m(:,15:24,:,:), 1)
u_h1 = dim_avg_n_Wrap(u_h(:,0:9,:,:), 1)
u_h2 = dim_avg_n_Wrap(u_h(:,15:24,:,:), 1)

v_m1 = dim_avg_n_Wrap(v_m(:,0:9,:,:), 1)
v_m2 = dim_avg_n_Wrap(v_m(:,15:24,:,:), 1)
v_h1 = dim_avg_n_Wrap(v_h(:,0:9,:,:), 1)
v_h2 = dim_avg_n_Wrap(v_h(:,15:24,:,:), 1)
printVarSummary(u_m1)

; p_m = pp2 - pp3
; p_h = pp1 - pp2
; copy_VarCoords(pp1, p_m)
; copy_VarCoords(pp1, p_h)

; u_m = uu2 - uu3
; u_h = uu1 - uu2
; copy_VarCoords(uu1, u_m)
; copy_VarCoords(uu1, u_h)

; v_m = vv2 - vv3
; v_h = vv1 - vv2
; copy_VarCoords(vv1, v_m)
; copy_VarCoords(vv1, v_h)

tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191116/"
figname = "check_model_diff_B10_A10"

lonstr = 0 ;20
lonend = 360 ;220
latstr = 0
latend = 90

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 16
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

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
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

do i = 0, dimsizes(Season)-1

res@gsnRightString = Season(i)
res@gsnLeftString = "Thermal"
plot(2*i) = gsn_csm_contour_map_ce(wks, p_h1(i,:,:), res)
plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i) = gsn_csm_vector(wks, u_h1(i,:,:), v_h1(i,:,:), vcres)
overlay(plot(2*i), plot_topo(2*i))
overlay(plot(2*i), plotv(2*i))

res@gsnLeftString = "Dynamical"
plot(2*i+1) = gsn_csm_contour_map_ce(wks, p_m1(i,:,:), res)
plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i+1) = gsn_csm_vector(wks, u_m1(i,:,:), v_m1(i,:,:), vcres)
overlay(plot(2*i+1), plot_topo(2*i+1))
overlay(plot(2*i+1), plotv(2*i+1))

res@gsnRightString = Season(i)
res@gsnLeftString = "Thermal"
plot(2*i+8) = gsn_csm_contour_map_ce(wks, p_h2(i,:,:), res)
plot_topo(2*i+8) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i+8) = gsn_csm_vector(wks, u_h2(i,:,:), v_h2(i,:,:), vcres)
overlay(plot(2*i+8), plot_topo(2*i+8))
overlay(plot(2*i+8), plotv(2*i+8))

res@gsnLeftString = "Dynamical"
plot(2*i+9) = gsn_csm_contour_map_ce(wks, p_m2(i,:,:), res)
plot_topo(2*i+9) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2*i+9) = gsn_csm_vector(wks, u_m2(i,:,:), v_m2(i,:,:), vcres)
overlay(plot(2*i+9), plot_topo(2*i+9))
overlay(plot(2*i+9), plotv(2*i+9))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "First 10 years"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot(0:7), (/4,2/), pres)


pres@gsnPanelMainString = "Last 10 years"
gsn_panel(wks, plot(8:15), (/4,2/), pres)

end