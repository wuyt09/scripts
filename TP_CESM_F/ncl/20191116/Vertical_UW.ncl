
begin

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

filiu1 = "F2000_CTRL.plevel.U.mon.0601-3012.nc"
filiu2 = "F2000_NOQ1.plevel.U.mon.0601-3012.nc"
filiu3 = "F2000_NOTP_NOQ1.plevel.U.mon.0601-3012.nc"

filiv1 = "F2000_CTRL.plevel.OMEGA.mon.0601-3012.nc"
filiv2 = "F2000_NOQ1.plevel.OMEGA.mon.0601-3012.nc"
filiv3 = "F2000_NOTP_NOQ1.plevel.OMEGA.mon.0601-3012.nc"

filit1 = "F2000_CTRL.plevel.T.mon.0601-3012.nc"
filit2 = "F2000_NOQ1.plevel.T.mon.0601-3012.nc"
filit3 = "F2000_NOTP_NOQ1.plevel.T.mon.0601-3012.nc"

fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")

ft1 = addfile(diri1+filit1, "r")
ft2 = addfile(diri2+filit2, "r")
ft3 = addfile(diri3+filit3, "r")

u1 = fu1->U
u2 = fu2->U
u3 = fu3->U
v1 = fv1->OMEGA
v2 = fv2->OMEGA
v3 = fv3->OMEGA
t1 = ft1->T
t2 = ft2->T
t3 = ft3->T
printVarSummary(u1)

;================================================================
;p1 = v1
;p2 = v2
;p3 = v3
;p1 = v1*(-1.)
;p2 = v2*(-1.)
;p3 = v3*(-1.) ; upward is positive value
p1 = t1
p2 = t2
p3 = t3


v1 = v1*(-50.)
v2 = v2*(-50.)
v3 = v3*(-50.)
;==================================================================
w = u1(0,0,{25:40},0)
w = 1
u1x = dim_avg_wgt_n_Wrap(u1(:,:,{25:40},:), w, 0, 2)
u2x = dim_avg_wgt_n_Wrap(u2(:,:,{25:40},:), w, 0, 2)
u3x = dim_avg_wgt_n_Wrap(u3(:,:,{25:40},:), w, 0, 2)
v1x = dim_avg_wgt_n_Wrap(v1(:,:,{25:40},:), w, 0, 2)
v2x = dim_avg_wgt_n_Wrap(v2(:,:,{25:40},:), w, 0, 2)
v3x = dim_avg_wgt_n_Wrap(v3(:,:,{25:40},:), w, 0, 2)
p1x = dim_avg_wgt_n_Wrap(p1(:,:,{25:40},:), w, 0, 2)
p2x = dim_avg_wgt_n_Wrap(p2(:,:,{25:40},:), w, 0, 2)
p3x = dim_avg_wgt_n_Wrap(p3(:,:,{25:40},:), w, 0, 2)

Season = (/"MAM","JJA","SON","DJF"/)
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

pathmask = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmask = "USGS-gtopo30_1.9x2.5_remap_c050602.nc"
notpmask = "USGS-gtopo30_1.9x2.5_remap_C20190410_TP1000m_wuyt.nc"
fh1 = addfile(pathmask+tpmask, "r")
fh2 = addfile(pathmask+notpmask, "r")

HGT1 = fh1->PHIS
HGT2 = fh2->PHIS
hgt1 = HGT1
hgt2 = HGT2
hgt1 = HGT1/9.8
hgt2 = HGT2/9.8
printVarSummary(hgt1)

hgt1x = dim_avg_n_Wrap(hgt1({25:40},:), 0)
hgt2x = dim_avg_n_Wrap(hgt2({25:40},:), 0)
printVarSummary(hgt1x)

tp1x = hgt1x
tp1x = 1013.25*(1-hgt1x*0.0065/288.15)^5.25145
tp2x = hgt2x
tp2x = 1013.25*(1-hgt2x*0.0065/288.15)^5.25145
printVarSummary(tp1x)


;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191116/"
figname = "Vertical_UW"

lonstr = 40
lonend = 180
latstr = -10
latend = 70

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 8
plot = new(nplot, graphic)
plotv = plot
plot_topo = plot

 xyres  = True
 xyres@gsnDraw  = False
 xyres@gsnFrame = False
 xyres@gsnYRefLine = 1000.0
 xyres@gsnBelowYRefLineColor = "black"
 ;xyres@xyLineColors = "gray"
 ;xyres@xyLineThicknesses = 1.
 xyres@trYMinF = 100.
 xyres@trYMaxF = 1000.
 ;xyres@trYReverse = True


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

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

;------------------------------------------
 res@trYMinF   = 100
 res@trYMaxF   = 1000
  
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

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
  res@vcRefAnnoFontHeightF = 0.03

;------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.0
res@cnMaxLevelValF = 5.0
res@cnMinLevelValF = -5.0
res@cnFillPalette = "BlueWhiteOrangeRed"
res@vcRefMagnitudeF         = 2.             ; define vector ref mag
res@lbLabelBarOn = False
res@gsnPresHgtHeightLabelOn = False

do i = 0, dimsizes(Season)-1

res@gsnRightString = Season(i)
res@gsnLeftString = "Thermal"
;plot(i) = gsn_csm_pres_hgt_vector(wks, p_h(i,:,{lonstr:lonend}), u_h(i,:,{lonstr:lonend}), v_h(i,:,{lonstr:lonend}), res)
plot(2*i) = gsn_csm_pres_hgt_vector(wks, p_h(i,:,{lonstr:lonend}), u_h(i,:,{lonstr:lonend}), v_h(i,:,{lonstr:lonend}), res)
;plot_topo(2*i) = gsn_csm_xy(wks, tp1x&lon, tp1x, xyres)
;overlay(plot(2*i), plot_topo(2*i))

res@gsnLeftString = "Dynamical"
;plot(i+4) = gsn_csm_pres_hgt_vector(wks, p_m(i,:,{lonstr:lonend}), u_m(i,:,{lonstr:lonend}), v_m(i,:,{lonstr:lonend}), res)
plot(2*i+1) = gsn_csm_pres_hgt_vector(wks, p_m(i,:,{lonstr:lonend}), u_m(i,:,{lonstr:lonend}), v_m(i,:,{lonstr:lonend}), res)
;plot_topo(2*i+1) = gsn_csm_xy(wks, tp1x&lon, tp1x, xyres)
;overlay(plot(2*i+1), plot_topo(2*i+1))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Vertical Cell & T (25-40N)"
pres@gsnPanelLabelBar  = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/4,2/), pres)
; ;================================================================================

; res@lbLabelBarOn = False

; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 10
; res@cnMaxLevelValF = 300
; res@cnMinLevelValF = 200
; res@cnFillPalette = "BlueWhiteOrangeRed"

; res@vcRefMagnitudeF         = 5.0             ; define vector ref mag

; do i = 0, dimsizes(Season)-1
; res@gsnRightString = Season(i)
; res@gsnLeftString = "CTRL"
; plot(i) = gsn_csm_pres_hgt_vector(wks, pp1(i,:,{lonstr:lonend}), uu1(i,:,{lonstr:lonend}), vv1(i,:,{lonstr:lonend}), res)
; ;plot_topo(i) = gsn_csm_xy(wks, tp1x&lon, tp1x, xyres)
; ;overlay(plot(i), plot_topo(i))


; res@gsnLeftString = "NOQ1"
; plot(i+4) = gsn_csm_pres_hgt_vector(wks, pp2(i,:,{lonstr:lonend}), uu2(i,:,{lonstr:lonend}), vv2(i,:,{lonstr:lonend}), res)
; ;plot_topo(i+4) = gsn_csm_xy(wks, tp1x&lon, tp1x, xyres)
; ;overlay(plot(i+4), plot_topo(i+4))


; res@gsnLeftString = "NOTP_NOQ1"
; plot(i+8) = gsn_csm_pres_hgt_vector(wks, pp3(i,:,{lonstr:lonend}), uu3(i,:,{lonstr:lonend}), vv3(i,:,{lonstr:lonend}), res)
; ;plot_topo(i+8) = gsn_csm_xy(wks, tp1x&lon, tp1x, xyres)
; ;overlay(plot(i+8), plot_topo(i+8))

; end do

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelMainString = "UW (25-40N)"
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = -0.01

; gsn_panel(wks, plot, (/3,4/), pres)
end