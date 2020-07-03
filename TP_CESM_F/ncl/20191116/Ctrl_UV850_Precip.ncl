
begin
LEV = 850

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

filiu1 = "F2000_CTRL.plevel.U.mon.0601-3012.nc"
filiu2 = "F2000_NOQ1.plevel.U.mon.0601-3012.nc"
filiu3 = "F2000_NOTP_NOQ1.plevel.U.mon.0601-3012.nc"

filiv1 = "F2000_CTRL.plevel.V.mon.0601-3012.nc"
filiv2 = "F2000_NOQ1.plevel.V.mon.0601-3012.nc"
filiv3 = "F2000_NOTP_NOQ1.plevel.V.mon.0601-3012.nc"

filic1 = "F2000_CTRL.PRECC.mon.0601-3012.nc"
filic2 = "F2000_NOQ1.PRECC.mon.0601-3012.nc"
filic3 = "F2000_NOTP_NOQ1.PRECC.mon.0601-3012.nc"

filil1 = "F2000_CTRL.PRECL.mon.0601-3012.nc"
filil2 = "F2000_NOQ1.PRECL.mon.0601-3012.nc"
filil3 = "F2000_NOTP_NOQ1.PRECL.mon.0601-3012.nc"

fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")

U1 = fu1->U(:,{LEV},:,:)
; u2 = fu2->U
; u3 = fu3->U
V1 = fv1->V(:,{LEV},:,:)
; v2 = fv2->V
; v3 = fv3->V
uzonal = dim_avg_n_Wrap(U1, 1)
u1 = U1
u1 = U1 - conform_dims(dimsizes(U1), uzonal, (/0,2/))
vzonal = dim_avg_n_Wrap(V1, 1)
v1 = V1
v1 = V1 - conform_dims(dimsizes(V1), vzonal, (/0,2/))
printVarSummary(u1)
;==================================================================

Season = (/"MAM","JJA","SON","DJF"/)

uu11 = month_to_seasonN(u1, Season)
uu1 = dim_avg_n_Wrap(uu11, 1)
; uu22 = month_to_seasonN(u2, Season)
; uu2 = dim_avg_n_Wrap(uu22, 1)
; uu33 = month_to_seasonN(u3, Season)
; uu3 = dim_avg_n_Wrap(uu33, 1)

vv11 = month_to_seasonN(v1, Season)
vv1 = dim_avg_n_Wrap(vv11, 1)
; vv22 = month_to_seasonN(v2, Season)
; vv2 = dim_avg_n_Wrap(vv22, 1)
; vv33 = month_to_seasonN(v3, Season)
; vv3 = dim_avg_n_Wrap(vv33, 1)
printVarSummary(uu1)

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

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191116/"
figname = "Ctrl_UV850_Precip"

lonstr = 50
lonend = 130
latstr = 0
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
 res_topo@cnLineColor            = "blue"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

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

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF =  -0.6 ;-1.0             ; move ref vector up

  ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.02
  vcres@vcRefLengthF            = 0.025            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.015
;------------------------------------------
vcres@vcRefMagnitudeF         = 10.0  ; define vector ref mag
res@lbLabelBarOn = False
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 0.5
; res@cnMaxLevelValF = 5.
; res@cnMinLevelValF = -5.
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/-10,-5,-3,-2,-1,-0.5,0.5,1,2,3,5,10/)
res@cnFillPalette = "CBR_drywet"

; do i = 0, dimsizes(Season)-1

; res@gsnRightString = Season(i)
; res@gsnLeftString = "Thermal"
; plot(2*i) = gsn_csm_contour_map(wks, p_h(i,:,:), res)
; plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
; plotv(2*i) = gsn_csm_vector(wks, u_h(i,{850},:,:), v_h(i,{850},:,:), vcres)
; overlay(plot(2*i), plot_topo(2*i))
; overlay(plot(2*i), plotv(2*i))

; res@gsnLeftString = "Dynamical"
; plot(2*i+1) = gsn_csm_contour_map(wks, p_m(i,:,:), res)
; plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
; plotv(2*i+1) = gsn_csm_vector(wks, u_m(i,{850},:,:), v_m(i,{850},:,:), vcres)
; overlay(plot(2*i+1), plot_topo(2*i+1))
; overlay(plot(2*i+1), plotv(2*i+1))

; end do

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelMainString = "UV850_Precip"
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = -0.01
; pres@pmLabelBarWidthF = 0.6
; pres@lbLabelFontHeightF = 0.015

; gsn_panel(wks, plot, (/4,2/), pres)

plot(0) = gsn_csm_map_ce(wks, res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(0), plot_topo(0))

vcres@gsnLeftString = "JJA"
;plotv(0) = gsn_csm_streamline(wks, uu1(0,:,:), vv1(0,:,:), vcres)
plotv(0) = gsn_csm_vector(wks, uu1(1,:,:), vv1(1,:,:), vcres)
overlay(plot(0), plotv(0))
draw(plot(0))
frame(wks)

plot(1) = gsn_csm_map_ce(wks, res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(1), plot_topo(1))

vcres@gsnLeftString = "DJF"
plotv(1) = gsn_csm_vector(wks, uu1(3,:,:), vv1(3,:,:), vcres)
overlay(plot(1), plotv(1))
draw(plot(1))
end