
begin

runs = (/"CTRL", "NOSH","NOTP_NOSH"/)
nm = dimsizes(runs)
uu = new((/nm,96,144/), float)
uu!0 = "run"
uu&run = runs
vv = uu
pp = uu
; pp = new((/nm,37,96,144/), float)
; pp!0 = "run"
; pp&run = runs

do i = 0, nm-1

  diri = "/home/yangsong3/data-model/wuyt/TP_"+runs(i)+"_B2000_F19G16_CAM4/pro/"

  filiu = "TP_"+runs(i)+".plevel.U.mon.26001-30712.nc"
  filiv = "TP_"+runs(i)+".plevel.V.mon.26001-30712.nc"
  filiz = "TP_"+runs(i)+".plevel.Z3.mon.26001-30712.nc"
  ; filip = "TP_"+runs(i)+".PS.mon.26001-30712.nc"
  ; filih = "TP_"+runs(i)+".plevel.Z3.mon.26001-30712.nc"

  fu = addfile(diri+filiu, "r")
  fv = addfile(diri+filiv, "r")
  fc = addfile(diri+filiz, "r")
  u = fu->U(:,{500},:,:)
  v = fv->V(:,{500},:,:)
  p = fc->Z3(:,{500},:,:)
  ;================================================================

  Season = "JJA"
  pp1 = month_to_season(p, Season)
  pp(i,:,:) = dim_avg_n_Wrap(pp1, 0)

  uu1 = month_to_season(u, Season)
  uu(i,:,:) = dim_avg_n_Wrap(uu1, 0)

  vv1 = month_to_season(v, Season)
  vv(i,:,:) = dim_avg_n_Wrap(vv1, 0)
end do
printVarSummary(uu)

p_m = pp(1,:,:) - pp(2,:,:)
p_h = pp(0,:,:) - pp(1,:,:)
copy_VarCoords(pp(0,:,:), p_m)
copy_VarCoords(pp(0,:,:), p_h)

u_m = uu(1,:,:) - uu(2,:,:)
u_h = uu(0,:,:) - uu(1,:,:)
copy_VarCoords(uu(0,:,:), u_m)
copy_VarCoords(uu(0,:,:), u_h)

v_m = vv(1,:,:) - vv(2,:,:)
v_h = vv(0,:,:) - vv(1,:,:)
copy_VarCoords(vv(0,:,:), v_m)
copy_VarCoords(vv(0,:,:), v_h)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_B/model_need_input_data/"
tpmaskname = "TPmask_C20190304.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/20190924_CTRL_NOSH_NOTPNOSH/"
figname = "draw_clim_"+Season+"_CTRL_S1_S2_UV500_Hgt500"

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


res@gsnLeftString = runs(0)
plot(2) = gsn_csm_contour_map(wks, pp(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(2) = gsn_csm_vector(wks, uu1, vv1, vcres)
overlay(plot(2), plot_topo(2))
overlay(plot(2), plotv(2))

res@gsnLeftString = runs(1)
plot(3) = gsn_csm_contour_map(wks, pp(1,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
plotv(3) = gsn_csm_vector(wks, uu2, vv2, vcres)
overlay(plot(3), plot_topo(3))
overlay(plot(3), plotv(3))

res@gsnLeftString = runs(2)
plot(4) = gsn_csm_contour_map(wks, pp(2,:,:), res)
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
