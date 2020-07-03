
begin

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

filiu1 = "F2000_CTRL.TS.mon.0601-3012.nc"
filiu2 = "F2000_NOQ1.TS.mon.0601-3012.nc"
filiu3 = "F2000_NOTP_NOQ1.TS.mon.0601-3012.nc"


fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

u1 = fu1->TS
u2 = fu2->TS
u3 = fu3->TS
printVarSummary(u1)


;Season = (/"DJF","MAM","JJA","SON"/)
Season = "JJA"
uu11 = month_to_season(u1, Season)
uu1 = dim_avg_n_Wrap(uu11, 0)
uu22 = month_to_season(u2, Season)
uu2 = dim_avg_n_Wrap(uu22, 0)
uu33 = month_to_season(u3, Season)
uu3 = dim_avg_n_Wrap(uu33, 0)
printVarSummary(uu1)

u_m = uu2 - uu3
u_h = uu1 - uu2
copy_VarCoords(uu1, u_m)
copy_VarCoords(uu1, u_h)

;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/JJA/"
figname = "draw_clim_season_TS"

lonstr = 40
lonend = 160
latstr = -10
latend = 70

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
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
res@mpGeophysicalLineColor = "grey25" ;"grey25"
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
  vcres@vcRefLengthF            = 0.025            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------
vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 5.
res@cnMaxLevelValF = 320.
res@cnMinLevelValF = 240.
res@cnFillPalette = "temp_19lev"



res@gsnLeftString = "CTRL"
plot(0) = gsn_csm_contour_map(wks, uu1, res)
plot_topo(0) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "NOQ1"
plot(1) = gsn_csm_contour_map(wks, uu2, res)
plot_topo(1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(1), plot_topo(1))

res@gsnLeftString = "NOTP_NOQ1"
plot(2) = gsn_csm_contour_map(wks, uu3, res)
plot_topo(2) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2), plot_topo(2))


pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = Season+"_TS_climatology"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot(0:2), (/1,3/), pres)


;------------------------------------------
vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
;vcres@vcRefLengthF            = 0.035            ; define length of vec ref
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 10.
res@cnMinLevelValF = -10.
res@cnFillPalette = "temp_diff_18lev"


res@gsnLeftString = "thermal"
plot(3) = gsn_csm_contour_map(wks, u_h, res)
plot_topo(3) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(3), plot_topo(3))

res@gsnLeftString = "mechanical"
plot(4) = gsn_csm_contour_map(wks, u_m, res)
plot_topo(4) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(4), plot_topo(4))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = Season+"_TS_difference"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot(3:4), (/1,2/), pres)

end