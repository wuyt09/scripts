
begin

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

filic1 = "F2000_CTRL.plevel.T.mon.0601-3012.nc"
filic2 = "F2000_NOQ1.plevel.T.mon.0601-3012.nc"
filic3 = "F2000_NOTP_NOQ1.plevel.T.mon.0601-3012.nc"

fc1 = addfile(diri1+filic1, "r")
fc2 = addfile(diri2+filic2, "r")
fc3 = addfile(diri3+filic3, "r")

c1 = fc1->T
c2 = fc2->T
c3 = fc3->T
printVarSummary(c1)

lev1 = c1&lev*100
lev1!0 = "lev"
lev1&lev = c1&lev
lev1@units = "Pa"
psfc = max(lev1)
ptop = min(lev1)
;lev1 = c1&lev
;lev1 = c1&lev

dP = dpres_plevel_Wrap(lev1, psfc, ptop, 0)
dp = conform_dims(dimsizes(c1), dP, 1)
c11 = c1*dp
c22 = c2*dp
c33 = c3*dp
copy_VarCoords(c1, c11)
copy_VarCoords(c2, c22)
copy_VarCoords(c3, c33)
p1 = dim_sum_n_Wrap(c11(:,{500:200},:,:), 1)/sum(dP({500:200}))
p2 = dim_sum_n_Wrap(c22(:,{500:200},:,:), 1)/sum(dP({500:200}))
p3 = dim_sum_n_Wrap(c33(:,{500:200},:,:), 1)/sum(dP({500:200}))
copy_VarCoords(c1(:,0,:,:), p1)
copy_VarCoords(c2(:,0,:,:), p2)
copy_VarCoords(c3(:,0,:,:), p3)
printVarSummary(p1)
;==================================================================

Season = (/"DJF","MAM","JJA","SON"/)
p11 = month_to_seasonN(p1, Season)
p22 = month_to_seasonN(p2, Season)
p33 = month_to_seasonN(p3, Season)

pp1 = dim_avg_n_Wrap(p11, 1)
pp2 = dim_avg_n_Wrap(p22, 1)
pp3 = dim_avg_n_Wrap(p33, 1)
printVarSummary(pp1)

p_h = pp1 - pp2
p_m = pp2 - pp3
copy_VarCoords(pp1, p_m)
copy_VarCoords(pp1, p_h)


tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_season_TropoT"

lonstr = 40 ;20
lonend = 160 ;220
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
res@mpGeophysicalLineColor = "grey45" ;"red"
res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = 180  ;(lonstr+lonend)/2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
;------------------------------------------
;------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.
res@cnFillPalette = "BlueWhiteOrangeRed"
res@lbLabelBarOn = False

do i = 0, dimsizes(Season)-1

res@gsnRightString = Season(i)
res@gsnLeftString = "CTRL - TP&NOQ1"
plot(2*i) = gsn_csm_contour_map(wks, p_h(i,:,:), res)
plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2*i), plot_topo(2*i))

res@gsnLeftString = "TP&NOQ1 - NOTP&NOQ1"
plot(2*i+1) = gsn_csm_contour_map(wks, p_m(i,:,:), res)
plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2*i+1), plot_topo(2*i+1))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "TropoT_difference"
pres@gsnPanelLabelBar  = True
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/4,2/), pres)
;draw(plot(6))
;frame(wks)
;================================================================================
delete([/plot,plotv,plot_topo/])
nplot = 12
plot = new(nplot, graphic)
plotv = plot
plot_topo = plot


res@lbLabelBarOn = False

;res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 5.
res@cnMaxLevelValF = 250.
res@cnMinLevelValF = 210.
res@cnFillPalette = "BlueWhiteOrangeRed"


do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)
res@gsnLeftString = "CTRL"
plot(i) = gsn_csm_contour_map(wks, pp1(i,:,:), res)
plot_topo(i) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(i), plot_topo(i))

res@gsnLeftString = "NOTP"
plot(i+4) = gsn_csm_contour_map(wks, pp2(i,:,:), res)
plot_topo(i+4) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(i+4), plot_topo(i+4))

res@gsnLeftString = "NOQ1"
plot(i+8) = gsn_csm_contour_map(wks, pp3(i,:,:), res)
plot_topo(i+8) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(i+8), plot_topo(i+8))
end do
pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "TropoT_climatology"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)
end