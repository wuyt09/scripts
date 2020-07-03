
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


Season = (/"MAM","JJA","SON","DJF"/)

uu11 = month_to_seasonN(u1, Season)
uu1 = dim_avg_n_Wrap(uu11, 1)
uu22 = month_to_seasonN(u2, Season)
uu2 = dim_avg_n_Wrap(uu22, 1)
uu33 = month_to_seasonN(u3, Season)
uu3 = dim_avg_n_Wrap(uu33, 1)
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

u_m = mask(u_m, tpmask.ge.0.9, False)
u_h = mask(u_h, tpmask.ge.0.9, False)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191116/"
figname = "TS"

lonstr = 0 ;20
lonend = 360 ;220
latstr = 10
latend = 80

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

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.
res@cnFillPalette = "temp_diff_18lev"

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "Thermal"
plot(2*i) = gsn_csm_contour_map(wks, u_h(i,:,:), res)
plot_topo(2*i) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2*i), plot_topo(2*i))

res@gsnLeftString = "Dynamical"
plot(2*i+1) = gsn_csm_contour_map(wks, u_m(i,:,:), res)
plot_topo(2*i+1) = gsn_csm_contour(wks, tpmask, res_topo)
overlay(plot(2*i+1), plot_topo(2*i+1))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "TS"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/4,2/), pres)

end