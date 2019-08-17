
begin
Cp = 1004.07 ; J/(Kg*K)
Season = (/"MAM","JJA","SON","DJF"/)
lonstr = 40
lonend = 180
latstr = -10
latend = 70

;---------------------------------------------------------
diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
fili = "q1.mon.mean.1979.2016.ERAI.nc" ; ERAI 

f = addfile(diri+fili, "r")
q1 = f->q1(:,:,{latstr:latend},{lonstr:lonend})    ; units K/s
printVarSummary(q1)

q1 = q1 * Cp ; units W/Kg

;----------------integration-----------------------------
plev = q1&lev*100.
psfc = max(plev)
ptop = min(plev)
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
printVarSummary(dP)

dp = conform_dims(dimsizes(q1), dP, 1)
q1dp = q1 * dp
copy_VarCoords(q1, q1dp)
printVarSummary(q1dp)
Q1 = dim_sum_n_Wrap(q1dp(:,{1000:1},:,:), 1)
copy_VarCoords(q1(:,0,:,:), Q1)
printVarSummary(Q1)
Q1 = Q1/9.8  ; units W/m2
;---------------------------------------------------------

q1Sea = month_to_seasonN(Q1, Season)
printVarSummary(q1Sea)
q1sea = dim_avg_n_Wrap(q1Sea, 1)
;q1sea = q1Sea(:,22,:,:)
printVarSummary(q1sea)
qsea_obs = q1sea    
delete([/q1,plev,psfc,ptop,dP,dp,q1dp,Q1,q1Sea,q1sea/])

;=====================================================================

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
filiu1 = "F2000_CTRL.plevel.DTCOND.mon.0601-3012.nc"
filiv1 = "F2000_CTRL.plevel.DTV.mon.0601-3012.nc"
filis1 = "F2000_CTRL.plevel.QRS.mon.0601-3012.nc"
filil1 = "F2000_CTRL.plevel.QRL.mon.0601-3012.nc"

fu1 = addfile(diri1+filiu1, "r")
fv1 = addfile(diri1+filiv1, "r")
fs1 = addfile(diri1+filis1, "r")
fl1 = addfile(diri1+filil1, "r")

u1 = fu1->DTCOND(:,:,{latstr:latend},{lonstr:lonend})
v1 = fv1->DTV(:,:,{latstr:latend},{lonstr:lonend})
s1 = fs1->QRS(:,:,{latstr:latend},{lonstr:lonend})
l1 = fl1->QRL(:,:,{latstr:latend},{lonstr:lonend})

q1 = u1
q1 = u1+v1+s1+l1

q1 = q1 * Cp ; units W/Kg

;----------------integration-----------------------------
plev = q1&lev*100.
psfc = max(plev)
ptop = min(plev)
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
printVarSummary(dP)

dp = conform_dims(dimsizes(q1), dP, 1)
q1dp = q1 * dp
copy_VarCoords(q1, q1dp)
printVarSummary(q1dp)
Q1 = dim_sum_n_Wrap(q1dp(:,{1000:1},:,:), 1)
copy_VarCoords(q1(:,0,:,:), Q1)
printVarSummary(Q1)
Q1 = Q1/9.8  ; units W/m2
;---------------------------------------------------------

q1Sea = month_to_seasonN(Q1, Season)
printVarSummary(q1Sea)
q1sea = dim_avg_n_Wrap(q1Sea, 1)
;q1sea = q1Sea(:,22,:,:)
printVarSummary(q1sea)
qsea_ctl = q1sea    
delete([/q1,plev,psfc,ptop,dP,dp,q1dp,Q1,q1Sea,q1sea/])

;=====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_season_Q1_OBS_CTRL"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 8
plot = new(nplot, graphic)
plot_topo = plot

 ;------------------------------------------------------------
 ; TOPO
 ;------------------------------------------------------------
 
 f2 = "~/wuyt/data/TOPO_1500.nc"
 topo_in = addfile(f2,"r")
 topo = topo_in->TP

 xlat = topo_in->LATIXY(:,0)
 xlon = topo_in->LONGXY(0,:)
 topo!0 = "lsmlat"
 topo!1 = "lsmlon"
 topo&lsmlat = xlat
 topo&lsmlon = xlon
 
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
 res_topo@cnLevels               = (/1500/)


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnAddCyclic = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2

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

res@lbLabelBarOn = False
res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 30.
res@cnMaxLevelValF = 300.
res@cnMinLevelValF = -300.

res@gsnRightString = "W/m~S~2~N~"
do i = 0, dimsizes(Season)-1
    res@gsnRightString = Season(i)
    res@gsnLeftString = "ERAI"
    ;plot(2*i) = gsn_csm_contour_map(wks, qsea_obs(i,:,:), res)
    ;plot_topo(2*i) = gsn_csm_contour(wks, topo, res_topo)
    ;overlay(plot(2*i), plot_topo(2*i))
    ;res@gsnLeftString = "CTRL"
    ;plot(2*i+1) = gsn_csm_contour_map(wks, qsea_ctl(i,:,:), res)
    ;plot_topo(2*i+1) = gsn_csm_contour(wks, topo, res_topo)
    ;overlay(plot(2*i+1), plot_topo(2*i+1))

    res@gsnLeftString = "ERAI"
    plot(i) = gsn_csm_contour_map(wks, qsea_obs(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    res@gsnLeftString = "CTRL"
    plot(i+4) = gsn_csm_contour_map(wks, qsea_ctl(i,:,:), res)
    plot_topo(i+4) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i+4), plot_topo(i+4))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1 climatology"
pres@gsnStringFontHeightF = 0.035
pres@gsnPanelLabelBar = True
pres@pmLabelBarWidthF = 0.6
pres@pmLabelBarOrthogonalPosF = -0.005
pres@pmLabelBarParallelPosF = 0.01
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,4/), pres)

end