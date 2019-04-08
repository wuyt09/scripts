
begin

diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
fili = "q1.mon.mean.1979.2016.ERAI.nc" ; ERAI 

f = addfile(diri+fili, "r")
q1 = f->q1    ; units K/s
printVarSummary(q1)

Cp = 1004.07 ; J/(Kg*K)
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

Season = (/"MAM","JJA","SON","DJF"/)

q1Sea = month_to_seasonN(Q1, Season)
printVarSummary(q1Sea)
q1sea = dim_avg_n_Wrap(q1Sea, 1)
;q1sea = q1Sea(:,22,:,:)
printVarSummary(q1sea)    

;========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/climatology/"
figname = "draw_clim_season_Q1_ERAI_wuyt"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4
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
 res_topo@cnLineColor            = "gray20"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)


res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@lbLabelBarOn = False
res@cnFillPalette = "cmp_b2r"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 50.
res@cnMaxLevelValF = 400.
res@cnMinLevelValF = -400.

res@gsnRightString = "W/m~S~2~N~"
do i = 0, dimsizes(Season)-1
    res@gsnLeftString = Season(i)
    plot(i) = gsn_csm_contour_map(wks, q1sea(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1 climatology"
pres@gsnStringFontHeightF = 0.035
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/2,2/), pres)

latstr = 20
latend = 50
lonstr = 60
lonend = 120
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

do i = 0, dimsizes(Season)-1
    res@gsnLeftString = Season(i)
    plot(i) = gsn_csm_contour_map(wks, q1sea(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do
gsn_panel(wks, plot, (/2,2/), pres)

end