
begin

year = ispan(1979, 2013, 1)
ny = dimsizes(year)
Season = (/"MAM","JJA","SON","DJF"/)
ns = dimsizes(Season)

;---------------------------------------------
path = "./"
name = "Q1_TP80E_dtrend_standardize.txt"

x = asciiread(path+name, (/ny,ns/), "float")
x@_FillValue = -32767
x!0 = "year"
x&year = year
x!1 = "season"
x&season = Season
printVarSummary(x)

index = x(:,0) ;MAM
indlow = ind(index .le. -0.8)
indhig = ind(index .ge. 0.8)
print(indlow)
numhig = dimsizes(indhig)
numlow = dimsizes(indlow)
yearhig = year(indhig)
yearlow = year(indlow)
print(yearhig)
print(yearlow)
;---------------------------------------------

;===================Composite========================
diri = "/home/yangsong3/wuyt/data/"
fili = "q1.mon.mean.1979.2013.nc" ; NCEP_DOE 
f = addfile(diri+fili, "r")
q1 = f->q1(:,:,::-1,:)
printVarSummary(q1)

Cp = 1004.07 ; J/(Kg*K)
q1 = q1 * Cp ; units W/Kg

;----------------integration------------------
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
;----------------------------------------------
Q1clm = clmMonTLL(Q1)
Q1ano = calcMonAnomTLL(Q1, Q1clm)
q1Sea = month_to_season(Q1ano, "MAM")
printVarSummary(q1Sea)
q1Sea!0 = "year"
q1Sea&year = year
printVarSummary(q1Sea)

q1hig = q1Sea(indhig,:,:)
q1low = q1Sea(indlow,:,:)
printVarSummary(q1hig)

higvar = dim_variance_n_Wrap(q1hig,0)
lowvar = dim_variance_n_Wrap(q1low, 0)
higave = dim_avg_n_Wrap(q1hig, 0)
lowave = dim_avg_n_Wrap(q1low, 0)
prob   = ttest(higave, higvar, numhig, lowave, lowvar, numlow, False, False)

diff  = higave - lowave
copy_VarCoords(q1Sea(0,:,:), diff)
copy_VarCoords(q1Sea(0,:,:), prob)
printVarSummary(prob)
printVarSummary(diff)

;========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/climatology/"
figname = "draw_Comp_MAM_Q1"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 3
plot = new(nplot, graphic)
plot_topo = plot
plot_prob = plot

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
;------------------------------------------------------------

;----------------------------------------------------
; draw prob
;----------------------------------------------------
;------------------------
; draw prob
;------------------------
res3  = True
res3@gsnDraw   = False
res3@gsnFrame  = False

res3@cnFillOn  = False 
res3@cnLinesOn = False
res3@cnLineLabelsOn = False
res3@cnInfoLabelOn  = False
res3@cnFillColor = "green4"

res3@gsnRightString = ""

res3@cnLevelSelectionMode = "ManualLevels"
res3@cnMinLevelValF = 0.00
res3@cnMaxLevelValF = 0.10
res3@cnLevelSpacingF = 0.01

opt  = True
opt@gsnShadeFillType = "pattern"
opt@gsnShadeLow = 17
opt@gsnShadeFillScaleF = 0.8
opt@gsnShadeFillDotSizeF = 0.003

;----------------------------------------------------

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
res@gsnRightString = "W/m~S~2~N~"

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 10.
res@cnMaxLevelValF = 100.
res@cnMinLevelValF = -100.
res@gsnLeftString = "HIGH"
plot(0) = gsn_csm_contour_map(wks, higave, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
plot_prob(0) = gsn_csm_contour(wks, prob, res3)
plot_prob(0) = gsn_contour_shade(plot_prob(0), 0.05, -999., opt)
overlay(plot(0), plot_prob(0))

res@gsnLeftString = "LOW"
plot(1) = gsn_csm_contour_map(wks, lowave, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
plot_prob(1) = gsn_csm_contour(wks, prob, res3)
plot_prob(1) = gsn_contour_shade(plot_prob(1), 0.05, -999., opt)
overlay(plot(1), plot_prob(1))

res@gsnLeftString = "DIFF"
plot(2) = gsn_csm_contour_map(wks, diff, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot_prob(2) = gsn_csm_contour(wks, prob, res3)
plot_prob(2) = gsn_contour_shade(plot_prob(2), 0.05, -999., opt)
overlay(plot(2), plot_prob(2))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Composite MAM Q1 anomaly"
pres@gsnStringFontHeightF = 0.035
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/2,2/), pres)
;=================================================================
latstr = 20
latend = 50
lonstr = 60
lonend = 120
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 10.
res@cnMaxLevelValF = 100.
res@cnMinLevelValF = -100.
res@gsnLeftString = "HIGH"
plot(0) = gsn_csm_contour_map(wks, higave, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
plot_prob(0) = gsn_csm_contour(wks, prob, res3)
plot_prob(0) = gsn_contour_shade(plot_prob(0), 0.05, -999., opt)
overlay(plot(0), plot_prob(0))

res@gsnLeftString = "LOW"
plot(1) = gsn_csm_contour_map(wks, lowave, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
plot_prob(1) = gsn_csm_contour(wks, prob, res3)
plot_prob(1) = gsn_contour_shade(plot_prob(1), 0.05, -999., opt)
overlay(plot(1), plot_prob(1))

res@gsnLeftString = "DIFF"
plot(2) = gsn_csm_contour_map(wks, diff, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot_prob(2) = gsn_csm_contour(wks, prob, res3)
plot_prob(2) = gsn_contour_shade(plot_prob(2), 0.05, -999., opt)
overlay(plot(2), plot_prob(2))

gsn_panel(wks, plot, (/2,2/), pres)

end