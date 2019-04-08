
begin

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)
ny = dimsizes(year)
Season = (/"MAM"/)
ns = dimsizes(Season)

;---------------------------------------------
path = "./"
name = "TTT_TP_dtrend_standardize_ERAI.txt"

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
numhig = dimsizes(indhig)
numlow = dimsizes(indlow)
yearhig = year(indhig)
yearlow = year(indlow)
print(yearhig)
print(yearlow)
;---------------------------------------------

;===================Composite========================
diri = "/home/yangsong3/data-observation/GPCP/"
fili = "precip.mon.mean.nc" 
f = addfile(diri+fili, "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

p = f->precip(iYYYY,:,:)
printVarSummary(p)

;----------------------------------------------
pclm = clmMonTLL(p)
pano = calcMonAnomTLL(p, pclm)
pSea = pano(4::12,:,:)
printVarSummary(pSea)
pSea!0 = "year"
pSea&year = year
printVarSummary(pSea)

phig = pSea(indhig,:,:)
plow = pSea(indlow,:,:)
printVarSummary(phig)

higvar = dim_variance_n_Wrap(phig,0)
lowvar = dim_variance_n_Wrap(plow, 0)
higave = dim_avg_n_Wrap(phig, 0)
lowave = dim_avg_n_Wrap(plow, 0)
prob   = ttest(higave, higvar, numhig, lowave, lowvar, numlow, False, False)

diff  = higave - lowave
copy_VarCoords(pSea(0,:,:), diff)
copy_VarCoords(pSea(0,:,:), prob)
printVarSummary(prob)
printVarSummary(diff)

;========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/observation/TTT/"
figname = "draw_Comp_May_precip_ERAI"

latstr = 0
latend = 70
lonstr = -30
lonend = 150

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

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@lbLabelBarOn = False
res@cnFillPalette = "CBR_drywet"
res@gsnRightString = "mm/day"

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.3
res@cnMaxLevelValF = 3.
res@cnMinLevelValF = -3.
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
pres@gsnPanelMainString = "Composite May precip anomaly"
pres@gsnStringFontHeightF = 0.035
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/2,2/), pres)
;=================================================================
end