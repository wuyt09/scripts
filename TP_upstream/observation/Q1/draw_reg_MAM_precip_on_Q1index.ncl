
begin

yrStrt = 1979
yrLast = 2013
year = ispan(yrStrt, yrLast, 1)
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
;---------------------------------------------
;===================Correlation========================
diri = "/home/yangsong3/data-observation/GPCP/"
fili = "precip.mon.mean.nc" 
f = addfile(diri+fili, "r")

TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

p = f->precip(iYYYY,:,:)
printVarSummary(p)

pSea = month_to_season(p, "MAM")
printVarSummary(pSea)
pSea!0 = "year"
pSea&year = year
printVarSummary(pSea)

psea = dtrend_msg_n(year, pSea, False, False, 0)
printVarSummary(psea)
;-------------------------------------------------------
Nr = ny

rc = regCoef_n(index, psea, 0, 0)
cc = escorc_n(index, psea, 0, 0)
prob = cc
prob = rtest(cc, Nr, 0)
copy_VarCoords(p(0,:,:), rc)
copy_VarCoords(p(0,:,:), prob)
printVarSummary(rc)

;========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/climatology/"
figname = "draw_reg_MAM_precip_on_Q1index"

latstr = 0
latend = 70
lonstr = -30
lonend = 150

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 1
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
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2

res@mpMaxLonF = lonend
res@mpMinLonF = lonstr
res@mpMaxLatF = latend
res@mpMinLatF = latstr

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

;res@lbLabelBarOn = False
res@cnFillPalette = "BlueWhiteOrangeRed"
;res@gsnRightString = "95%"

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.005
res@cnMaxLevelValF = 0.02
res@cnMinLevelValF = -0.02

res@tiMainString = "Reg MAM Precip on TP Q1 timeseries"
plot(0) = gsn_csm_contour_map(wks, rc, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
plot_prob(0) = gsn_csm_contour(wks, prob, res3)
plot_prob(0) = gsn_contour_shade(plot_prob(0), 0.05, -999., opt)
overlay(plot(0), plot_prob(0))

draw(plot)
end