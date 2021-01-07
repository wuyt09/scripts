
begin

yrStrt = 1979
yrLast = 2019
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = (/"MAM"/)

fx = "time_series_TP_skt_DJF_MAM_ERA5.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(:,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))

; hig_ind = hig_ind+1
; low_ind = low_ind+1

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "snowmelt_mon_1x1_1979_2019.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
var = short2flt(f->smlt(iYYYY,:,:))
var = var*1000.
var@units = "mm of water equivalent"
printVarSummary(var)

var_sea = month_to_season(var, season)
printVarSummary(var_sea)

var_hig = var_sea(hig_ind,:,:)
var_low = var_sea(low_ind,:,:)

;=============ttest====================================

var_hig_ave = dim_avg_n_Wrap(var_hig, 0)
var_low_ave = dim_avg_n_Wrap(var_low, 0)

var_hig_var = dim_variance_n_Wrap(var_hig, 0)
var_low_var = dim_variance_n_Wrap(var_low, 0)
printVarSummary(var_hig_ave)

prob = var_hig_ave
prob = ttest(var_hig_ave, var_hig_var, dimsizes(hig_ind), var_low_ave, var_low_var, dimsizes(low_ind), False, 0)
printVarSummary(prob)

diff_var = var_hig_ave
diff_var = var_hig_ave - var_low_ave
printVarSummary(diff_var)

;===================================================================
; plot
;===================================================================
lonstr = 60
lonend = 110
latstr = 15
latend = 45

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Comp_"+season+"_SnowMelt_by_TP_skt_DJF"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)
plot_topo = plot
plot_prob = new(1, graphic)

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey40"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02


;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------


res@cnFillPalette = "sunshine_9lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 1.
res@cnMaxLevelValF = 5.
res@cnLevelSpacingF = 0.5
; res@lbLabelBarOn = False

res@gsnLeftString = "(a) High DJF TP skt"
plot(0) = gsn_csm_contour_map(wks, var_hig_ave, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "(b) Low DJF TP skt"
plot(1) = gsn_csm_contour_map(wks, var_low_ave, res)
plot(1) = gsn_csm_contour_map(wks, var_hig_ave, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

res@cnFillPalette = "BlueDarkRed18"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1.
res@cnMaxLevelValF = 1.
res@cnLevelSpacingF = 0.1
res@lbLabelBarOn = True
res@gsnLeftString = "(c) a-b"
plot(2) = gsn_csm_contour_map(wks, diff_var, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

plot_prob(0) = gsn_csm_contour(wks,prob,res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),0.05,6)
overlay(plot(2), plot_prob(0))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Sonw Melt in "+season
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = -0.015
gsn_panel(wks, plot, (/1,3/), pres)


end