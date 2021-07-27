begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

yrStrt = 1951
yrLast = 2007

year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"


diri = "/home/yangsong3/wuyt/sysu/scripts/model_precip/APHRO/"
fili = "APHRO_MA_025deg_V1003R1.monmean.1951-2007.nc"
f = addfile(diri+fili, "r")

precip = f->precip
printVarSummary(precip)

prec = month_to_season(precip, season)
printVarSummary(prec)

rc = regCoef_n(year, prec, 0, 0)
copy_VarCoords(prec(0,:,:), rc)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

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

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "trend_pr_JJA_APHRO_MA"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnAddCyclic = False
res@gsnMaximize = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red"
; res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

; res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.01
res@cnMaxLevelValF = 0.1
res@cnMinLevelValF = -0.1
; res@cnLevelSpacingF = 1.
; res@cnMaxLevelValF = 10
; res@cnMinLevelValF = -10
res@cnFillPalette = "MPL_BrBG"
res@gsnRightString = "mm/day/yr"


res@gsnLeftString = "APHRO"
plot(0) = gsn_csm_contour_map(wks, rc, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

draw(plot(0))
end