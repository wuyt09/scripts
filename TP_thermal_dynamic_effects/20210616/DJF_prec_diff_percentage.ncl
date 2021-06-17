
begin

nlat = 96
nlon = 144
season = "DJF"

lonstr = -20
lonend = 160
latstr = 0
latend = 70

casename = (/"CTRL","NOQ1", "NOTP","NOTP_NOQ1"/) ;
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)

z = new((/ncase,nyear,nlat,nlon/), float)

do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
  zname = "F2000_"+casename(i)+".PREC.mon.0601-3012.nc"
  fz = addfile(diri+zname, "r")
  Z = fz->PREC
  Z = Z*1000*86400
  Z@units = "mm/day"

  z(i,:,:,:) = month_to_season(Z, season)
end do

z!0 = "case"
z&case = casename
z&time = year
printVarSummary(z)

Clim = dim_avg_n_Wrap(z(0,:,:,:), 0)

;=========================ttest================================================
var = dim_variance_n_Wrap(z, 1)
ave = dim_avg_n_Wrap(z, 1)

sig = 0.05
prob_h = ttest(ave(0,:,:), var(0,:,:), nyear, ave(1,:,:), var(1,:,:), nyear, False, False)
copy_VarCoords(ave(0,:,:), prob_h)
prob_d = ttest(ave(0,:,:), var(0,:,:), nyear, ave(2,:,:), var(2,:,:), nyear, False, False)
copy_VarCoords(ave(0,:,:), prob_d)
prob_t = ttest(ave(0,:,:), var(0,:,:), nyear, ave(3,:,:), var(3,:,:), nyear, False, False)
copy_VarCoords(ave(0,:,:), prob_t)

h = ave(0,:,:) - ave(1,:,:)
copy_VarCoords(ave(0,:,:), h)
d = ave(0,:,:) - ave(2,:,:)
copy_VarCoords(ave(0,:,:), d)
hd = h+d
copy_VarCoords(ave(0,:,:), hd)
t = ave(0,:,:) - ave(3,:,:)
copy_VarCoords(ave(0,:,:), t)
diff = hd - t
copy_VarCoords(ave(0,:,:), diff)
printVarSummary(h)

h_percent = h/ave(1,:,:)*100
d_percent = d/ave(2,:,:)*100
t_percent = t/ave(3,:,:)*100
copy_VarCoords(h, h_percent)
copy_VarCoords(d, d_percent)
copy_VarCoords(t, t_percent)
printVarSummary(h_percent)


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)

;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_thermal_dynamic_effects/20210616/"
figname = season+"_prec_diff_percentage"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(9, graphic)
plot_topo = plot
plot_prob = new(4, graphic)

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "forestgreen"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

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
 res3@cnFillColor  =  "purple"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.002
 res3@cnFillScaleF  = 0.8
 ;----------------------------------------------
 ;----------------------------------------------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 3
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

; res@lbLabelBarOn = False

res@pmLabelBarHeightF = 0.1
res@pmLabelBarWidthF = 0.7
res@pmLabelBarOrthogonalPosF = 0.2

res@gsnRightString = "mm/day"

res@cnFillPalette = "CBR_wet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 1
res@cnMaxLevelValF = 10
res@cnLevelSpacingF = 1.
res@gsnLeftString = "(a) CTRL Clim"
plot(0) = gsn_csm_contour_map(wks, Clim, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))


res@cnFillPalette = "CBR_drywet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 0.5
res@gsnLeftString = "(b) Thermal effects"
plot(1) = gsn_csm_contour_map(wks, h, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
plot_prob(0) = gsn_csm_contour(wks,prob_h,res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),sig,17)
overlay(plot(1), plot_prob(0))

res@gsnLeftString = "(c) Dynamical effects"
plot(2) = gsn_csm_contour_map(wks, d, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))
plot_prob(1) = gsn_csm_contour(wks,prob_d,res3)
plot_prob(1) = ShadeLtContour(plot_prob(1),sig,17)
overlay(plot(2), plot_prob(1))

res@gsnLeftString = "(d) b+c"
plot(3) = gsn_csm_contour_map(wks, hd, res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

res@gsnLeftString = "(e) Total effects"
plot(4) = gsn_csm_contour_map(wks, t, res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))
plot_prob(2) = gsn_csm_contour(wks,prob_t,res3)
plot_prob(2) = ShadeLtContour(plot_prob(2),sig,17)
overlay(plot(4), plot_prob(2))

res@gsnLeftString = "(f) e-d"
plot(5) = gsn_csm_contour_map(wks, diff, res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))


res@cnFillPalette = "precip4_diff_19lev"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -100
res@cnMaxLevelValF = 100
res@cnLevelSpacingF = 10
res@gsnRightString = "%"
res@gsnLeftString ="(g) Thermal effects percentage"
plot(6) = gsn_csm_contour_map(wks, h_percent, res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))

res@gsnLeftString = "(h) Dynamical  effects percentage"
plot(7) = gsn_csm_contour_map(wks, d_percent, res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7), plot_topo(7))

res@gsnLeftString = "(i) Total effects percentage"
plot(8) = gsn_csm_contour_map(wks, t_percent, res)
plot_topo(8) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(8), plot_topo(8))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = season+" mean prec"
gsn_panel(wks, plot, (/3,3/), pres)

end
