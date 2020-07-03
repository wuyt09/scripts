
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601
nyear = 25

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"
SF = new((/ncase,25,96,144/), float)
SC = SF
OMEGA = SF
TS = SF

LEV = 500

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.PRECSC.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.PRECSL.1979-2005.nc"
    fili3 = casename(i)+".clm2.h0.FSNO.1979-2005.nc"
    fili4 = casename(i)+".cam.h0.OMEGA.1979-2005.nc"
    fili5 = casename(i)+".cam.h0.TREFHT.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")
    f5 = addfile(diri+fili5, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    SF1 = f1->PRECSC(iYYYY,:,:)
    SF2 = f2->PRECSL(iYYYY,:,:)
    SFX = SF1 
    SFX = SF1 + SF2
    SFX = SFX * 1000 * 86400.
    SFX@units = "mm/day"
    SF(i,:,:,:) = month_to_season(SFX, season)

    SCX = f3->FSNO(iYYYY,:,:)
    SC(i,:,:,:) = month_to_season(SCX, season)

    OMEGAX = f4->OMEGA(iYYYY,{LEV},:,:)
    OMEGA(i,:,:,:) = month_to_season(OMEGAX, season)

    TSX = f5->TREFHT(iYYYY,:,:)
    TSX = TSX - 273.15
    TS(i,:,:,:) = month_to_season(TSX, season)


delete([/diri,fili1,fili2,fili3,fili4,fili5,f1,f2,f3,f4,f5,SFX,SCX,OMEGAX,TSX/])
end do
printVarSummary(SF)

;=======================================================================
; ttest
;=======================================================================

sf_ave = dim_avg_n_Wrap(SF, 1)
sf_var = dim_variance_n_Wrap(SF, 1)
probsf = ttest(sf_ave(0,:,:), sf_var(0,:,:), nyear, sf_ave(1,:,:), sf_var(1,:,:), nyear, False, False)
copy_VarCoords(sf_ave(0,:,:), probsf)
printVarSummary(probsf)

sc_ave = dim_avg_n_Wrap(SC, 1)
sc_var = dim_variance_n_Wrap(SC, 1)
probsc = ttest(sc_ave(0,:,:), sc_var(0,:,:), nyear, sc_ave(1,:,:), sc_var(1,:,:), nyear, False, False)
copy_VarCoords(sc_ave(0,:,:), probsc)
printVarSummary(probsc)

w_ave = dim_avg_n_Wrap(OMEGA, 1)
w_var = dim_variance_n_Wrap(OMEGA, 1)
probw = ttest(w_ave(0,:,:), w_var(0,:,:), nyear, w_ave(1,:,:), w_var(1,:,:), nyear, False, False)
copy_VarCoords(w_ave(0,:,:), probw)
printVarSummary(probw)

ts_ave = dim_avg_n_Wrap(TS, 1)
ts_var = dim_variance_n_Wrap(TS, 1)
probts = ttest(ts_ave(0,:,:), ts_var(0,:,:), nyear, ts_ave(1,:,:), ts_var(1,:,:), nyear, False, False)
copy_VarCoords(ts_ave(0,:,:), probts)
printVarSummary(probts)


diffsf = sf_ave(0,:,:)
diffsf = sf_ave(1,:,:) - sf_ave(0,:,:)

diffsc = sc_ave(0,:,:)
diffsc = sc_ave(1,:,:) - sc_ave(0,:,:)

diffw = w_ave(0,:,:)
diffw = w_ave(1,:,:) - w_ave(0,:,:)

diffts = ts_ave(0,:,:)
diffts = ts_ave(1,:,:) - ts_ave(0,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_snowcover_snowfall_omega500_ts"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)
plot_topo = plot
plot_prob = new(4, graphic)

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
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
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035
;res@gsnRightString       = "K"

res@pmLabelBarOrthogonalPosF = 0.12
res@pmLabelBarWidthF = 0.7
res@lbLabelFontHeightF = 0.02
;----------------------------------------------------
res@cnFillPalette  = "CBR_wet"
res@gsnLeftString  = "Snowfall"
res@gsnRightString = "Ctrl"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.
res@cnMaxLevelValF = 5.
res@cnLevelSpacingF = 0.5
plot(0) = gsn_csm_contour_map(wks, sf_ave(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@cnFillPalette  = "precip_diff_12lev"
res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.05
plot(1) = gsn_csm_contour_map(wks, diffsf, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

plot_prob(0) = gsn_csm_contour(wks,probsf,res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),0.055,17)
overlay(plot(1), plot_prob(0))
;----------------------------------------------------
res@cnFillPalette  = "CBR_wet"
res@gsnLeftString  = "Snow Cover"
res@gsnRightString = "Ctrl"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.
res@cnMaxLevelValF = 1.
res@cnLevelSpacingF = 0.1
plot(2) = gsn_csm_contour_map(wks, sc_ave(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@cnFillPalette = "precip_diff_12lev"
res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
plot(3) = gsn_csm_contour_map(wks, diffsc, res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))

plot_prob(1) = gsn_csm_contour(wks,probsc,res3)
plot_prob(1) = ShadeLtContour(plot_prob(1),0.055,17)
overlay(plot(3), plot_prob(1))

;----------------------------------------------------
res@cnFillPalette  = "MPL_RdBu"
res@gsnLeftString  = "OMEGA"+LEV
res@gsnRightString = "Ctrl"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.2
res@cnMaxLevelValF = 0.2
res@cnLevelSpacingF = 0.05
plot(4) = gsn_csm_contour_map(wks, w_ave(0,:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))


res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.05
res@cnMaxLevelValF = 0.05
res@cnLevelSpacingF = 0.005
plot(5) = gsn_csm_contour_map(wks, diffw, res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

plot_prob(2) = gsn_csm_contour(wks,probw,res3)
plot_prob(2) = ShadeLtContour(plot_prob(2),0.055,17)
overlay(plot(5), plot_prob(2))
;----------------------------------------------------
res@cnFillPalette  = "WhiteBlue"
res@gsnSpreadColorEnd = 2
res@gsnSpreadColorStart = 253
res@gsnLeftString  = "TS"
res@gsnRightString = "Ctrl"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -30.
res@cnMaxLevelValF = 0.
res@cnLevelSpacingF = 3
plot(6) = gsn_csm_contour_map(wks, ts_ave(0,:,:), res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))

res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3.
res@cnMaxLevelValF = 0.
res@cnLevelSpacingF = 0.3
plot(7) = gsn_csm_contour_map(wks, diffts, res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7), plot_topo(7))


plot_prob(3) = gsn_csm_contour(wks,probts,res3)
plot_prob(3) = ShadeLtContour(plot_prob(3),0.055,17)
overlay(plot(7), plot_prob(3))
;----------------------------------------------------

plots = plot
do i = 0,3
    plots(i) = plot(2*i)
    plots(i+4) = plot(2*i+1)
end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean"
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

; gsn_panel(wks, plot, (/4,2/), pres)
gsn_panel(wks, plots, (/2,4/), pres)
end