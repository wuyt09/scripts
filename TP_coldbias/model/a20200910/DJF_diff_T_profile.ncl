
begin

lonstr = 60
lonend = 150
latstr = 15
latend = 50


yrStrt = 198102
yrLast = 200601
nyear = 25

casename = (/"Ctrl_FAMIP","Hist_FAMIP_CMIP6","Hist_PO"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)

season = "DJF"
nlat = 96
nlon = 144
nlev = 17

LAT = 35


t = new((/ncase,nyear,nlev,nlon/), float) ;,nlat
t!0 = "case"
t&case = casename

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"

    fili1 = casename(i)+".cam.h0.T.1979-2005.nc"
    ; fili2 = casename(i)+".cam.h0.U.1979-2005.nc"
    f1 = addfile(diri+fili1, "r")
    ; f2 = addfile(diri+fili2, "r")
    ; f3 = addfile(diri+fili3, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    T = f1->T(iYYYY,:,:,:)
    tx = month_to_season(T, season)
    ; tx = mask(tx, TOPO.ge.1500, True)
    ; t(i,:,:) = wgt_areaave_Wrap(tx(:,:,{25:40},{75:105}), 1, 1, 0)
    t(i,:,:,:) = tx(:,:,{LAT},:)
end do

printVarSummary(t)

t_ave = dim_avg_n_Wrap(t, 1)
t_var = dim_variance_n_Wrap(t, 1)
t_diff = t(0:1,0,:,:)
prob = t(0:1,0,:,:)
do i = 0, 1
    t_diff(i,:,:) = t_ave(i+1,:,:) - t_ave(0,:,:)
    prob(i,:,:) = ttest(t_ave(0,:,:), t_var(0,:,:), nyear, t_ave(i+1,:,:), t_var(i+1,:,:), nyear, False, False)
end do
printVarSummary(prob)


;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/"
figname = season+"_diff_T_profile"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(2, graphic)
plot_prob = plot
plot_topo = plot

title = (/"All Basin", "Pacific Basin" /)
;=====================================================
    res_topo                   = True

    res_topo@gsnDraw  = False           ;先别直接画出来
    res_topo@gsnFrame = False           ;也不要翻页

    res_topo@gsnLeftString  = " "              ;不要标示出变量名和单位等信息
    res_topo@gsnRightString = " "

    ;res_topo@gsnMaximize        = True
    ;res_topo@trYMaxF            = 1000.0
    res_topo@gsnYRefLine        = 1000 ;(/100,1000/)          ; create a reference line 
    res_topo@gsnBelowYRefLineColor = "black" ;(/"black","black"/) ; above ref line fill grey
    ;res_topo@gsnXYBelowFillColors = "goldenrod2"
    res_topo@gsnAboveYRefLineColor = "blue"
    res_topo@xyLineColors = "black"
    res_topo@xyLineThicknesses = 4.
;===================================================================
;===================================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.7

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@gsnStringFontHeightF = 0.03
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@pmLabelBarOrthogonalPosF = -0.06
    ; res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@tiYAxisString = ""

    res@tmYLMode      = "Explicit"
    res@tmYLValues    = (/"1000","850","750","500","300","250"/) ;,"150"
    res@tmYLLabels    = (/"1000","850","750","500","300","250"/) ;,"150"

    res@trYMinF   = 250
    res@trYMaxF   = 1000
    res@gsnPresHgtHeightLabelOn = False
    res@tmYROn = False

;=====================================================

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

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.2
res@cnMaxLevelValF = 2.
res@cnMinLevelValF = -2.
res@cnFillPalette = "temp_diff_18lev" ; "MPL_RdBu"

; do i = 0, 1
;     res@gsnLeftString = title(i)+" run - Ctrl run"
;     plot(i) = gsn_csm_pres_hgt(wks, t_diff(i,:,{lonstr:lonend}), res)
;     plot_prob(i) = gsn_csm_contour(wks,prob(i,:,{lonstr:lonend}),res3)
;     plot_prob(i) = ShadeLtContour(plot_prob(i),0.055,17)
;     overlay(plot(i), plot_prob(i))
; end do

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelLabelBar  = True
; pres@gsnPanelMainString = season+" mean: T profiles over "+LAT+"~S~o~N~N"
; gsn_panel(wks, plot, (/1,2/), pres)

res@gsnLeftString = "Air Temp. diff"
res@gsnRightString = "lat: "+LAT+"~S~o~N~N"
plot(0) = gsn_csm_pres_hgt(wks, t_diff(0,:,{lonstr:lonend}), res)
plot_prob(0) = gsn_csm_contour(wks,prob(0,:,{lonstr:lonend}),res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),0.055,17)
overlay(plot(0), plot_prob(0))
draw(plot(0))

end