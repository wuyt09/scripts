
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
LEV = 500

fx = "time_series_TP_skt_DJF_MAM_ERA5.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(:,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))

; hig_ind = hig_ind+1
; low_ind = low_ind+1

LATs = 24
LATe = 30
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

topo_z = dim_avg_n_Wrap(topo({LATs:LATe},:), 0)
topo_zonal = 1013.25*(1-topo_z*0.0065/288.15)^5.25145
copy_VarCoords(topo({LATs},:), topo_zonal)
;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili = "w_mon_1x1_1979_2019.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
varx = short2flt(f->w(iYYYY,:,:,:))
varx = varx*100.
varx@units = "10^-2 Pa/s"
; var@units = "W/m~S~2~N~"
; var = var * (-1.)
printVarSummary(varx)

var = dim_avg_n_Wrap(varx(:,:,{LATs:LATe},:), 2)
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
lonstr = 70
lonend = 130
; latstr = 10
; latend = 50

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Comp_MAM_Omega_profile_by_TP_skt_DJF"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)
plot_topo = plot
plot_prob = new(1, graphic)

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
    res@tmYLValues    = (/"1000","850","750","500","300","250","150","100"/) ;
    res@tmYLLabels    = (/"1000","850","750","500","300","250","150","100"/) ;

    res@trYMinF   = 100
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

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 5.
res@cnMaxLevelValF = 50.
res@cnMinLevelValF = -50.
res@cnFillPalette = "temp_diff_18lev" ; "MPL_RdBu"


res@gsnLeftString = "(a) High DJF TP t2m"
plot(0) = gsn_csm_pres_hgt(wks, var_hig_ave(:,{lonstr:lonend}), res)
plot_topo(0) = gsn_csm_xy(wks, topo&lon, topo_zonal, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "(b) Low DJF TP t2m"
plot(1) = gsn_csm_pres_hgt(wks, var_low_ave(:,{lonstr:lonend}), res)
plot_topo(1) = gsn_csm_xy(wks, topo&lon, topo_zonal, res_topo)
overlay(plot(1), plot_topo(1))

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 10.
res@cnMinLevelValF = -10.
res@gsnLeftString = "(c) a-b"
plot(2) = gsn_csm_pres_hgt(wks, diff_var(:,{lonstr:lonend}), res)
plot_topo(2) = gsn_csm_xy(wks, topo&lon, topo_zonal, res_topo)
overlay(plot(2), plot_topo(2))
plot_prob(0) = gsn_csm_contour(wks,prob,res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),0.05,6)
overlay(plot(2), plot_prob(0))

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Omega profile ("+LATs+"-"+LATe+"N) in MAM"
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = -0.015
gsn_panel(wks, plot, (/1,3/), pres)


end