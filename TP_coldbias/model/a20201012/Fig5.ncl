
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

year = ispan(1981, 2005, 1)
nyear = dimsizes(year)
nlat = 96
nlon = 144

Num = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"

ncase = dimsizes(casename)
season = "DJF"
; season = "JJA"
var = (/"UpR","DownR","netR","DownS","UpS","netS","SH","LH","G0"/) ;"G0=netS-netR-LH-SH" "DownL"
nvar = dimsizes(var)

TS = new((/ncase,nvar,nyear,nlat,nlon/), float) ;, "1e+36"
do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"

    fili1 = casename(i)+".cam.h0.FSDS.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.FSNS.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.FLDS.1979-2005.nc"
    fili4 = casename(i)+".cam.h0.FLNS.1979-2005.nc"
    fili5 = casename(i)+".cam.h0.SHFLX.1979-2005.nc"
    fili6 = casename(i)+".cam.h0.LHFLX.1979-2005.nc"
    fili7 = casename(i)+".clm2.h0.FGR.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")
    f5 = addfile(diri+fili5, "r")
    f6 = addfile(diri+fili6, "r")
    f7 = addfile(diri+fili7, "r")
    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    x1 = f1->FSDS(iYYYY,:,:)
    x2 = f2->FSNS(iYYYY,:,:)
    x3 = f3->FLDS(iYYYY,:,:)
    x4 = f4->FLNS(iYYYY,:,:)
    x5 = f5->SHFLX(iYYYY,:,:)
    x6 = f6->LHFLX(iYYYY,:,:)
    ; printVarSummary(x5)
    ; printVarSummary(x6)
    ; exit
    x7 = x6
    ; x7 = (/f7->FGR(iYYYY,:,:)/)
    x7 = x2-x4-x5-x6

    us = x2
    us = x1 - x2
    ; us = -1.*us
    ur = x4
    ur = x3 + x4
    ; ur = -1.*ur
    ; x5 = -1.*x5
    ; x6 = -1.*x6
    x4 = -1.*x4

    TS(i,0,:,:,:) = month_to_season(ur, season)
    TS(i,1,:,:,:) = month_to_season(x3, season)
    TS(i,2,:,:,:) = month_to_season(x4, season)
    TS(i,3,:,:,:) = month_to_season(x1, season)
    TS(i,4,:,:,:) = month_to_season(us, season)
    TS(i,5,:,:,:) = month_to_season(x2, season)
    TS(i,6,:,:,:) = month_to_season(x5, season)
    TS(i,7,:,:,:) = month_to_season(x6, season)
    TS(i,8,:,:,:) = month_to_season(x7, season)

    ; TS(i,0,:,:) = dim_avg_n_Wrap(X1, 0)
    ; TS(i,1,:,:) = dim_avg_n_Wrap(US, 0)
    ; TS(i,2,:,:) = dim_avg_n_Wrap(X2, 0)
    ; TS(i,3,:,:) = dim_avg_n_Wrap(X3, 0)
    ; TS(i,4,:,:) = dim_avg_n_Wrap(UR, 0)
    ; TS(i,5,:,:) = dim_avg_n_Wrap(X4, 0)
    ; TS(i,6,:,:) = dim_avg_n_Wrap(X5, 0)
    ; TS(i,7,:,:) = dim_avg_n_Wrap(X6, 0)
    ; TS(i,8,:,:) = dim_avg_n_Wrap(X7, 0)
end do
TS!0 = "case"
TS&case = casename
TS!1 = "var"
TS&var = var
printVarSummary(TS)

t_var = dim_variance_n_Wrap(TS, 2)
t_ave = dim_avg_n_Wrap(TS, 2)
prob = t_ave(0,:,:,:)
do i = 0, nvar-1
    prob(i,:,:) = ttest(t_ave(0,i,:,:), t_var(0,i,:,:), nyear, t_ave(1,i,:,:), t_var(1,i,:,:), nyear, False, False)
end do
printVarSummary(prob)

ts_diff = t_ave(0,:,:,:)
ts_diff = t_ave(1,:,:,:) - t_ave(0,:,:,:)
copy_VarCoords(t_ave(0,0,:,:), ts_diff(0,:,:))
printVarSummary(ts_diff)

;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
topo = smth9_Wrap(topo, 0.5, 0.25, True)
;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/"
figname = "Fig5"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(9, graphic)
plot_topo = plot
plot_prob = plot

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
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

;------------------------
; draw prob
;------------------------
sig = 0.95

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

res3@cnFillPattern = 17 ;(/17/)
res3@cnFillColor  =  "black" ;"green4"
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1 ;1.2

;--------------------------------------------------

res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
res@cnFillPalette         = "CBR_coldhot"
res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035


res@gsnRightString = "W/m~S~2~N~"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -10
res@cnMaxLevelValF = 10
res@cnLevelSpacingF = 1.
do i = 0,5
    res@gsnLeftString = Num(i)+" "+var(i)
    plot(i) = gsn_csm_contour_map(wks, ts_diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    plot_prob(i) = gsn_csm_contour(wks, prob(i,:,:), res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),(1-sig),6)
    overlay(plot(i), plot_prob(i))
end do

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 0.5

do i = 6,8
    res@gsnLeftString = Num(i)+" "+var(i)
    plot(i) = gsn_csm_contour_map(wks, ts_diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    plot_prob(i) = gsn_csm_contour(wks, prob(i,:,:), res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),(1-sig),6)
    overlay(plot(i), plot_prob(i))
end do

; pres = True
; pres@gsnMaximize = True
; ;pres@gsnPanelLabelBar = True
; ; pres@gsnPanelMainString = season+" mean change: downward is positive" ;(500hPa to 200hPa)
; ; pres@pmLabelBarOrthogonalPosF = -0.02
; ; pres@lbLabelFontHeightF = 0.015
; gsn_panel(wks, plot, (/3,3/), pres)


pres = True
; pres@gsnMaximize = True
pres@gsnDraw  = False
pres@gsnFrame = False

height = 0.6
pres@gsnDraw = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
pres@gsnPanelTop = 0.9
pres@gsnPanelBottom = pres@gsnPanelTop - height
gsn_panel(wks, plot(0:5), (/2,3/), pres)


pres@gsnPanelLabelBar = True
pres@gsnPanelTop = 0.43 ;pres@gsnPanelBottom
pres@gsnPanelBottom = pres@gsnPanelTop - 0.33
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot(6:8), (/1,3/), pres)
end
