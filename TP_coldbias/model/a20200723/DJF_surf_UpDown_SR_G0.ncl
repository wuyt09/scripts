
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601


casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"

ncase = dimsizes(casename)
season = "DJF"
; season = "JJA"
var = (/"DownS","UpS","netS","DownR","UpR","netR","SH","LH","G0"/) ;"G0=netS-netR-LH-SH" "DownL"
nvar = dimsizes(var)

TS = new((/ncase,nvar,96,144/), float) ;, "1e+36"
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
    x7 = (/f7->FGR(iYYYY,:,:)/)

    us = x2
    us = x1 - x2
    us = -1.*us
    ur = x4
    ur = x3 + x4
    ur = -1.*ur
    x5 = -1.*x5
    x6 = -1.*x6
    x4 = -1.*x4

    X1 = month_to_season(x1, season)
    TS(i,0,:,:) = dim_avg_n_Wrap(X1, 0)
    US = month_to_season(us, season)
    TS(i,1,:,:) = dim_avg_n_Wrap(US, 0)
    X2 = month_to_season(x2, season)
    TS(i,2,:,:) = dim_avg_n_Wrap(X2, 0)
    X3 = month_to_season(x3, season)
    TS(i,3,:,:) = dim_avg_n_Wrap(X3, 0)
    UR = month_to_season(ur, season)
    TS(i,4,:,:) = dim_avg_n_Wrap(UR, 0)
    X4 = month_to_season(x4, season)
    TS(i,5,:,:) = dim_avg_n_Wrap(X4, 0)

    X5 = month_to_season(x5, season)
    TS(i,6,:,:) = dim_avg_n_Wrap(X5, 0)
    X6 = month_to_season(x6, season)
    TS(i,7,:,:) = dim_avg_n_Wrap(X6, 0)

    X7 = month_to_season(x7, season)
    TS(i,8,:,:) = dim_avg_n_Wrap(X7, 0)
end do
TS!0 = "case"
TS&case = casename
TS!1 = "var"
TS&var = var
printVarSummary(TS)

ts_diff = TS(0,:,:,:)
ts_diff = TS(1,:,:,:) - TS(0,:,:,:)
copy_VarCoords(TS(0,0,:,:), ts_diff(0,:,:))
printVarSummary(ts_diff)

;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot
;=========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200723/"
figname = season+"_surf_UpDown_SR_G0"


wks = gsn_open_wks("pdf", figpath+figname)
plot = new(9, graphic)
plot_topo = plot

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

res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
res@cnFillPalette         = "CBR_coldhot"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035

;res@gsnRightString = "W/m~S~2~C~"

; res@gsnRightString = "Ctrl"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -500
; res@cnMaxLevelValF = 500
; res@cnLevelSpacingF = 50

; do i = 0, 8
;     res@gsnLeftString = var(i)
;     plot(i) = gsn_csm_contour_map(wks, TS(0,i,:,:), res)
;     plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
;     overlay(plot(i), plot_topo(i))
;     ; write loop content
; end do

res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -10
res@cnMaxLevelValF = 10
res@cnLevelSpacingF = 1.
do i = 0,5
    res@gsnLeftString = var(i)
    plot(i) = gsn_csm_contour_map(wks, ts_diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    ; res@gsnLeftString = var(i+3)
    ; plot(i+3) = gsn_csm_contour_map(wks, ts_diff(i+3,:,:), res)
    ; plot_topo(i+3) = gsn_csm_contour(wks, topo, res_topo)
    ; overlay(plot(i+3), plot_topo(i+3))
end do

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -2
res@cnMaxLevelValF = 2
res@cnLevelSpacingF = 0.2

; res@gsnLeftString = var(2)
; plot(2) = gsn_csm_contour_map(wks, ts_diff(2,:,:), res)
; plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(2), plot_topo(2))

; res@gsnLeftString = var(5)
; plot(5) = gsn_csm_contour_map(wks, ts_diff(5,:,:), res)
; plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(5), plot_topo(5))

do i = 6,8
    res@gsnLeftString = var(i)
    plot(i) = gsn_csm_contour_map(wks, ts_diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean change: downward is positive" ;(500hPa to 200hPa)
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/3,3/), pres)
end
