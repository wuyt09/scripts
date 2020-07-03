
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "Ann"
SF = new((/ncase,96,144/), float)
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
    SF0 = month_to_annual(SFX, 1)
    SF(i,:,:) = dim_avg_n_Wrap(SF0, 0)

    SCX = f3->FSNO(iYYYY,:,:)
    SC0 = month_to_annual(SCX, 1)
    SC(i,:,:) = dim_avg_n_Wrap(SC0, 0)

    OMEGAX = f4->OMEGA(iYYYY,{LEV},:,:)
    OMEGA0 = month_to_annual(OMEGAX, 1)
    OMEGA(i,:,:) = dim_avg_n_Wrap(OMEGA0, 0)

    TSX = f5->TREFHT(iYYYY,:,:)
    TS0 = month_to_annual(TSX, 1)
    TS(i,:,:) = dim_avg_n_Wrap(TS0, 0)

delete([/diri,fili1,fili2,fili3,fili4,fili5,f1,f2,f3,f4,f5,SFX,SCX,OMEGAX,TSX/])
end do
printVarSummary(SF)

sf = SF
sf(1,:,:) = SF(1,:,:) - SF(0,:,:)
sc = SC
sc(1,:,:) = SC(1,:,:) - SC(0,:,:)
omega = OMEGA
omega(1,:,:) = OMEGA(1,:,:) - OMEGA(0,:,:)
;omega = omega*100
ts = TS
ts(0,:,:) = ts(0,:,:) - 273.15
ts(1,:,:) = TS(1,:,:) - TS(0,:,:)
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
res@gsnRightString = "Clim"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.
res@cnMaxLevelValF = 5.
res@cnLevelSpacingF = 0.5
plot(0) = gsn_csm_contour_map(wks, sf(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@cnFillPalette  = "precip_diff_12lev"
res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.05
plot(1) = gsn_csm_contour_map(wks, sf(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
;----------------------------------------------------
res@cnFillPalette  = "CBR_wet"
res@gsnLeftString  = "Snow Cover"
res@gsnRightString = "Clim"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0.1
res@cnMaxLevelValF = 1.
res@cnLevelSpacingF = 0.1
plot(2) = gsn_csm_contour_map(wks, sc(0,:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@cnFillPalette = "precip_diff_12lev"
res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
plot(3) = gsn_csm_contour_map(wks, sc(1,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))
;----------------------------------------------------
res@cnFillPalette  = "MPL_RdBu"
res@gsnLeftString  = "OMEGA"+LEV
res@gsnRightString = "Clim"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.2
res@cnMaxLevelValF = 0.2
res@cnLevelSpacingF = 0.05
plot(4) = gsn_csm_contour_map(wks, omega(0,:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.05
res@cnMaxLevelValF = 0.05
res@cnLevelSpacingF = 0.005
plot(5) = gsn_csm_contour_map(wks, omega(1,:,:), res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))
;----------------------------------------------------
res@cnFillPalette  = "WhiteBlue"
res@gsnSpreadColorEnd = 2
res@gsnSpreadColorStart = 253
res@gsnLeftString  = "TS"
res@gsnRightString = "Clim"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -30.
res@cnMaxLevelValF = 0.
res@cnLevelSpacingF = 3
plot(6) = gsn_csm_contour_map(wks, ts(0,:,:), res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))

res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3.
res@cnMaxLevelValF = 0.
res@cnLevelSpacingF = 0.3
plot(7) = gsn_csm_contour_map(wks, ts(1,:,:), res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7), plot_topo(7))
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