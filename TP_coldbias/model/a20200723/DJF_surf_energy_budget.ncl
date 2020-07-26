
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

TS = new((/ncase,6,96,144/), float, "1e+36")
do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.SHFLX.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.LHFLX.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.FSNS.1979-2005.nc"
    fili4 = casename(i)+".cam.h0.FLNS.1979-2005.nc"
    fili5 = casename(i)+".clm2.h0.FGR.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")
    f5 = addfile(diri+fili5, "r")
    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    x1 = f1->SHFLX(iYYYY,:,:)
    x2 = f2->LHFLX(iYYYY,:,:)
    x3 = f3->FSNS(iYYYY,:,:)
    x4 = f4->FLNS(iYYYY,:,:)
    ; x5 = f5->FLDS(iYYYY,:,:)
    x5 = x4
    x5 = (/f5->FGR(iYYYY,:,:)/)
    ; printVarSummary(x3)
    ; printVarSummary(x4)
    ; printVarSummary(x5)
    ; exit

    ; TS(i,0,:,:,:) = clmMonTLL(x1)
    ; TS(i,1,:,:,:) = clmMonTLL(x2)
    ; TS(i,2,:,:,:) = clmMonTLL(x3)
    ; TS(i,3,:,:,:) = clmMonTLL(x4)
    ; TS(i,4,:,:,:) = clmMonTLL(x5)

    X1 = month_to_season(x1, season)
    TS(i,0,:,:) = dim_avg_n_Wrap(X1, 0)
    X2 = month_to_season(x2, season)
    TS(i,1,:,:) = dim_avg_n_Wrap(X2, 0)
    X3 = month_to_season(x3, season)
    TS(i,2,:,:) = dim_avg_n_Wrap(X3, 0)
    X4 = month_to_season(x4, season)
    TS(i,3,:,:) = dim_avg_n_Wrap(X4, 0)
    X5 = month_to_season(x5, season)
    TS(i,4,:,:) = dim_avg_n_Wrap(X5, 0)

    x0 = x1
    x0 = x3 - x1 - x2 - x4
    X0 = month_to_season(x0, season)
    TS(i,5,:,:) = dim_avg_n_Wrap(X0, 0)

end do
printVarSummary(TS)
; ts = TS(0:ncase-2,:,:,:)
; do i = 0, ncase-2
;     ts(i,:,:,:) = TS(i+1,:,:,:) - TS(0,:,:,:)
; end do
; printVarSummary(ts)

ts_diff = (/TS(0,:,:,:)/)
ts_diff = (/TS(1,:,:,:) - TS(0,:,:,:)/)
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
figname = season+"_surf_energy_budget"

leftsti = (/"SH","LH","netS","netR","G0","netS-netR-LH-SH"/)

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6, graphic)
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

do i = 0,3
    res@gsnLeftString = leftsti(i)
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -100
    res@cnMaxLevelValF = 100
    res@cnLevelSpacingF = 10.
    res@gsnRightString = "Ctrl"
    plot(i) = gsn_csm_contour_map(wks, TS(0,i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

do i = 4, 5
    res@gsnLeftString = leftsti(i)
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -10
    res@cnMaxLevelValF = 10
    res@cnLevelSpacingF = 1.
    res@gsnRightString = "Ctrl"
    plot(i) = gsn_csm_contour_map(wks, TS(0,i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean surface energy budget" ;(500hPa to 200hPa)
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,3/), pres)


end
