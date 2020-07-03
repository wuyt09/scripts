
begin
lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

nlat = 96
nlon = 144
nyear = 25
;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)
season = "DJF"

TS = new((/ncase,nyear,5,nlat,nlon/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.SHFLX.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.LHFLX.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.FSNS.1979-2005.nc"
    fili4 = casename(i)+".cam.h0.FLNS.1979-2005.nc"
    fili5 = casename(i)+".clm2.h0.FSR.1979-2005.nc"

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
    x5 = f5->FSR(iYYYY,:,:)
    ; x6 = x1+x2+x4;+x5

    TSX = new((/300,5,nlat,nlon/), float)
    TSX(:,0,:,:) = x1
    TSX(:,1,:,:) = x2
    TSX(:,2,:,:) = x3
    TSX(:,3,:,:) = x4
    TSX(:,4,:,:) = x5
    ; TSX(:,3,:,:) = (/x6/)
    TS(i,:,:,:,:) = month_to_season(TSX, season)
end do
printVarSummary(TS)


vt_ave = dim_avg_n_Wrap(TS, 1)
vt_var = dim_variance_n_Wrap(TS, 1)
prob = ttest(vt_ave(0,:,:,:), vt_var(0,:,:,:), nyear, vt_ave(1,:,:,:), vt_var(1,:,:,:), nyear, False, False)
copy_VarCoords(vt_ave(0,:,:,:), prob)
printVarSummary(prob)

vt_diff = vt_ave(0,:,:,:)
vt_diff = vt_ave(1,:,:,:) - vt_ave(0,:,:,:)
printVarSummary(vt_diff)

; print(vt_ave(:,{38},{98}))
; print(vt_diff({38},{98}))


;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, TS&lon, TS&lat, 0)

; tdiff = mask(tdiff, TOPO.ge.1500, True)

; txdiff = wgt_areaave_Wrap(tdiff(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
; leftsti = (/"Total","IndianOnly","PacificOnly","AtlanticOnly",\
;            "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)
; print(leftsti+" "+txdiff)
leftsti = (/"SH","LH","netS","netR","G0"/)
; leftsti = (/"SH","LH","netS","netR","SH+LH+netR"/)
;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = "DJF_surf_flux"

; lonstr = 60
; lonend = 120
; latstr = 20
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(10, graphic)
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
res@cnFillPalette         = "hotcold_18lev"
; res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@gsnStringFontHeightF = 0.035
res@gsnRightString       = "W/m~S~2~N~"


do i = 0, 1
    res@gsnLeftString = leftsti(i)
    res@cnMinLevelValF = -50
    res@cnMaxLevelValF = 50
    res@cnLevelSpacingF = 5
    res@gsnRightString = "Ctrl"
    plot(i) = gsn_csm_contour_map(wks, vt_ave(0,i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    delete(res@cnLevelSpacingF)

    res@cnMinLevelValF = -5
    res@cnMaxLevelValF = 5
    res@cnLevelSpacingF = 0.5
    res@gsnRightString = "Diff"
    plot(i+5) = gsn_csm_contour_map(wks, vt_diff(i,:,:), res)
    plot_topo(i+5) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i+5), plot_topo(i+5))

end do


do i = 2, 4
    res@gsnLeftString = leftsti(i)
    res@cnMinLevelValF = -300
    res@cnMaxLevelValF = 300
    res@cnLevelSpacingF = 50
    res@gsnRightString = "Ctrl"
    plot(i) = gsn_csm_contour_map(wks, vt_ave(0,i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    res@cnMinLevelValF = -10
    res@cnMaxLevelValF = 10
    res@cnLevelSpacingF = 1
    res@gsnRightString = "Diff"
    plot(i+5) = gsn_csm_contour_map(wks, vt_diff(i,:,:), res)
    plot_topo(i+5) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i+5), plot_topo(i+5))
end do


pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean: All Flux"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,5/), pres)

end
