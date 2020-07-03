
begin
lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

;U = new((/ncase,324,96,144/), float)
TS = new((/ncase,25,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili = casename(i)+".cam.h0.FSNS.1979-2005.nc"
    f = addfile(diri+fili, "r")
    TIME = f->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    TSX = f->FSNS(iYYYY,:,:)
    ;tsx(i,:,:,:) = month_to_season(TS, "DJF")
    TS(i,:,:,:) = month_to_annual(TSX, 1)
delete([/diri,fili,f,TSX/])
end do
printVarSummary(TS)

ts = dim_avg_n_Wrap(TS, 1)
; ;ts = month_to_season(TS, "DJF")
; ts0 = reshape(TS, (/5,27,12,96,144/))
; ts1 = dim_avg_n_Wrap(ts0, 1)

tdiff = ts(0:ncase-2,:,:)
do i = 0, ncase-2
    tdiff(i,:,:) = ts(i+1,:,:) - ts(0,:,:)
end do
printVarSummary(tdiff)

;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, TS&lon, TS&lat, 0)
TS = mask(TS, TOPO.ge.1500, True)

tmaskdiff = mask(tdiff, TOPO.ge.1500, True)

txdiff = wgt_areaave_Wrap(tmaskdiff(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
leftsti = (/"Total","IndianOnly","PacificOnly","AtlanticOnly",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)
print(leftsti+" "+txdiff)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20191206/"
figname = "Ann_Hist_Basin_Ctrl_netSR"

lonstr = 60
lonend = 120
latstr = 10
latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(ncase, graphic)
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
res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3
res@cnMaxLevelValF = 3
res@cnLevelSpacingF = 0.3
res@gsnStringFontHeightF = 0.035
res@gsnRightString       = "K"

do i = 0, ncase-2
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_contour_map(wks, tdiff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Annual mean: netSR"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,4/), pres)

end
