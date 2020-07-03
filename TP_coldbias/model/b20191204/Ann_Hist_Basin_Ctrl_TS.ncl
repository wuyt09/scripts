
begin

casename = (/"Ctrl","Hist","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

;U = new((/ncase,324,96,144/), float)
TS = new((/ncase,27,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias_"+casename(i)+"_FAMIP/"
    fili = casename(i)+".esm.cam.h0.TREFHT.1979-2005.nc"
    f = addfile(diri+fili, "r")
    TSX = f->TREFHT
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

txdiff = wgt_areaave_Wrap(tdiff(:,{25:40},{80:105}), 1, 1, 0)
leftsti = (/"Total","IndianOnly","PacificOnly","AtlanticOnly",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)
print(leftsti+" "+txdiff)

;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = "Ann_Hist_Basin_Ctrl_TAS"

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
res@cnMinLevelValF = -1.5
res@cnMaxLevelValF = 1.5
res@cnLevelSpacingF = 0.1
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
pres@gsnPanelMainString = "Annual mean: TS"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,4/), pres)

end
