begin

casename = (/"Ctrl","Hist","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

PREC = new((/ncase,27,96,144/), float)


do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias_"+casename(i)+"_FAMIP/"
    filipc = casename(i)+".esm.cam.h0.PRECC.1979-2005.nc"
    filipl = casename(i)+".esm.cam.h0.PRECL.1979-2005.nc"
    fc = addfile(diri+filipc, "r")
    fl = addfile(diri+filipl, "r")

    PRECC = fc->PRECC
    PRECL = fl->PRECL
    PRECX = PRECC
    PRECX = (PRECC + PRECL)*86400*100.
    PREC(i,:,:,:) = month_to_annual(PRECX, 1)
delete([/diri,filipc,filipl,fc,fl,PRECC,PRECL,PRECX/])
end do

printVarSummary(PREC)
prec = dim_avg_n_Wrap(PREC, 1)
;print(prec(0,{20:40},{110:120}))
;print(prec(0,{20:40},{75:105}))
precdiff = prec(0:ncase-2,:,:)

do i = 0, ncase-2
    precdiff(i,:,:) = prec(i+1,:,:) - prec(0,:,:)
end do
printVarSummary(precdiff)

txdiff = wgt_areaave_Wrap(precdiff(:,{25:40},{80:105}), 1, 1, 0)
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
figname = "Ann_Hist_Basin_Ctrl_Precip"

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


res = True
res@gsnMaximize = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpFillOn = True
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpMinLatF = latstr
res@mpMaxLatF = latend

res@cnFillOn = True
res@cnFillPalette = "precip4_diff_19lev"
res@lbLabelBarOn = False
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False


res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
; res@cnMinLevelValF       = -0.1               ; set min contour level
; res@cnMaxLevelValF       = 0.1                 ; set max contour level
; res@cnLevelSpacingF      = 0.01              ; set contour spacing
res@cnMinLevelValF       = -0.2               ; set min contour level
res@cnMaxLevelValF       = 0.2                 ; set max contour level
res@cnLevelSpacingF      = 0.02              ; set contour spacing
res@gsnRightString       = "mm/day"
res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

do i = 0, ncase-2
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_contour_map(wks, precdiff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do


pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Annual mean: Precip"
pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,4/), pres)
end 