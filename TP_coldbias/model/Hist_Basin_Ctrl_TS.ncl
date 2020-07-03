
begin

casename = (/"Ctrl","Hist","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO"/) ;"Hist_AO",

ncase = dimsizes(casename)

;TS = new((/ncase,324,96,144/), float)
tsx = new((/ncase,27,96,144/), float)
do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias_"+casename(i)+"_FAMIP/"
    fili = casename(i)+".esm.cam.h0.TS.1979-2005.nc"
    f = addfile(diri+fili, "r")
    TS = f->TS
    ;tsx(i,:,:,:) = month_to_season(TS, "DJF")
    tsx(i,:,:,:) = month_to_annual(TS, 1)
delete([/diri,fili,f/])
end do
printVarSummary(TS)

ts = dim_avg_n_Wrap(tsx, 1)
; ;ts = month_to_season(TS, "DJF")
; ts0 = reshape(TS, (/5,27,12,96,144/))
; ts1 = dim_avg_n_Wrap(ts0, 1)

TO = ts(0,:,:)
TO = ts(1,:,:) - ts(0,:,:)
IO = ts(0,:,:)
IO = ts(2,:,:) - ts(0,:,:)
PO = ts(0,:,:)
PO = ts(3,:,:) - ts(0,:,:)
AO = ts(0,:,:)
AO = ts(4,:,:) - ts(0,:,:)
PIO = ts(0,:,:)
PIO = ts(5,:,:) - ts(0,:,:)
PAO = ts(0,:,:)
PAO = ts(6,:,:) - ts(0,:,:)

to = wgt_areaave_Wrap(TO({25:40},{80:105}), 1, 1, 0)
io = wgt_areaave_Wrap(IO({25:40},{80:105}), 1, 1, 0)
po = wgt_areaave_Wrap(PO({25:40},{80:105}), 1, 1, 0)
ao = wgt_areaave_Wrap(AO({25:40},{80:105}), 1, 1, 0)
pio = wgt_areaave_Wrap(PIO({25:40},{80:105}), 1, 1, 0)
pao = wgt_areaave_Wrap(PAO({25:40},{80:105}), 1, 1, 0)
print("total "+to)
print("IO "+io)
print("PO "+po)
print("AO "+ao)
print("PIO "+pio)
print("PAO "+pao)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = "Hist_Basin_Ctrl_TS"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(dimsizes(casename)-1, graphic)


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
;res@cnLevelSpacingF       = 0.5      ; contour spacing
res@cnFillPalette         = "hotcold_18lev"
res@lbLabelBarOn           = False

res@mpMinLatF            = -30      ; range to zoom in on
res@mpMaxLatF            =  80.
;res@mpMinLonF            =  30.
;res@mpMaxLonF            = 120.
res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -2
res@cnMaxLevelValF = 2
res@cnLevelSpacingF = 0.2
res@gsnStringFontHeightF = 0.025

res@gsnLeftString = "Total"
plot(0) = gsn_csm_contour_map(wks, TO, res)

res@gsnLeftString = "Indian Ocean"
plot(1) = gsn_csm_contour_map(wks, IO, res)

res@gsnLeftString = "Pacific Ocean"
plot(2) = gsn_csm_contour_map(wks, PO, res)

res@gsnLeftString = "Atlantic Ocean"
plot(3) = gsn_csm_contour_map(wks, AO, res)

res@gsnLeftString = "Pacific&Indian Ocean"
plot(4) = gsn_csm_contour_map(wks, PIO, res)

res@gsnLeftString = "Pacific&Atlantic Ocean"
plot(5) = gsn_csm_contour_map(wks, PAO, res)

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "TS"

gsn_panel(wks, plot, (/3,2/), pres)

end
