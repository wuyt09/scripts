
begin

yrStrt = 1979
yrLast = 2005


casename = (/"Ctrl","Hist","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

TS = new((/ncase,324,96,144/), float)

filictrl = "../data/sst_HadOIBl_bc_1.9x2.5_1850_2013_c140701.nc"
filihist = "../data/sst_cmip5_bc_1.9x2.5_1950_1999.nc"
filiio   = "../data/IO_sst_cmip5_bc_1.9x2.5_c20191112.nc"
filipo   = "../data/PO_sst_cmip5_bc_1.9x2.5_c20191112.nc"
filiao   = "../data/AO_sst_cmip5_bc_1.9x2.5_c20191112.nc"
filipio   = "../data/PIO_sst_cmip5_bc_1.9x2.5_c20191117.nc"
filipao   = "../data/PAO_sst_cmip5_bc_1.9x2.5_c20191119.nc"
filiiao   = "../data/IAO_sst_cmip5_bc_1.9x2.5_c20191121.nc"
fctrl = addfile(filictrl, "r")
fhist = addfile(filihist, "r")
fio   = addfile(filiio, "r")
fpo   = addfile(filipo, "r")
fao   = addfile(filiao, "r")
fpio   = addfile(filipio, "r")
fpao   = addfile(filipao, "r")
fiao   = addfile(filiiao, "r")

Time = fctrl->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

TS(0,:,:,:) = fctrl->SST_cpl(iYYYY,:,:)
TS(1,:,:,:) = fhist->SST_cpl(iYYYY,:,:)
TS(2,:,:,:) = fio->SST_cpl(iYYYY,:,:)
TS(3,:,:,:) = fpo->SST_cpl(iYYYY,:,:)
TS(4,:,:,:) = fao->SST_cpl(iYYYY,:,:)
TS(5,:,:,:) = fpio->SST_cpl(iYYYY,:,:)
TS(6,:,:,:) = fpao->SST_cpl(iYYYY,:,:)
TS(7,:,:,:) = fiao->SST_cpl(iYYYY,:,:)


ts = dim_avg_n_Wrap(TS, 1)
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
IAO = ts(0,:,:)
IAO = ts(7,:,:) - ts(0,:,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/"
figname = "Hist_Basin_Ctrl_SST"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
;res@cnLevelSpacingF       = 0.5      ; contour spacing
res@cnFillPalette         = "hotcold_18lev"
res@lbLabelBarOn          = False
res@cnInfoLabelOn         = False
res@cnLineLabelsOn        = False

res@mpMinLatF            = -30      ; range to zoom in on
res@mpMaxLatF            =  60.
;res@mpMinLonF            =  30.
;res@mpMaxLonF            = 120.
res@mpCenterLonF         = 180.

res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
res@gsnStringFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -2
res@cnMaxLevelValF = 2
res@cnLevelSpacingF = 0.2

leftsti = (/"Total","IndianOnly","PacificOnly","AtlanticOnly",\
           "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)

res@gsnLeftString = leftsti(0)
plot(0) = gsn_csm_contour_map(wks, TO, res)

res@gsnLeftString = leftsti(1)
plot(1) = gsn_csm_contour_map(wks, IO, res)

res@gsnLeftString = leftsti(2)
plot(2) = gsn_csm_contour_map(wks, PO, res)

res@gsnLeftString = leftsti(3)
plot(3) = gsn_csm_contour_map(wks, AO, res)

res@gsnLeftString = leftsti(4)
plot(4) = gsn_csm_contour_map(wks, PIO, res)

res@gsnLeftString = leftsti(5)
plot(5) = gsn_csm_contour_map(wks, PAO, res)

res@gsnLeftString = leftsti(6)
plot(6) = gsn_csm_contour_map(wks, IAO, res)

; res@cnMinLevelValF = 25
; res@cnMaxLevelValF = 30
; res@cnLevelSpacingF = 0.5

; res@gsnLeftString = "Ctrl"
; plot(7) = gsn_csm_contour_map(wks, ts(0,:,:), res)

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01
pres@gsnPanelMainString = "SST difference"

gsn_panel(wks, plot, (/4,2/), pres)

end