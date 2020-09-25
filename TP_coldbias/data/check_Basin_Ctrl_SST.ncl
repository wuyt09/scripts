
begin

yrStrt = 1979
yrLast = 2005


casename = (/"Ctrl","Hist","Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"/)

ncase = dimsizes(casename)

TS = new((/ncase,324,96,144/), float)

filictrl = "./sst_HadOIBl_bc_1.9x2.5_1850_2012_c130411.nc"
filihist = "./Hist_bias_sst_cmip5_bc_1.9x2.5_c20191204.nc"
filiio   = "./IO_sst_cmip5_bc_1.9x2.5_c20191204.nc"
filipo   = "./PO_sst_cmip5_bc_1.9x2.5_c20191204.nc"
filiao   = "./AO_sst_cmip5_bc_1.9x2.5_c20191204.nc"
filipio   = "./PIO_sst_cmip5_bc_1.9x2.5_c20191225.nc"
filipao   = "./PAO_sst_cmip5_bc_1.9x2.5_c20191225.nc"
filiiao   = "./IAO_sst_cmip5_bc_1.9x2.5_c20191204.nc"
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
tdiff = ts(0:ncase-2,:,:)
do i = 0, ncase-2
    tdiff(i,:,:) = ts(i+1,:,:) - ts(0,:,:)
end do
printVarSummary(tdiff)

;=========================================================
; plot 
;=========================================================
figpath = "./"
figname = "check_Basin_Ctrl_SST"

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
res@mpMaxLatF            =  30.
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

do i = 0, ncase-2
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_contour_map(wks, tdiff(i,:,:), res)
end do

; res@cnMinLevelValF = 25
; res@cnMaxLevelValF = 30
; res@cnLevelSpacingF = 0.5
res@lbLabelBarOn = True
res@pmLabelBarOrthogonalPosF = 0.4
res@gsnLeftString = "All Basin"
res@gsnRightString = "K"
plot(0) = gsn_csm_contour_map(wks, tdiff(0,:,:), res)

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarOrthogonalPosF = -0.01
; pres@gsnPanelMainString = "SST difference"

; gsn_panel(wks, plot, (/2,1/), pres)

draw(plot(0))
frame(wks)
end