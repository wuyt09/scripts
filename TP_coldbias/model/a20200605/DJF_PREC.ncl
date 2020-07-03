
begin


lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"
ncase = dimsizes(casename)
nlat = 96
nlon = 144
season = "DJF"
TS = new((/ncase,25,3,nlat,nlon/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.PREC.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.PRECS.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    SC = f1->PREC(iYYYY,:,:)
    SL = f2->PRECS(iYYYY,:,:)

    TSX = new((/300,3,nlat,nlon/), float)

    TSX(:,0,:,:) = (/SC/)
    TSX(:,1,:,:) = (/SL/) 
    TSX(:,2,:,:) = TSX(:,0,:,:) - TSX(:,1,:,:)
    copy_VarCoords(SC, TSX(:,0,:,:))
    TSX = TSX * 1000 * 86400.
    TSX@units = "mm/day"
    TS(i,:,:,:,:) = month_to_season(TSX, season)
delete([/diri,fili1,fili2,f1,f2,TSX/])
end do
printVarSummary(TS)

varlist = (/"Precipitation","Snowfall","Rainfall"/)

TS!0 = "case"
TS&case = casename
TS!1 = "year"
TS&year = ispan(1981, 2005, 1)
TS!2 = "var"
TS&var = varlist

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
; TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, TS&lon, TS&lat, 0)
; TS = mask(TS, TOPO.ge.1500, True)

ts = dim_avg_n_Wrap(TS, 1)

tdiff = ts(0:ncase-2,:,:,:)
do i = 0, ncase-2
    tdiff(i,:,:,:) = ts(i+1,:,:,:) - ts(0,:,:,:)
end do
printVarSummary(tdiff)


; txdiff = wgt_areaave_Wrap(tdiff(:,{25:40},{75:105}), 1, 1, 0)
; leftsti = (/"Total","IndianOnly","PacificOnly","AtlanticOnly",\
;            "Pacific&Indian","Pacific&Atlantic","Indian&Atlantic"/)
; print(sprintf("%4.2f", txdiff)+" "+leftsti)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = season+"_PREC"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)
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
res@cnFillPalette         = "CBR_drywet"
res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.5
res@cnMaxLevelValF = 0.5
res@cnLevelSpacingF = 0.1
res@gsnStringFontHeightF = 0.035
;res@gsnRightString       = "K"

do i = 0, 2
    res@gsnLeftString = varlist(i)
    plot(i) = gsn_csm_contour_map(wks, tdiff(0,i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean diff"
; pres@pmLabelBarOrthogonalPosF = -0.02
pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/3,1/), pres)


end

