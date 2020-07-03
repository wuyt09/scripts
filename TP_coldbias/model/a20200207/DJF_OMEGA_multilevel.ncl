
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)
LEV = (/500,300,200,100/)
nlev = dimsizes(LEV)
printVarSummary(LEV)
season = "DJF"
OMEGA = new((/ncase,nlev,96,144/), float)
; OMEGA_X = new((/ncase,17,144/), float)
; OMEGA_Y = new((/ncase,17,96/), float)


do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili = casename(i)+".cam.h0.OMEGA.1979-2005.nc"
    
    f1 = addfile(diri+fili, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    OMEGAX = f1->OMEGA(iYYYY,{LEV},:,:)
    OMEGA0 = month_to_season(OMEGAX, season)
    OMEGA(i,:,:,:) = dim_avg_n_Wrap(OMEGA0, 0)
    OMEGA(i,:,:,:) = OMEGA(i,:,:,:)*100.

    ; OMEGA_X(i,:,:) = dim_avg_n_Wrap(OMEGA(:,{latstr:latend},:), 1)
    ; OMEGA_Y(i,:,:) = dim_avg_n_Wrap(OMEGA(:,:,{lonstr:lonend}), 2)
    ; OMEGA_X(i,:,:) = OMEGA(:,{34},:)
    ; OMEGA_Y(i,:,:) = OMEGA(:,:,{90})

end do

OMEGA@units = "10^-2 Pa/s"
printVarSummary(OMEGA)

diff = OMEGA(0,:,:,:)
diff = OMEGA(1,:,:,:) - OMEGA(0,:,:,:)
; printMinMax(OMEGA(0,:,:,:), False)
; printMinMax(diff, False)
; exit

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_OMEGA_multilevel"

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

res@cnFillPalette  = "MPL_RdBu"

do i = 0,3
res@gsnLeftString  = "OMEGA "+LEV(i)+"hPa"
res@gsnRightString = "Clim"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -20
res@cnMaxLevelValF = 20
res@cnLevelSpacingF = 2.
plot(i) = gsn_csm_contour_map(wks, OMEGA(0,i,:,:), res)
plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(i), plot_topo(i))

res@gsnRightString = "Diff"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 0.5
plot(i+4) = gsn_csm_contour_map(wks, diff(i,:,:), res)
plot_topo(i+4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(i+4), plot_topo(i+4))
end do


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar  = True
gsn_panel(wks, plot, (/2,4/), pres)

end