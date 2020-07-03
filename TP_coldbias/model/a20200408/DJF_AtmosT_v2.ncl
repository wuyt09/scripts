
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)

season = "DJF"

vt = new((/ncase,17,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.T.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.PS.1979-2005.nc"
    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    
    TIME = f2->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    VT = f1->T(iYYYY,:,:,:)
    vtx = month_to_season(VT, season)
    vt(i,:,:,:) = dim_avg_n_Wrap(vtx, 0)
end do
printVarSummary(vt)

plev = vt&level
plev = plev*100
psfc = max(plev)
ptop = min(plev)
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dp = conform_dims(dimsizes(vt), dP, 1)
vtdp = vt*dp

vtt = dim_sum_n_Wrap(vtdp, 1)/sum(dP)
copy_VarCoords(vt(:,0,:,:), vtt)

vt_diff = vtt(0,:,:)
vt_diff = vtt(1,:,:) - vtt(0,:,:)
print(vtt(0,{35},{80}))
print(vtt(1,{35},{80}))
print(vt_diff({35},{80}))
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_AtmosT_v2"

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

res@gsnRightString = "K"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = 30
; res@cnMaxLevelValF = 80
; res@cnLevelSpacingF = 5.
res@gsnLeftString = "Ctrl"
plot(0) = gsn_csm_contour_map(wks, vtt(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "Sen"
plot(1) = gsn_csm_contour_map(wks, vtt(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.1
res@gsnLeftString = "Diff"
plot(2) = gsn_csm_contour_map(wks, vt_diff, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean AtmosT (1000hPa-300hPa) " ;(500hPa to 200hPa)
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/1,3/), pres)

end