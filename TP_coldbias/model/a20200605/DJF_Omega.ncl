
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

vt = new((/ncase,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.OMEGA.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.PS.1979-2005.nc"
    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    
    TIME = f2->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    VT = f1->OMEGA(iYYYY,:,:,:)
    ; VT = VTX(:,1:7,:,:) ; 1000hPa-300hPa
    
    ; plev = VT&level
    ; plev!0 = "lev"
    ; plev&lev = VT&level
    ; plev = plev*100
    ; psfc = max(plev)
    ; ptop = min(plev)
    ; dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    ; copy_VarCoords(plev, dP)
    ; dp = conform_dims(dimsizes(VT), dP, 1)
    ; vtdp = VT*dp
    ; copy_VarCoords(VT, vtdp)
    ; vt0 = dim_sum_n_Wrap(vtdp(:,{500:300},:,:), 1)/sum(dP({500:300}))
    ; copy_VarCoords(VT(:,0,:,:), vt0)

    PS = f2->PS(iYYYY,:,:)

    lev = VT&level
    plev = lev
    plev = plev *100
    plev!0 = "lev"
    plev&lev = lev
    plev@units = "Pa"

    ptop = min(plev)
    psfc = PS
    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    ;dP(:,{600},{30},{90}) = dP@_FillValue
    printVarSummary(dP)
    dP = where(ismissing(VT), dP@_FillValue, dP)
    ; exit

    VTdP = VT*dP
    copy_VarCoords(VT, VTdP)
    print(VTdP(20,{600},{29},{84.5}))
    print(dP(20,{600},{29},{84.5}))
    print(VTdP(20,{850},{29},{84.5}))
    print(dP(20,{850},{29},{84.5}))
    vt0 = dim_sum_n_Wrap(VTdP(:,{600:300},:,:), 1)/dim_sum_n_Wrap(dP(:,{600:300},:,:), 1)
    copy_VarCoords(VTdP(:,0,:,:), vt0)

    vtx = month_to_season(vt0, season)
    vt(i,:,:) = dim_avg_n_Wrap(vtx, 0)

    ; vtx = month_to_season(VT(:,{750},:,:), season)
    ; vt(i,:,:) = dim_avg_n_Wrap(vtx, 0)
end do
printVarSummary(vt)
; vt = vt*1.e-4
vt_diff = vt(0,:,:)
vt_diff = vt(1,:,:) - vt(0,:,:)
print(vt(0,{29},{84.5}))
print(vt(1,{29},{84.5}))
print(vt_diff({29},{84.5}))
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = season+"_OMEGA"

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
res@cnFillPalette         = "MPL_RdBu"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035

res@gsnRightString = "m/s"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.2
res@cnMaxLevelValF = 0.2
res@cnLevelSpacingF = 0.02
res@gsnLeftString = "Ctrl"
plot(0) = gsn_csm_contour_map(wks, vt(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "Sen"
plot(1) = gsn_csm_contour_map(wks, vt(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.04
res@cnMaxLevelValF = 0.04
res@cnLevelSpacingF = 0.005
res@gsnLeftString = "Diff"
plot(2) = gsn_csm_contour_map(wks, vt_diff, res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean OMEGA (600hPa-300hPa) " ;(500hPa to 200hPa)
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/1,3/), pres)

end