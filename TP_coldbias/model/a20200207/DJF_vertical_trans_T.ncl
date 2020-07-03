
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
diri1 = "/home/yangsong3/data-model/wuyt/TPbias/TEMP/"

do i = 0, ncase-1
    diri2 = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    ;fili1 = casename(i)+".VgradT.1981-2005.mon.mean.nc"
    fili1 = casename(i)+".Vertical.1981-2005.mon.mean.nc"
    fili2 = casename(i)+".cam.h0.PS.1979-2005.nc"
    f1 = addfile(diri1+fili1, "r")
    f2 = addfile(diri2+fili2, "r")
    
    TIME = f2->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    VT = f1->P3
    PS = f2->PS(iYYYY,:,:)

    lev = VT&lev
    plev = lev
    plev = plev *100
    plev!0 = "lev"
    plev&lev = lev
    plev@units = "Pa"

    ptop = min(plev)
    psfc = PS
    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    ; printVarSummary(dP)
    ; exit

    VTdP = VT*dP
    copy_VarCoords(VT, VTdP)
    ; print(VTdP(20,{850},{30},{90}))
    ; exit
    ; ; vt0 = dim_avg_n_Wrap(VTdP(:,{1000:200},:,:), 1)/dim_sum_n_Wrap(dP(:,{1000:200},:,:), 1)
    ; vt0 = dim_avg_n_Wrap(VTdP(:,{500:200},:,:), 1)/dim_sum_n_Wrap(dP(:,{500:200},:,:), 1)
    vt0 = dim_avg_n_Wrap(VTdP, 1)/dim_sum_n_Wrap(dP, 1)
    copy_VarCoords(VTdP(:,0,:,:), vt0)
    vtx = month_to_season(vt0, season)
    vt(i,:,:) = dim_avg_n_Wrap(vtx, 0)
end do
vt = vt*86400
vt@units = "K/day"
printVarSummary(vt)

vt_diff = vt(0,:,:)
vt_diff = vt(1,:,:) - vt(0,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200207/"
figname = season+"_vertical_trans_T"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(2, graphic)
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

res@gsnRightString = "K/day"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1.
res@cnMaxLevelValF = 1.
res@cnLevelSpacingF = 0.1
res@gsnLeftString = "Ctrl Clim"
plot(0) = gsn_csm_contour_map(wks, vt(0,:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -0.1
res@cnMaxLevelValF = 0.1
res@cnLevelSpacingF = 0.01
res@gsnLeftString = "Diff"
plot(1) = gsn_csm_contour_map(wks, vt_diff, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean vertical temp trans"
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,1/), pres)

end