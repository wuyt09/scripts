
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 1979
yrLast = 2015

season = "DJF"


diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
fili1 = "air.mon.mean.nc"
fili2 = "Surface_pressure.mon.mean.nc"
f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")

TIME = f2->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

VT = short2flt(f1->t(iYYYY,:,:,:))
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

PS = short2flt(f2->sp(iYYYY,:,:))

lev = VT&level
plev = lev
plev = plev *100
plev!0 = "lev"
plev&lev = lev
plev@units = "Pa"

ptop = min(plev)
psfc = PS
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
printVarSummary(dP)
; exit

VTdP = VT*dP
copy_VarCoords(VT, VTdP)
print(VTdP(20,{600},{30},{90}))
print(dP(20,{600},{30},{90}))
exit

vt0 = dim_sum_n_Wrap(VTdP(:,{1000:300},:,:), 1)/dim_sum_n_Wrap(dP(:,{1000:300},:,:), 1)
copy_VarCoords(VTdP(:,0,:,:), vt0)

vtx = month_to_season(vt0, season)
vt = dim_avg_n_Wrap(vtx, 0)
printVarSummary(vt)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_AtmosT_ERAI"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
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
plot(0) = gsn_csm_contour_map(wks, vt, res) ;(0,:,:)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
draw(plot(0))
exit
res@gsnLeftString = "Sen"
plot(1) = gsn_csm_contour_map(wks, vt(1,:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -1
; res@cnMaxLevelValF = 1
; res@cnLevelSpacingF = 0.1
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