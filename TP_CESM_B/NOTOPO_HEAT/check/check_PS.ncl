
begin

diri = "./"
fili = "TP_hybrid_HEAT.cam.h0.0250-*.nc"
filis = systemfunc("ls "+diri+fili)
print(filis)

f = addfiles(filis, "r")

ps = f[:]->PS
printVarSummary(ps)

;=================================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/NOTOPO_HEAT/check/"
figname = "check_PS"

wks = gsn_open_wks("ps", figpath+figname)

;==================================================================
; spatial distribution at specfic levels
;==================================================================

nplot = 12
plot = new(nplot,graphic)
plot_topo = plot
gsn_define_colormap(wks,"cmp_flux")

 ;------------------------------------------------------------
 ; TOPO
 ;------------------------------------------------------------
 
 ff = "~/wuyt/data/TOPO_1500.nc"
 topo_in = addfile(ff,"r")
 topo = topo_in->TP

 xlat = topo_in->LATIXY(:,0)
 xlon = topo_in->LONGXY(0,:)
 topo!0 = "lsmlat"
 topo!1 = "lsmlon"
 topo&lsmlat = xlat
 topo&lsmlon = xlon
 
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 2
 res_topo@cnLineColor            = "grey4"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)

;---------------------------------------
latstr = 20
latend = 50
lonstr = 60
lonend = 110

res                 =  True
res@gsnDraw         = False
res@gsnFrame        = False
res@gsnAddCyclic    = True
;res@gsnSpreadColors = True

res@mpOutlineOn  = True
res@mpFillOn     = True
;res@mpCenterLonF = 0.
res@mpMinLatF    = latstr
res@mpMaxLatF    = latend
res@mpMinLonF    = lonstr
res@mpMaxLonF    = lonend
;res@mpCenterLonF = 180
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn     = True
res@cnLinesOn    = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn  = False
res@lbLabelBarOn   = False

res@cnLineLabelsOn         = False
res@cnLineThicknessF       = 5
res@cnLineColor            = "black"
res@cnInfoLabelOn          = False
;res@cnLevelSelectionMode   ="ExplicitLevels"
;res@cnLevels               = (/0.9/)
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 88000 ;100000
res@cnMinLevelValF = 84000 ;60000
res@cnLevelSpacingF = 500 ;4000

res@gsnStringFontHeightF = 0.035
do i = 0, 11
    plot(i) = gsn_csm_contour_map(wks, ps(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "surface pressure"
gsn_panel(wks, plot, (/3,4/), pres)


end

