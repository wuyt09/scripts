
begin

diri = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_F2000_F19G16_CAM4/pro/"
fili = "heat_float_input.nc"

f = addfile(diri+fili, "r")
heat = f->heat
heat = heat*86400
printVarSummary(heat)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/check_heating/"
figname = "4_check_input_heating_tpmask"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 3
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
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "black"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)

;---------------------------------------

res                 =  True
res@gsnDraw         = False
res@gsnFrame        = False
res@gsnAddCyclic    = True
;res@gsnSpreadColors = True

res@mpOutlineOn  = True
res@mpFillOn     = True
;res@mpCenterLonF = 0.
;res@mpMinLatF    = latstr
;res@mpMaxLatF    = latend
;res@mpMinLonF    = lonstr
;res@mpMaxLonF    = lonend
;res@mpCenterLonF = 180
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn     = True
res@cnLinesOn    = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn  = False
res@lbLabelBarOn   = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF  = 5  ;5e-5
res@cnMinLevelValF  = -5  ;-5e-5
res@cnLevelSpacingF = 0.5  ;5e-6

res@gsnStringFontHeightF = 0.035

res@gsnLeftString = "850hPa"
plot(0) = gsn_csm_contour_map(wks, heat(2,{850},:,:), res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
res@gsnLeftString = "500hPa"
plot(1) = gsn_csm_contour_map(wks, heat(2,{500},:,:), res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
res@gsnLeftString = "200hPa"
plot(2) = gsn_csm_contour_map(wks, heat(2,{200},:,:), res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))


pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/3,1/), pres)

end