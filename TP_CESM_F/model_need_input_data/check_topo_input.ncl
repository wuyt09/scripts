
begin

diri = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
fili = "USGS-gtopo30_1.9x2.5_remap_C20190410_TP1000m_wuyt.nc"

f = addfile(diri+fili, "r")

hgt = f->PHIS
printVarSummary(hgt)

;=================================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/model_need_input_data/"
figname = "check_topo_input"

wks = gsn_open_wks("pdf", figpath+figname)

;==================================================================
; spatial distribution at specfic levels
;==================================================================

nplot = 1
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
;res@lbLabelBarOn   = False

res@cnLineLabelsOn         = False
res@cnLineThicknessF       = 5
res@cnLineColor            = "black"
res@cnInfoLabelOn          = False
;res@cnLevelSelectionMode   ="ExplicitLevels"
;res@cnLevels               = (/0.9/)

res@gsnStringFontHeightF = 0.035

plot = gsn_csm_contour_map(wks, hgt, res)

draw(plot)
frame(wks)


end
