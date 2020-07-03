
begin

lonstr = 50
lonend = 150
latstr = 10
latend = 60

diri = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
; fili1 = "USGS-gtopo30_1.9x2.5_remap_C20190525_TP1000m_wuyt.nc"
; fili2 = "USGS-gtopo30_1.9x2.5_remap_c050602.nc"
; fili3 = "TPmask_1000m_C20190525.nc"
fili1 = "USGS-gtopo30_1.9x2.5_remap_C20190410_TP1000m_wuyt.nc"
fili2 = "USGS-gtopo30_1.9x2.5_remap_c050602.nc"
fili3 = "TPmask_1000m_C20190410.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")

hgt1 = f1->PHIS
hgt1 = hgt1/9.8
hgt2 = f2->PHIS
hgt2 = hgt2/9.8
tpmask = f3->tp
printVarSummary(hgt1)

dirix = "/home/yangsong3/wuyt/data/"
filix = "TOPO.nc"
fz = addfile(dirix+filix, "r")
zx = fz->TP
printVarSummary(zx)
zx&lsmlon = fz->LONGXY(0,:)
zx&lsmlat = fz->LATIXY(:,0)
printVarSummary(zx)
;=================================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/model_need_input_data/"
figname = "check_topo_input"

wks = gsn_open_wks("pdf", figpath+figname)

;==================================================================
; spatial distribution at specfic levels
;==================================================================

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
res@mpMinLatF    = latstr
res@mpMaxLatF    = latend
res@mpMinLonF    = lonstr
res@mpMaxLonF    = lonend
;res@mpCenterLonF = 180
res@mpGeophysicalLineThicknessF = 1

res@cnFillOn     = False ;True
res@cnLinesOn    = True
res@cnLineLabelsOn         = True
res@cnInfoLabelOn          = False
res@cnLabelDrawOrder = "PostDraw"
res@cnLineLabelFontHeightF = 0.025
res@cnLineLabelFontColor = "red"

res@cnLineThicknessF      = 3
;res@cnLineColor            = "black"
res@cnLineDashPatterns = (/0,2,4/)
res@cnLineColor = "green4"
;res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.035
res@gsnRightString = "m"

res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/500,1000,1500,2000/)

res@gsnLeftString = "CTRL"
plot(0) = gsn_csm_contour_map(wks, hgt2, res)

res@gsnLeftString = "NOTP"
plot(1) = gsn_csm_contour_map(wks, zx, res)

delete(res@cnLevels)
res@cnLevelSelectionMode   ="ExplicitLevels"
res@cnLevels               = (/0.9/)
res@gsnLeftString = "Mask data"
plot(2) = gsn_csm_contour_map(wks, tpmask, res)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/3,1/), pres)


end
