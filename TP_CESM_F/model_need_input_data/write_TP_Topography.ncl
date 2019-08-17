
begin

diri = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
;fili = "USGS-gtopo30_1.9x2.5_remap_C20190525_TP1000m_wuyt.nc"
fili = "USGS-gtopo30_1.9x2.5_remap_c050602.nc"
filitp = "TPmask_1000m_C20190525.nc"

f = addfile(diri+fili, "r")
ftp = addfile(diri+filitp, "r")
phis = f->PHIS
tpmask = ftp->tp

printVarSummary(phis)


hgt = phis  ;! it should not set to zero

lat = phis&lat
lon = phis&lon
printVarSummary(lat)

outline = phis
outline = 0

x1d = ndtooned(tpmask)
i1d = ind(x1d.ge.0.9)
printVarSummary(i1d)
indx = ind_resolve(i1d, dimsizes(tpmask))
printVarSummary(indx)

print(indx(:,1))
print(indx(:,0))
print(" "+lat(indx(:,0))+"  "+lon(indx(:,1)))

strt = min(indx(:,0))
last = max(indx(:,0))
count = 0
size = 2*(last - strt +1)
xx = new(size, integer)
do i = strt,last,1
    j = num(indx(:,0).eq.i)
    xx(2*(i-strt)) = count
    count = count+j
    xx(2*(i-strt)+1) = count - 1
end do
print(xx)

strt = min(indx(:,1))
last = max(indx(:,1))
count = 0
size = 2*(last - strt +1)
yy = new(size, integer)
do i = strt,last,1
    j = num(indx(:,1).eq.i)
    yy(2*(i-strt)) = count
    count = count+j
    yy(2*(i-strt)+1) = count - 1
end do
print(yy)

outline(indx(xx,0),indx(xx,1)) = 1
outline(indx(yy,0),indx(yy,1)) = 1
printVarSummary(outline)

;============test====================
y1d = ndtooned(outline)
j1d = ind(y1d.ge.0.9)
indy = ind_resolve(j1d, dimsizes(tpmask))
print(" "+lat(indy(:,0))+"  "+lon(indy(:,1)))
;exit
;=================================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/model_need_input_data/"
figname = "write_TP_Topography"

wks = gsn_open_wks("ps", figpath+figname)

;==================================================================
; spatial distribution at specfic levels
;==================================================================

nplot = 2
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
res@mpFillOn     = False
;res@mpCenterLonF = 0.
;res@mpMinLatF    = latstr
;res@mpMaxLatF    = latend
;res@mpMinLonF    = lonstr
;res@mpMaxLonF    = lonend
;res@mpCenterLonF = 180
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn     = True
res@cnFillPalette = "WhiteBlueGreenYellowRed"
res@cnLinesOn    = True
res@cnLineLabelsOn = True
res@cnInfoLabelOn  = False
;res@lbLabelBarOn   = False

res@cnLineThicknessF       = 2
res@cnLineColor            = "green4"
res@cnInfoLabelOn          = False
res@cnLevelSelectionMode   ="ExplicitLevels"
res@cnLevels               = (/0.9/)

res@gsnStringFontHeightF = 0.035

res@gsnLeftString = "TP outline"
plot(0) = gsn_csm_contour_map(wks, outline, res)

res@cnLineColor            = "red"
res@cnLineDashPattern = 2
plot(1) = gsn_csm_contour(wks, tpmask, res)
overlay(plot(0), plot(1))

plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

delete(plot(1))
delete(res@cnLevels)
res@cnLevels               = (/(1000*9.8)/)
res@gsnLeftString = "hgt"
plot(1) = gsn_csm_contour_map(wks, hgt, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/2,1/), pres)


end