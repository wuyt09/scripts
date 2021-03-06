
begin

diri = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
fili = "USGS-gtopo30_1.9x2.5_remap_c050602.nc"

f = addfile(diri+fili, "w")

phis = f->PHIS

printVarSummary(phis)

tpmask = phis
tpmask = 0

hgt = phis

lat = phis&lat
lon = phis&lon
printVarSummary(lat)
print(lon)

do j = 0,dimsizes(lat)-1
    do i = 0, 34  ; 85E
        case1 = lon(i).ge.50 .and. lon(i).le.110
        case2 = lat(j).ge.20 .and. lat(j).le.50
        case3 = phis(j,i).ge.(1500*9.8)
        if case1 .and. case2 .and. case3
            tpmask(j,i) = 1
        else 
            tpmask(j,i) = 0
        end if
    end do
end do

do j = 0,dimsizes(lat)-1
    do i = 35, dimsizes(lon)-1
        case1 = lon(i).ge.50 .and. lon(i).le.110
        case2 = lat(j).ge.20 .and. lat(j).le.40
        case3 = phis(j,i).ge.(1500*9.8)
        if case1 .and. case2 .and. case3
            tpmask(j,i) = 1
        else 
            tpmask(j,i) = 0
        end if
    end do
end do
; do j = 0, dimsizes(lat)-1
;    do i = 0,dimsizes(lon)-1
;        case1 = lon(i).ge.62 .and. lon(i).le.105
;        case2 = lat(j).ge.23 .and. lat(j).le.45
;        case3 = phis(j,i).ge.(1500*9.8)
;        if case1 .and. case2 .and. case3
;            tpmask(j,i) = 1
;        else 
;            tpmask(j,i) = 0
;        end if
;    end do
; end do
x = tpmask
tpmask = smth9_Wrap(x, 0.5, 0.25, True)

tpmask@long_name = "TP mask, TP area(1500m)=1, others=0"
printVarSummary(tpmask)
printVarSummary(hgt)

;===============================================================
outpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
maskname = "TPmask_1500m_C20200713.nc"
;hgtname = "USGS-gtopo30_1.9x2.5_remap_C20190304_wuyt.nc"

system("rm -f "+outpath+maskname)
;system("rm -f "+outpath+hgtname)

fm = addfile(outpath+maskname, "c")
;fh = addfile(outpath+hgtname, "c")

fm->tp = tpmask

;=================================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/model_need_input_data/"
figname = "write_TPmask"

wks = gsn_open_wks("ps", figpath+figname)
lonstr = 60
lonend = 110
latstr = 20
latend = 50
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
;res@lbLabelBarOn   = False

res@cnLineLabelsOn         = False
res@cnLineThicknessF       = 5
res@cnLineColor            = "black"
res@cnInfoLabelOn          = False
res@cnLevelSelectionMode   ="ExplicitLevels"
res@cnLevels               = (/0.5/)

res@gsnStringFontHeightF = 0.035

res@gsnLeftString = "TP mask"
plot(0) = gsn_csm_contour_map(wks, tpmask, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

delete(res@cnLevels)
res@cnLevels               = (/(1500*9.8)/)
res@gsnLeftString = "hgt"
plot(1) = gsn_csm_contour_map(wks, hgt, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))

pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/2,1/), pres)


end
