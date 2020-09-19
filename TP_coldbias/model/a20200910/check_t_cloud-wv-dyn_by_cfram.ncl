
begin

  lonstr = 60 ;10
  lonend = 110 ;150
  latstr = 20 ;-10
  latend = 50 ;70

yrStrt = 198102
yrLast = 200601

title = (/"baseline","t_cloud","t_wv","t_dyn"/)
ntitle = dimsizes(title)

;================read metadata=============================================
diri  = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_Ctrl_FAMIP/a20191206/" 
file1 = "Ctrl_FAMIP.cam.h0.T.NoInterpolation.1979-2005.nc"
f1 = addfile(diri+file1,"r")
lon = tofloat(f1->lon)
lat = tofloat(f1->lat)
lev = tofloat(f1->lev)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)
delete(f1)

;===========================================================================
;================read metadata=============================================

nvname = (/"solar_in","swdn_sfc","swup_sfc","tsurf","huss","ps","o3","ta","hus","cld_amt","liq_wat","ice_wat"/)
nv_size = dimsizes(nvname)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram_decomposeT/" ;cfram20200805


UNDEF = -999.
x = new((/ntitle,6+6*nlev,nlat,nlon/),float,UNDEF)
x!0 = "var"
x!1 = "lev"
x!2 = "lat"
x!3 = "lon"
x&var = title
x&lat = lat
x&lon = lon
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"

do ll = 0,ntitle-1
    fili = diri2+title(ll)+"_input_1.dat"
    print(fili)
    do i = 0,5
        x(ll,i,:,:) = fbindirread(fili, i, (/nlat,nlon/), "float")
    end do
    do i = 6, nv_size-1
        do j = 0, nlev-1
            irec = 6+(i-6)*nlev+j
            x(ll,irec,:,:) = fbindirread(fili, irec, (/nlat,nlon/), "float")
        end do
    end do
end do

printVarSummary(x)

x_diff = x(0:ntitle-2,:,:,:)
do i = 0, ntitle-2
    x_diff(i,:,:,:) = x(i+1,:,:,:)-x(0,:,:,:)
end do
printVarSummary(x_diff)


;=================================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/" ;diri2
figname = "check_t_cloud-cloud-wv-dyn_by_cfram"
wks  = gsn_open_wks("pdf",figpath+figname)
gsn_define_colormap(wks,"temp_diff_18lev")

nplot = nv_size
plot = new(nplot, graphic)
plot_topo = plot


;=======================================
; this is for outline of TP
;=======================================

ff = "~/wuyt/data/TOPO_1500.nc"
topo_in = addfile(ff,"r")
topo = topo_in ->TP

mlat = topo_in->LATIXY(:,0)
mlon = topo_in->LONGXY(0,:)
topo!0 = "lsmlat"
topo!1 = "lsmlon"
topo&lsmlat = mlat
topo&lsmlon = mlon

;----------------------
res_topo            = True
res_topo@cnFillOn   = False
res_topo@cnLinesOn  = True
res_topo@gsnFrame   = False
res_topo@gsnDraw    = False
res_topo@cnLineLabelsOn         = False
res_topo@cnLineThicknessF       = 2
res_topo@cnLineColor            = "grey32"
;res_topo@cnLineDashPattern      = "1"
res_topo@cnInfoLabelOn          = False
res_topo@cnLevelSelectionMode   ="ExplicitLevels"
res_topo@cnLevels               = (/1500/)
;==========================================================

res = True
res@gsnDraw   = False
res@gsnFrame  = False
;res@gsnSpreadColors = True
res@gsnAddCyclic = False

res@mpFillOn  = True
;res@mpMinLonF = 180.
;res@mpMaxLonF = 360.
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
;res@mpCenterLonF = 85.
res@mpGeophysicalLineThicknessF = 2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False

;res@lbLabelAutoStride = True
;res@lbLabelFontHeightF = 0.01
;res@pmLabelBarOrthogonalPosF = 0.16  ; 水平位置
;res@pmLabelBarHeightF = 0.02    ;色标高度
;res@pmLabelBarWidthF = 0.7
;res@pmLabelBarHeightF = 0.15
;res@lbLabelFontHeightF = 0.016
res@lbLabelBarOn         = False 
res@lbLabelsOn           = False

res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
res@tmXTOn = False
res@tmYROn = False

;res@gsnLeftStringFontHeightF  = 0.02
;res@gsnRightStringFontHeightF  = 0.02
res@gsnStringFontHeightF = 0.035
; res@gsnRightString = "W/m~S~2~N~"

;-------------------------------------------------------

; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1,1.5,2,2.5,3/)
; res@cnFillColors = (/4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMaxLevelValF = 10
; res@cnMinLevelValF = -10
; res@cnLevelSpacingF = 1

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06


do ll = 0, ntitle-2
    do i = 0, 5
        res@gsnLeftString = nvname(i)
        plot(i) = gsn_csm_contour_map(wks, x_diff(ll,i,:,:), res)
        plot_topo(i) = gsn_csm_contour_map(wks, topo, res_topo)
        overlay(plot(i), plot_topo(i))
    end do
    do i = 6,nv_size-1
        res@gsnLeftString = nvname(i)
        jj = 6+(i-6)*nlev+20 ; 20th level
        plot(i) = gsn_csm_contour_map(wks, x_diff(ll,j,:,:), res)
        plot_topo(i) = gsn_csm_contour_map(wks, topo, res_topo)
        overlay(plot(i), plot_topo(i))
    end do

pres@gsnPanelMainString = title(ll+1)+" - baseline"
gsn_panel(wks,plot,(/3,4/),pres)
end do






end