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

do j = 0,dimsizes(lat)-1
    do i = 0, 36  ; 90E
        case1 = lon(i).ge.70 .and. lon(i).le.105
        case2 = lat(j).ge.23 .and. lat(j).le.40
        case3 = phis(j,i).ge.(1000*9.8)
        if case1 .and. case2 .and. case3
            tpmask(j,i) = 1
        else 
            tpmask(j,i) = 0
        end if
    end do
end do

do j = 0, dimsizes(lat)-1
    do i = 37,dimsizes(lon)-1
        case1 = lon(i).ge.70 .and. lon(i).le.105
        case2 = lat(j).ge.23 .and. lat(j).le.40
        case3 = phis(j,i).ge.(1000*9.8)
        if case1 .and. case2 .and. case3
            tpmask(j,i) = 1
        else 
            tpmask(j,i) = 0
        end if
    end do
end do

do j = 0,dimsizes(lat)-1
    do i = 0, dimsizes(lon)-1
        case1 = lon(i).ge.60 .and. lon(i).le.90
        case2 = lat(j).ge.5 .and. lat(j).le.35
        case3 = phis(j,i).le.(1000*9.8)
        if case1 .and. case2 .and. case3
            tpmask(j,i) = 1
        ;else 
        ;    tpmask(j,i) = 0
        end if
    end do
end do

x = tpmask
;do j = 0, dimsizes(lat)-1
;    do i = 0, dimsizes(lon)-1
        tpmask = smth9_Wrap(x, 0.5, 0.25, True)
        ; write loop content
;    end do
    ; write loop content
;end do
x = tpmask
;do j = 0, dimsizes(lat)-1
;    do i = 0, dimsizes(lon)-1
        tpmask = smth9_Wrap(x, 0.5, 0.25, True)
        ; write loop content
;    end do
    ; write loop content
;end do
x = tpmask
;do j = 0, dimsizes(lat)-1
;    do i = 0, dimsizes(lon)-1
        tpmask = smth9_Wrap(x, 0.5, 0.25, True)
        ; write loop content
;    end do
    ; write loop content
;end do


;;=====================================================================
filiq = "./Q1_JJA.nc"
fq = addfile(filiq, "r")
Q1 = fq->q1
q1sea = dim_avg_n_Wrap(Q1, 0)
; Cp = 1004.07 ; J/(Kg*K)
; Season = "JJA"
 lonstr = 40
 lonend = 180
 latstr = -10
 latend = 70

; diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
; filiu1 = "F2000_CTRL.plevel.DTCOND.mon.0601-3012.nc"
; filiv1 = "F2000_CTRL.plevel.DTV.mon.0601-3012.nc"
; filis1 = "F2000_CTRL.plevel.QRS.mon.0601-3012.nc"
; filil1 = "F2000_CTRL.plevel.QRL.mon.0601-3012.nc"
; filisp = "F2000_CTRL.PS.mon.0601-3012.nc"

; fu1 = addfile(diri1+filiu1, "r")
; fv1 = addfile(diri1+filiv1, "r")
; fs1 = addfile(diri1+filis1, "r")
; fl1 = addfile(diri1+filil1, "r")
; fsp = addfile(diri1+filisp, "r")

; ;u1 = fu1->DTCOND(:,:,{latstr:latend},{lonstr:lonend})
; ;v1 = fv1->DTV(:,:,{latstr:latend},{lonstr:lonend})
; ;s1 = fs1->QRS(:,:,{latstr:latend},{lonstr:lonend})
; ;l1 = fl1->QRL(:,:,{latstr:latend},{lonstr:lonend})
; ;sp = fsp->PS(:,{latstr:latend},{lonstr:lonend})

; u1 = fu1->DTCOND;(:,:,{latstr:latend},{lonstr:lonend})
; v1 = fv1->DTV;(:,:,{latstr:latend},{lonstr:lonend})
; s1 = fs1->QRS;(:,:,{latstr:latend},{lonstr:lonend})
; l1 = fl1->QRL;(:,:,{latstr:latend},{lonstr:lonend})
; sp = fsp->PS;(:,{latstr:latend},{lonstr:lonend})

; q1 = u1
; q1 = u1+v1+s1+l1

; q1 = q1 * Cp ; units W/Kg

; ;----------------integration-----------------------------
; plev = q1&lev*100.
; psfc = sp
; ptop = min(plev)
; dp = dpres_plevel_Wrap(plev, sp, ptop, 0)
; ;printVarSummary(dP)

; ;dp = conform_dims(dimsizes(q1), dP, 1)
; q1dp = q1 * dp
; copy_VarCoords(q1, q1dp)
; printVarSummary(q1dp)
; Q1 = dim_sum_n_Wrap(q1dp, 1)
; copy_VarCoords(q1(:,0,:,:), Q1)
; printVarSummary(Q1)
; Q1 = Q1/9.8  ; units W/m2
; ;---------------------------------------------------------

; q1Sea = month_to_season(Q1, Season)
; printVarSummary(q1Sea)
; q1sea = dim_avg_n_Wrap(q1Sea, 0)
; ;q1sea = q1Sea(:,22,:,:)
; printVarSummary(q1sea)

; system("rm -rf Q1_JJA.nc")
; filiout = addfile("./Q1_JJA.nc", "c")
; filiout->q1 = q1Sea

do j = 0,dimsizes(lat)-1
    do i = 0, dimsizes(lon)-1
        case1 = lon(i).ge.70 .and. lon(i).le.100
        case2 = lat(j).ge.25 .and. lat(j).le.35
        case3 = q1sea(j,i).ge.250.
        if case1 .and. case2 .and. case3
            tpmask(j,i) = 0
        ;else 
        ;    tpmask(j,i) = 0
        end if
    end do
end do

;===============================================================
outpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
maskname = "Slope_mask_1000m_C20190710.nc"
;hgtname = "USGS-gtopo30_1.9x2.5_remap_C20190304_wuyt.nc"

system("rm -f "+outpath+maskname)
;system("rm -f "+outpath+hgtname)

fm = addfile(outpath+maskname, "c")
;fh = addfile(outpath+hgtname, "c")

fm->slope = tpmask

;=================================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/model_need_input_data/"
figname = "write_southslope_mask"

wks = gsn_open_wks("pdf", figpath+figname)

;==================================================================
; spatial distribution at specfic levels
;==================================================================

nplot = 2
plot = new(nplot,graphic)
plot_topo = plot
gsn_define_colormap(wks,"BlueWhiteOrangeRed")

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

res@cnFillOn     = False
res@cnLinesOn    = True
res@cnLineLabelsOn = False
res@cnInfoLabelOn  = False
;res@lbLabelBarOn   = False
res@cnFillPalette = "BlueWhiteOrangeRed"

res@cnLineLabelsOn         = False
res@cnLineThicknessF       = 5
res@cnLineColor            = "black"
res@cnInfoLabelOn          = False
res@cnLevelSelectionMode   ="ExplicitLevels"
res@cnLevels               = (/0.2,0.4,0.6,0.8,0.9/)

res@gsnStringFontHeightF = 0.035

res@gsnLeftString = "slope mask"
plot(0) = gsn_csm_contour_map(wks, tpmask, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

draw(plot)
frame(wks)
end