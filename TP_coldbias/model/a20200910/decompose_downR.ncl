
begin

  lonstr = 60 ;10
  lonend = 110 ;150
  latstr = 20 ;-10
  latend = 50 ;70

yrStrt = 198102
yrLast = 200601

nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)

;================read metadata=============================================

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/) ;,"Hist_IO","Hist_PO","Hist_AO","Hist_PIO","Hist_PAO","Hist_IAO"

ncase = dimsizes(casename)
season = "DJF"

ssrd = new((/ncase,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili = casename(i)+".cam.h0.FLDS.1979-2005.nc"
    f = addfile(diri+fili, "r")

    TIME = f->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    SSRD = f->FLDS(iYYYY,:,:)
    ssrdx = month_to_season(SSRD, season)
    ssrd(i,:,:) = dim_avg_n_Wrap(ssrdx, 0)
end do
printVarSummary(ssrd)

lon = tofloat(ssrd&lon)
lat = tofloat(ssrd&lat)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
delete(f)

;===========================================================================

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram20200917/" ;cfram20200805
fili2 = diri2+"baseline_radsfc_ranc_1.grd"
fili3 = diri2+"cloud_radsfc_ranc_1.grd"
fili4 = diri2+"wv_radsfc_ranc_1.grd"
fili5 = diri2+"t_radsfc_ranc_1.grd"
fili6 = diri2+"warm_radsfc_ranc_1.grd"


UNDEF = -999.
x = new((/nvars,nlat,nlon/),float,UNDEF)
x!0 = "var"
x!1 = "lat"
x!2 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x@long_name = "forcing"
x@units     = "W/m2"

cld = x
wv = x
t = x
y = x

do nv = 0,nvars-1
    x(nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
    cld(nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
    wv(nv,:,:) = fbindirread(fili4,nv,(/nlat,nlon/),"float")
    t(nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
    y(nv,:,:) = fbindirread(fili6,nv,(/nlat,nlon/),"float")
end do

; if (any(isnan_ieee(x))) then   ; if x got a element is Not A Number, then the Fill_value
;   value = -999.
;   replace_ieeenan (x,value,0)
;   x@_FillValue = value
; end if
printVarSummary(x)
printMinMax(x(0,:,:),1)


nx = 2
title = (/"(a) Cloud","(b) WV","(c) Temp","(d) a+b+c","(e) model total","(f) error: e-d"/)
x_diff = new((/dimsizes(title),nlat,nlon/), float)
x_diff!0 = "term"
x_diff&term = title

x_diff(0,:,:) = cld(nx,:,:) - x(nx,:,:)
x_diff(1,:,:) = wv(nx,:,:) - x(nx,:,:)
x_diff(2,:,:) = t(nx,:,:) - x(nx,:,:)
x_diff(3,:,:) = x_diff(0,:,:) + x_diff(1,:,:) + x_diff(2,:,:)
x_diff(4,:,:) = ssrd(1,:,:) - ssrd(0,:,:) ;y(nx,:,:) - x(nx,:,:)
x_diff(5,:,:) = x_diff(4,:,:) - x_diff(3,:,:)
copy_VarCoords(x(0,:,:), x_diff(0,:,:))
printVarSummary(x_diff)

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

TOPO = linint2_Wrap(mlon, mlat, topo, True, lon, lat, 0)

; p_cloud = x_diff(0,{latstr:latend},{lonstr:lonend})
; p_error = x_diff(5,{latstr:latend},{lonstr:lonend})
; p_cloud = mask(p_cloud, TOPO({latstr:latend},{lonstr:lonend}).ge.1500, True)
; p_error = mask(p_error, TOPO({latstr:latend},{lonstr:lonend}).ge.1500, True)
p_cloud = x_diff(0,:,:)
p_error = x_diff(5,:,:)
p_cloud = mask(p_cloud, TOPO.ge.1500, True)
p_error = mask(p_error, TOPO.ge.1500, True)
pa_corr = pattern_cor(p_cloud({latstr:latend},{lonstr:lonend}), p_error({latstr:latend},{lonstr:lonend}), 1., 1)
print(pa_corr)

x_med = x_diff(0,:,:)
x_diff(0,:,:) = (1+pa_corr)*x_med
x_diff(5,:,:) = x_diff(5,:,:) - pa_corr*x_med
x_diff(3,:,:) = x_diff(0,:,:) + x_diff(1,:,:) + x_diff(2,:,:)

;==========================================================
;****************************************************
; plot
;****************************************************

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/" ;diri2
figname = "decompose_downR"
wks  = gsn_open_wks("pdf",figpath+figname)
gsn_define_colormap(wks,"temp_diff_18lev")

nplot = dimsizes(title)
plot = new(nplot, graphic)
plot_topo = plot

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
res@gsnRightString = "W/m~S~2~N~"

;-------------------------------------------------------

; res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1,1.5,2,2.5,3/)
; res@cnFillColors = (/4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10
res@cnLevelSpacingF = 1

do i = 0, dimsizes(title)-1
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_contour_map(wks, x_diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

; i = dimsizes(title)-1
; res@gsnLeftString = title(i)
; plot(i) = gsn_csm_contour_map(wks, x_diff(i,:,:), res)
; plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(i), plot_topo(i))

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06

gsn_panel(wks,plot,(/2,3/),pres)

end 