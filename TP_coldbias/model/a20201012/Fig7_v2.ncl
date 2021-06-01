
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

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram_decomposeT/" ;cfram20200805
fili2 = diri2+"baseline_radsfc_ranc_1.grd"
fili3 = diri2+"t_cloud_radsfc_ranc_1.grd"
fili4 = diri2+"t_wv_radsfc_ranc_1.grd"
fili5 = diri2+"t_dyn_radsfc_ranc_1.grd"

diri3 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram20200917/"
fili6 = diri3+"t_radsfc_ranc_1.grd"


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
dyn = x
y = x

do nv = 0,nvars-1
    x(nv,:,:) = fbindirread(fili2,nv,(/nlat,nlon/),"float")
    cld(nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
    wv(nv,:,:) = fbindirread(fili4,nv,(/nlat,nlon/),"float")
    dyn(nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
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
title = (/"(a) Cloud","(b) WV","(c) Atm. dyn.","(d) a+b+c","(e) total by Air Temp.","(f) error: e-d"/)
x_diff = new((/dimsizes(title),nlat,nlon/), float)
x_diff!0 = "term"
x_diff&term = title

x_diff(0,:,:) = cld(nx,:,:) - x(nx,:,:)
x_diff(1,:,:) = wv(nx,:,:) - x(nx,:,:)
x_diff(2,:,:) = dyn(nx,:,:) - x(nx,:,:)
x_diff(3,:,:) = x_diff(0,:,:) + x_diff(1,:,:) + x_diff(2,:,:)
x_diff(4,:,:) = y(nx,:,:) - x(nx,:,:)
x_diff(5,:,:) = x_diff(4,:,:) - x_diff(3,:,:)
copy_VarCoords(x(0,:,:), x_diff(0,:,:))
printVarSummary(x_diff)

;=======================================
; this is for outline of TP
;=======================================

ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, lon, lat, 0)

topo = smth9_Wrap(topo, 0.5, 0.25, True)

;===================================================================================
; correction
;===================================================================================

clat = cos(0.01745329*lat)
error = x_diff(5,:,:)
cld_corr = x_diff(2,:,:)
error = mask(error, TOPO.ge.1500, True)
cld_corr = mask(cld_corr, TOPO.ge.1500, True)

; error = x_diff(5,{latstr:latend},{lonstr:lonend})
; cld_corr = x_diff(0,{latstr:latend},{lonstr:lonend})

upper_term = error*cld_corr
uterm = wgt_areaave_Wrap(upper_term, 1.0, 1.0, 0)
lower_term = cld_corr*cld_corr
lterm = wgt_areaave_Wrap(lower_term, 1.0, 1.0, 0)
alpha = uterm/lterm
print("alpha = "+alpha)

x_diff(2,:,:) = x_diff(2,:,:)*(1+alpha)
x_diff(3,:,:) = x_diff(0,:,:) + x_diff(1,:,:) + x_diff(2,:,:)
x_diff(5,:,:) = x_diff(5,:,:)-alpha*x_diff(0,:,:)
;==========================================================
;****************************************************
; plot
;****************************************************

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20201012/" ;diri2
figname = "Fig7_v2"
wks  = gsn_open_wks("pdf",figpath+figname)
gsn_define_colormap(wks,"BlueWhiteOrangeRed") ;BlueYellowRed

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
res_topo@cnLineThicknessF       = 3
res_topo@cnLineColor            = "grey40"
;res_topo@cnLineDashPattern      = "1"
res_topo@cnInfoLabelOn          = False
res_topo@cnLevelSelectionMode   ="ExplicitLevels"
res_topo@cnLevels               = (/1500/)
res_topo@gsnLeftString = ""
res_topo@gsnRightString = ""
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

res@cnLevelSelectionMode = "ExplicitLevels"
; res@cnLevels = (/-12,-9,-7,-5,-3,-2,-1,-0.5,0.5,1,2,3,5,7,9,12/) ;
; res@cnFillColors = (/10,26,35,47,53,66,74,90,122,154,170,177,189,200,218,234,250/)
res@cnLevels = (/-7,-5,-3,-2,-1,-0.5,0.5,1,2,3,5,7/) ;
res@cnFillColors = (/10,35,47,66,74,90,122,154,170,177,200,218,250/)
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMaxLevelValF = 5
; res@cnMinLevelValF = -5
; res@cnLevelSpacingF = 0.5

do i = 0, dimsizes(title)-1
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_contour_map(wks, x_diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

cnres = res
delete(cnres@mpMinLonF)
delete(cnres@mpMaxLonF)
delete(cnres@mpMaxLatF)
delete(cnres@mpMinLatF)
delete(cnres@mpFillOn)
delete(cnres@mpGeophysicalLineThicknessF)
cnres@gsnLeftString = ""
cnres@gsnRightString = ""
cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLineLabelFontColor = "red"
cnres@cnLineColor = "purple"
cnres@cnLineThicknessF = 4
cnres@cnInfoLabelOn = True
cnres@cnInfoLabelOrthogonalPosF = -0.1
cnres@cnInfoLabelFontHeightF = 0.02
cnres@cnLabelDrawOrder = "PostDraw"
; cnres@gsnContourNegLineDashPattern = 16
; cnres@gsnContourPosLineDashPattern = 0

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMaxLevelValF = 0
cnres@cnMinLevelValF = -5
cnres@cnLevelSpacingF = 1.5

plot(4) = gsn_csm_contour(wks, x_diff(4,:,:), cnres)
overlay(plot(3),plot(4))

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True

pres@lbLabelFontHeightF = 0.01
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06

gsn_panel(wks,plot,(/2,2/),pres)

end 