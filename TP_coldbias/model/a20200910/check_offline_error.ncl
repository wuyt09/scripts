
begin

  lonstr = 60 ;10
  lonend = 110 ;150
  latstr = 20 ;-10
  latend = 50 ;70


;================read metadata=============================================
casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)
nlat = 96
nlon = 144
yrStrt = 198102
yrLast = 200601
season = "DJF"

ssrd = new((/ncase,nlat,nlon/), float)
ssr = ssrd

do i = 0, ncase-1
    diri  = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/" 
    file1 = casename(i)+".cam.h0.FLDS.1979-2005.nc"
    ; file2 = casename(i)+".cam.h0.FLNS.1979-2005.nc"
    f1 = addfile(diri+file1,"r")
    ; f2 = addfile(diri+file2,"r")

    Time = f1->time
    YYYY = cd_calendar(Time, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    SSRD = f1->FLDS(iYYYY,:,:)
    ; SSR  = f2->FLNS(iYYYY,:,:)

    SSRDx = month_to_season(SSRD, season)
    ssrd(i,:,:) = dim_avg_n_Wrap(SSRDx, 0)
    ; ssr(i,:,:)  = month_to_season(SSR, season)
end do
printVarSummary(ssrd)
    
; ; level = tofloat(f1->lev)
; lon = tofloat(f1->lon)
; lat = tofloat(f1->lat); nlev = 27
; delete(f1)
lat = ssrd&lat
lon = ssrd&lon

;===========================================================================

nvname = (/"rad_conv","lwu","lwd","swu","swd","lwu_toa","lwd_toa","swu_toa","swd_toa"/)  ; ,"error"
nvars = dimsizes(nvname)

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/cfram/" ;cfram20200805
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
    ; cld(nv,:,:) = fbindirread(fili3,nv,(/nlat,nlon/),"float")
    ; wv(nv,:,:) = fbindirread(fili4,nv,(/nlat,nlon/),"float")
    ; t(nv,:,:) = fbindirread(fili5,nv,(/nlat,nlon/),"float")
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

pp = new((/4,nlat,nlon/), float)
pp(0,:,:) = ssrd(0,:,:)
pp(1,:,:) = x(nx,:,:)
pp(2,:,:) = (/ssrd(1,:,:)/)
pp(3,:,:) = y(nx,:,:)
copy_VarCoords(ssrd(0,:,:), pp(0,:,:))

diff = pp
diff(0,:,:) = pp(1,:,:) - pp(0,:,:)
diff(1,:,:) = pp(3,:,:) - pp(2,:,:)
diff(2,:,:) = pp(2,:,:) - pp(0,:,:)
diff(3,:,:) = pp(3,:,:) - pp(1,:,:)
copy_VarCoords(pp, diff)

;==========================================================
;****************************************************
; plot
;****************************************************

figpath = "./" ;diri2
figname = "check_offline_error"

wks  = gsn_open_wks("pdf",figpath+figname)
gsn_define_colormap(wks,"temp_diff_18lev")

title = (/"(a) Ctrl_DownR","(b) Base_DownR","(c) Hist_DownR","(d) Warm_DownR",\
         "(d) b-a","(e) d-c","(f) c-a","(g) d-b"/)

nplot = dimsizes(title)
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
; res@lbLabelBarOn         = False 
; res@lbLabelsOn           = False

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
res@cnMaxLevelValF = 400
res@cnMinLevelValF = 0
res@cnLevelSpacingF = 20

do i = 0,3  
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_contour_map(wks, pp(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do



res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 0
res@cnMinLevelValF = -200
res@cnLevelSpacingF = 10
res@gsnLeftString = title(4)
plot(4) = gsn_csm_contour_map(wks, diff(0,:,:), res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@gsnLeftString = title(5)
plot(5) = gsn_csm_contour_map(wks, diff(1,:,:), res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10
res@cnLevelSpacingF = 1
res@gsnLeftString = title(6)
plot(6) = gsn_csm_contour_map(wks, diff(2,:,:), res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6), plot_topo(6))

res@gsnLeftString = title(7)
plot(7) = gsn_csm_contour_map(wks, diff(3,:,:), res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7), plot_topo(7))

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
; pres@gsnPanelLabelBar = True

; pres@lbLabelFontHeightF = 0.015
; pres@pmLabelBarWidthF = 0.8    ;plot manager 
; pres@pmLabelBarHeightF = 0.06

gsn_panel(wks,plot,(/2,4/),pres)



end

