;*****************************************
; compare with xiaoming's paper figures
;*****************************************

begin
  lonstr = 60 ;10
  lonend = 110 ;150
  latstr = 20 ;-10
  latend = 50 ;70
  ; lonstr = 10
  ; lonend = 150
  ; latstr = -10
  ; latend = 70

nvname = (/"albedo","wv","cloud","cloud_sw","cloud_lw","co2","o3","solar",\
           "dyn","atm_dyn","ocean_dyn","shflx","lhflx"/)  ; ,"error"

;================read metadata=============================================
diri1 = "./"
fili1 = "Ctrl_FAMIP.cam.h0.2005-12.nc"
f1 = addfile(diri1+fili1,"r")
lon = tofloat(f1->lon)
lat = tofloat(f1->lat)
level = tofloat(f1->lev)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = 27
nvars = dimsizes(nvname)
delete(f1)
;===========================================================================

diri2 = "../../"
; file2 = diri2+"20200416/forcing_1.grd"
file2 = diri2+"20200416/partial_T_1.grd"
file3 = diri2+"data/t2m_warm.dat"
file4 = diri2+"data/t2m_base.dat"



lev = new(nlev,float,"No_FillValue")
lev(0:nlev-2) = level
lev(nlev-1) = 1001

setfileoption("bin","ReadByteOrder","Native")

UNDEF = -999.
x = new((/nvars,nlev,nlat,nlon/),float,UNDEF)
x!0 = "var"
x!1 = "lev"
x!2 = "lat"
x!3 = "lon"
x&var = nvname
x&lat = lat
x&lon = lon
x&lev = lev
x&lat@units = "degrees_north"
x&lon@units = "degrees_east"
x&lev@units = "hPa"
x@long_name = "Delta T"
x@units     = "K"

do nv = 0,nvars-1
   do nl = 0,nlev-1
      x(nv,nl,:,:) = fbindirread(file2,nv*nlev+nl,(/nlat,nlon/),"float")
   end do
end do

if (any(isnan_ieee(x))) then   ; if x got a element is Not A Number, then the Fill_value
  value = -999.
  replace_ieeenan (x,value,0)
  x@_FillValue = value
end if
printVarSummary(x)
printMinMax(x(0,:,:,:),1)
print(x(:,26,{50},{120}))

;===============================================================================================
obs = new((/2,nlat,nlon/),float,UNDEF)
obs(0,:,:) = fbindirread(file3,0,(/nlat,nlon/),"float")
obs(1,:,:) = fbindirread(file4,0,(/nlat,nlon/),"float")
copy_VarMeta(x(0,0,:,:),obs(0,:,:))
printVarSummary(obs)
print("end of data read")
d_obs = obs(0,:,:) - obs(1,:,:)
copy_VarMeta(obs(0,:,:), d_obs)

total_x = dim_sum_n(x(0:2,:,:,:),0)+dim_sum_n(x(5:7,:,:,:),0)+dim_sum_n(x(9:12,:,:,:),0)
copy_VarMeta(x(0,:,:,:), total_x)
printVarSummary(total_x)

error = d_obs - total_x(26,:,:)
copy_VarMeta(total_x(26,:,:), error)
printVarSummary(error)
;print(error({35:50},{100:120}))


; total_x(17,:,:) = where(abs(error).lt.1.5,total_x(17,:,:),0)
; x(9,17,:,:)= where(abs(error).lt.1.5,x(9,17,:,:),0) ; atm_dyn
; x(10,17,:,:)= where(abs(error).lt.1.5,x(10,17,:,:),0) ; atm_dyn
; error = where(abs(error).gt.1.5,0,error)
; copy_VarMeta(x(0,37,:,:),error)

x_inn = x(:,26,:,:)

;==========================================================
;****************************************************
; plot
;****************************************************

figpath = "./"
wks  = gsn_open_wks("pdf",figpath+"draw_process_spatial_cfram")
gsn_define_colormap(wks,"temp_diff_18lev")

plot = new(10, graphic)
plot_topo  = new(10,graphic)
  dum0 = plot
  dum1 = plot
  dum2 = plot
  dum3 = plot

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
res_topo@cnLineThicknessF       = 5
res_topo@cnLineColor            = "black"
res_topo@cnInfoLabelOn          = False
res_topo@cnLevelSelectionMode   ="ExplicitLevels"
res_topo@cnLevels               = (/1500/)

;*****************************************************
; this is for CFRAM results
;*****************************************************

;==========================================================
; this is for the box of my focus area
;==========================================================

 gres = True
 gres@gsLineColor = "blue"
 gres@gsLineThicknessF = 3

 glatw = (/25,40,40,25,25/)
 glonw = (/40,40,60,60,40/)
 glate = (/25,40,40,25,25/)
 glone = (/110,110,120,120,110/)
 glatc = (/25,40,40,25,25/)
 glonc = (/75,75,105,105,75/)
 glatn = (/45,55,55,45,45/)
 glonn = (/75,75,105,105,75/)


;==========================================================

res = True
res@gsnDraw   = False
res@gsnFrame  = False
;res@gsnSpreadColors = True

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
res@gsnRightString = "K"

;-------------------------------------------------------

res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1,1.5,2,2.5,3/)
res@cnFillColors = (/4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)
;res@cnLevelSelectionMode = "ManualLevels"
;res@cnMaxLevelValF = 3
;res@cnMinLevelValF = -3
;res@cnLevelSpacingF = 0.2


res@gsnLeftString = "(a) Heat Storage"
plot(0) = gsn_csm_contour_map(wks,x_inn(10,:,:),res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0),plot_topo(0))
    dum0(0) = gsn_add_polyline(wks, plot(0), glonw, glatw, gres)
    dum1(0) = gsn_add_polyline(wks, plot(0), glone, glate, gres)
    ;dum2(0) = gsn_add_polyline(wks, plot(0), glonc, glatc, gres)
    ;dum3(0) = gsn_add_polyline(wks, plot(0), glonn, glatn, gres)


res@gsnLeftString = "(b) Surf.LH"
plot(1) = gsn_csm_contour_map(wks,x_inn(12,:,:),res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1),plot_topo(1))
    dum0(1) = gsn_add_polyline(wks, plot(1), glonw, glatw, gres)
    dum1(1) = gsn_add_polyline(wks, plot(1), glone, glate, gres)
    ;dum2(1) = gsn_add_polyline(wks, plot(1), glonc, glatc, gres)
    ;dum3(1) = gsn_add_polyline(wks, plot(1), glonn, glatn, gres)

res@gsnLeftString = "(c) Surf.SH"
plot(2) = gsn_csm_contour_map(wks,x_inn(11,:,:),res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2),plot_topo(2))
    dum0(2) = gsn_add_polyline(wks, plot(2), glonw, glatw, gres)
    dum1(2) = gsn_add_polyline(wks, plot(2), glone, glate, gres)
    ;dum2(2) = gsn_add_polyline(wks, plot(2), glonc, glatc, gres)
    ;dum3(2) = gsn_add_polyline(wks, plot(2), glonn, glatn, gres)


;**************************************************
;delete(res@cnMaxLevelValF)
;delete(res@cnMinLevelValF)
;res@cnMinLevelValF = -1.5
;res@cnMaxLevelValF = 1.5
;res@cnLevelSpacingF = 0.1

res@gsnLeftString = "(d) Albedo"
plot(3) = gsn_csm_contour_map(wks,x_inn(0,:,:),res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3),plot_topo(3))
    dum0(3) = gsn_add_polyline(wks, plot(3), glonw, glatw, gres)
    dum1(3) = gsn_add_polyline(wks, plot(3), glone, glate, gres)
    ;dum2(3) = gsn_add_polyline(wks, plot(3), glonc, glatc, gres)
    ;dum3(3) = gsn_add_polyline(wks, plot(3), glonn, glatn, gres)


res@gsnLeftString = "(e) Cloud"
plot(4) = gsn_csm_contour_map(wks,x_inn(2,:,:),res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4),plot_topo(4))
    dum0(4) = gsn_add_polyline(wks, plot(4), glonw, glatw, gres)
    dum1(4) = gsn_add_polyline(wks, plot(4), glone, glate, gres)
    ;dum2(4) = gsn_add_polyline(wks, plot(4), glonc, glatc, gres)
    ;dum3(4) = gsn_add_polyline(wks, plot(4), glonn, glatn, gres)


res@gsnLeftString = "(f) Atm_Dyn"
plot(5) = gsn_csm_contour_map(wks,x_inn(9,:,:),res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5),plot_topo(5))
    dum0(5) = gsn_add_polyline(wks, plot(5), glonw, glatw, gres)
    dum1(5) = gsn_add_polyline(wks, plot(5), glone, glate, gres)
    ;dum2(5) = gsn_add_polyline(wks, plot(5), glonc, glatc, gres)
    ;dum3(5) = gsn_add_polyline(wks, plot(5), glonn, glatn, gres)


;***************************************************

;res@cnMinLevelValF = -1
;res@cnMaxLevelValF = 1
;res@cnLevelSpacingF = 0.05

res@gsnLeftString = "(g) WV"
plot(6) = gsn_csm_contour_map(wks,x_inn(1,:,:),res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6),plot_topo(6))
    dum0(6) = gsn_add_polyline(wks, plot(6), glonw, glatw, gres)
    dum1(6) = gsn_add_polyline(wks, plot(6), glone, glate, gres)
    ;dum2(6) = gsn_add_polyline(wks, plot(6), glonc, glatc, gres)
    ;dum3(6) = gsn_add_polyline(wks, plot(6), glonn, glatn, gres)


res@gsnLeftString = "(h) error"
plot(7) = gsn_csm_contour_map(wks,error,res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7),plot_topo(7))
    dum0(7) = gsn_add_polyline(wks, plot(7), glonw, glatw, gres)
    dum1(7) = gsn_add_polyline(wks, plot(7), glone, glate, gres)
    ;dum2(7) = gsn_add_polyline(wks, plot(7), glonc, glatc, gres)
    ;dum3(7) = gsn_add_polyline(wks, plot(7), glonn, glatn, gres)


res@gsnLeftString = "(i) CO2"
plot(8) = gsn_csm_contour_map(wks,x_inn(5,:,:),res)
plot_topo(8) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(8),plot_topo(8))
    dum0(8) = gsn_add_polyline(wks, plot(8), glonw, glatw, gres)
    dum1(8) = gsn_add_polyline(wks, plot(8), glone, glate, gres)
    ;dum2(8) = gsn_add_polyline(wks, plot(8), glonc, glatc, gres)
    ;dum3(8) = gsn_add_polyline(wks, plot(8), glonn, glatn, gres)


res@gsnLeftString = "(j) SOLAR"
plot(9) = gsn_csm_contour_map(wks,x_inn(7,:,:),res)
plot_topo(9) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(9),plot_topo(9))
    dum0(9) = gsn_add_polyline(wks, plot(9), glonw, glatw, gres)
    dum1(9) = gsn_add_polyline(wks, plot(9), glone, glate, gres)
    ;dum2(8) = gsn_add_polyline(wks, plot(8), glonc, glatc, gres)
    ;dum3(8) = gsn_add_polyline(wks, plot(8), glonn, glatn, gres)
;------------------------------------------------------
; the following use for diff color tables
;------------------------------------------------------
;-------------------------------------------------------

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True ;False
pres@gsnDraw     = True ;False
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "DJF"  ;(MERRA data)

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06

; pplot = gsn_panel_return(wks,plot(0:2),(/1,3/),pres)
; pplot2 = gsn_panel_return(wks,plot(3:8),(/2,3/),pres)
; ;pplot3 = gsn_panel_return(wks,plot(6:8),(/1,3/),pres)

; bb = NhlGetBB(pplot)
; height = max(bb(0,:)-min(bb(:,1)))
; pres@gsnDraw = True

; pres@gsnPanelTop = 0.9
; pres@gsnPanelBottom = pres@gsnPanelTop - height
; gsn_panel(wks,plot(0:2),(/1,3/),pres)
; pres@gsnPanelMainString = ""
; pres@gsnPanelTop = pres@gsnPanelBottom
; pres@gsnPanelBottom = pres@gsnPanelTop - height*2
; gsn_panel(wks,plot(3:8),(/2,3/),pres)
; ;pres@gsnPanelTop = pres@gsnPanelBottom
; ;pres@gsnPanelBottom = pres@gsnPanelTop - height
; ;gsn_panel(wks,plot(6:8),(/1,3/),pres)

gsn_panel(wks,plot(0:7),(/3,3/),pres)


end


