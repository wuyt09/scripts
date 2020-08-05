;*****************************************
; compare with xiaoming's paper figures
;*****************************************

begin
  lonstr = 60 ;10
  lonend = 110 ;150
  latstr = 20 ;-10
  latend = 50 ;70

nvname = (/"albedo","wv","cloud","cloud_sw","cloud_lw","co2","o3","solar",\
           "dyn","atm_dyn","ocean_dyn","shflx","lhflx","SR_Q"/)  ; ,"error"

;================read metadata=============================================
diri  = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_Ctrl_FAMIP/a20191206/" 
file1 = "Ctrl_FAMIP.cam.h0.T.NoInterpolation.1979-2005.nc"
f1 = addfile(diri+file1,"r")
lon = tofloat(f1->lon)
lat = tofloat(f1->lat)
level = tofloat(f1->lev)

nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = 27
nvars = dimsizes(nvname)
delete(f1)
;===========================================================================

diri2 = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200723/cfram/"
file2 = diri2+"forcing_SR.grd"
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
x@long_name = "forcing"
x@units     = "W/m2"

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


; ;===============================================================================================
; obs = new((/2,nlat,nlon/),float,UNDEF)
; obs(0,:,:) = fbindirread(file3,0,(/nlat,nlon/),"float")
; obs(1,:,:) = fbindirread(file4,0,(/nlat,nlon/),"float")
; copy_VarMeta(x(0,0,:,:),obs(0,:,:))
; printVarSummary(obs)
; print("end of data read")
; d_obs = obs(0,:,:) - obs(1,:,:)
; copy_VarMeta(obs(0,:,:), d_obs)

total_x1 = dim_sum_n(x(0:2,:,:,:),0)+dim_sum_n(x(5:7,:,:,:),0)+dim_sum_n(x(9:12,:,:,:),0)
total_x2 = dim_sum_n(x(0:2,:,:,:),0)+dim_sum_n(x(5:7,:,:,:),0)+dim_sum_n(x(11:13,:,:,:),0)+x(9,:,:,:)
copy_VarMeta(x(0,:,:,:), total_x1)
copy_VarMeta(x(0,:,:,:), total_x2)
; printVarSummary(total_x1)

; error1 = d_obs - total_x1(26,:,:)
; error2 = d_obs - total_x2(26,:,:)
; copy_VarMeta(total_x1(26,:,:), error1)
; copy_VarMeta(total_x2(26,:,:), error2)
; printVarSummary(error1)
;print(error({35:50},{100:120}))


x_inn = x(:,26,:,:)

;==========================================================
;****************************************************
; plot
;****************************************************

figpath = diri2+"post-pro/"
figname = "draw_process_spatial_forcing_SR"
wks  = gsn_open_wks("pdf",figpath+figname)
gsn_define_colormap(wks,"temp_diff_18lev")

nplot = 12
plot = new(nplot, graphic)
plot_topo  = new(nplot,graphic)
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
res_topo@cnLineThicknessF       = 2
res_topo@cnLineColor            = "grey32"
;res_topo@cnLineDashPattern      = "1"
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
res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1,1.5,2,2.5,3/)
res@cnFillColors = (/4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/)

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMaxLevelValF = 5
; res@cnMinLevelValF = -5
; res@cnLevelSpacingF = 0.5


res@gsnLeftString = "(a) Heat Storage"
plot(0) = gsn_csm_contour_map(wks,x_inn(10,:,:),res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0),plot_topo(0))

res@gsnLeftString = "(b) Surf.LH"
plot(1) = gsn_csm_contour_map(wks,x_inn(12,:,:),res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1),plot_topo(1))

res@gsnLeftString = "(c) Surf.SH"
plot(2) = gsn_csm_contour_map(wks,x_inn(11,:,:),res)
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2),plot_topo(2))


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


res@gsnLeftString = "(e) Cloud"
plot(4) = gsn_csm_contour_map(wks,x_inn(2,:,:),res)
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4),plot_topo(4))

res@gsnLeftString = "(f) Atm_Dyn"
plot(5) = gsn_csm_contour_map(wks,x_inn(9,:,:),res)
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5),plot_topo(5))

;***************************************************

;res@cnMinLevelValF = -1
;res@cnMaxLevelValF = 1
;res@cnLevelSpacingF = 0.05

res@gsnLeftString = "(g) WV"
plot(6) = gsn_csm_contour_map(wks,x_inn(1,:,:),res)
plot_topo(6) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(6),plot_topo(6))

res@gsnLeftString = "(h) SR_Q"
plot(7) = gsn_csm_contour_map(wks,x_inn(13,:,:),res)
plot_topo(7) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(7),plot_topo(7))

res@gsnLeftString = "(h) sum of (a)-(g)"
plot(8) = gsn_csm_contour_map(wks,total_x1(26,:,:),res)
plot_topo(8) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(8),plot_topo(8))

; res@gsnLeftString = "(j) error: OBS - (i)"
; plot(9) = gsn_csm_contour_map(wks,error1,res)
; plot_topo(9) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(9),plot_topo(9))

res@gsnLeftString = "(j) sum of (b)-(h)"
plot(9) = gsn_csm_contour_map(wks,total_x2(26,:,:),res)
plot_topo(9) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(9),plot_topo(9))

; res@gsnLeftString = "(l) error: OBS - (k)"
; plot(11) = gsn_csm_contour_map(wks,error2,res)
; plot_topo(11) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(11),plot_topo(11))

; res@gsnLeftString = "(j) SOLAR"
; plot(9) = gsn_csm_contour_map(wks,x_inn(7,:,:),res)
; plot_topo(9) = gsn_csm_contour(wks, topo, res_topo)
; overlay(plot(9),plot_topo(9))


;------------------------------------------------------
; the following use for diff color tables
;------------------------------------------------------
;-------------------------------------------------------

pres = True
pres@gsnMaximize = True
pres@gsnFrame    = True
pres@gsnDraw     = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "DJF"  ;(MERRA data)

pres@lbLabelFontHeightF = 0.015
pres@pmLabelBarWidthF = 0.8    ;plot manager 
pres@pmLabelBarHeightF = 0.06

;pplot = gsn_panel_return(wks,plot(0:2),(/1,3/),pres)
;pplot2 = gsn_panel_return(wks,plot(3:8),(/2,3/),pres)

;bb = NhlGetBB(pplot)
;height = max(bb(0,:)-min(bb(:,1)))
;pres@gsnDraw = True

;pres@gsnPanelTop = 0.9
;pres@gsnPanelBottom = pres@gsnPanelTop - height
;gsn_panel(wks,plot(0:2),(/1,3/),pres)
;pres@gsnPanelMainString = ""
;pres@gsnPanelTop = pres@gsnPanelBottom
;pres@gsnPanelBottom = pres@gsnPanelTop - height*2
;gsn_panel(wks,plot(3:8),(/2,3/),pres)


gsn_panel(wks,plot,(/3,4/),pres)


end


