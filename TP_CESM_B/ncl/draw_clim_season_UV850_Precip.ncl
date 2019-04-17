
begin

yrStat = 260
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOSH_B2000_F19G16_CAM4/pro/"

filiu1 = "B2000_f19g16_CP_CTRL.plevel.U.mon.260-304.nc"
filiu2 = "TP_NOTOPO.plevel.U.mon.260-304.nc"
filiu3 = "TP_NOSH.plevel.U.mon.260-304.nc"

filiv1 = "B2000_f19g16_CP_CTRL.plevel.V.mon.260-304.nc"
filiv2 = "TP_NOTOPO.plevel.V.mon.260-304.nc"
filiv3 = "TP_NOSH.plevel.V.mon.260-304.nc"

filic1 = "B2000_f19g16_CP_CTRL.cam.h0.PRECC.260-304.nc"
filic2 = "TP_NOTOPO.PRECC.mon.26001-30508.nc"
filic3 = "TP_NOSH.PRECC.mon.26001-30712.nc"

filil1 = "B2000_f19g16_CP_CTRL.cam.h0.PRECL.260-304.nc"
filil2 = "TP_NOTOPO.PRECL.mon.26001-30508.nc"
filil3 = "TP_NOSH.PRECL.mon.26001-30712.nc"

fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")

u1 = fu1->U
u2 = fu2->U
u3 = fu3->U
v1 = fv1->V
v2 = fv2->V
v3 = fv3->V
printVarSummary(u1)

;================================================================
fc1 = addfile(diri1+filic1, "r")
fc2 = addfile(diri2+filic2, "r")
fc3 = addfile(diri3+filic3, "r")
fl1 = addfile(diri1+filil1, "r")
fl2 = addfile(diri2+filil2, "r")
fl3 = addfile(diri3+filil3, "r")

time1 = fc1->time
time2 = fc2->time
time3 = fc3->time
TIME1 = cd_calendar(time1, -2)/100
TIME2 = cd_calendar(time2, -2)/100
TIME3 = cd_calendar(time3, -2)/100
;print(TIME1(0:40))
;print(TIME2)

iYYYY1 = ind(TIME1.ge.26002 .and. TIME1.le.30501)
iYYYY2 = ind(TIME2.ge.26002 .and. TIME2.le.30501)
iYYYY3 = ind(TIME3.ge.26002 .and. TIME3.le.30501)

c1 = fc1->PRECC(iYYYY1,:,:)
c2 = fc2->PRECC(iYYYY2,:,:)
c3 = fc3->PRECC(iYYYY3,:,:)
l1 = fl1->PRECL(iYYYY1,:,:)
l2 = fl2->PRECL(iYYYY2,:,:)
l3 = fl3->PRECL(iYYYY3,:,:)
printVarSummary(c1)
;==================================================================

p1 = c1+l1
p2 = c2+l2
p3 = c3+l3
p1 = p1*1000*86400. ; convert units from m/s to mm/day
p2 = p2*1000*86400.
p3 = p3*1000*86400.
copy_VarCoords(c1, p1)
copy_VarCoords(c2, p2)
copy_VarCoords(c3, p3)
printVarSummary(p1)

Season = (/"DJF","MAM","JJA","SON"/)
pp11 = month_to_seasonN(p1, Season)
pp1 = dim_avg_n_Wrap(pp11, 1)
pp22 = month_to_seasonN(p2, Season)
pp2 = dim_avg_n_Wrap(pp22, 1)
pp33 = month_to_seasonN(p3, Season)
pp3 = dim_avg_n_Wrap(pp33, 1)

uu11 = month_to_seasonN(u1, Season)
uu1 = dim_avg_n_Wrap(uu11, 1)
uu22 = month_to_seasonN(u2, Season)
uu2 = dim_avg_n_Wrap(uu22, 1)
uu33 = month_to_seasonN(u3, Season)
uu3 = dim_avg_n_Wrap(uu33, 1)

vv11 = month_to_seasonN(v1, Season)
vv1 = dim_avg_n_Wrap(vv11, 1)
vv22 = month_to_seasonN(v2, Season)
vv2 = dim_avg_n_Wrap(vv22, 1)
vv33 = month_to_seasonN(v3, Season)
vv3 = dim_avg_n_Wrap(vv33, 1)
printVarSummary(uu1)

p_m = pp1 - pp3
p_t = pp1 - pp2
p_h = p_t - p_m
copy_VarCoords(pp1, p_m)
copy_VarCoords(pp1, p_t)
copy_VarCoords(pp1, p_h)

u_m = uu1 - uu3
u_t = uu1 - uu2
u_h = u_t - u_m
copy_VarCoords(uu1, u_m)
copy_VarCoords(uu1, u_t)
copy_VarCoords(uu1, u_h)

v_m = vv1 - vv3
v_t = vv1 - vv2
v_h = v_t - v_m
copy_VarCoords(vv1, v_m)
copy_VarCoords(vv1, v_t)
copy_VarCoords(vv1, v_h)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_season_UV850_Precip"

lonstr = 40
lonend = 180
latstr = -20
latend = 60

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
plot = new(nplot, graphic)
plotv = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "red" ;"grey25"
res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up

  ;vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.025            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------
vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 2.
res@cnMaxLevelValF = 15.
res@cnMinLevelValF = 1.
res@cnFillPalette = "CBR_wet"

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "CTRL"
plot(i) = gsn_csm_contour_map(wks, pp1(i,:,:), res)
plotv(i) = gsn_csm_vector(wks, uu1(i,{850},:,:), vv1(i,{850},:,:), vcres)
overlay(plot(i), plotv(i))

res@gsnLeftString = "NOTP"
plot(i+4) = gsn_csm_contour_map(wks, pp2(i,:,:), res)
plotv(i+4) = gsn_csm_vector(wks, uu2(i,{850},:,:), vv2(i,{850},:,:), vcres)
overlay(plot(i+4), plotv(i+4))

res@gsnLeftString = "NOSH"
plot(i+8) = gsn_csm_contour_map(wks, pp3(i,:,:), res)
plotv(i+8) = gsn_csm_vector(wks, uu3(i,{850},:,:), vv3(i,{850},:,:), vcres)
overlay(plot(i+8), plotv(i+8))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV850_Precip "
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)


;------------------------------------------
vcres@vcRefMagnitudeF         = 1.0  ; define vector ref mag
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 4.
res@cnMinLevelValF = -4.
res@cnFillPalette = "CBR_drywet"

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "CTRL - NOTP"
plot(i) = gsn_csm_contour_map(wks, p_t(i,:,:), res)
;plot(1) = gsn_csm_vector_scalar_map(wks, u_t({850},:,:), v_t({850},:,:), p_t, res)
plotv(i) = gsn_csm_vector(wks, u_t(i,{850},:,:), v_t(i,{850},:,:), vcres)
overlay(plot(i), plotv(i))

res@gsnLeftString = "CTRL - NOSH"
plot(i+4) = gsn_csm_contour_map(wks, p_m(i,:,:), res)
;plot(3) = gsn_csm_vector_scalar_map(wks, u_m({850},:,:), v_m({850},:,:), p_m, res)
plotv(i+4) = gsn_csm_vector(wks, u_m(i,{850},:,:), v_m(i,{850},:,:), vcres)
overlay(plot(i+4), plotv(i+4))

res@gsnLeftString = "NOSH - NOTP"
plot(i+8) = gsn_csm_contour_map(wks, p_h(i,:,:), res)
;plot(5) = gsn_csm_vector_scalar_map(wks, u_h({850},:,:), v_h({850},:,:), p_h, res)
plotv(i+8) = gsn_csm_vector(wks, u_h(i,{850},:,:), v_h(i,{850},:,:), vcres)
overlay(plot(i+8), plotv(i+8))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "UV850_Precip "
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)

end