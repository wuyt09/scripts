
begin

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
;diri2 = "/home/yangsong3/data-model/wuyt/TP_NOTP_F2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

filiu1 = "F2000_CTRL.plevel.DTCOND.mon.0601-3012.nc"
;filiu2 = "F2000_NOTP.plevel.DTCOND.mon.0601-3012.nc"
filiu2 = "F2000_NOQ1.plevel.DTCOND.mon.0601-3012.nc"
filiu3 = "F2000_NOTP_NOQ1.plevel.DTCOND.mon.0601-3012.nc"

filiv1 = "F2000_CTRL.plevel.DTV.mon.0601-3012.nc"
;filiv2 = "F2000_NOTP.plevel.DTV.mon.0601-3012.nc"
filiv2 = "F2000_NOQ1.plevel.DTV.mon.0601-3012.nc"
filiv3 = "F2000_NOTP_NOQ1.plevel.DTV.mon.0601-3012.nc"

filis1 = "F2000_CTRL.plevel.QRS.mon.0601-3012.nc"
;filis2 = "F2000_NOTP.plevel.QRS.mon.0601-3012.nc"
filis2 = "F2000_NOQ1.plevel.QRS.mon.0601-3012.nc"
filis3 = "F2000_NOTP_NOQ1.plevel.QRS.mon.0601-3012.nc"

filil1 = "F2000_CTRL.plevel.QRL.mon.0601-3012.nc"
;filil2 = "F2000_NOTP.plevel.QRL.mon.0601-3012.nc"
filil2 = "F2000_NOQ1.plevel.QRL.mon.0601-3012.nc"
filil3 = "F2000_NOTP_NOQ1.plevel.QRL.mon.0601-3012.nc"


fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")

fs1 = addfile(diri1+filis1, "r")
fs2 = addfile(diri2+filis2, "r")
fs3 = addfile(diri3+filis3, "r")

fl1 = addfile(diri1+filil1, "r")
fl2 = addfile(diri2+filil2, "r")
fl3 = addfile(diri3+filil3, "r")

u1 = fu1->DTCOND
u2 = fu2->DTCOND
u3 = fu3->DTCOND

v1 = fv1->DTV
v2 = fv2->DTV
v3 = fv3->DTV

s1 = fs1->QRS
s2 = fs1->QRS
s3 = fs3->QRS

l1 = fl1->QRL
l2 = fl2->QRL
l3 = fl3->QRL
printVarSummary(u1)

Cpair = 1004.07
G = 9.8 

q1 = u1
q1 = (u1+v1+s1+l1)*Cpair/G

q2 = u2
q2 = (u2+v2+s2+l2)*Cpair/G

q3 = u3
q3 = (u3+v3+s3+l3)*Cpair/G

plev = u1&lev*100
plev!0 = "lev"
plev&lev = u1&lev
plev@units = "Pa"
psfc = max(plev)
ptop = min(plev)
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dp = conform_dims(dimsizes(q1), dP, 1)

q1dp = q1
q1dp = q1*dp
q2dp = q2
q2dp = q2*dp
q3dp = q3
q3dp = q3*dp

q11 = dim_sum_n_Wrap(q1dp, 1);/sum(dP)
q22 = dim_sum_n_Wrap(q2dp, 1);/sum(dP)
q33 = dim_sum_n_Wrap(q3dp, 1);/sum(dP)
copy_VarCoords(q1(:,0,:,:), q11)
copy_VarCoords(q2(:,0,:,:), q22)
copy_VarCoords(q3(:,0,:,:), q33)

Season = (/"MAM","JJA","SON","DJF"/)

uu11 = month_to_seasonN(q11, Season)
printVarSummary(uu11)
uu1 = dim_avg_n_Wrap(uu11, 1)
uu22 = month_to_seasonN(q22, Season)
uu2 = dim_avg_n_Wrap(uu22, 1)
uu33 = month_to_seasonN(q33, Season)
uu3 = dim_avg_n_Wrap(uu33, 1)
printVarSummary(uu1)

u_m = uu1 - uu3
u_t = uu1 - uu2
u_h = u_t - u_m
copy_VarCoords(uu1, u_m)
copy_VarCoords(uu1, u_t)
copy_VarCoords(uu1, u_h)

;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_season_Q1"

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
res@mpGeophysicalLineColor = "grey45" ;"grey25"
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
res@cnLevelSpacingF = 30.
res@cnMaxLevelValF = 300.
res@cnMinLevelValF = -300.
res@cnFillPalette = "BlueWhiteOrangeRed"

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "CTRL"
plot(i) = gsn_csm_contour_map(wks, uu1(i,:,:), res)

res@gsnLeftString = "NOQ1" ;NOTP
plot(i+4) = gsn_csm_contour_map(wks, uu2(i,:,:), res)

;res@gsnLeftString = "NOQ1"
res@gsnLeftString = "NOTP_NOQ1"
plot(i+8) = gsn_csm_contour_map(wks, uu3(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1_climatology"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)


;------------------------------------------
vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
;vcres@vcRefLengthF            = 0.035            ; define length of vec ref
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 1.
res@cnMinLevelValF = -1.
res@cnFillPalette = "BlueWhiteOrangeRed"

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "CTRL - NOTP"
plot(i) = gsn_csm_contour_map(wks, u_t(i,:,:), res)

res@gsnLeftString = "CTRL - NOQ1"
plot(i+4) = gsn_csm_contour_map(wks, u_m(i,:,:), res)

res@gsnLeftString = "NOQ1 - NOTP"
plot(i+8) = gsn_csm_contour_map(wks, u_h(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1_difference"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)

end