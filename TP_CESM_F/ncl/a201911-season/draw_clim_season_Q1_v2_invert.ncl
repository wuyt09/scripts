
begin

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/daily/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOTP_F2000_F19G16_CAM4/pro/daily/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/daily/"

fili1 = systemfunc("ls "+diri1+"q1.2*.nc")
fili2 = systemfunc("ls "+diri2+"q1.2*.nc")
fili3 = systemfunc("ls "+diri3+"q1.2*.nc")
print(fili1)

f1 = addfiles(fili1, "r")
f2 = addfiles(fili2, "r")
f3 = addfiles(fili3, "r")

Q11 = f1[:]->q1
printVarSummary(Q11)
nlat = dimsizes(Q11&lat)
nlon = dimsizes(Q11&lon)
nlev = dimsizes(Q11&lev)
Q1 = reshape(Q11, (/4,365,nlev,nlat,nlon/))
q1 = dim_avg_n_Wrap(Q1, 0)
copy_VarCoords(Q11(0:364,:,:,:), q1)
printVarSummary(q1)

Q22 = f2[:]->q1
Q33 = f3[:]->q1
Q2 = reshape(Q22, (/4,365,nlev,nlat,nlon/))
q2 = dim_avg_n_Wrap(Q2, 0)
copy_VarCoords(Q22(0:364,:,:,:), q2)
Q3 = reshape(Q33, (/4,365,nlev,nlat,nlon/))
q3 = dim_avg_n_Wrap(Q3, 0)
copy_VarCoords(Q33(0:364,:,:,:), q3)

;q1 = q1*1004.07  ;Cp J/Kg/K
q1 = q1*86400
q2 = q2*86400
q3 = q3*86400

plev = Q11&lev*100
plev!0 = "lev"
plev&lev = Q11&lev
psfc = max(plev)
ptop = min(plev)
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dp = conform_dims(dimsizes(q1), dP, 1)
q1dp = q1*dp
q2dp = q2*dp
q3dp = q3*dp
copy_VarCoords(q1, q1dp)

q11 = dim_sum_n_Wrap(q1dp, 1)/sum(dP)
copy_VarCoords(q1(:,0,:,:), q11)
;q11 = q11/9.8
printVarSummary(q11)
q22 = dim_sum_n_Wrap(q2dp, 1)/sum(dP)
copy_VarCoords(q2(:,0,:,:), q22)
q33 = dim_sum_n_Wrap(q3dp, 1)/sum(dP)
copy_VarCoords(q3(:,0,:,:), q33)

calday = (/31,28,31,30,31,30,31,31,30,31,30,31/)
Q1mon = q11(0:11,:,:)
Q2mon = q22(0:11,:,:)
Q3mon = q33(0:11,:,:)
days = 0
do i = 0, dimsizes(calday)-1
  daye = days+calday(i)-1
  Q1mon(i,:,:) = dim_avg_n_Wrap(q11(days:daye,:,:), 0)
  Q2mon(i,:,:) = dim_avg_n_Wrap(q22(days:daye,:,:), 0)
  Q3mon(i,:,:) = dim_avg_n_Wrap(q33(days:daye,:,:), 0)
  days = daye+1
  print(days)
end do
printVarSummary(Q1mon)

Season = (/"DJF","MAM","JJA","SON"/)
Q1sea = month_to_seasonN(Q1mon, Season)
q1sea = dim_avg_n_Wrap(Q1sea, 1)
printVarSummary(q1sea)
Q2sea = month_to_seasonN(Q2mon, Season)
q2sea = dim_avg_n_Wrap(Q2sea, 1)
Q3sea = month_to_seasonN(Q3mon, Season)
q3sea = dim_avg_n_Wrap(Q3sea, 1)


;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_season_Q1_v2_invert"

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

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.
res@cnFillPalette = "BlueWhiteOrangeRed"

do i = 0, dimsizes(Season)-1
res@gsnRightString = Season(i)

res@gsnLeftString = "CTRL"
plot(i) = gsn_csm_contour_map(wks, q1sea(i,:,:), res)

res@gsnLeftString = "NOTP"
plot(i+4) = gsn_csm_contour_map(wks, q2sea(i,:,:), res)

res@gsnLeftString = "NOQ1"
plot(i+8) = gsn_csm_contour_map(wks, q3sea(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "Q1_climatology"
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,4/), pres)


end