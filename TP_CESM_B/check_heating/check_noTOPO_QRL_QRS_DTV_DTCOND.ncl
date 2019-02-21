
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40


diri  = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
fili1 = "TP_noTOPO.cam.h0.QRL.250-349.nc"
;fili1 = "TP_noTOPO.cam.h0.T.250-349.nc"
fili2 = "TP_noTOPO.cam.h0.QRS.250-349.nc"
fili3 = "TP_noTOPO.cam.h0.DTCOND.250-349.nc"
fili4 = "TP_noTOPO.cam.h0.DTV.250-349.nc"

filia = "TP_noTOPO.cam.h0.hyam.250-349.nc"
filib = "TP_noTOPO.cam.h0.hybm.250-349.nc"
filip = "TP_noTOPO.cam.h0.PS.250-349.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")
f4 = addfile(diri+fili4, "r")
fa = addfile(diri+filia, "r")
fb = addfile(diri+filib, "r")
fp = addfile(diri+filip, "r")

x01 = f1->T(0:1187,:,:,:)
x02 = f2->QRS(0:1187,:,:,:)
x03 = f3->DTCOND(0:1187,:,:,:)
x04 = f4->DTV(0:1187,:,:,:)
hyam = fa->hyam(0:1187)
hybm = fb->hybm(0:1187)
ps   = fp->PS(0:1187,:,:)
printVarSummary(x01)

plevo = (/1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700,\
         650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175,\
         150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1/)
;plevo = (/1000,850,750,500,300,200,100/)
plevo@units = "hPa"
p0  = 1000.
x1 = vinth2p(x01, hyam, hybm, plevo, ps, 1, p0, 1, False)
x2 = vinth2p(x02, hyam, hybm, plevo, ps, 1, p0, 1, False)
x3 = vinth2p(x03, hyam, hybm, plevo, ps, 1, p0, 1, False)
x4 = vinth2p(x04, hyam, hybm, plevo, ps, 1, p0, 1, False)
printVarSummary(x1)
x1!1 = "lev"
x1&lev = plevo
copy_VarMeta(x1, x2)
copy_VarMeta(x1, x3)
copy_VarMeta(x1, x4)
printVarSummary(x1)
print(x1&lev)


x1_JJA = month_to_season(x1, "JJA")
x2_JJA = month_to_season(x2, "JJA")
x3_JJA = month_to_season(x3, "JJA")
x4_JJA = month_to_season(x4, "JJA")

x1_DJF = month_to_season(x1, "DJF")
x2_DJF = month_to_season(x2, "DJF")
x3_DJF = month_to_season(x3, "DJF")
x4_DJF = month_to_season(x4, "DJF")
printVarSummary(x1_JJA)

x11_JJA = dim_avg_n_Wrap(x1_JJA, 0)
x22_JJA = dim_avg_n_Wrap(x2_JJA, 0)
x33_JJA = dim_avg_n_Wrap(x3_JJA, 0)
x44_JJA = dim_avg_n_Wrap(x4_JJA, 0)

x11_DJF = dim_avg_n_Wrap(x1_DJF, 0)
x22_DJF = dim_avg_n_Wrap(x2_DJF, 0)
x33_DJF = dim_avg_n_Wrap(x3_DJF, 0)
x44_DJF = dim_avg_n_Wrap(x4_DJF, 0)
printVarSummary(x11_JJA)

JJA = new((/4,dimsizes(plevo),dimsizes(x01&lat),dimsizes(x01&lon)/), float)
DJF = new((/4,dimsizes(plevo),dimsizes(x01&lat),dimsizes(x01&lon)/), float)
JJA(0,:,:,:) = x11_JJA
JJA(1,:,:,:) = x22_JJA
JJA(2,:,:,:) = x33_JJA
JJA(3,:,:,:) = x44_JJA
DJF(0,:,:,:) = x11_DJF
DJF(1,:,:,:) = x22_DJF
DJF(2,:,:,:) = x33_DJF
DJF(3,:,:,:) = x44_DJF
printVarSummary(JJA)


;jja = dim_avg_n_Wrap(JJA(:,:,{latstr:latend},{lonstr:lonend}), (/2,3/))
;djf = dim_avg_n_Wrap(DJF(:,:,{latstr:latend},{lonstr:lonend}), (/2,3/))
;printVarSummary(jja)

;jja = dim_avg_n_Wrap(JJA(:,:,{latstr:latend},:), 2)
;djf = dim_avg_n_Wrap(DJF(:,:,{latstr:latend},:), 2)
;printVarSummary(jja)

jja = JJA(:,:,{35},{90})
djf = DJF(:,:,{35},{90})

;==============================================================
 figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/check_heating/"
  figname = "check_noTOPO_QRL_QRS_DTV_DTCOND"
 wks = gsn_open_wks("pdf", "check_noTOPO_QRL_QRS_DTV_DTCOND")
 plot = new(2, graphic)
 dum1 = new(4, graphic)
 dum2 = dum1
 
 Colors = (/"blue","red","green4","gold"/)
 Texts = (/"QRL","QRS","DTCOND","DTV"/)

 res = True
 res@gsnDraw  = False
 res@gsnFrame = False

 res@gsnXRefLine = 0
 res@gsnXRefLineDashPattern = 2
 res@gsnXRefLineColor = "grey"
 res@gsnXRefLineThicknessF = 3
 res@trXMinF = -6.
 res@trXMaxF = 6
 res@trYMinF   = 10
 res@trYMaxF   = 1000
 
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)
 res@xyLineColors = Colors
 res@xyLineThicknessF = 3
 res@xyDashPattern = 0

 res@trYReverse = True
 res@trYLog = True
 
 res@gsnRightString = "K"
 res@gsnLeftString = "noTOPO_JJA"
 plot(0) = gsn_csm_xy(wks, jja(:,{1000:50})*86400, plevo(0:28), res)  ;
 res@gsnLeftString = "noTOPO_DJF"
 plot(1) = gsn_csm_xy(wks, djf(:,{1000:50})*86400, plevo(0:28), res)  ;

 res_lines = True
 res_lines@gsLineDashPattern = 0
 res_lines@gsLineThicknessF = 3

 res_text  = True
 res_text@txFontHeightF = 0.018
 res_text@txJust = "CenterLeft"
 xx = (/-5.5,-4.5/)
 text_pos = -4.3

 do i = 0,3
    yy = (/300+i*50,300+i*50/)
    res_lines@gsLineColor = Colors(i)
    dum1(i) = gsn_add_polyline(wks, plot(0), xx, yy, res_lines)
    dum2(i) = gsn_add_text(wks, plot(0), Texts(i), text_pos, yy(0), res_text)
end do

;res = True
;res@gsnDraw  = False
;res@gsnFrame = False

;res@cnFillOn = True
;res@cnLineLabelsOn = False
;res@cnInfoLabelOn  = False
;res@cnFillPalette  = "precip2_17lev"
;res@lbLabelAutoStride = True

;res@tmYLMode      = "Explicit"
;res@tmYLValues    = (/"1000","850","750","500","300","100"/)
;res@tmYLLabels    = (/"1000","850","750","500","300","100"/)

;res@gsnRightString = "lat :25-40N"
;res@cnLevelSelectionMode = "ManualLevels"
;res@cnMaxLevelValF  = 5e-5
;res@cnMinLevelValF  = 5e-6
;res@cnLevelSpacingF = 5e-6

;res@gsnLeftString = "CTRL_JJA"
;plot(0) = gsn_csm_pres_hgt(wks, jja(2,{1000:100},:), res)
;res@gsnLeftString = "CTRL_DJF"
;plot(1) = gsn_csm_pres_hgt(wks, djf(2,{1000:100},:), res)

pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/1,2/), pres)
end
