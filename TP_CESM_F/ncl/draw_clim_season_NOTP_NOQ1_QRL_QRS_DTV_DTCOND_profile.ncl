
begin

lonstr = 65
lonend = 110
latstr = 20
latend = 45
;lonstr = 80
;lonend = 100
;latstr = 27.5
;latend = 37.5


diri  = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"
fili1 = "F2000_NOTP_NOQ1.plevel.QRL.mon.0601-3012.nc"
fili2 = "F2000_NOTP_NOQ1.plevel.QRS.mon.0601-3012.nc"
fili3 = "F2000_NOTP_NOQ1.plevel.DTCOND.mon.0601-3012.nc"
fili4 = "F2000_NOTP_NOQ1.plevel.DTV.mon.0601-3012.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")
f4 = addfile(diri+fili4, "r")

x1 = f1->QRL(:,:,{latstr:latend},{lonstr:lonend})
x2 = f2->QRS(:,:,{latstr:latend},{lonstr:lonend})
x3 = f3->DTCOND(:,:,{latstr:latend},{lonstr:lonend})
x4 = f4->DTV(:,:,{latstr:latend},{lonstr:lonend})
printVarSummary(x1)
lev = x1&lev
lev!0 = "level"
lev&level = lev

Season = (/"DJF","MAM","JJA","SON"/)
x1_JJA = month_to_seasonN(x1, Season)
x2_JJA = month_to_seasonN(x2, Season)
x3_JJA = month_to_seasonN(x3, Season)
x4_JJA = month_to_seasonN(x4, Season)
printVarSummary(x1_JJA)

x11_JJA = dim_avg_n_Wrap(x1_JJA, 1)
x22_JJA = dim_avg_n_Wrap(x2_JJA, 1)
x33_JJA = dim_avg_n_Wrap(x3_JJA, 1)
x44_JJA = dim_avg_n_Wrap(x4_JJA, 1)

printVarSummary(x11_JJA)

JJA = new((/4,dimsizes(Season),dimsizes(x1&lev),dimsizes(x1&lat),dimsizes(x1&lon)/), float)
JJA(0,:,:,:,:) = x11_JJA
JJA(1,:,:,:,:) = x22_JJA
JJA(2,:,:,:,:) = x33_JJA
JJA(3,:,:,:,:) = x44_JJA
printVarSummary(JJA)

JJA2 = new((/4,dimsizes(Season),dimsizes((x1&lev)),dimsizes(x1&lat),dimsizes(x1&lon)/), float)
JJA2(2:3,:,:,:,:) = JJA(2:3,:,:,:,:)
JJA2(1,:,:,:,:) = JJA(0,:,:,:,:)+JJA(1,:,:,:,:) ;radiation cooling
JJA2(0,:,:,:,:) = JJA2(1,:,:,:,:) + JJA2(2,:,:,:,:) + JJA2(3,:,:,:,:)

;jja = wgt_areaave_Wrap(JJA(:,:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)
;jja2 = wgt_areaave_Wrap(JJA2(:,:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)

;=======================================================
maskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
maskname = "TPmask_1000m_C20190410.nc"
fmask = addfile(maskpath+maskname, "r")
tp = fmask->tp({latstr:latend},{lonstr:lonend})
tp0 = conform_dims(dimsizes(JJA2), tp, (/3,4/))
printVarSummary(tp0)
JJA2 = mask(JJA2, tp0.ge.0.9, True)
;print(JJA2(0,0,4,{35},{90}))
;print(JJA2(0,0,4,{20},{90}))

jja = wgt_areaave_Wrap(JJA(:,:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)
jja2 = wgt_areaave_Wrap(JJA2(:,:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)
printVarSummary(jja)

;==============================================================
 figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
 figname = "draw_clim_season_NOTP_NOQ1_QRL_QRS_DTV_DTCOND_profile"
 wks = gsn_open_wks("pdf", figpath+figname)
 plot = new(4, graphic)
 dum1 = new(4, graphic)
 dum2 = dum1
 
 ;Colors = (/"green","red","green4","gold"/)
 ;Texts = (/"QRL","QRS","DTCOND","DTV"/)

 Colors = (/"black","goldenrod","red","green4"/)
 Texts = (/"total","QRL+QRS","DTCOND","DTV"/)

 res = True
 res@gsnDraw  = False
 res@gsnFrame = False

 res@gsnXRefLine = 0
 res@gsnXRefLineDashPattern = 2
 res@gsnXRefLineColor = "grey"
 res@gsnXRefLineThicknessF = 3
 res@trXMinF = -6.
 res@trXMaxF = 6
 res@trYMinF   = 100
 res@trYMaxF   = 1000
 res@trYLog = True
 
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)
 res@xyLineColors = Colors
 res@xyLineThicknessF = 3
 res@xyDashPattern = 0

 res@trYReverse = True
 
 res@gsnRightString = "K"
 res@gsnStringFontHeightF = 0.025

 do j = 0,dimsizes(Season)-1
 res@gsnLeftString = Season(j)
 plot(j) = gsn_csm_xy(wks, jja2(:,j,{1000:50})*86400, lev({1000:50}), res)  ;
 
 ;res@gsnLeftString = "CTRL_DJF"
 ;plot(1) = gsn_csm_xy(wks, jja2(:,i,{1000:50})*86400, lev({1000:50}), res)  ;
 end do

 res_lines = True
 res_lines@gsLineDashPattern = 0
 res_lines@gsLineThicknessF = 3

 res_text  = True
 res_text@txFontHeightF = 0.018
 res_text@txJust = "CenterLeft"
 xx = (/-5.5,-4.5/)
 text_pos = -4.3

 do i = 0,3
    yy = (/200+i*100,200+i*100/)
    res_lines@gsLineColor = Colors(i)
    dum1(i) = gsn_add_polyline(wks, plot(0), xx, yy, res_lines)
    dum2(i) = gsn_add_text(wks, plot(0), Texts(i), text_pos, yy(0), res_text)
end do


 pres = True
 pres@gsnMaximize = True
 pres@gsnPanelMainString = "NOTP_NOQ1"
gsn_panel(wks, plot, (/2,2/), pres)
end
