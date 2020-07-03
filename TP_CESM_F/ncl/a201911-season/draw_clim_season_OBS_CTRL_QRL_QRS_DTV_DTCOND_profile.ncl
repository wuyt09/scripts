
begin

lonstr = 65
lonend = 110
latstr = 20
latend = 45
;lonstr = 80
;lonend = 100
;latstr = 27.5
;latend = 37.5
Season = (/"MAM","JJA","SON","DJF"/)
;===========================================================================
diri = "/home/yangsong3/wuyt/data/ERA-Interim/monthly/"
fili = "q1.mon.mean.1979.2016.ERAI.nc" ; ERAI 

f = addfile(diri+fili, "r")
Q1 = f->q1(:,:,{latstr:latend},{lonstr:lonend})    ; units K/s

;===========================================================================
diri  = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
fili1 = "F2000_CTRL.plevel.QRL.mon.0601-3012.nc"
fili2 = "F2000_CTRL.plevel.QRS.mon.0601-3012.nc"
fili3 = "F2000_CTRL.plevel.DTCOND.mon.0601-3012.nc"
fili4 = "F2000_CTRL.plevel.DTV.mon.0601-3012.nc"

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

;=============================================
q1 = linint2_Wrap(Q1&lon, Q1&lat, Q1, False, x1&lon, x1&lat, 0)
printVarSummary(q1)

Q1sea = month_to_seasonN(q1, Season)
q1sea = dim_avg_n_Wrap(Q1sea, 1)
printVarSummary(q1sea)
;=============================================
JJA = new((/4,dimsizes(Season),dimsizes(x1&lev),dimsizes(x1&lat),dimsizes(x1&lon)/), float)
JJA(0,:,:,:,:) = x11_JJA
JJA(1,:,:,:,:) = x22_JJA
JJA(2,:,:,:,:) = x33_JJA
JJA(3,:,:,:,:) = x44_JJA
printVarSummary(JJA)

JJA2 = new((/5,dimsizes(Season),dimsizes((x1&lev)),dimsizes(x1&lat),dimsizes(x1&lon)/), float)
JJA2(2:3,:,:,:,:) = JJA(2:3,:,:,:,:)
JJA2(1,:,:,:,:) = JJA(0,:,:,:,:)+JJA(1,:,:,:,:) ;radiation cooling
JJA2(0,:,:,:,:) = JJA2(1,:,:,:,:) + JJA2(2,:,:,:,:) + JJA2(3,:,:,:,:)
JJA2(4,:,:,:,:) = q1sea

;=======================================================
maskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
maskname = "TPmask_1000m_C20190410.nc"
fmask = addfile(maskpath+maskname, "r")
tp = fmask->tp({latstr:latend},{lonstr:lonend})
tp0 = conform_dims(dimsizes(JJA2), tp, (/3,4/))
printVarSummary(tp0)
;JJA2 = mask(JJA2, tp0.ge.0.9, True)
;print(JJA2(0,0,4,{35},{90}))
;print(JJA2(0,0,4,{20},{90}))

wgty = JJA2(0,0,0,{25:40},0)
wgtx = JJA2(0,0,0,0,{75:105})
;wgty = JJA2(0,0,0,{25:40},0)
;wgtx = JJA2(0,0,0,0,{75:105})
;wgty = JJA2(0,0,0,{30:35},0)
;wgtx = JJA2(0,0,0,0,{90:100})
wgtx = 1
wgty = 1
;jja = wgt_areaave_Wrap(JJA(:,:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 1)
jja2 = wgt_areaave_Wrap(JJA2(:,:,:,{25:40},{75:105}), wgty, wgtx, 1)
;jja2 = wgt_areaave_Wrap(JJA2(:,:,:,{25:40},{75:105}), wgty, wgtx, 1)
;jja2 = wgt_areaave_Wrap(JJA2(:,:,:,{30:35},{90:100}), 1, 1, 1)
printVarSummary(jja2)

;==============================================================
 figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
 figname = "draw_clim_season_OBS_CTRL_QRL_QRS_DTV_DTCOND_profile"
 wks = gsn_open_wks("pdf", figpath+figname)
 plot = new(4, graphic)
 dum1 = new(5, graphic)
 dum2 = dum1
 
 ;Colors = (/"green","red","green4","gold"/)
 ;Texts = (/"QRL","QRS","DTCOND","DTV"/)

 Colors = (/"black","goldenrod","red","green4","black"/)
 Texts = (/"CTRL_total","QRL+QRS","DTCOND","DTV","ERA_total"/)

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
 res@xyLineThicknessF = 5
 res@xyDashPatterns = (/0,0,0,0,2/)

 res@trYReverse = True
 
 res@gsnRightString = "K/day"
 res@gsnStringFontHeightF = 0.025

 do j = 0,dimsizes(Season)-1
 res@gsnLeftString = Season(j)
 plot(j) = gsn_csm_xy(wks, jja2(:,j,{1000:50})*86400, lev({1000:50}), res)  ;
 
 ;res@gsnLeftString = "CTRL_DJF"
 ;plot(1) = gsn_csm_xy(wks, jja2(:,i,{1000:50})*86400, lev({1000:50}), res)  ;
 end do

 res_lines = True
 res_lines@gsLineDashPattern = 0
 res_lines@gsLineThicknessF = 4

 res_text  = True
 res_text@txFontHeightF = 0.018
 res_text@txJust = "CenterLeft"
 xx = (/-5.5,-4.5/)
 text_pos = -4.3

 do i = 0,4
    yy = (/150+i*100,150+i*100/)
    res_lines@gsLineColor = Colors(i)
    if (i.eq.4) then
        res_lines@gsLineDashPattern = 2
    end if
    dum1(i) = gsn_add_polyline(wks, plot(0), xx, yy, res_lines)
    dum2(i) = gsn_add_text(wks, plot(0), Texts(i), text_pos, yy(0), res_text)
end do


 pres = True
 pres@gsnMaximize = True
 ;pres@gsnPanelMainString = "CTRL"
gsn_panel(wks, plot, (/2,2/), pres)
end
