
begin

yrStrt = 260
yrLast = 304

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"

filiu1 = "B2000_f19g16_CP_CTRL.plevel.DTCOND.mon.260-304.nc"
filiu2 = "TP_NOTOPO.plevel.DTCOND.mon.260-304.nc"

filiv1 = "B2000_f19g16_CP_CTRL.plevel.DTV.mon.260-304.nc"
filiv2 = "TP_NOTOPO.plevel.DTV.mon.260-304.nc"

fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")

u1 = fu1->DTCOND(:,:,{35},{90})
u2 = fu2->DTCOND(:,:,{35},{90})
v1 = fv1->DTV(:,:,{35},{90})
v2 = fv2->DTV(:,:,{35},{90})

lev = u1&lev
lev!0 = "level"
lev&level = lev

x1 = u1
x1 = u1+v1
x1 = x1*86400
x1@long_name = "Sum of heating"
x2 = u2
x2 = u2+v2
x2 = x2*86400
x2@long_name = "Sum of heating"
printVarSummary(x1)

;diff = x2(0,:) - x1(0,:)
;copy_VarCoords(x2(0,:), diff)
;printVarSummary(diff)
;print(diff)

diff = x2 - x1
copy_VarCoords(x2, diff)
printVarSummary(diff)


;x = x&lev

pp = new((/3,12,dimsizes(lev)/), "float")
pp(0,:,:) = x1(0:11,:)
pp(1,:,:) = x2(0:11,:)
pp(2,:,:) = diff(0:11,:)
copy_VarCoords(x1(0:11,:), pp(0,:,:))
printVarSummary(pp)
;=====================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "cala_fitting_CTRL_heating_to_NOTP_heating"

 Colors = (/"green4","red","goldenrod"/)
 Texts = (/"CTRL","NOTP","Diff"/)

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 12
plot = new(nplot, graphic)
dum1 = new(3, graphic)
dum2 = dum1

 res = True
 res@gsnDraw  = False
 res@gsnFrame = False
 ;res@gsnMaximize = True

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
 
 res@gsnRightString = "K/day"
 res@gsnStringFontHeightF = 0.025
 
 do j = 0, 11
    res@gsnLeftString = "Month = "+(j+1)
    plot(j) = gsn_csm_xy(wks, pp(:,j,{1000:100}), lev({1000:100}), res)
 end do

 res_lines = True
 res_lines@gsLineDashPattern = 0
 res_lines@gsLineThicknessF = 3

 res_text  = True
 res_text@txFontHeightF = 0.018
 res_text@txJust = "CenterLeft"
 xx = (/-5.5,-4.5/)
 text_pos = -4.3

 do i = 0,2
    yy = (/300+i*100,300+i*100/)
    res_lines@gsLineColor = Colors(i)
    dum1(i) = gsn_add_polyline(wks, plot(0), xx, yy, res_lines)
    dum2(i) = gsn_add_text(wks, plot(0), Texts(i), text_pos, yy(0), res_text)
end do

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/3,4/), pres)
end