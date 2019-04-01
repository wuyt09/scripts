
begin

yrStat = 260
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_B2000_F19G16_CAM4/pro/"

filiu1 = "B2000_f19g16_CP_CTRL.plevel.DTCOND.mon.260-304.nc"
filiu2 = "TP_NOTOPO.plevel.DTCOND.mon.260-304.nc"
filiu3 = "TP_HEAT.plevel.DTCOND.mon.260-304.nc"

filiv1 = "B2000_f19g16_CP_CTRL.plevel.DTV.mon.260-304.nc"
filiv2 = "TP_NOTOPO.plevel.DTV.mon.260-304.nc"
filiv3 = "TP_HEAT.plevel.DTV.mon.260-304.nc"


fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")

fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")

u1 = fu1->DTCOND
u2 = fu2->DTCOND
u3 = fu3->DTCOND
v1 = fv1->DTV
v2 = fv2->DTV
v3 = fv3->DTV
printVarSummary(u1)

x1 = u1
x1 = u1+v1
x1 = x1*86400
x1@long_name = "Sum of heating"
x2 = u2
x2 = u2+v2
x2 = x2*86400
x2@long_name = "Sum of heating"
x3 = u3
x3 = u3+v3
x3 = x3*86400
x3@long_name = "Sum of heating"

;print(x1(0,{750},{35},{90}))
;print(x2(0,{750},{35},{90}))
;print(x3(0,{750},{35},{90}))

uu1 = clmMonTLLL(x1)
uu2 = clmMonTLLL(x2)
uu3 = clmMonTLLL(x3)
;print(uu1(0,{750},{35},{90}))
;print(uu2(0,{750},{35},{90}))
;print(uu3(0,{750},{35},{90}))
u_13 = uu1 - uu3
u_23 = uu2 - uu3
copy_VarCoords(uu1, u_13)
copy_VarCoords(uu2, u_23)
;print(u_13(0,:,{35},{90}))
;print(u_23(0,:,{35},{90}))
;exit
;=======================================================

maskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_B/model_need_input_data/"
maskname = "TPmask_C20190304.nc"
fmask = addfile(maskpath+maskname, "r")
tp = fmask->tp
tp0 = conform_dims(dimsizes(uu1), tp, (/2,3/))
printVarSummary(tp0)
u_13 = mask(u_13, tp0.ge.0.9, True)
u_23 = mask(u_23, tp0.ge.0.9, True)
print(u_23(0,4,{35},{90}))
print(u_23(0,4,{20},{90}))

lev = u1&lev
lev!0 = "level"
lev&level = lev

lat = u1&lat
lon = u1&lon
p = new((/2,12,dimsizes(lev),dimsizes(lat),dimsizes(lon)/), float)
p(0,:,:,:,:) = u_13
p(1,:,:,:,:) = u_23

pp = wgt_areaave_Wrap(p(:,:,:,{20:45},{60:110}), 1, 1, 1)
;pp = p(:,:,:,{35},{90})
printVarSummary(pp)

Texts = (/"CTRL - HEAT", "NOTP - HEAT"/)
Colors = (/"red","green4"/)
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/ncl/"
figname = "draw_clim_mon_diff_CTRL_NOTP_HEAT_profile"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 12
plot = new(nplot, graphic)
dum5 = new(2,graphic)
dum6 = dum5

 res = True
 res@gsnDraw  = False
 res@gsnFrame = False

 res@gsnXRefLine = 0
 res@gsnXRefLineDashPattern = 2
 res@gsnXRefLineColor = "grey"
 res@gsnXRefLineThicknessF = 3
 res@trXMinF = -2.
 res@trXMaxF = 2.
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
 res@gsnStringFontHeightF = 0.035

 do j = 0,11
    res@gsnLeftString = "Month = "+(j+1)
    plot(j) = gsn_csm_xy(wks, pp(:,j,{1000:50}), lev({1000:50}), res)  ;
 end do

    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineDashPattern = 0                     ; solid line
    res_lines@gsLineThicknessF  = 3                     ; line thicker

    ;res_lines@gsMarkerColor     = "blue"
    ;res_lines@gsMarkerIndex     = 16
    ;res_lines@gsMarkerSizeF     = 0.01
    
    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.025                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

        ;yy = (/100,100/)
    xx = (/0.5,1.0/)
    text_pos = 1.1
    do i = 0,dimsizes(Texts)-1
        yy = (/150+i*50,150+i*50/)
        res_lines@gsLineColor       = Colors(i)
        dum5(i) = gsn_add_polyline(wks, plot(0), xx, yy, res_lines)
        dum6(i) = gsn_add_text(wks, plot(0), Texts(i), text_pos, yy(0), res_text)
    end do
pres = True
pres@gsnMaximize = True
;pres@gsnPanelMainString = "CTRL (Sum of heating entire TP)"

gsn_panel(wks, plot, (/3,4/), pres)
end