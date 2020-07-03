
begin

lonW1 = 75
lonE1 = 105
lonW2 = 105
lonE2 = 120

casename = (/"CTRL","NOQ1","NOTP_NOQ1"/) ;,"NOTP"
ncase = dimsizes(casename)
; month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; nm = dimsizes(month)
season = (/"JJA"/)

u = new((/ncase,37,96,144/), float)
v = u
w = u
t = u

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"

    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    fu = addfile(diri+filiu, "r")
    U1 = fu->U
    u1 = month_to_season(U1, season)
    u(i,:,:,:) = dim_avg_n_Wrap(u1, 0)
end do

; ux = u(0:1,:,:,:)
; ux(0,:,:,:) = u(0,:,:,:) - u(1,:,:,:)
; ux(1,:,:,:) = u(1,:,:,:) - u(2,:,:,:)
; ux_south = dim_avg_n_Wrap(ux(:,:,:,{lonW1:lonE1}), 3)
; ux_east  = dim_avg_n_Wrap(ux(:,:,:,{lonW2:lonE2}), 3)

u_south = dim_avg_n_Wrap(u(:,:,:,{lonW1:lonE1}), 3)
u_east  = dim_avg_n_Wrap(u(:,:,:,{lonW2:lonE2}), 3)

ux_south = u_south(0:1,:,:)
ux_south(0,:,:) = u_south(0,:,:) - u_south(1,:,:)
ux_south(1,:,:) = u_south(1,:,:) - u_south(2,:,:)

ux_east = u_east(0:1,:,:)
ux_east(0,:,:) = u_east(0,:,:) - u_east(1,:,:)
ux_east(1,:,:) = u_east(1,:,:) - u_east(2,:,:)

printVarSummary(ux_east)

;----------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200102/"

figname = "JJA_U_Shear"

wks = gsn_open_wks("eps", figpath+figname)
nplot = 6
plot = new(nplot, graphic)
plotc = plot
plot_topo = plot


;=====================================================
   res_topo            = True
   res_topo@cnFillOn   = False
   res_topo@cnLinesOn  = True
   res_topo@gsnFrame   = False
   res_topo@gsnDraw    = False
   res_topo@cnLineLabelsOn         = False
   res_topo@cnLineThicknessF       = 5
   res_topo@cnLineDashPattern      = 2
   res_topo@cnLineColor            = "gray40"
   res_topo@cnInfoLabelOn          = False
   res_topo@cnLevelSelectionMode   ="ExplicitLevels"
   res_topo@cnLevels               = (/0.9/)
   res_topo@gsnLeftString = ""
   res_topo@gsnRightString = ""
;=====================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.7

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@gsnStringFontHeightF = 0.03
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    ;res@pmLabelBarOrthogonalPosF = 0.12
    res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@tiYAxisString = ""

;=====================================================

    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "green4"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = False
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnPresHgtHeightLabelOn = False

;=====================================================

res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10
res@cnLevelSpacingF = 2

cnres@cnLevelSelectionMode = "ExplicitLevels"
cnres@cnLevels = (/-30,-25,-20,-15,-10,-5,5,10,15,20,25,30/)

cnres@gsnContourPosLineDashPattern = 0
cnres@gsnContourNegLineDashPattern = 2

res@gsnLeftString = "(a) Thermal"
res@gsnRightString = lonW1+"~S~o~N~E-"+lonE1+"~S~o~N~E"
plot(0) = gsn_csm_pres_hgt(wks, ux_south(0,{1000:100},{-20:80}), res)
plotc(0) = gsn_csm_pres_hgt(wks, u_south(0,{1000:100},{-20:80}), cnres)
overlay(plot(0), plotc(0))

res@gsnLeftString = "(b) Dynamical"
plot(1) = gsn_csm_pres_hgt(wks, ux_south(1,{1000:100},{-20:80}), res)
plotc(1) = gsn_csm_pres_hgt(wks, u_south(1,{1000:100},{-20:80}), cnres)
overlay(plot(1), plotc(1))

res@gsnLeftString = "(c) Thermal"
res@gsnRightString = lonW2+"~S~o~N~E-"+lonE2+"~S~o~N~E"
plot(2) = gsn_csm_pres_hgt(wks, ux_east(0,{1000:100},{-20:80}), res)
plotc(2) = gsn_csm_pres_hgt(wks, u_east(0,{1000:100},{-20:80}), cnres)
overlay(plot(2), plotc(2))

res@gsnLeftString = "(d) Thermal"
plot(3) = gsn_csm_pres_hgt(wks, ux_east(1,{1000:100},{-20:80}), res)
plotc(3) = gsn_csm_pres_hgt(wks, u_east(1,{1000:100},{-20:80}), cnres)
overlay(plot(3), plotc(3))

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar  = True
gsn_panel(wks, plot, (/2,2/), pres)

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMaxLevelValF = 20
; res@cnMinLevelValF = -20
; res@cnLevelSpacingF = 4
; res@lbLabelBarOn = False
; do i = 0,2
;     res@gsnLeftString = casename(i)
;     res@gsnRightString = lonW1+"~S~o~N~E-"+lonE1+"~S~o~N~E"
;     plot(2*i) = gsn_csm_pres_hgt(wks, u_south(i,{1000:100},{-20:80}), res)
;     res@gsnRightString = lonW2+"~S~o~N~E-"+lonE2+"~S~o~N~E"
;     plot(2*i+1) = gsn_csm_pres_hgt(wks, u_east(i,{1000:100},{-20:80}), res)
; end do

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; gsn_panel(wks, plot, (/3,2/), pres)


end

