begin

lonstr = 80
lonend = 105
latstr = 30
latend = 40
yrStrt = 198102
yrLast = 200601

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)
nyear = 25
LEV = (/200,500,850/)
nlev = dimsizes(LEV)
U = new((/ncase,nyear,nlev,96,144/), float)
V = U
SPEED = U
Z = U
T = new((/ncase,nyear,96,144/), float)

season = "DJF"

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    filiu = casename(i)+".cam.h0.U.1979-2005.nc"
    filiv = casename(i)+".cam.h0.V.1979-2005.nc"    
    filiz = casename(i)+".cam.h0.Z3.1979-2005.nc"
    filit = casename(i)+".cam.h0.TS.1979-2005.nc"

    fu = addfile(diri+filiu, "r")
    fv = addfile(diri+filiv, "r")
    fz = addfile(diri+filiz, "r")
    ft = addfile(diri+filit, "r")

    TIME = fu->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    UX = fu->U(iYYYY,{LEV},:,:)
    VX = fv->V(iYYYY,{LEV},:,:)
    SPEEDX = wind_speed(UX, VX)
    ZX = fz->Z3(iYYYY,{LEV},:,:)
    TX = ft->TS(iYYYY,:,:)

    ; U(i,:,:,:) = month_to_annual(UX, 1)
    ; V(i,:,:,:) = month_to_annual(VX, 1)
    ; SPEED(i,:,:,:) = month_to_annual(SPEEDX, 1)
    U(i,:,:,:,:) = month_to_season(UX, season)
    V(i,:,:,:,:) = month_to_season(VX, season)
    SPEED(i,:,:,:,:) = month_to_season(SPEEDX, season)
    Z(i,:,:,:,:) = month_to_season(ZX, season)
    T(i,:,:,:) = month_to_season(TX, season)

delete([/diri,filiu,filiv,filiz,filit,fu,fv,fz,ft,UX,VX,SPEEDX,ZX,TX/])
end do
printVarSummary(U)

;=======================================================================
; ttest
;=======================================================================

u_ave = dim_avg_n_Wrap(U, 1)
u_var = dim_variance_n_Wrap(U, 1)
probu = ttest(u_ave(0,:,:,:), u_var(0,:,:,:), nyear, u_ave(1,:,:,:), u_var(1,:,:,:), nyear, False, False)
copy_VarCoords(u_ave(0,:,:,:), probu)
printVarSummary(probu)

v_ave = dim_avg_n_Wrap(V, 1)
v_var = dim_variance_n_Wrap(V, 1)
probv = ttest(v_ave(0,:,:,:), v_var(0,:,:,:), nyear, v_ave(1,:,:,:), v_var(1,:,:,:), nyear, False, False)
copy_VarCoords(v_ave(0,:,:,:), probv)
printVarSummary(probv)

z_ave = dim_avg_n_Wrap(Z, 1)
z_var = dim_variance_n_Wrap(Z, 1)
probz = ttest(z_ave(0,:,:,:), z_var(0,:,:,:), nyear, z_ave(1,:,:,:), z_var(1,:,:,:), nyear, False, False)
copy_VarCoords(z_ave(0,:,:,:), probz)
printVarSummary(probz)

t_ave = dim_avg_n_Wrap(T, 1)
t_var = dim_variance_n_Wrap(T, 1)
probt = ttest(t_ave(0,:,:), t_var(0,:,:), nyear, t_ave(1,:,:), t_var(1,:,:), nyear, False, False)
copy_VarCoords(t_ave(0,:,:), probt)
printVarSummary(probt)


diff_u = u_ave(0,:,:,:)
diff_u = u_ave(1,:,:,:) - u_ave(0,:,:,:)
diff_v = v_ave(0,:,:,:)
diff_v = v_ave(1,:,:,:) - v_ave(0,:,:,:)
diff_z = z_ave(0,:,:,:)
diff_z = z_ave(1,:,:,:) - z_ave(0,:,:,:)
diff_t = t_ave(0,:,:)
diff_t = t_ave(1,:,:) - t_ave(0,:,:)
;=========================================================
; outline of tp 
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200408/"
figname = season+"_200hPa_500hPa_850hPa_circulation"

; lonstr = 0
; lonend = 360
; latstr = -90
; latend = 90
lonstr = 0
lonend = 360
latstr = -30
latend = 60

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(6, graphic)
plotv = plot
plot_topo = plot
plot_prob = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

res = True
res@gsnMaximize = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpFillOn = True
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpMinLatF = latstr
res@mpMaxLatF = latend

res@cnFillOn = True
res@cnFillPalette = "BlueWhiteOrangeRed";WhiteBlueGreenYellowRed
;res@lbLabelBarOn = False
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
res@pmLabelBarOrthogonalPosF = 0.2


res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

vcres = True
vcres@gsnDraw  = False
vcres@gsnFrame = False
vcres@vcMinDistanceF            = 0.025
vcres@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
vcres@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
;vcres@vcRefAnnoArrowUseVecColor = False
vcres@vcRefAnnoString2On = False
vcres@vcRefAnnoFontHeightF = 0.02

vcres@gsnLeftString = ""
vcres@gsnRightString = ""

vcres@vcGlyphStyle            = "LineArrow"     ; turn on curly vectors CurlyVector
;vcres@vcLineArrowColor        = "white"           ; change vector color
vcres@vcLineArrowThicknessF   = 3.0               ; change vector thickness
;vcres@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last


res@gsnLeftString = "UV200 & U200 (Jetstream)"
res@gsnRightString = "Ctrl"
vcres@vcRefMagnitudeF           = 40.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
; res@cnMinLevelValF       = -0.1               ; set min contour level
; res@cnMaxLevelValF       = 0.1                 ; set max contour level
; res@cnLevelSpacingF      = 0.01              ; set contour spacing
res@cnMinLevelValF       = -40               ; set min contour level
res@cnMaxLevelValF       = 40                 ; set max contour level
res@cnLevelSpacingF      = 10.              ; set contour spacing
plot(0) = gsn_csm_contour_map(wks, u_ave(0,0,:,:), res)
plotv(0) = gsn_csm_vector(wks, u_ave(0,0,:,:), v_ave(0,0,:,:), vcres)
overlay(plot(0), plotv(0))
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))

res@gsnLeftString = "UV200 & U200 (Jetstream)"
res@gsnRightString = "diff"
vcres@vcRefMagnitudeF           = 5.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
; res@cnMinLevelValF       = -0.1               ; set min contour level
; res@cnMaxLevelValF       = 0.1                 ; set max contour level
; res@cnLevelSpacingF      = 0.01              ; set contour spacing
res@cnMinLevelValF       = -5               ; set min contour level
res@cnMaxLevelValF       = 5                 ; set max contour level
res@cnLevelSpacingF      = 0.5              ; set contour spacing
plot(1) = gsn_csm_contour_map(wks, diff_u(0,:,:), res)
plotv(1) = gsn_csm_vector(wks, diff_u(0,:,:), diff_v(0,:,:), vcres)
overlay(plot(1), plotv(1))
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


res@gsnLeftString = "UV500 & Z500"
res@gsnRightString = "Ctrl"
vcres@vcRefMagnitudeF           = 20.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
; res@cnMinLevelValF       = -0.1               ; set min contour level
; res@cnMaxLevelValF       = 0.1                 ; set max contour level
; res@cnLevelSpacingF      = 0.01              ; set contour spacing
res@cnMinLevelValF       = 5000               ; set min contour level
res@cnMaxLevelValF       = 6000                 ; set max contour level
res@cnLevelSpacingF      = 100.              ; set contour spacing
plot(2) = gsn_csm_contour_map(wks, z_ave(0,1,:,:), res)
plotv(2) = gsn_csm_vector(wks, u_ave(0,1,:,:), v_ave(0,1,:,:), vcres)
overlay(plot(2), plotv(2))
plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2), plot_topo(2))

res@gsnLeftString = "UV500 & Z500"
res@gsnRightString = "diff"
vcres@vcRefMagnitudeF           = 2.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
res@cnMinLevelValF       = -30               ; set min contour level
res@cnMaxLevelValF       = 30                 ; set max contour level
res@cnLevelSpacingF      = 5.              ; set contour spacing
plot(3) = gsn_csm_contour_map(wks, diff_z(1,:,:), res)
plotv(3) = gsn_csm_vector(wks, diff_u(1,:,:), diff_v(1,:,:), vcres)
overlay(plot(3), plotv(3))
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))



res@gsnLeftString = "UV800 & Ts"
res@gsnRightString = "Ctrl"
vcres@vcRefMagnitudeF           = 10.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
; res@cnMinLevelValF       = -0.1               ; set min contour level
; res@cnMaxLevelValF       = 0.1                 ; set max contour level
; res@cnLevelSpacingF      = 0.01              ; set contour spacing
res@cnMinLevelValF       = 240               ; set min contour level
res@cnMaxLevelValF       = 300                 ; set max contour level
res@cnLevelSpacingF      = 5.              ; set contour spacing
plot(4) = gsn_csm_contour_map(wks, t_ave(0,:,:), res)
plotv(4) = gsn_csm_vector(wks, u_ave(0,2,:,:), v_ave(0,2,:,:), vcres)
overlay(plot(4), plotv(4))
plot_topo(4) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(4), plot_topo(4))

res@gsnLeftString = "UV850 & Ts"
res@gsnRightString = "diff"
vcres@vcRefMagnitudeF           = 2.0             ; define vector ref mag
vcres@vcRefLengthF              = 0.045           ; define length of vec ref
res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
res@cnMinLevelValF       = -2               ; set min contour level
res@cnMaxLevelValF       = 2                 ; set max contour level
res@cnLevelSpacingF      = 0.2              ; set contour spacing
plot(5) = gsn_csm_contour_map(wks, diff_t, res)
plotv(5) = gsn_csm_vector(wks, diff_u(2,:,:), diff_v(2,:,:), vcres)
overlay(plot(5), plotv(5))
plot_topo(5) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(5), plot_topo(5))

;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "green4"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.005
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

 plot_prob(0) = gsn_csm_contour(wks,probu(0,:,:),res3)
 plot_prob(0) = ShadeLtContour(plot_prob(0),0.055,17)
 overlay(plot(1), plot_prob(0))

 plot_prob(1) = gsn_csm_contour(wks,probz(1,:,:),res3)
 plot_prob(1) = ShadeLtContour(plot_prob(1),0.055,17)
 overlay(plot(3), plot_prob(1))

 plot_prob(2) = gsn_csm_contour(wks,probt,res3)
 plot_prob(2) = ShadeLtContour(plot_prob(2),0.055,17)
 overlay(plot(5), plot_prob(2))

 diff_u = mask(diff_u, probu.le.0.055 .and. probv.le.0.055, True)
 diff_v = mask(diff_v, probu.le.0.055 .and. probv.le.0.055, True)
 vcres@vcLineArrowColor = "red"
 vcres@vcRefMagnitudeF           = 5.0             ; define vector ref mag
 plot_prob(3) = gsn_csm_vector(wks, diff_u(0,:,:), diff_v(0,:,:), vcres)
 overlay(plot(1), plot_prob(3))
 
 vcres@vcRefMagnitudeF           = 2.0             ; define vector ref mag
 plot_prob(4) = gsn_csm_vector(wks, diff_u(1,:,:), diff_v(1,:,:), vcres)
 overlay(plot(3), plot_prob(4))

 vcres@vcRefMagnitudeF           = 2.0             ; define vector ref mag
 plot_prob(5) = gsn_csm_vector(wks, diff_u(2,:,:), diff_v(2,:,:), vcres)
 overlay(plot(5), plot_prob(5))






pres = True
pres@gsnDraw = True
pres@gsnFrame = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" mean: Circulation"

gsn_panel(wks, plot, (/3,2/), pres)

frame(wks)
pres@gsnPanelMainString = ""
gsn_panel(wks, plot(1::2), (/3,1/), pres)
end 