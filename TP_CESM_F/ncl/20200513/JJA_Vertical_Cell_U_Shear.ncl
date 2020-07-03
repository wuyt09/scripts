
begin

latS   =  0
latN   =  60

lonW1 = 70
lonE1 = 90
lonW2 = 110
lonE2 = 120

casename = (/"CTRL","NOQ1"/) ;,"NOTP","NOTP_NOQ1"
ncase = dimsizes(casename)
; month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; nm = dimsizes(month)
season = (/"JJA"/)
nyear= 25

u_lon = new((/ncase,nyear,2,37,96/), float)
u_lon!0 = "case"
u_lon&case = casename
u_lon!2 = "cell"
; u_lon&cell = (/lonW1+"~S~o~N~E-"+lonE1+"~S~o~N~E",lonW2+"~S~o~N~E-"+lonE2+"~S~o~N~E"/)
u_lon&cell = (/"70-90E","110-120E"/)
v_lon = u_lon
w_lon = u_lon

do i = 0, dimsizes(casename)-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"

    filiu = "F2000_"+casename(i)+".plevel.U.mon.0601-3012.nc"
    fu = addfile(diri+filiu, "r")
    U1 = fu->U
    U_lon1 = dim_avg_n_Wrap(U1(:,:,:,{lonW1:lonE1}), 3)
    U_lon2 = dim_avg_n_Wrap(U1(:,:,:,{lonW2:lonE2}), 3)
    u_lon1 = month_to_season(U_lon1, season)
    u_lon2 = month_to_season(U_lon2, season)
    u_lon(i,:,0,:,:) = u_lon1
    u_lon(i,:,1,:,:) = u_lon2


    filiv = "F2000_"+casename(i)+".plevel.V.mon.0601-3012.nc"
    fv = addfile(diri+filiv, "r")
    V1 = fv->V
    V_lon1 = dim_avg_n_Wrap(V1(:,:,:,{lonW1:lonE1}), 3)
    V_lon2 = dim_avg_n_Wrap(V1(:,:,:,{lonW2:lonE2}), 3)
    v_lon1 = month_to_season(V_lon1, season)
    v_lon2 = month_to_season(V_lon2, season)
    v_lon(i,:,0,:,:) = v_lon1
    v_lon(i,:,1,:,:) = v_lon2

    filiw = "F2000_"+casename(i)+".plevel.OMEGA.mon.0601-3012.nc"
    fw = addfile(diri+filiw, "r")
    W1 = fw->OMEGA
    W1 = W1*(-1.)
    W_lon1 = dim_avg_n_Wrap(W1(:,:,:,{lonW1:lonE1}), 3)
    W_lon2 = dim_avg_n_Wrap(W1(:,:,:,{lonW2:lonE2}), 3)
    w_lon1 = month_to_season(W_lon1, season)
    w_lon2 = month_to_season(W_lon2, season)
    w_lon(i,:,0,:,:) = w_lon1
    w_lon(i,:,1,:,:) = w_lon2
end do
printVarSummary(w_lon)
;============================================================
; ttest
;============================================================
iflag = False
sig = 0.05

u_lon_ave = dim_avg_n_Wrap(u_lon, 1)
u_lon_var = dim_variance_n_Wrap(u_lon, 1)

v_lon_ave = dim_avg_n_Wrap(v_lon, 1)
v_lon_var = dim_variance_n_Wrap(v_lon, 1)

w_lon_ave = dim_avg_n_Wrap(w_lon, 1)
w_lon_var = dim_variance_n_Wrap(w_lon, 1)
printVarSummary(w_lon_ave)

probu = ttest(u_lon_ave(0,:,:,:), u_lon_var(0,:,:,:), nyear, u_lon_ave(1,:,:,:), u_lon_var(1,:,:,:), nyear, iflag, False)
probv = ttest(v_lon_ave(0,:,:,:), v_lon_var(0,:,:,:), nyear, v_lon_ave(1,:,:,:), v_lon_var(1,:,:,:), nyear, iflag, False)
probw = ttest(w_lon_ave(0,:,:,:), w_lon_var(0,:,:,:), nyear, w_lon_ave(1,:,:,:), w_lon_var(1,:,:,:), nyear, iflag, False)
copy_VarCoords(u_lon_ave(0,:,:,:), probu)
copy_VarCoords(v_lon_ave(0,:,:,:), probv)
copy_VarCoords(w_lon_ave(0,:,:,:), probw)
printVarSummary(probw)
;============================================================

u_lon_diff = u_lon_ave(0,:,:,:) - u_lon_ave(1,:,:,:)
v_lon_diff = v_lon_ave(0,:,:,:) - v_lon_ave(1,:,:,:)
w_lon_diff = w_lon_ave(0,:,:,:) - w_lon_ave(1,:,:,:)
w_lon_diff = w_lon_diff * 100.
copy_VarCoords(u_lon_ave(0,:,:,:), u_lon_diff)
copy_VarCoords(v_lon_ave(0,:,:,:), v_lon_diff)
copy_VarCoords(w_lon_ave(0,:,:,:), w_lon_diff)

; u_lon_diff = mask(u_lon_diff, probu.le.sig, True)
v_lon_diff_sig = mask(v_lon_diff, probv.le.sig, True)
w_lon_diff_sig = mask(w_lon_diff, probw.le.sig, True)
v_lon_diff_nosig = mask(v_lon_diff, probv.gt.sig, True)
w_lon_diff_nosig = mask(w_lon_diff, probw.gt.sig, True)
copy_VarCoords(v_lon_diff, v_lon_diff_nosig)
copy_VarCoords(v_lon_diff, v_lon_diff_sig)
copy_VarCoords(w_lon_diff, w_lon_diff_nosig)
copy_VarCoords(w_lon_diff, w_lon_diff_sig)


;----------------------------------------------------------------
; topo file
;----------------------------------------------------------------

ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
topo_lon1 = dim_avg_n_Wrap(topo(:,{lonW1:lonE1}), 1)
topo_lon2 = dim_avg_n_Wrap(topo(:,{lonW2:lonE2}), 1)
elev_lon1 = 1013.25*(1-topo_lon1*0.0065/288.15)^5.25145
elev_lon2 = 1013.25*(1-topo_lon2*0.0065/288.15)^5.25145
copy_VarCoords(topo_lon1, elev_lon1)
copy_VarCoords(topo_lon2, elev_lon2)
;----------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"

figname = "JJA_Vertival_Cell_U_Shear"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 2
plot = new(nplot, graphic)
plotc = plot
plot_topo = plot
plot_prob = plot
plot_sig = plot


;=====================================================
   res_topo            = True
   res_topo@gsnFrame   = False
   res_topo@gsnDraw    = False
   res_topo@gsnLeftString = ""
   res_topo@gsnRightString = ""
   ; res_topo@gsnMaximize         = True
   ; res_topo@trYMaxF          = 1000.
   res_topo@xyLineThicknessF = 5
   res_topo@gsnYRefLine      = (/100,1000/)
   res_topo@gsnBelowYRefLineColor            = (/"black","black"/)
;=====================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.8

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@gsnStringFontHeightF = 0.025
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    ;res@pmLabelBarOrthogonalPosF = 0.12
    res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@tiYAxisString = ""


    res@vcRefLengthF            = 0.04            ;the length of the reference magnitude
    res@vcRefMagnitudeF         = 5.0
;     res@vcRefAnnoFontHeightF    = 0.0095
    res@vcMinDistanceF          = 0.018            ;the minimum distance of neighboring vectors
    res@vcGlyphStyle            = "CurlyVector"   ;the style of vector
;     res@vcLineArrowHeadMaxSizeF = 0.011           ;change the vector line arrows
;     res@vcLineArrowHeadMinSizeF = 0.007
    res@vcLineArrowThicknessF   = 3.0
;     res@vcMapDirection          = False           ;vectors to be pointing in the right direction
;     res@vcRefAnnoOn             = True           ;turn on vector reference box
;     res@vcRefAnnoPerimOn        = False
;     res@vcRefAnnoString2        = " "
    ; res@vpClipOn               = False      ;turn off vector clipping
    res@vcRefAnnoOrthogonalPosF = -0.4          ;move ref vector down
;     res@vcRefAnnoParallelPosF   = 1.0
    res@vcRefAnnoString2On = False
    res@vcRefAnnoFontHeightF = 0.02
    res@vcLineArrowColor = "black" ;grey50

    res@tmYLMode = "Explicit"
    res@tmYLValues = (/1000,700,500,300,200,100/)
    res@tmYLLabels = (/1000,700,500,300,200,100/)

;=====================================================

    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "green4"
    cnres@cnLineThicknessF = 8
    cnres@cnInfoLabelOn = False
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnPresHgtHeightLabelOn = False

;=====================================================

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
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------
 ; sig vector
    resv          = True
    resv@tmYLMode = "Explicit" 
    resv@tmYLValues = (/1000,700,500,300,200,100/)
    resv@tmYLLabels = (/1000,700,500,300,200,100/)
    resv@gsnDraw                 = False
    resv@gsnFrame                = False      ;don't advance frame
    resv@vcRefLengthF            = 0.07            ;the length of the reference magnitude
    resv@vcRefAnnoFontHeightF    = 0.0095
    resv@vcMinDistanceF          = 0.015            ;the minimum distance of neighboring vectors
;   res@vcGlyphStyle            = "CurlyVector"   ;the style of vector
    resv@vcLineArrowHeadMaxSizeF = 0.011           ;change the vector line arrows
    resv@vcLineArrowHeadMinSizeF = 0.007
    resv@vcLineArrowThicknessF   = 1.0
    resv@vcMapDirection          = False           ;vectors to be pointing in the right direction
    resv@vcRefAnnoOn             = True           ;turn on vector reference box
    resv@vcRefAnnoPerimOn        = False
    resv@vcRefAnnoString2        = " "
    resv@vpClipOn                = False      ;turn off vector clipping
;    resv@vcVectorDrawOrder       = "PostDraw"      ;draw vectors last
    resv@vcRefAnnoOrthogonalPosF = -0.1          ;move ref vector down
    resv@vcRefAnnoParallelPosF   = 1.0

    resv@vcRefMagnitudeF       = 4.0               ;the magnitude used as the reference magnitude
    resv@vcRefAnnoString1      = "4"
    resv@vcLineArrowColor      = "black"           ; white vectors
    resv@gsnLeftString         = " "
    resv@gsnRightString        = " "


    resv          = True
    resv@gsnDraw                 = False
    resv@gsnFrame                = False      ;don't advance frame
    resv@vcRefLengthF            = 0.04            ;the length of the reference magnitude
    resv@vcRefMagnitudeF         = 5.0
    resv@vcMinDistanceF          = 0.018            ;the minimum distance of neighboring vectors
    resv@vcGlyphStyle            = "CurlyVector"   ;the style of vector
    resv@vcLineArrowThicknessF   = 3.0
    resv@vcRefAnnoOrthogonalPosF = -0.4          ;move ref vector down
    resv@vcRefAnnoString2On = False
    resv@vcRefAnnoFontHeightF = 0.02
    resv@vcLineArrowColor = "black"

    resv@tmYLMode = "Explicit"
    resv@tmYLValues = (/1000,700,500,300,200,100/)
    resv@tmYLLabels = (/1000,700,500,300,200,100/)
 ;----------------------------------------------
 ;----------------------------------------------

; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMaxLevelValF = 10
; res@cnMinLevelValF = -10
; res@cnLevelSpacingF = 2

cnres@cnLevelSelectionMode = "ExplicitLevels"
cnres@cnLevels = (/-30,-25,-20,-15,-10,-5,5,10,15,20,25,30/)

cnres@gsnContourPosLineDashPattern = 0
cnres@gsnContourNegLineDashPattern = 2

res@gsnLeftString = "(a) Diff in U and Vertival Cell"
res@gsnRightString = lonW1+"~S~o~N~E-"+lonE1+"~S~o~N~E"
plot(0) = gsn_csm_pres_hgt_vector(wks, u_lon_diff(0,{1000:100},{latS:latN}), v_lon_diff(0,{1000:100},{latS:latN}), w_lon_diff(0,{1000:100},{latS:latN}), res)

; res@vcLineArrowColor = "black"
; plot_sig(0) = gsn_csm_pres_hgt_vector(wks, u_lon_diff(0,{1000:100},{latS:latN}), v_lon_diff_sig(0,{1000:100},{latS:latN}), w_lon_diff_sig(0,{1000:100},{latS:latN}), res)
; plot_sig(0) = gsn_csm_vector(wks, v_lon_diff_sig(0,{1000:100},{latS:latN}), w_lon_diff_sig(0,{1000:100},{latS:latN}), resv)
; overlay(plot(0), plot_sig(0))

plotc(0) = gsn_csm_pres_hgt(wks, u_lon_ave(0,0,{1000:100},{latS:latN}), cnres)
overlay(plot(0), plotc(0))

plot_prob(0) = gsn_csm_contour(wks, probu(0,:,:), res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),0.055,17)
overlay(plot(0), plot_prob(0))

plot_topo(0) =  gsn_csm_xy(wks,elev_lon1&lat,elev_lon1,res_topo)
overlay(plot(0),plot_topo(0))



res@gsnLeftString = "(b) Diff in U and Vertival Cell"
res@gsnRightString = lonW2+"~S~o~N~E-"+lonE2+"~S~o~N~E"
plot(1) = gsn_csm_pres_hgt_vector(wks, u_lon_diff(1,{1000:100},{latS:latN}), v_lon_diff(1,{1000:100},{latS:latN}), w_lon_diff(1,{1000:100},{latS:latN}), res)

; res@vcLineArrowColor = "black"
; plot_sig(0) = gsn_csm_pres_hgt_vector(wks, u_lon_diff(0,{1000:100},{latS:latN}), v_lon_diff_sig(0,{1000:100},{latS:latN}), w_lon_diff_sig(0,{1000:100},{latS:latN}), res)
; plot_sig(0) = gsn_csm_vector(wks, v_lon_diff_sig(0,{1000:100},{latS:latN}), w_lon_diff_sig(0,{1000:100},{latS:latN}), resv)
; overlay(plot(0), plot_sig(0))

plotc(1) = gsn_csm_pres_hgt(wks, u_lon_ave(0,1,{1000:100},{latS:latN}), cnres)
overlay(plot(1), plotc(1))

plot_prob(1) = gsn_csm_contour(wks, probu(1,:,:), res3)
plot_prob(1) = ShadeLtContour(plot_prob(1),0.055,17)
overlay(plot(1), plot_prob(1))

plot_topo(1) =  gsn_csm_xy(wks,elev_lon2&lat,elev_lon2,res_topo)
overlay(plot(1),plot_topo(1))

; plotc(1) = gsn_csm_pres_hgt(wks, u_lon_ave(0,1,{1000:100},{latS:latN}), cnres)
; overlay(plot(1), plotc(1))

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar  = True
gsn_panel(wks, plot, (/2,1/), pres)

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

