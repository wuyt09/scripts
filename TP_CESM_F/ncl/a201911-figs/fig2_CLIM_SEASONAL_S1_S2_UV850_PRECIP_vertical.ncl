begin

;=================================================================
; NOQ1, NOTP_NOQ1
;=================================================================
lonstr = 40 ;20
lonend = 180 ;220
latstr = -10
latend = 70
season = (/"MAM","JJA","SON","DJF"/)
;==================================================================
;==================================================================
dirinq1 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"

nq1p1name = "F2000_NOQ1.PRECC.mon.0601-3012.nc"
nq1p2name = "F2000_NOQ1.PRECL.mon.0601-3012.nc"
nq1tname = "F2000_NOQ1.TS.mon.0601-3012.nc"
nq1uname = "F2000_NOQ1.plevel.U.mon.0601-3012.nc"
nq1vname = "F2000_NOQ1.plevel.V.mon.0601-3012.nc"
nq1zname = "F2000_NOQ1.PHIS.mon.0601-3012.nc"

fnq1p1 = addfile(dirinq1+nq1p1name, "r")
fnq1p2 = addfile(dirinq1+nq1p2name, "r")
fnq1t  = addfile(dirinq1+nq1tname, "r")
fnq1u  = addfile(dirinq1+nq1uname, "r")
fnq1v  = addfile(dirinq1+nq1vname, "r")
fnq1z  = addfile(dirinq1+nq1zname, "r")

P1_nq1 = fnq1p1->PRECC(:,{latstr:latend},{lonstr:lonend})
P2_nq1 = fnq1p2->PRECL(:,{latstr:latend},{lonstr:lonend})
P_nq1 = P1_nq1
P_nq1 = (P1_nq1 + P2_nq1)*1000*86400.
T_nq1 = fnq1t->TS(:,{latstr:latend},{lonstr:lonend})
U_nq1 = fnq1u->U(:,{850},{latstr:latend},{lonstr:lonend})
V_nq1 = fnq1v->V(:,{850},{latstr:latend},{lonstr:lonend})
Z_nq1 = fnq1z->PHIS(:,{latstr:latend},{lonstr:lonend})
Z_nq1 = Z_nq1/9.8

p_nq1 = month_to_seasonN(P_nq1, season)
t_nq1 = month_to_seasonN(T_nq1, season)
u_nq1 = month_to_seasonN(U_nq1, season)
v_nq1 = month_to_seasonN(V_nq1, season)
z_nq1 = month_to_seasonN(Z_nq1, season)
nq1_p = dim_avg_n_Wrap(p_nq1, 1)
nq1_t = dim_avg_n_Wrap(t_nq1, 1)
nq1_u = dim_avg_n_Wrap(u_nq1, 1)
nq1_v = dim_avg_n_Wrap(v_nq1, 1)
nq1_z = dim_avg_n_Wrap(z_nq1, 1)

print("end of NOQ1")
;==================================================================
dirintp = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

ntpp1name = "F2000_NOTP_NOQ1.PRECC.mon.0601-3012.nc"
ntpp2name = "F2000_NOTP_NOQ1.PRECL.mon.0601-3012.nc"
ntptname = "F2000_NOTP_NOQ1.TS.mon.0601-3012.nc"
ntpuname = "F2000_NOTP_NOQ1.plevel.U.mon.0601-3012.nc"
ntpvname = "F2000_NOTP_NOQ1.plevel.V.mon.0601-3012.nc"
ntpzname = "F2000_NOTP_NOQ1.PHIS.mon.0601-3012.nc"

fntpp1 = addfile(dirintp+ntpp1name, "r")
fntpp2 = addfile(dirintp+ntpp2name, "r")
fntpt  = addfile(dirintp+ntptname, "r")
fntpu  = addfile(dirintp+ntpuname, "r")
fntpv  = addfile(dirintp+ntpvname, "r")
fntpz  = addfile(dirintp+ntpzname, "r")

P1_ntp = fntpp1->PRECC(:,{latstr:latend},{lonstr:lonend})
P2_ntp = fntpp2->PRECL(:,{latstr:latend},{lonstr:lonend})
P_ntp = P1_ntp
P_ntp = (P1_ntp + P2_ntp)*1000*86400.
T_ntp = fntpt->TS(:,{latstr:latend},{lonstr:lonend})
U_ntp = fntpu->U(:,{850},{latstr:latend},{lonstr:lonend})
V_ntp = fntpv->V(:,{850},{latstr:latend},{lonstr:lonend})
Z_ntp = fntpz->PHIS(:,{latstr:latend},{lonstr:lonend})
Z_ntp = Z_ntp/9.8

p_ntp = month_to_seasonN(P_ntp, season)
t_ntp = month_to_seasonN(T_ntp, season)
u_ntp = month_to_seasonN(U_ntp, season)
v_ntp = month_to_seasonN(V_ntp, season)
z_ntp = month_to_seasonN(Z_ntp, season)
ntp_p = dim_avg_n_Wrap(p_ntp, 1)
ntp_t = dim_avg_n_Wrap(t_ntp, 1)
ntp_u = dim_avg_n_Wrap(u_ntp, 1)
ntp_v = dim_avg_n_Wrap(v_ntp, 1)
ntp_z = dim_avg_n_Wrap(z_ntp, 1)

print("end of NOTP_NOQ1")
;=================================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "fig2_CLIM_SEASONAL_S1_S2_UV850_PRECIP_vertical"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 8
plot = new(nplot, graphic)
plotv = plot
plotc = plot
plot_topo = plot
plot_map  = plot

;----------------------------------------------------
mpres = True
mpres@gsnFrame = False
mpres@gsnDraw  = False
;mpres@gsnAddCyclic = False

mpres@mpOutlineOn = True
mpres@mpFillOn = True
mpres@mpGeophysicalLineThicknessF = 2
mpres@mpGeophysicalLineColor = "grey45" ;"grey25"
mpres@mpFillDrawOrder = "PreDraw"
mpres@mpOutlineDrawOrder = "PostDraw"
;mpres@mpOceanFillColor       = "lightskyblue1"
;mpres@mpLandFillColor        = "gray"

mpres@mpMinLatF = latstr
mpres@mpMaxLatF = latend
mpres@mpMinLonF = lonstr
mpres@mpMaxLonF = lonend
mpres@mpCenterLonF = (lonstr+lonend)/2

mpres@gsnStringFontHeightF = 0.035
mpres@tmXBLabelFontHeightF = 0.025
mpres@tmYLLabelFontHeightF = 0.025

Alefttitle = (/"(a) ","(b) ","(c) ","(d) "/)
Blefttitle = (/"(e) ","(f) ","(g) ","(h) "/)

do i = 0, dimsizes(season)-1
mpres@gsnLeftString = Alefttitle(i)+"NO_HEAT"
mpres@gsnRightString = season(i)
plot_map(2*i) = gsn_csm_map(wks, mpres)

mpres@gsnLeftString = Blefttitle(i)+"NO_TP&NO_HEAT"
mpres@gsnRightString = season(i)
plot_map(2*i+1) = gsn_csm_map(wks, mpres)
end do
;----------------------------------------------------
res_topo = True
res_topo@gsnFrame = False
res_topo@gsnDraw  = False
res_topo@gsnAddCyclic = False
res_topo@gsnLeftString = ""
res_topo@gsnRightString = ""

res_topo@cnFillOn = False
res_topo@cnLinesOn = True
res_topo@cnLineLabelsOn = False
res_topo@cnInfoLabelOn = False
res_topo@cnLineThicknessF = 3

res_topo@cnLevelSelectionMode = "ExplicitLevels"
res_topo@cnLevels = (/1500/)
res_topo@cnLineColor = "purple"
res_topo@cnLineDashPattern = 1

do i = 0, dimsizes(season)-1
plot_topo(2*i) = gsn_csm_contour(wks, nq1_z(i,:,:), res_topo)
plot_topo(2*i+1) = gsn_csm_contour(wks, ntp_z(i,:,:), res_topo)
end do
;----------------------------------------------------
;----------------------------------------------------

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnAddCyclic = False
res@gsnLeftString = ""
res@gsnRightString = ""

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
;res@pmLabelBarOrthogonalPosF = 0.15

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 2.
res@cnMaxLevelValF = 16.
res@cnMinLevelValF = 0.
res@cnFillPalette = "CBR_wet"

do i = 0, dimsizes(season)-1
;res@gsnRightString = season(i)
;res@gsnLeftString = "NO_HEAT"
plot(2*i) = gsn_csm_contour(wks, nq1_p(i,:,:), res)

;res@gsnLeftString = "NO_TP + NO_HEAT"
plot(2*i+1) = gsn_csm_contour(wks, ntp_p(i,:,:), res)
end do

;----------------------------------------------------
  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@gsnLeftString = ""
  vcres@gsnRightString = ""
  
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
  vcres@vcRefAnnoOrthogonalPosF =  -1.0             ; move ref vector up
  vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag

  do i = 0, dimsizes(season)-1
  plotv(2*i) = gsn_csm_vector(wks, nq1_u(i,:,:), nq1_v(i,:,:), vcres)
  
  plotv(2*i+1) = gsn_csm_vector(wks, ntp_u(i,:,:), ntp_v(i,:,:), vcres)
  end do
;----------------------------------------------------

    cnres = True
    cnres@gsnDraw  = False
    cnres@gsnFrame = False
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnInfoLabelOn = False
    cnres@cnLineThicknessF = 3
    cnres@cnLineColor = "red"
    cnres@cnLineLabelFontHeightF = 0.02
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLabelDrawOrder = "PostDraw"
    cnres@cnLineLabelPerimColor = "green4"
    ;cnres@cnLinePalette = "CBR_coldhot"

    cnres@cnLevelSelectionMode = "ManualLevels"
    cnres@cnMinLevelValF = 270
    cnres@cnMaxLevelValF = 300
    cnres@cnLevelSpacingF = 10

    
    do i = 0, dimsizes(season)-1
    plotc(2*i) = gsn_csm_contour(wks, nq1_t(i,:,:), cnres)

    plotc(2*i+1) = gsn_csm_contour(wks, ntp_t(i,:,:), cnres)
    end do

;----------------------------------------------------
do i = 0, dimsizes(season)-1
overlay(plot_map(2*i), plot(2*i))
overlay(plot_map(2*i), plot_topo(2*i))
overlay(plot_map(2*i), plotc(2*i))
overlay(plot_map(2*i), plotv(2*i))

overlay(plot_map(2*i+1), plot(2*i+1))
overlay(plot_map(2*i+1), plot_topo(2*i+1))
overlay(plot_map(2*i+1), plotc(2*i+1))
overlay(plot_map(2*i+1), plotv(2*i+1))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarWidthF = 0.55

gsn_panel(wks, plot_map, (/4,2/), pres)

end