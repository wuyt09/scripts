begin

;=================================================================
; OBS, CTRL, NOQ1, NOTP_NOQ1
;=================================================================
lonstr = 40 ;20
lonend = 180 ;220
latstr = -10
latend = 70

;==================================================================
yrStrt = 1979
yrLast = 2015

diriobs = "/home/yangsong3/data-observation/ERA-interim/monthly/"
diriobs_p = "/home/yangsong3/data-observation/GPCP/"

obspname = "precip.mon.mean.nc"
obstname = "2_metre_temperature.mon.mean.nc"
obsuname = "uwnd.mon.mean.nc"
obsvname = "vwnd.mon.mean.nc"
obssname = "Surface_pressure.mon.mean.nc"

fobsp = addfile(diriobs_p+obspname, "r")
TIMEp = fobsp->time
YYYYp = cd_calendar(TIMEp, -1)/100
iYYYYp = ind(YYYYp.ge.yrStrt .and. YYYYp.le.yrLast)
P_obs = fobsp->precip(iYYYYp,{latstr:latend},{lonstr:lonend})
printVarSummary(P_obs)

fobst = addfile(diriobs+obstname, "r")
TIMEt = fobst->time
YYYYt = cd_calendar(TIMEt, -1)/100
iYYYYt = ind(YYYYt.ge.yrStrt .and. YYYYt.le.yrLast)
T_obs = short2flt(fobst->t2m(iYYYYt,{latstr:latend},{lonstr:lonend}))
printVarSummary(T_obs)

fobsu = addfile(diriobs+obsuname, "r")
TIMEu = fobsu->time
YYYYu = cd_calendar(TIMEu, -1)/100
iYYYYu = ind(YYYYu.ge.yrStrt .and. YYYYu.le.yrLast)
U_obs = short2flt(fobsu->u(iYYYYu,{850},{latstr:latend},{lonstr:lonend}))
printVarSummary(U_obs)

fobsv = addfile(diriobs+obsvname, "r")
TIMEv = fobsv->time
YYYYv = cd_calendar(TIMEv, -1)/100
iYYYYv = ind(YYYYv.ge.yrStrt .and. YYYYv.le.yrLast)
V_obs = short2flt(fobsv->v(iYYYYv,{850},{latstr:latend},{lonstr:lonend}))
printVarSummary(V_obs)

fobss = addfile(diriobs+obssname, "r")
TIMEs = fobss->time
YYYYs = cd_calendar(TIMEs, -1)/100
iYYYYs = ind(YYYYs.ge.yrStrt .and. YYYYs.le.yrLast)
S_obs = short2flt(fobss->sp(iYYYYs,{latstr:latend},{lonstr:lonend}))
printVarSummary(S_obs)

delete([/TIMEp,TIMEt,TIMEu,TIMEv,TIMEs/])
delete([/YYYYp,YYYYt,YYYYu,YYYYv,YYYYs/])
delete([/iYYYYp,iYYYYt,iYYYYu,iYYYYv,iYYYYs/])

U_obs = mask(U_obs, S_obs.ge.85000., True)
V_obs = mask(V_obs, S_obs.ge.85000., True)

obs_p = dim_avg_n_Wrap(P_obs, 0)
obs_t = dim_avg_n_Wrap(T_obs, 0)
obs_u = dim_avg_n_Wrap(U_obs, 0)
obs_v = dim_avg_n_Wrap(V_obs, 0)

print("end of observation")

;==================================================================
dirictl = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"

ctlp1name = "F2000_CTRL.PRECC.mon.0601-3012.nc"
ctlp2name = "F2000_CTRL.PRECL.mon.0601-3012.nc"
ctltname = "F2000_CTRL.TS.mon.0601-3012.nc"
ctluname = "F2000_CTRL.U.mon.0601-3012.nc"
ctlvname = "F2000_CTRL.V.mon.0601-3012.nc"

fctlp1 = addfile(dirictl+ctlp1name, "r")
fctlp2 = addfile(dirictl+ctlp2name, "r")
fctlt  = addfile(dirictl+ctltname, "r")
fctlu  = addfile(dirictl+ctluname, "r")
fctlv  = addfile(dirictl+ctlvname, "r")

P1_ctl = fctlp1->PRECC(:,{latstr:latend},{lonstr:lonend})
P2_ctl = fctlp2->PRECL(:,{latstr:latend},{lonstr:lonend})
P_ctl = P1_ctl
P_ctl = (P1_ctl + P2_ctl)*1000*86400.  ; convert units from m/s to mm/day
T_ctl = fctlt->TS(:,{latstr:latend},{lonstr:lonend})
U_ctl = fctlu->U(:,{850},{latstr:latend},{lonstr:lonend})
V_ctl = fctlv->V(:,{850},{latstr:latend},{lonstr:lonend})

ctl_p = dim_avg_n_Wrap(P_ctl, 0)
ctl_t = dim_avg_n_Wrap(T_ctl, 0)
ctl_u = dim_avg_n_Wrap(U_ctl, 0)
ctl_v = dim_avg_n_Wrap(V_ctl, 0)

print("end of CTRL")

;==================================================================
dirinq1 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"

nq1p1name = "F2000_NOQ1.PRECC.mon.0601-3012.nc"
nq1p2name = "F2000_NOQ1.PRECL.mon.0601-3012.nc"
nq1tname = "F2000_NOQ1.TS.mon.0601-3012.nc"
nq1uname = "F2000_NOQ1.U.mon.0601-3012.nc"
nq1vname = "F2000_NOQ1.V.mon.0601-3012.nc"

fnq1p1 = addfile(dirinq1+nq1p1name, "r")
fnq1p2 = addfile(dirinq1+nq1p2name, "r")
fnq1t  = addfile(dirinq1+nq1tname, "r")
fnq1u  = addfile(dirinq1+nq1uname, "r")
fnq1v  = addfile(dirinq1+nq1vname, "r")

P1_nq1 = fnq1p1->PRECC(:,{latstr:latend},{lonstr:lonend})
P2_nq1 = fnq1p2->PRECL(:,{latstr:latend},{lonstr:lonend})
P_nq1 = P1_nq1
P_nq1 = (P1_nq1 + P2_nq1)*1000*86400.
T_nq1 = fnq1t->TS(:,{latstr:latend},{lonstr:lonend})
U_nq1 = fnq1u->U(:,{850},{latstr:latend},{lonstr:lonend})
V_nq1 = fnq1v->V(:,{850},{latstr:latend},{lonstr:lonend})

nq1_p = dim_avg_n_Wrap(P_nq1, 0)
nq1_t = dim_avg_n_Wrap(T_nq1, 0)
nq1_u = dim_avg_n_Wrap(U_nq1, 0)
nq1_v = dim_avg_n_Wrap(V_nq1, 0)

print("end of NOQ1")
;==================================================================
dirintp = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

ntpp1name = "F2000_NOTP_NOQ1.PRECC.mon.0601-3012.nc"
ntpp2name = "F2000_NOTP_NOQ1.PRECL.mon.0601-3012.nc"
ntptname = "F2000_NOTP_NOQ1.TS.mon.0601-3012.nc"
ntpuname = "F2000_NOTP_NOQ1.U.mon.0601-3012.nc"
ntpvname = "F2000_NOTP_NOQ1.V.mon.0601-3012.nc"

fntpp1 = addfile(dirintp+ntpp1name, "r")
fntpp2 = addfile(dirintp+ntpp2name, "r")
fntpt  = addfile(dirintp+ntptname, "r")
fntpu  = addfile(dirintp+ntpuname, "r")
fntpv  = addfile(dirintp+ntpvname, "r")

P1_ntp = fntpp1->PRECC(:,{latstr:latend},{lonstr:lonend})
P2_ntp = fntpp2->PRECL(:,{latstr:latend},{lonstr:lonend})
P_ntp = P1_ntp
P_ntp = (P1_ntp + P2_ntp)*1000*86400.
T_ntp = fntpt->TS(:,{latstr:latend},{lonstr:lonend})
U_ntp = fntpu->U(:,{850},{latstr:latend},{lonstr:lonend})
V_ntp = fntpv->V(:,{850},{latstr:latend},{lonstr:lonend})

ntp_p = dim_avg_n_Wrap(P_ntp, 0)
ntp_t = dim_avg_n_Wrap(T_ntp, 0)
ntp_u = dim_avg_n_Wrap(U_ntp, 0)
ntp_v = dim_avg_n_Wrap(V_ntp, 0)

print("end of NOTP_NOQ1")
;=================================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "fig2_CLIM_ANN_OBS_CTRL_THERMAL_DYNAMIC"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4
plot = new(nplot, graphic)
plotv = plot
plotc = plot

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnAddCyclic = False
res@gsnLeftString = ""
res@gsnRightString = ""

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 3
res@mpGeophysicalLineColor = "grey45" ;"grey25"
res@mpFillDrawOrder = "PreDraw"
res@mpOutlineDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2

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
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 10.
res@cnMinLevelValF = 0.
res@cnFillPalette = "CBR_wet"

res@gsnLeftString = "OBS"
plot(0) = gsn_csm_contour_map(wks, obs_p, res)

res@gsnLeftString = "CTRL"
plot(1) = gsn_csm_contour_map(wks, ctl_p, res)

res@gsnLeftString = "NOQ1"
plot(2) = gsn_csm_contour_map(wks, nq1_p, res)

res@gsnLeftString = "NOTP_NOQ1"
plot(3) = gsn_csm_contour_map(wks, ntp_p, res)

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
  vcres@vcRefAnnoOrthogonalPosF =  -0.65 ;-1.0             ; move ref vector up
  vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag

  plotv(0) = gsn_csm_vector(wks, obs_u, obs_v, vcres)

  plotv(1) = gsn_csm_vector(wks, ctl_u, ctl_v, vcres)
  
  plotv(2) = gsn_csm_vector(wks, nq1_u, nq1_v, vcres)
  
  plotv(3) = gsn_csm_vector(wks, ntp_u, ntp_v, vcres)
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
    ;cnres@cnLinePalette = "CBR_coldhot"

    cnres@cnLevelSelectionMode = "ManualLevels"
    cnres@cnMinLevelValF = 270
    cnres@cnMaxLevelValF = 300
    cnres@cnLevelSpacingF = 5

    plotc(0) = gsn_csm_contour(wks, obs_t, cnres)

    plotc(1) = gsn_csm_contour(wks, ctl_t, cnres)

    plotc(2) = gsn_csm_contour(wks, nq1_t, cnres)

    plotc(3) = gsn_csm_contour(wks, ntp_t, cnres)

;----------------------------------------------------

overlay(plot(0), plotv(0))
overlay(plot(0), plotc(0))

overlay(plot(1), plotv(1))
overlay(plot(1), plotc(1))

overlay(plot(2), plotv(2))
overlay(plot(2), plotc(2))

overlay(plot(3), plotv(3))
overlay(plot(3), plotc(3))


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/2,2/), pres)


end