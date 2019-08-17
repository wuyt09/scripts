
begin

;=================================================================
; Wind850+Pr, Wind500+Hgt500, cross-section WE, cross-section SN
;=================================================================

yrStrt = 1979
yrLast = 2016

lonstr = 40 ;20
lonend = 180 ;220
latstr = -10
latend = 70

;==================================================================
dirip = "/home/yangsong3/data-observation/GPCP/"
filip = "precip.mon.mean.nc"

fp = addfile(dirip+filip, "r")
TIME = fp->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
precip = fp->precip(iYYYY,:,:)
printVarSummary(precip)

;==================================================================
dirim = "/home/yangsong3/wuyt/data/"
filim = "topo.sfc.1x1.nc"
fm = addfile(dirim+filim, "r")
topo = fm->topo
printVarSummary(topo)

dirin = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
filin = "TPmask_1000m_C20190410.nc"
fn = addfile(dirin+filin, "r")
tpmask = fn->tp
printVarSummary(tpmask)

;==================================================================
diri = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"
filiw = "omega.mon.mean.nc"
filit = "air.mon.mean.nc"
filiz = "hgt.mon.mean.nc"

fu = addfile(diri+filiu, "r")
TIMEu = fu->time
YYYYu = cd_calendar(TIMEu, -1)/100
iYYYYu = ind(YYYYu.ge.yrStrt .and. YYYYu.le.yrLast)
uwnd = short2flt(fu->uwnd(iYYYYu,:,:,:))
printVarSummary(uwnd)

fv = addfile(diri+filiv, "r")
TIMEv = fv->time
YYYYv = cd_calendar(TIMEv, -1)/100
iYYYYv = ind(YYYYv.ge.yrStrt .and. YYYYv.le.yrLast)
vwnd = short2flt(fv->vwnd(iYYYYv,:,:,:))
printVarSummary(vwnd)

fw = addfile(diri+filiw, "r")
TIMEw = fw->time
YYYYw = cd_calendar(TIMEw, -1)/100
iYYYYw = ind(YYYYw.ge.yrStrt .and. YYYYw.le.yrLast)
wwnd = short2flt(fw->omega(iYYYYw,:,:,:))
printVarSummary(wwnd)

ft = addfile(diri+filit, "r")
TIMEt = ft->time
YYYYt = cd_calendar(TIMEt, -1)/100
iYYYYt = ind(YYYYt.ge.yrStrt .and. YYYYt.le.yrLast)
T = short2flt(ft->air(iYYYYt,:,:,:))
printVarSummary(T)

fz = addfile(diri+filiz, "r")
TIMEz = fz->time
YYYYz = cd_calendar(TIMEz, -1)/100
iYYYYz = ind(YYYYz.ge.yrStrt .and. YYYYz.le.yrLast)
Z = short2flt(fz->hgt(iYYYYz,{500},:,:))
printVarSummary(Z)

delete([/TIME,TIMEu,TIMEv,TIMEw,TIMEt,TIMEz/])
delete([/YYYY,YYYYu,YYYYv,YYYYw,YYYYt,YYYYz/])
delete([/iYYYY,iYYYYu,iYYYYv,iYYYYw,iYYYYt,iYYYYz/])
;==================================================================

prec = dim_avg_n_Wrap(precip, 0)
u = dim_avg_n_Wrap(uwnd, 0)
v = dim_avg_n_Wrap(vwnd, 0)
w = dim_avg_n_Wrap(wwnd, 0)
t = dim_avg_n_Wrap(T, 0)
z = dim_avg_n_Wrap(Z, 0)

;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "fig1_CLIM_ANN_OBS_TP_EFFECTS_NCEP2"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4
plot = new(nplot, graphic)


;----------------------------------------------------
res_topo = True
res_topo@gsnFrame = False
res_topo@gsnDraw  = False

res_topo@mpMaxLatF = 50
res_topo@mpMinLatF = 20
res_topo@mpMaxLonF = 115
res_topo@mpMinLonF = 55
res_topo@gsnLeftString = ""
res_topo@gsnRightString = ""

res_topo@cnFillOn = False
res_topo@cnLinesOn = True
res_topo@cnLineLabelsOn = True
res_topo@cnInfoLabelOn = False
res_topo@cnLineThicknessF = 2

res_topo@cnLevelSelectionMode = "ExplicitLevels"
res_topo@cnLevels = (/500,1000,1500/)
res_topo@cnLineColor = "green4"
;----------------------------------------------------
optt = True
optt@gsnDraw  = False
optt@gsnFrame = False
optt@gsnLeftString = ""
optt@gsnRightString = ""

optt@cnFillOn = True
optt@cnLinesOn = False
optt@cnLineLabelsOn = False
optt@cnInfoLabelOn = False

optt@cnLevelSelectionMode = "ExplicitLevels"
optt@cnLevels = (/0.9/)
optt@cnFillColor = "black"
optt@lbLabelBarOn = False

opt = True
opt@gsnShadeFillType = "pattern"
opt@gsnShadeHigh = 14
opt@gsnShadeFillScaleF = 1.2
opt@gsnShadeFillDotSizeF = 0.015

;----------------------------------------------------
plot(0) = gsn_csm_contour_map(wks, topo, res_topo)
plot_tp = gsn_csm_contour(wks, tpmask, optt)
plot_tp = gsn_contour_shade(plot_tp, -999., 0.9, opt)
overlay(plot(0), plot_tp)
;----------------------------------------------------

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnLeftString = ""
res@gsnRightString = ""

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 3
res@mpGeophysicalLineColor = "grey45" ;"grey25"
res@mpFillDrawOrder = "PreDraw"
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

;----------------------------------------------------
vcres@vcRefAnnoOrthogonalPosF =  -0.55 ;-1.0             ; move ref vector up
vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
res@lbLabelBarOn = True
res@pmLabelBarOrthogonalPosF = -0.01
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 10.
res@cnMinLevelValF = 0.
res@cnFillPalette = "CBR_wet"

plot(1) = gsn_csm_contour_map(wks, prec, res)
plotv1 = gsn_csm_vector(wks, u({850},:,:), v({850},:,:), vcres)
overlay(plot(1), plotv1)

draw(plot(1))
frame(wks)
end 