
begin

;=================================================================
; Wind850+Pr, Wind500+Hgt500, cross-section WE, cross-section SN
;=================================================================

yrStrt = 1979
yrLast = 2016
season = (/"MAM","JJA","SON","DJF"/)

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
precip = fp->precip(iYYYY,{latstr:latend},{lonstr:lonend})
printVarSummary(precip)

;==================================================================
diritopo = "/home/yangsong3/wuyt/data/"
filitopo = "topo.sfc.1x1.nc"
ftopo = addfile(diritopo+filitopo, "r")
topo = ftopo->topo({latstr:latend},{lonstr:lonend})
printVarSummary(topo)

;==================================================================
diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"
filis = "Surface_pressure.mon.mean.nc"
filit = "2_metre_temperature.mon.mean.nc"

fu = addfile(diri+filiu, "r")
TIMEu = fu->time
YYYYu = cd_calendar(TIMEu, -1)/100
iYYYYu = ind(YYYYu.ge.yrStrt .and. YYYYu.le.yrLast)
uwnd = short2flt(fu->u(iYYYYu,{850},{latstr:latend},{lonstr:lonend}))
printVarSummary(uwnd)

fv = addfile(diri+filiv, "r")
TIMEv = fv->time
YYYYv = cd_calendar(TIMEv, -1)/100
iYYYYv = ind(YYYYv.ge.yrStrt .and. YYYYv.le.yrLast)
vwnd = short2flt(fv->v(iYYYYv,{850},{latstr:latend},{lonstr:lonend}))
printVarSummary(vwnd)

fs = addfile(diri+filis, "r")
TIMEs = fs->time
YYYYs = cd_calendar(TIMEs, -1)/100
iYYYYs = ind(YYYYs.ge.yrStrt .and. YYYYs.le.yrLast-1)
SP = short2flt(fs->sp(iYYYYs,{latstr:latend},{lonstr:lonend}))
printVarSummary(SP)

ft = addfile(diri+filit, "r")
TIMEt = ft->time
YYYYt = cd_calendar(TIMEt, -1)/100
iYYYYt = ind(YYYYt.ge.yrStrt .and. YYYYt.le.yrLast-1)
air = short2flt(ft->t2m(iYYYYt,{latstr:latend},{lonstr:lonend}))
printVarSummary(air)

delete([/filip,filiu,filiv,filis,filit/])
delete([/fp,fu,fv,fs,ft/])
delete([/TIME,TIMEu,TIMEv,TIMEs,TIMEt/])
delete([/YYYY,YYYYu,YYYYv,YYYYs,YYYYt/])
delete([/iYYYY,iYYYYu,iYYYYv,iYYYYs,iYYYYt/])
;==================================================================
PRECIP = month_to_seasonN(precip, season)
U = month_to_seasonN(uwnd, season)
V = month_to_seasonN(vwnd, season)
SP0 = month_to_seasonN(SP, season)
T2M = month_to_seasonN(air, season)


prec = dim_avg_n_Wrap(PRECIP, 1)
u = dim_avg_n_Wrap(U, 1)
v = dim_avg_n_Wrap(V, 1)
sp = dim_avg_n_Wrap(SP0, 1)
t2m = dim_avg_n_Wrap(T2M, 1)
printVarSummary(u)

u = mask(u, sp.ge.85000., True)
v = mask(v, sp.ge.85000., True)

;topo1 = dim_avg_n_Wrap(topo(:,{25:40},:), 1)
;topo2 = dim_avg_n_Wrap(topo(:,:,{75:105}),2)
;tp1 = topo1
;tp2 = topo2
;tp1 = 1013.25*(1-topo1*0.0065/288.15)^5.25145 ; pressure-height
;tp2 = 1013.25*(1-topo2*0.0065/288.15)^5.25145
;printVarSummary(tp1)
delete([/PRECIP,U,V,SP0,T2M/])

;==================================================================
dirim = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"

filiu = "F2000_CTRL.U.mon.0601-3012.nc"
filiv = "F2000_CTRL.V.mon.0601-3012.nc"
filic = "F2000_CTRL.PRECC.mon.0601-3012.nc"
filil = "F2000_CTRL.PRECL.mon.0601-3012.nc"
filip = "F2000_CTRL.PS.mon.0601-3012.nc"
filit = "F2000_CTRL.TS.mon.0601-3012.nc"
;filiz = "USGS-gtopo30_1.9x2.5_remap_c050602.nc"
filiz = "TOPO.nc"

fu = addfile(dirim+filiu, "r")
fv = addfile(dirim+filiv, "r")
fc = addfile(dirim+filic, "r")
fl = addfile(dirim+filil, "r")
fp = addfile(dirim+filip, "r")
ft = addfile(dirim+filit, "r")
fz = addfile(dirim+filiz, "r")

UX = fu->U(:,{850},{latstr:latend},{lonstr:lonend})
VX = fv->V(:,{850},{latstr:latend},{lonstr:lonend})
TX = ft->TS(:,{latstr:latend},{lonstr:lonend})
CX = fc->PRECC(:,{latstr:latend},{lonstr:lonend})
LX = fl->PRECL(:,{latstr:latend},{lonstr:lonend})
HX = fp->PS(:,{latstr:latend},{lonstr:lonend})
PX = CX
PX = (CX+LX)*1000.*86400
printVarSummary(PX)

;zx = fz->PHIS({latstr:latend},{lonstr:lonend})
;zx = zx/10 ;9.8
;zx@units = "meters"
zx = fz->TP
printVarSummary(zx)
zx&lsmlon = fz->LONGXY(0,:)
zx&lsmlat = fz->LATIXY(:,0)
printVarSummary(zx)

UX0 = month_to_seasonN(UX, season)
VX0 = month_to_seasonN(VX, season)
TX0 = month_to_seasonN(TX, season)
PX0 = month_to_seasonN(PX, season)
HX0 = month_to_seasonN(HX, season)

ux = dim_avg_n_Wrap(UX0, 1)
vx = dim_avg_n_Wrap(VX0, 1)
tx = dim_avg_n_Wrap(TX0, 1)
px = dim_avg_n_Wrap(PX0, 1)
hx = dim_avg_n_Wrap(HX0, 1)

ux = mask(ux, hx.ge.85000., True)
vx = mask(vx, hx.ge.85000., True)

;==================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "fig1_CLIM_SEASONAL_OBS_CTRL_UV850_PRECIP_vertical"

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 8
plot = new(nplot, graphic)
plotv = plot
plotc = plot
plot_topo = plot

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
res_topo@cnLineDashPattern = 14
;----------------------------------------------------
;----------------------------------------------------

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
res@pmLabelBarOrthogonalPosF = 0.15
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
;----------------------------------------------------
;----------------------------------------------------

vcres@vcRefAnnoOrthogonalPosF =  -1.0             ; move ref vector up
vcres@vcRefMagnitudeF         = 3.0  ; define vector ref mag
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 2.
res@cnMaxLevelValF = 16.
res@cnMinLevelValF = 0.
res@cnFillPalette = "CBR_wet"

Alefttitle = (/"(a) ","(b) ","(c) ","(d) "/)
Blefttitle = (/"(e) ","(f) ","(g) ","(h) "/)

do i = 0, dimsizes(season)-1
res@gsnLeftString = Alefttitle(i)+"OBS"
res@gsnRightString = season(i)
plot(2*i) = gsn_csm_contour_map(wks, prec(i,:,:), res)
plotv(2*i) = gsn_csm_vector(wks, u(i,:,:), v(i,:,:), vcres)
plotc(2*i) = gsn_csm_contour(wks, t2m(i,:,:), cnres)
plot_topo(2*i) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(2*i), plot_topo(2*i))
overlay(plot(2*i), plotc(2*i))
overlay(plot(2*i), plotv(2*i))

res@gsnLeftString = Blefttitle(i)+"CTRL"
res@gsnRightString = season(i)
plot(2*i+1) = gsn_csm_contour_map(wks, px(i,:,:), res)
plotv(2*i+1) = gsn_csm_vector(wks, ux(i,:,:), vx(i,:,:), vcres)
plotc(2*i+1) = gsn_csm_contour(wks, tx(i,:,:), cnres)
plot_topo(2*i+1) = gsn_csm_contour(wks, zx, res_topo)
overlay(plot(2*i+1), plot_topo(2*i+1))
overlay(plot(2*i+1), plotc(2*i+1))
overlay(plot(2*i+1), plotv(2*i+1))

end do
;----------------------------------------------------
;----------------------------------------------------

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarWidthF = 0.55
gsn_panel(wks, plot, (/4,2/), pres)
end 