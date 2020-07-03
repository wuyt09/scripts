
begin

yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
;season = (/"MAM","SON"/)
season = (/"MAM"/) ;,"JJA","SON"

;================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;================================================================

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
fili = "air.mon.mean.nc"

f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f->t(iYYYY,:,::-1,:))
printVarSummary(T)
printMinMax(T, True)

filip = "Surface_pressure.mon.mean.nc"
fp = addfile(diri+filip, "r")
PS = short2flt(fp->sp(iYYYY,::-1,:))

print(T&latitude)
nlat = dimsizes(T&latitude)
nlon = dimsizes(T&longitude)
;===========integration========================================
levs = 600
levt = 300

lev = T&level
plev = lev*100.
psfc = PS
ptop = min(plev)

dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dP!1 = "lev"
dP&lev = lev
printVarSummary(dP)

Tdp = T*dP
copy_VarCoords(T, Tdp)
TTT = dim_sum_n_Wrap(Tdp(:,{levs:levt},:,:), 1)/dim_sum_n_Wrap(dP(:,{levs:levt},:,:),1)
copy_VarCoords(T(:,0,:,:), TTT)
printVarSummary(TTT)
printMinMax(TTT, True)

;===============================================================

; Tclm = clmMonTLL(TTT)
; Tano = calcMonAnomTLL(TTT, Tclm)
; ;print(Tano(20,{25:40},{75:105}))
; printVarSummary(Tano)


; t = new((/3,nyear,nlat,nlon/), float)
; t_Mar = TTT(2::12,:,:)
; t_Apr = TTT(3::12,:,:)
; t_MarApr = (t_Mar+t_Apr)/2
; t(0,:,:,:) = t_Mar
; t(1,:,:,:) = t_Apr
; t(2,:,:,:) = t_MarApr
; copy_VarCoords(TTT(0,:,:), t(0,0,:,:))
t = month_to_season(TTT, season)
printVarSummary(t)
t = mask(t, topo.ge.1500, True)

ttt = wgt_areaave_Wrap(t(:,{25:40},{75:105}), 1, 1, 0)
printVarSummary(ttt)

; t_dtr = dtrend_msg_n(year, ttt, False, False, 1)
; t_std = dim_standardize_n_Wrap(t_dtr, 1, 1)

;===================================================================
; ERSST v5 data
;===================================================================
dirisst = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
filisst = "sst.mnmean.v5.nc"
fsst = addfile(dirisst+filisst, "r")

TIME := fsst->time
YYYY := cd_calendar(TIME, -1)/100
iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

Time = fsst->time(iYYYY)
time = ut_calendar(Time, -1)
lat = fsst->lat(::-1)
lon = fsst->lon

sst = fsst->sst(iYYYY,::-1,:)
sst&time = time
printVarSummary(sst)
; sst = dtrend_msg_n(sst&time, sst, False, False, 0)
SST = month_to_season(sst, season)
SST&time = ispan(yrStrt, yrLast, 1)
printVarSummary(SST)


;===================================================================
; regression / correlation 
;===================================================================

rc = regCoef_n(ttt, SST, 0, 0)
copy_VarCoords(sst(0,:,:), rc)
printVarSummary(rc)

cc = escorc_n(ttt, SST, 0, 0)
copy_VarCoords(sst(0,:,:), cc)
printVarSummary(cc)

Nr = dimsizes(SST&time)
prob = cc
prob = rtest(cc, Nr, 0)
printVarSummary(prob)
prob = mask(prob, cc.eq.0, False)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/obs/"
filiname = "reg_TTT_on_SST"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)
plot_prob = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
res@pmLabelBarOrthogonalPosF = 0.1

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnFillPalette = "CBR_coldhot"
res@gsnLeftString = season+" mean"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.8
res@cnMinLevelValF = -0.8
res@tiMainString = "regression coefficient: ERAI TP Tropo Temp & ERSST.v5 "

sig = 0.95
res@gsnRightString = sig*100+"%"

plot(0) = gsn_csm_contour_map(wks, rc, res)

;------------------------
; draw prob
;------------------------
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

;res3@cnFillPattern = 17 ;(/17/)
res3@cnFillColor  =  "black" ;"green4"
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1 ;1.2

plot_prob(0) = gsn_csm_contour(wks, prob, res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),(1-sig),6)
overlay(plot(0), plot_prob(0))
draw(plot)
end