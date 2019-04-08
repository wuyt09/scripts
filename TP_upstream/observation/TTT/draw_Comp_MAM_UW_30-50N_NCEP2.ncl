
begin

yrStrt = 1979
yrLast = 2016
year = ispan(yrStrt, yrLast, 1)
ny = dimsizes(year)
Season = (/"MAM"/)
ns = dimsizes(Season)

;---------------------------------------------
path = "./"
name = "TTT_TP_dtrend_standardize_NCEP2.txt"

x = asciiread(path+name, (/ny,ns/), "float")
x@_FillValue = -32767
x!0 = "year"
x&year = year
x!1 = "season"
x&season = Season
printVarSummary(x)

index = x(:,0) ;MAM
indlow = ind(index .le. -0.8)
indhig = ind(index .ge. 0.8)
numhig = dimsizes(indhig)
numlow = dimsizes(indlow)
yearhig = year(indhig)
yearlow = year(indlow)
print(yearhig)
print(yearlow)
;---------------------------------------------

diri = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
filu = "uwnd.mon.mean.nc"
filw = "omega.mon.mean.nc"

fu = addfile(diri+filu, "r")
fw = addfile(diri+filw, "r")

TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

U = short2flt(fu->uwnd(iYYYY,:,{30:50},:))
W = short2flt(fw->omega(iYYYY,:,{30:50},:))
W = W*(-100.)
U = lonFlip(U)
W = lonFlip(W)

printVarSummary(U)

uclm = clmMonTLLL(U)
printVarSummary(uclm)
u0 = calcMonAnomTLLL(U, uclm)
wclm = clmMonTLLL(W)
w0 = calcMonAnomTLLL(W, wclm)
printVarSummary(u0)

;u = dim_avg_n_Wrap(u0(:,:,{30:50},:), 2)
;w = dim_avg_n_Wrap(w0(:,:,{30:50},:), 2)
wgt = U&lat
wgt = 1
u = dim_avg_wgt_n_Wrap(u0(:,:,:,:), wgt, 1, 2)
w = dim_avg_wgt_n_Wrap(w0(:,:,:,:), wgt, 1, 2)
printVarSummary(u)

u_spr = month_to_season(u, "MAM")
w_spr = month_to_season(w, "MAM")
u_spr!0 = "year"
u_spr&year = year
copy_VarCoords(u_spr, w_spr)
printVarSummary(u_spr)

uhig = u_spr(indhig,:,:)
ulow = u_spr(indlow,:,:)
whig = w_spr(indhig,:,:)
wlow = w_spr(indlow,:,:)
printVarSummary(uhig)

uhigvar = dim_variance_n_Wrap(uhig, 0)
ulowvar = dim_variance_n_Wrap(ulow, 0)
uhigave = dim_avg_n_Wrap(uhig, 0)
ulowave = dim_avg_n_Wrap(ulow, 0)
probu = ttest(uhigave, uhigvar, numhig, ulowave, ulowvar, numlow, False, False)

whigvar = dim_variance_n_Wrap(whig, 0)
wlowvar = dim_variance_n_Wrap(wlow, 0)
whigave = dim_avg_n_Wrap(whig, 0)
wlowave = dim_avg_n_Wrap(wlow, 0)
probw = ttest(whigave, whigvar, numhig, wlowave, wlowvar, numlow, False, False)

udiff = uhigave - ulowave
wdiff = whigave - wlowave
copy_VarCoords(u(0,:,:), udiff)
copy_VarCoords(w(0,:,:), wdiff)

;udiff = mask(udiff, (probu.lt.0.05).or.(probw.lt.0.05), True)
;wdiff = mask(wdiff, (probu.lt.0.05).or.(probw.lt.0.05), True)

;delete([/U,W,U0,W0,uclm,wclm,u,w/])
;===================speed of w======================================
whspeed = whigave
whspeed = whigave/(-100.)
wlspeed = wlowave
wlspeed = wlowave/(-100.)
wdspeed = wdiff
wdspeed = wdiff/(-100.)
printVarSummary(whspeed)
;===================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/observation/TTT/"
figname = "draw_Comp_MAM_UW_30-50N_NCEP2"

lonstr = -30
lonend = 150

latS = 25
latN = 40

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 3
plot = new(nplot, graphic)
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@vpHeightF = 0.3        ; Changes the aspect ratio
res@vpWidthF  = 0.55

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
res@cnFillDrawOrder = "PreDraw"

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

;------------------------------------------
 res@trYMinF   = 100
 res@trYMaxF   = 1000
  
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)

  res@vcRefAnnoOrthogonalPosF =  -0.4 ;-1.0             ; move ref vector up
  ;res@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;res@vcMinDistanceF          = 0.017
  res@vcGlyphStyle            = "LineArrow"
  res@vcMinDistanceF          = 0.03
  res@vcRefLengthF            = 0.025            ; define length of vec ref
  
  res@vcLineArrowThicknessF   = 3.0
  ;res@vcLineArrowHeadMinSizeF = 0.008
  ;res@vcLineArrowHeadMaxSizeF = 0.3
  res@vcVectorDrawOrder = "PostDraw" 
 
  res@vcRefAnnoString2On = False
  res@vcRefAnnoFontHeightF = 0.02

;------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.005
res@cnMaxLevelValF = 0.02
res@cnMinLevelValF = -0.02
res@cnFillPalette = "BlueWhiteOrangeRed"
res@vcRefMagnitudeF         = 1.             ; define vector ref mag
res@lbLabelBarOn = False


res@gsnRightString = " "
res@gsnLeftString = "HIGH"
plot(0) = gsn_csm_pres_hgt_vector(wks, whspeed(:,{lonstr:lonend}), uhigave(:,{lonstr:lonend}), whigave(:,{lonstr:lonend}), res)

res@gsnLeftString = "LOW"
plot(1) = gsn_csm_pres_hgt_vector(wks, wlspeed(:,{lonstr:lonend}), ulowave(:,{lonstr:lonend}), wlowave(:,{lonstr:lonend}), res)

res@gsnLeftString = "DIFF"
plot(2) = gsn_csm_pres_hgt_vector(wks, wdspeed(:,{lonstr:lonend}), udiff(:,{lonstr:lonend}), wdiff(:,{lonstr:lonend}), res)


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,1/), pres)


end 
