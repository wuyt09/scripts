
begin
yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
season = "May" ;"JJA"

; indexfili = "./TTT_TP_dtrend_standardize_Mar_Apr_ERAI.txt"
indexfili = "./TTT_TP_dtrend_standardize_MAM_JJA_SON_DJF_ERAI.txt"
inx = asciiread(indexfili, (/nyear,4/), "float")
inx!0 = "year"
inx&year = year

index = inx(:,1)
printVarSummary(index)

indhig = ind(index.ge.0.8)
indlow = ind(index.le.-0.8)
; indhig = ind(index.ge.1)
; indlow = ind(index.le.-1)
yearhig = year(indhig)
yearlow = year(indlow)
print(yearlow)
print(yearhig)
; exit
numhig = dimsizes(yearhig)
numlow = dimsizes(yearlow)

;=================================================================
filitp = "~/wuyt/data/topo.sfc.1x1.nc"
ftp = addfile(filitp, "r")
topo = ftp->topo(::-1,:)
printVarSummary(topo)
;=================================================================

;==================================================================
; Composite
;==================================================================
 lonstr = 50
 lonend = 180
 latstr = -10
 latend = 50
 LEV = 850

diri = "/home/yangsong3/data-observation/ERA-interim/monthly/"
filiu = "uwnd.mon.mean.nc"
filiv = "vwnd.mon.mean.nc"
fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")
TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

U = short2flt(fu->u(iYYYY,{LEV},::-1,:)) ;{latstr:latend},{lonstr:lonend}
V = short2flt(fv->v(iYYYY,{LEV},::-1,:)) ;{latstr:latend},{lonstr:lonend}

filip = "Surface_pressure.mon.mean.nc"
fp = addfile(diri+filip, "r")
PS = short2flt(fp->sp(iYYYY,::-1,:)) ;{latstr:latend},{lonstr:lonend}
U = mask(U, PS.ge.85000, True)
V = mask(V, PS.ge.85000, True)
nlat = dimsizes(U&latitude)
nlon = dimsizes(U&longitude)
U_SEA = U(4::12,:,:) ;month_to_season(U, season)
V_SEA = V(4::12,:,:) ;month_to_season(V, season)
copy_VarCoords(U(0,:,:), U_SEA(0,:,:))
U_SEA!0 = "time"
U_SEA&time = year
copy_VarCoords(U_SEA, V_SEA)
printVarSummary(U_SEA)
U_clm = dim_avg_n_Wrap(U_SEA, 0)
V_clm = dim_avg_n_Wrap(V_SEA, 0)
copy_VarCoords(U(0,:,:), U_clm)
copy_VarCoords(V(0,:,:), V_clm)
printVarSummary(U_clm)
; printVarSummary(U)
; printVarSummary(V)
; printVarSummary(PS)
filiz = "/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc"
fz = addfile(filiz, "r")
Zx = fz->precip(0:443,:,:) ;{latstr:latend},{lonstr:lonend}
printVarSummary(Zx)

Z = linint2_Wrap(Zx&lon, Zx&lat, Zx, False, U&longitude, U&latitude, 0)
Z_SEA = Z(4::12,:,:) ;month_to_season(Z, season)
copy_VarCoords(U_SEA, Z_SEA)
Z_clm = dim_avg_n_Wrap(Z_SEA, 0)
copy_VarCoords(Z(0,:,:), Z_clm)
;==================================================================

;==================================================================

U_hig = U_SEA(indhig,:,:)
U_low = U_SEA(indlow,:,:)
V_hig = V_SEA(indhig,:,:)
V_low = V_SEA(indlow,:,:)
Z_hig = Z_SEA(indhig,:,:)
Z_low = Z_SEA(indlow,:,:)
printVarSummary(U_hig)

Z_higx = dim_avg_n_Wrap(Z_hig, 0)
Z_lowx = dim_avg_n_Wrap(Z_low, 0)

;------------ttest--------------------------------------
uhig_var = dim_variance_n_Wrap(U_hig, 0)
ulow_var = dim_variance_n_Wrap(U_low, 0)
vhig_var = dim_variance_n_Wrap(V_hig, 0)
vlow_var = dim_variance_n_Wrap(V_low, 0)
zhig_var = dim_variance_n_Wrap(Z_hig, 0)
zlow_var = dim_variance_n_Wrap(Z_low, 0)

uhig_ave = dim_avg_n_Wrap(U_hig, 0)
ulow_ave = dim_avg_n_Wrap(U_low, 0)
vhig_ave = dim_avg_n_Wrap(V_hig, 0)
vlow_ave = dim_avg_n_Wrap(V_low, 0)
zhig_ave = dim_avg_n_Wrap(Z_hig, 0)
zlow_ave = dim_avg_n_Wrap(Z_low, 0)


probu = ttest(uhig_ave, uhig_var, numhig, ulow_ave, ulow_var, numlow, False, False)
probv = ttest(vhig_ave, vhig_var, numhig, vlow_ave, vlow_var, numlow, False, False)
probz = ttest(zhig_ave, zhig_var, numhig, zlow_ave, zlow_var, numlow, False, False)

u_diff = uhig_ave - ulow_ave
copy_VarCoords(uhig_ave, u_diff)
v_diff = vhig_ave - vlow_ave
copy_VarCoords(vhig_ave, v_diff)
z_diff = zhig_ave - zlow_ave
copy_VarCoords(zhig_ave, z_diff)

sig = 0.05
; u_diff = mask(u_diff, (probu.lt.sig).or.(probv.lt.sig), True)
; v_diff = mask(v_diff, (probu.lt.sig).or.(probv.lt.sig), True)
; z_diff = mask(z_diff, probz.lt.sig, True)
;==================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/obs/"
figname = "Comp_"+season+"_UV850_Precip"

wks = gsn_open_wks("pdf", figpath+figname)
plotv = new(2, graphic)
plotz = plotv
plot_topo = plotv

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "grey45"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res = True
res@gsnDraw  = False
res@gsnFrame = False
; res@gsnAddCyclic = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 3
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"
res@mpFillOn = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillOn = True
res@cnFillPalette = "CBR_wet"
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 20.
res@cnMinLevelValF = 0.
; res@lbLabelBarOn = False
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0    ; -0.35          ; move ref vector up

  vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.045            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------

;---------------------画到同一张图------------------------------------------------
    resP                     = True         ; modify the panel plot
    resP@gsnMaximize         = True         ; large format
    resP@gsnPanelLabelBar    = True         ; add common colorbar
    resP@lbLabelAutoStride   = True         ; auto stride on labels
    resP@lbLabelFontHeightF = 0.012
    resP@pmLabelBarOrthogonalPosF     = -0.005                  ; 调色标和图的间距
;--------------------------------------------------------------------------------  

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.2
res@cnMaxLevelValF = 2.
res@cnMinLevelValF = -2.
res@cnFillPalette = "CBR_drywet"
res@pmLabelBarOrthogonalPosF = 0.1

vcres@vcRefMagnitudeF         = 1.0 
res@gsnRightString = "Comp"


res@gsnLeftString = "UV850_Precip"
res@gsnRightString = season
plotz(0) = gsn_csm_contour_map_ce(wks, z_diff, res) ;
plotv(0) = gsn_csm_vector(wks, u_diff, v_diff, vcres)
overlay(plotz(0), plotv(0))
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plotz(0), plot_topo(0))
draw(plotz(0))

; resP@gsnPanelMainString = "UV850_Precip"
; ;gsn_panel(wks, plotz(4:7), (/1,4/), resP)
; gsn_panel(wks, plotz, (/2,3/), resP)

end