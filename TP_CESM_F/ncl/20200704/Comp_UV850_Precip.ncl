
begin
yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

indexfili = "./TTT_TP_dtrend_standardize_Mar_Apr_May_ERAI.txt"
inx = asciiread(indexfili, (/nyear,3/), "float")
inx!0 = "year"
inx&year = year

title = (/"Mar","Apr","May"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "Comp_TTT_UV850_Precip"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new((/4,3/), graphic)
plotz1 = plot
plotz2 = plot
plotv = plot
plot_topo = plotv


do ll = 0, 2
index = inx(:,ll)
printVarSummary(index)

indhig := ind(index.ge.1.)
indlow := ind(index.le.-1.)
yearhig := year(indhig)
yearlow := year(indlow)
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

U = short2flt(fu->u(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))
V = short2flt(fv->v(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))

filip = "Surface_pressure.mon.mean.nc"
fp = addfile(diri+filip, "r")
PS = short2flt(fp->sp(iYYYY,{latstr:latend},{lonstr:lonend}))
U = mask(U, PS.ge.85000, True)
V = mask(V, PS.ge.85000, True)
nlat = dimsizes(U&latitude)
nlon = dimsizes(U&longitude)
U_SEA = reshape(U, (/nyear,12,nlat,nlon/)) ;month_to_seasonN(U, season)
V_SEA = reshape(V, (/nyear,12,nlat,nlon/)) ;month_to_seasonN(V, season)
copy_VarCoords(U(0,:,:), U_SEA(0,0,:,:))
U_SEA!0 = "time"
U_SEA&time = year
copy_VarCoords(U_SEA, V_SEA)
printVarSummary(U_SEA)
U_clm = dim_avg_n_Wrap(U_SEA, 0)
V_clm = dim_avg_n_Wrap(V_SEA, 0)
copy_VarCoords(U(0:11,:,:), U_clm)
copy_VarCoords(V(0:11,:,:), V_clm)
printVarSummary(U_clm)

filiz = "hgt.mon.mean.nc"
fz = addfile(diri+filiz, "r")
Z = short2flt(fz->z(iYYYY,{500},{latstr:latend},{lonstr:lonend}))
Z = Z/9.8
Z_SEA = reshape(Z, (/37,12,nlat,nlon/)) ;month_to_seasonN(Z, season)
copy_VarCoords(U_SEA, Z_SEA)
Z_clm = dim_avg_n_Wrap(Z_SEA, 0)
copy_VarCoords(Z(0:11,:,:), Z_clm)


filip = "/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc"
fp = addfile(filip, "r")
Px = fp->precip(0:443,{latstr:latend},{lonstr:lonend})
P = linint2_Wrap(Px&lon, Px&lat, Px, False, U&longitude, U&latitude, 0)
P_SEA = reshape(P, (/37,12,nlat,nlon/)) ;month_to_seasonN(Z, season)
copy_VarCoords(U_SEA, P_SEA)
P_clm = dim_avg_n_Wrap(P_SEA, 0)
copy_VarCoords(P(0:11,:,:), P_clm)
;==================================================================

;==================================================================
monS = 3
monE = 6
U_hig := U_SEA(indhig,monS-1:monE-1,:,:)
U_low := U_SEA(indlow,monS-1:monE-1,:,:)
V_hig := V_SEA(indhig,monS-1:monE-1,:,:)
V_low := V_SEA(indlow,monS-1:monE-1,:,:)
Z_hig := Z_SEA(indhig,monS-1:monE-1,:,:)
Z_low := Z_SEA(indlow,monS-1:monE-1,:,:)
P_hig := P_SEA(indhig,monS-1:monE-1,:,:)
P_low := P_SEA(indlow,monS-1:monE-1,:,:)
printVarSummary(U_hig)

P_higx = dim_avg_n_Wrap(P_hig, 0)
P_lowx = dim_avg_n_Wrap(P_low, 0)

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
phig_ave = dim_avg_n_Wrap(P_hig, 0)
plow_ave = dim_avg_n_Wrap(P_low, 0)


probu = ttest(uhig_ave, uhig_var, numhig, ulow_ave, ulow_var, numlow, False, False)
probv = ttest(vhig_ave, vhig_var, numhig, vlow_ave, vlow_var, numlow, False, False)
; probz = ttest(zhig_ave, zhig_var, numhig, zlow_ave, zlow_var, numlow, False, False)
probp = ttest(zhig_ave, zhig_var, numhig, zlow_ave, zlow_var, numlow, False, False)

u_diff = uhig_ave - ulow_ave
copy_VarCoords(uhig_ave, u_diff)
v_diff = vhig_ave - vlow_ave
copy_VarCoords(vhig_ave, v_diff)
; z_diff = zhig_ave - zlow_ave
; copy_VarCoords(zhig_ave, z_diff)
p_diff = phig_ave - plow_ave
copy_VarCoords(phig_ave, p_diff)

sig = 0.05
; u_diff = mask(u_diff, (probu.lt.sig).or.(probv.lt.sig), True)
; v_diff = mask(v_diff, (probu.lt.sig).or.(probv.lt.sig), True)
; z_diff = mask(z_diff, probz.lt.sig, True)
;==================================================================



 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey45"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnAddCyclic = False

res@mpOutlineOn = True
res@mpGeophysicalLineThicknessF = 2
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
res@lbLabelBarOn = False
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
;=====================================================
    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "purple" ;goldenrod2 ;dodgerblue2
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = False
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""

;=====================================================

;---------------------画到同一张图------------------------------------------------
    resP                     = True         ; modify the panel plot
    resP@gsnMaximize         = True         ; large format
    resP@gsnPanelLabelBar    = True         ; add common colorbar
    resP@lbLabelAutoStride   = True         ; auto stride on labels
    resP@lbLabelFontHeightF = 0.012
    resP@pmLabelBarOrthogonalPosF     = -0.005                  ; 调色标和图的间距
;--------------------------------------------------------------------------------  

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.
res@cnFillPalette = "CBR_drywet"

vcres@vcRefMagnitudeF         = 1.0 
; res@gsnRightString = "Comp"

cnres@cnLevelSelectionMode = "ExplicitLevels"
cnres@cnLevels = (/5885/);12520,12560

do i = 0, monE - monS
  res@gsnLeftString = Month(i+monS-1)
  res@gsnRightString = title(ll)+" TTT Comp"
  plot(i,ll) = gsn_csm_contour_map_ce(wks, p_diff(i,:,:), res) ;
  plotv(i,ll) = gsn_csm_vector(wks, u_diff(i,:,:), v_diff(i,:,:), vcres)
  overlay(plot(i,ll), plotv(i,ll))

  cnres@cnLineColor = "red"
  plotz1(i,ll) = gsn_csm_contour(wks, zhig_ave(i,:,:), cnres)
  overlay(plot(i,ll), plotz1(i,ll))
  cnres@cnLineColor = "blue"
  plotz2(i,ll) = gsn_csm_contour(wks, zlow_ave(i,:,:), cnres)
  overlay(plot(i,ll), plotz2(i,ll))

  plot_topo(i,ll) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(i,ll), plot_topo(i,ll))
end do
end do

plotx = reshape(plot, 4*3)

resP@gsnPanelMainString = "UV850_Precip_WPSH"
;gsn_panel(wks, plotz(4:7), (/1,4/), resP)
gsn_panel(wks, plotx, (/4,3/), resP)

end