
begin
yrStrt = 1979
yrLast = 2015
year = ispan(yrStrt, yrLast, 1)
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

indexfili = "./TTT_TP_dtrend_standardize_MAM_SON_ERAI.txt"
index = asciiread(indexfili, (/37,2/), "float")
index!0 = "year"
index&year = year

index_spr = index(:,0)
index_aut = index(:,1)
printVarSummary(index_spr)

indhig1 = ind(index_spr.ge.0.8)
indlow1 = ind(index_spr.le.-0.8)
indhig2 = ind(index_aut.ge.0.8)
indlow2 = ind(index_aut.le.-0.8)

yearhig1 = year(indhig1)
yearlow1 = year(indlow1)
yearhig2 = year(indhig2)
yearlow2 = year(indlow2)
print(yearhig1)
print(yearlow1)
exit
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
 lonend = 160
 latstr = -10
 latend = 50
 season = (/"MAM","SON"/)
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
U_SEA = reshape(U, (/37,12,nlat,nlon/)) ;month_to_seasonN(U, season)
V_SEA = reshape(V, (/37,12,nlat,nlon/)) ;month_to_seasonN(V, season)
U_SEA!0 = "time"
U_SEA&time = year
copy_VarCoords(U_SEA, V_SEA)
printVarSummary(U_SEA)
U_clm = dim_avg_n_Wrap(U_SEA, 0)
V_clm = dim_avg_n_Wrap(V_SEA, 0)
copy_VarCoords(U(0:11,:,:), U_clm)
copy_VarCoords(V(0:11,:,:), V_clm)
printVarSummary(U_clm)

filiz = "/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc"
fz = addfile(filiz, "r")
Zx = fz->precip(0:443,{latstr:latend},{lonstr:lonend})
Z = linint2_Wrap(Zx&lon, Zx&lat, Zx, False, U&longitude, U&latitude, 0)
Z_SEA = reshape(Z, (/37,12,nlat,nlon/)) ;month_to_seasonN(Z, season)
copy_VarCoords(U_SEA, Z_SEA)
Z_clm = dim_avg_n_Wrap(Z_SEA, 0)
copy_VarCoords(Z(0:11,:,:), Z_clm)
;==================================================================

U_hig_spr = U_SEA(indhig1,2:5,:,:)
U_low_spr = U_SEA(indlow1,2:5,:,:)
U_hig_aut = U_SEA(indhig2,7:10,:,:)
U_low_aut = U_SEA(indlow2,7:10,:,:)
printVarSummary(U_hig_spr)
V_hig_spr = V_SEA(indhig1,2:5,:,:)
V_low_spr = V_SEA(indlow1,2:5,:,:)
V_hig_aut = V_SEA(indhig2,7:10,:,:)
V_low_aut = V_SEA(indlow2,7:10,:,:)

Z_hig_spr = Z_SEA(indhig1,2:5,:,:)
Z_low_spr = Z_SEA(indlow1,2:5,:,:)
Z_hig_aut = Z_SEA(indhig2,7:10,:,:)
Z_low_aut = Z_SEA(indlow2,7:10,:,:)

U_spr = dim_avg_n_Wrap(U_hig_spr, 0) - dim_avg_n_Wrap(U_low_spr, 0)
copy_VarCoords(U(0:3,:,:), U_spr)
V_spr = dim_avg_n_Wrap(V_hig_spr, 0) - dim_avg_n_Wrap(V_low_spr, 0)
copy_VarCoords(V(0:3,:,:), V_spr)
Z_spr = dim_avg_n_Wrap(Z_hig_spr, 0) - dim_avg_n_Wrap(Z_low_spr, 0)
copy_VarCoords(Z(0:3,:,:), Z_spr)

U_aut = dim_avg_n_Wrap(U_hig_aut, 0) - dim_avg_n_Wrap(U_low_aut, 0)
copy_VarCoords(U(0:3,:,:), U_aut)
V_aut = dim_avg_n_Wrap(V_hig_aut, 0) - dim_avg_n_Wrap(V_low_aut, 0)
copy_VarCoords(V(0:3,:,:), V_aut)
Z_aut = dim_avg_n_Wrap(Z_hig_aut, 0) - dim_avg_n_Wrap(Z_low_aut, 0)
copy_VarCoords(Z(0:3,:,:), Z_aut)

printMinMax(U_spr, True)
printMinMax(V_spr, True)
printMinMax(Z_spr, True)
printMinMax(U_clm, True)
printMinMax(V_clm, True)
printMinMax(Z_clm, True)
;==================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/seasonal/"
figname = "draw_Clim_Comp_UV850_Precip_MAM_SON"

wks = gsn_open_wks("pdf", figpath+figname)
plotv = new(8, graphic)
plotz = plotv
plot_topo = plotv

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

;---------------------画到同一张图------------------------------------------------
    resP                     = True         ; modify the panel plot
    resP@gsnMaximize         = True         ; large format
    resP@gsnPanelLabelBar    = True         ; add common colorbar
    resP@lbLabelAutoStride   = True         ; auto stride on labels
    resP@lbLabelFontHeightF = 0.012
    resP@pmLabelBarOrthogonalPosF     = -0.005                  ; 调色标和图的间距
;--------------------------------------------------------------------------------  

vcres@vcRefMagnitudeF         = 10.0  ; define vector ref mag
res2 = res 
res2@cnLevelSelectionMode = "ManualLevels"
res2@cnLevelSpacingF = 0.5
res2@cnMaxLevelValF = 5.
res2@cnMinLevelValF = -5.
res2@cnFillPalette = "CBR_drywet"
vcres2 = vcres
vcres2@vcRefMagnitudeF         = 1.0 
res@gsnRightString = "Clim"
res2@gsnRightString = "Comp"

do i = 0, 7
  plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
end do
do i = 0, 3
    res@gsnLeftString = Month(i+2)
    plotz(i) = gsn_csm_contour_map_ce(wks, Z_clm(i+2,:,:), res)
    plotv(i) = gsn_csm_vector(wks, U_clm(i+2,:,:), V_clm(i+2,:,:), vcres)
    overlay(plotz(i), plotv(i))
    overlay(plotz(i), plot_topo(i))
    res2@gsnLeftString = Month(i+2)
    plotz(i+4) = gsn_csm_contour_map_ce(wks, Z_spr(i,:,:), res2) ;
    plotv(i+4) = gsn_csm_vector(wks, U_spr(i,:,:), V_spr(i,:,:), vcres2)
    overlay(plotz(i+4), plotv(i+4))
    overlay(plotz(i+4), plot_topo(i+4))
end do

resP@gsnPanelMainString = "Spring"
gsn_panel(wks, plotz(4:7), (/1,4/), resP)
; gsn_panel(wks, plotz, (/2,4/), resP)

do i = 0, 7
  plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
end do
do i = 0, 3
    res@gsnLeftString = Month(i+7)
    plotz(i) = gsn_csm_contour_map_ce(wks, Z_clm(i+7,:,:), res)
    plotv(i) = gsn_csm_vector(wks, U_clm(i+7,:,:), V_clm(i+7,:,:), vcres)
    overlay(plotz(i), plotv(i))
    overlay(plotz(i), plot_topo(i))
    res2@gsnLeftString = Month(i+7)
    plotz(i+4) = gsn_csm_contour_map_ce(wks, Z_aut(i,:,:), res2) ;
    plotv(i+4) = gsn_csm_vector(wks, U_aut(i,:,:), V_aut(i,:,:), vcres2)
    overlay(plotz(i+4), plotv(i+4))
    overlay(plotz(i+4), plot_topo(i+4))
end do
resP@gsnPanelMainString = "Autumn"
gsn_panel(wks, plotz(4:7), (/1,4/), resP)
; gsn_panel(wks, plotz, (/2,4/), resP)
end