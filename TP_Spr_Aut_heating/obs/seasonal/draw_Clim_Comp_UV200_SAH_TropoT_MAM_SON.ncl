
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
print(yearhig2)
print(yearlow2)

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
 LEV = 200

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

filiz = "hgt.mon.mean.nc"
fz = addfile(diri+filiz, "r")
Z = short2flt(fz->z(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))
Z = Z/9.8
Z_SEA = reshape(Z, (/37,12,nlat,nlon/)) ;month_to_seasonN(Z, season)
copy_VarCoords(U_SEA, Z_SEA)
Z_clm = dim_avg_n_Wrap(Z_SEA, 0)
copy_VarCoords(Z(0:11,:,:), Z_clm)

filit = "air.mon.mean.nc"
ft = addfile(diri+filit, "r")
T = short2flt(ft->t(iYYYY,:,{latstr:latend},{lonstr:lonend}))
filip = "Surface_pressure.mon.mean.nc"
fp = addfile(diri+filip, "r")
PS = short2flt(fp->sp(iYYYY,{latstr:latend},{lonstr:lonend}))
lev = T&level
plev = lev*100
ptop = min(plev)
psfc = PS
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
printVarSummary(dP)
dP!1 = "lev"
dP&lev = lev
TdP = T
TdP = T*dP
t = dim_sum_n_Wrap(TdP(:,{500:200},:,:), 1)/dim_sum_n_Wrap(dP(:,{500:200},:,:),1)
copy_VarCoords(T(:,0,:,:), t)
printMinMax(t, True)

t_SEA = reshape(t, (/37,12,nlat,nlon/)) ;month_to_seasonN(Z, season)
copy_VarCoords(U_SEA, t_SEA)
t_clm = dim_avg_n_Wrap(t_SEA, 0)
copy_VarCoords(t(0:11,:,:), t_clm)

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

Z_hig_sprx = dim_avg_n_Wrap(Z_hig_spr, 0)
Z_low_sprx = dim_avg_n_Wrap(Z_low_spr, 0)
Z_hig_autx = dim_avg_n_Wrap(Z_hig_aut, 0)
Z_low_autx = dim_avg_n_Wrap(Z_low_aut, 0)

t_hig_spr = t_SEA(indhig1,2:5,:,:)
t_low_spr = t_SEA(indlow1,2:5,:,:)
t_hig_aut = t_SEA(indhig2,7:10,:,:)
t_low_aut = t_SEA(indlow2,7:10,:,:)

U_spr = dim_avg_n_Wrap(U_hig_spr, 0) - dim_avg_n_Wrap(U_low_spr, 0)
copy_VarCoords(U(0:3,:,:), U_spr)
V_spr = dim_avg_n_Wrap(V_hig_spr, 0) - dim_avg_n_Wrap(V_low_spr, 0)
copy_VarCoords(V(0:3,:,:), V_spr)
Z_spr = dim_avg_n_Wrap(Z_hig_spr, 0) - dim_avg_n_Wrap(Z_low_spr, 0)
copy_VarCoords(Z(0:3,:,:), Z_spr)
t_spr = dim_avg_n_Wrap(t_hig_spr, 0) - dim_avg_n_Wrap(t_low_spr, 0)
copy_VarCoords(t(0:3,:,:), t_spr)

U_aut = dim_avg_n_Wrap(U_hig_aut, 0) - dim_avg_n_Wrap(U_low_aut, 0)
copy_VarCoords(U(0:3,:,:), U_aut)
V_aut = dim_avg_n_Wrap(V_hig_aut, 0) - dim_avg_n_Wrap(V_low_aut, 0)
copy_VarCoords(V(0:3,:,:), V_aut)
Z_aut = dim_avg_n_Wrap(Z_hig_aut, 0) - dim_avg_n_Wrap(Z_low_aut, 0)
copy_VarCoords(Z(0:3,:,:), Z_aut)
t_aut = dim_avg_n_Wrap(t_hig_aut, 0) - dim_avg_n_Wrap(t_low_aut, 0)
copy_VarCoords(t(0:3,:,:), t_aut)

printMinMax(U_spr, True)
printMinMax(V_spr, True)
printMinMax(Z_spr, True)
printMinMax(t_spr, True)
printMinMax(U_clm, True)
printMinMax(V_clm, True)
printMinMax(Z_clm, True)
printMinMax(t_clm, True)

;==================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_Spr_Aut_heating/obs/seasonal/"
figname = "draw_Clim_Comp_UV200_SAH_TropoT_MAM_SON"

wks = gsn_open_wks("pdf", figpath+figname)
plotc = new(8, graphic)
plotz = plotc
plotv = plotc
plotz1 = plotc
plotz2 = plotc

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
res@cnFillPalette = "sunshine_9lev"
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 255.
res@cnMinLevelValF = 249.
res@lbLabelBarOn = False
;------------------------------------------
resZ = True
resZ@gsnDraw  = False
resZ@gsnFrame = False

resZ@cnFillOn = False
resZ@cnLinesOn = True
resZ@cnLineLabelsOn = True
resZ@cnLineLabelFontColor = "red"
resZ@cnLineColor = "green4"
resZ@cnLineThicknessF = 4
resZ@cnInfoLabelOn = True
resZ@cnInfoLabelOrthogonalPosF = -0.12
resZ@cnInfoLabelFontHeightF = 0.02
;res@cnLineDrawOrder = "PostDraw"
resZ@cnLabelDrawOrder = "PostDraw"
;res@cnLevelFlag = "LineAndLabel"
resZ@cnLineLabelInterval = 1
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
vcres@vcRefMagnitudeF         = 30.0  ; define vector ref mag
resZ@lbLabelBarOn = False
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 10.
; res@cnMaxLevelValF = 12550.
; res@cnMinLevelValF = 12400.
resZ@cnLevelSelectionMode = "ExplicitLevels" ; "ManualLevels"      ;use explicit levels
;resZ@cnLevels             = (/12500,12520,12540,12560,12575/) ;(/1500/) ;(/16750/) ;
resZ@cnLevels =(/12460,12480,12500,12520,12540/)

;---------------------画到同一张图------------------------------------------------
    resP                     = True         ; modify the panel plot
    resP@gsnDraw             = True
    resP@gsnFrame            = True
    resP@gsnMaximize         = True         ; large format
    resP@lbLabelAutoStride   = True         ; auto stride on labels
    resP@lbLabelFontHeightF = 0.012
    resP@pmLabelBarOrthogonalPosF     = -0.005  ; 调色标和图的间距
    resP@gsnPanelLabelBar = True
;--------------------------------------------------------------------------------  

res2 = res
res2@cnLevelSelectionMode = "ManualLevels"
res2@cnLevelSpacingF = 0.3
res2@cnMaxLevelValF = 3.
res2@cnMinLevelValF = -3.
res2@cnFillPalette = "CBR_coldhot"

resZ2 = resZ 
resZ2@cnLineThicknessF = 3
; resZ2@cnLevelSelectionMode = "ManualLevels"
; resZ2@cnLevelSpacingF = 30.
; resZ2@cnMaxLevelValF = 150.
; resZ2@cnMinLevelValF = -150.
resZ2@cnLevels =(/12420,12440,12460,12480,12500/)
vcres2 = vcres
vcres2@vcRefMagnitudeF         = 5.0 

res@gsnRightString = "Clim"
res2@gsnRightString = "Comp"
do i = 0, 3
    res@gsnLeftString = Month(i+2)
    plotc(i) = gsn_csm_contour_map_ce(wks, t_clm(i+2,:,:), res)
    plotz(i) = gsn_csm_contour(wks, Z_clm(i+2,:,:), resZ)
    plotv(i) = gsn_csm_vector(wks, U_clm(i+2,:,:), V_clm(i+2,:,:), vcres)
    overlay(plotc(i), plotz(i))
    overlay(plotc(i), plotv(i))
    res2@gsnLeftString = Month(i+2)
    plotc(i+4) = gsn_csm_contour_map_ce(wks, t_spr(i,:,:), res2)
    ;plotz(i+4) = gsn_csm_contour(wks, Z_spr(i,:,:), resZ2) ;
    resZ2@cnLineDashPattern = 0
    resZ2@cnLineColor = "red"
    plotz1(i+4) = gsn_csm_contour(wks, Z_hig_sprx(i,:,:), resZ2)
    resZ2@cnLineDashPattern = 16
    resZ2@cnLineColor = "blue"
    plotz2(i+4) = gsn_csm_contour(wks, Z_low_sprx(i,:,:), resZ2)
    plotv(i+4) = gsn_csm_vector(wks, U_spr(i,:,:), V_spr(i,:,:), vcres2)
    overlay(plotc(i+4), plotz1(i+4))
    overlay(plotc(i+4), plotz2(i+4))
    overlay(plotc(i+4), plotv(i+4))
end do
resP@gsnPanelMainString = "Spring"
gsn_panel(wks, plotc(4:7), (/1,4/), resP)
; gsn_panel(wks, plotc, (/2,4/), resP)
; ;resP@gsnPanelMainString = "Spring"
; pplot1 = gsn_panel_return(wks, plotc(0:3), (/1,4/), resP)
; pplot2 = gsn_panel_return(wks, plotc(4:7), (/1,4/), resP)
; bb = NhlGetBB(pplot1)
; height = max(bb(0,:)-min(bb(:,1)))

; resP@gsnDraw = True
                
; resP@gsnPanelMainString = "Spring"
; resP@gsnPanelTop = 1. ;0.94
; resP@gsnPanelBottom = resP@gsnPanelTop - height
; gsn_panel(wks, plotc(0:3), (/1,4/), resP)

; resP@gsnPanelLabelBar = True
; resP@gsnPanelMainString = ""
; resP@gsnPanelTop = resP@gsnPanelBottom
; resP@gsnPanelBottom = resP@gsnPanelTop - height
; gsn_panel(wks, plotc(4:7), (/1,4/), resP)
; frame(wks)

resZ2@cnLevels = resZ@cnLevels ;(/12500,12520,12540,12560,12580/)

do i = 0, 3
    res@gsnLeftString = Month(i+7)
    plotc(i) = gsn_csm_contour_map_ce(wks, t_clm(i+7,:,:), res)
    plotz(i) = gsn_csm_contour(wks, Z_clm(i+7,:,:), resZ)
    plotv(i) = gsn_csm_vector(wks, U_clm(i+7,:,:), V_clm(i+7,:,:), vcres)
    overlay(plotc(i), plotz(i))
    overlay(plotc(i), plotv(i))
    res2@gsnLeftString = Month(i+7)
    plotc(i+4) = gsn_csm_contour_map_ce(wks, t_aut(i,:,:), res2)
    ;plotz(i+4) = gsn_csm_contour(wks, Z_aut(i,:,:), resZ2) ;
    resZ2@cnLineDashPattern = 0
    resZ2@cnLineColor = "red"
    plotz1(i+4) = gsn_csm_contour(wks, Z_hig_autx(i,:,:), resZ2)
    resZ2@cnLineDashPattern = 16
    resZ2@cnLineColor = "blue"
    plotz2(i+4) = gsn_csm_contour(wks, Z_low_autx(i,:,:), resZ2)
    plotv(i+4) = gsn_csm_vector(wks, U_aut(i,:,:), V_aut(i,:,:), vcres2)
    overlay(plotc(i+4), plotz1(i+4))
    overlay(plotc(i+4), plotz2(i+4))
    overlay(plotc(i+4), plotv(i+4))
end do
resP@gsnPanelMainString = "Autumn"
gsn_panel(wks, plotc(4:7), (/1,4/), resP)
; gsn_panel(wks, plotc, (/2,4/), resP)
end