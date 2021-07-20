
begin

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
season = "JJA"

LEV = 200

lonstr = 0
lonend = 180
latstr = 0
latend = 60

;-----------------------------------------------------------
indexfili = "./time_series_TP_TT_Mar_Apr_May_ERA5.txt"
inx = asciiread(indexfili, (/nyear,4/), "float")
inx!0 = "year"
inx&year = year

index = inx(:,3)
printVarSummary(index)

indhig := ind(index.ge.1.)
indlow := ind(index.le.-1.)
yearhig := year(indhig)
yearlow := year(indlow)
print(yearlow)
print(yearhig)
numhig = dimsizes(yearhig)
numlow = dimsizes(yearlow)


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)


;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili1 = "t_mon_1x1_1979_2020.nc.nc"
fili2 = "u_mon_1x1_1979_2020.nc"
fili3 = "v_mon_1x1_1979_2020.nc"
fili4 = "hgt_mon_1x1_1979_2019.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")
f4 = addfile(diri+fili4, "r")
TIME = f1->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

T = short2flt(f1->t(iYYYY,{500},{latstr:latend},{lonstr:lonend}))
U = short2flt(f2->u(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))
V = short2flt(f3->v(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))
Z = short2flt(f4->z(iYYYY,{LEV},{latstr:latend},{lonstr:lonend}))
Z = Z/9.8
printVarSummary(U)

lat = U&latitude
lon = U&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)

u = new((/6,nyear,nlat,nlon/), float)
u!0 = "time"
title = (/"Mar","Apr","May","Jun","Jul","JJA"/)
u&time = title
v = u
t = u
z = u

do i = 0, 4
    i_strt = i+2
    ; print(i_strt)
    ; print(U(i_strt::12,20,20))
    u(i,:,:,:) = U(i_strt::12,:,:)
    v(i,:,:,:) = V(i_strt::12,:,:)
    z(i,:,:,:) = Z(i_strt::12,:,:)
    t(i,:,:,:) = T(i_strt::12,:,:)
end do
u(5,:,:,:) = month_to_season(U, season)
v(5,:,:,:) = month_to_season(V, season)
z(5,:,:,:) = month_to_season(Z, season)
t(5,:,:,:) = month_to_season(T, season)
printVarSummary(u)

u_hig = u(:,indhig,:,:)
u_low = u(:,indlow,:,:)
v_hig = v(:,indhig,:,:)
v_low = v(:,indlow,:,:)
z_hig = z(:,indhig,:,:)
z_low = z(:,indlow,:,:)
t_hig = t(:,indhig,:,:)
t_low = t(:,indlow,:,:)

uave_hig = dim_avg_n_Wrap(u_hig, 1)
uave_low = dim_avg_n_Wrap(u_low, 1)
vave_hig = dim_avg_n_Wrap(v_hig, 1)
vave_low = dim_avg_n_Wrap(v_low, 1)
zave_hig = dim_avg_n_Wrap(z_hig, 1)
zave_low = dim_avg_n_Wrap(z_low, 1)
tave_hig = dim_avg_n_Wrap(t_hig, 1)
tave_low = dim_avg_n_Wrap(t_low, 1)

uvar_hig = dim_variance_n_Wrap(u_hig, 1)
uvar_low = dim_variance_n_Wrap(u_low, 1)
vvar_hig = dim_variance_n_Wrap(v_hig, 1)
vvar_low = dim_variance_n_Wrap(v_low, 1)
tvar_hig = dim_variance_n_Wrap(t_hig, 1)
tvar_low = dim_variance_n_Wrap(t_low, 1)
zvar_hig = dim_variance_n_Wrap(z_hig, 1)
zvar_low = dim_variance_n_Wrap(z_low, 1)

prob_u = ttest(uave_hig, uvar_hig, numhig, uave_low, uvar_low, numlow, False, False)
prob_v = ttest(vave_hig, vvar_hig, numhig, vave_low, vvar_low, numlow, False, False)
prob_t = ttest(tave_hig, tvar_hig, numhig, tave_low, tvar_low, numlow, False, False)
prob_z = ttest(zave_hig, zvar_hig, numhig, zave_low, zvar_low, numlow, False, False)
printVarSummary(prob_u)

diff_u = uave_hig - uave_low
copy_VarCoords(uave_hig, diff_u)
diff_v = vave_hig - vave_low
copy_VarCoords(vave_hig, diff_v)
diff_t = tave_hig - tave_low
copy_VarCoords(tave_hig, diff_t)
diff_z = zave_hig - zave_low
copy_VarCoords(zave_hig, diff_z)
printVarSummary(diff_u)

sig = 0.05

;===============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210708_obs/TTT/"
figname = "Comp_May_TTT_UV200_T500_SAH_Apr_Aug_ERA5"

wks = gsn_open_wks("pdf", figpath+figname)

nplot = 6
plot = new(nplot, graphic)
plotv = plot
plotc1 = plot
plotc2 = plot
plot_topo = plot

;----------------------------------------------------
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 4
 res_topo@cnLineDashPattern      = 0
 res_topo@cnLineColor            = "green4"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = 2000 ;(/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;----------------------------------------------------

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True
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

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

; res@gsnRightString = ""
;-------------------------------------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -0.35          ; move ref vector up

  vcres@vcGlyphStyle            = "LineArrow" ; "CurlyVector"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.025
  vcres@vcRefLengthF            = 0.05            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
  vcres@gsnLeftString = ""
  vcres@gsnRightString = ""

;-------------------------------------------------------------------------
;===========================================================
cnres = True
cnres@gsnDraw  = False
cnres@gsnFrame = False

cnres@cnFillOn = False
cnres@cnLinesOn = True
cnres@cnLineLabelsOn = True
cnres@cnLineLabelFontColor = "red"
cnres@cnLineColor = "purple"
cnres@cnLineThicknessF = 5
cnres@cnInfoLabelOn = False
;cnres@cnLineDrawOrder = "PostDraw"
cnres@cnLabelDrawOrder = "PostDraw"
;cnres@cnLevelFlag = "LineAndLabel"
cnres@cnLineLabelInterval = 1

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnLevelSpacingF = 50.
cnres@cnMaxLevelValF = 12650.
cnres@cnMinLevelValF = 12450.
cnres@gsnLeftString = ""
cnres@gsnRightString = ""
;===========================================================

res@cnLevelSelectionMode = "ManualLevels"
res@cnFillPalette = "CBR_coldhot"
res@cnLevelSpacingF   = 0.2
res@cnMaxLevelValF    = 2.
res@cnMinLevelValF    = -2.
res@lbLabelBarOn = False

vcres@vcRefMagnitudeF         = 5.0 

do i = 0, 5
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_contour_map(wks, diff_t(i,:,:), res)
    plotv(i) = gsn_csm_vector(wks, diff_u(i,:,:), diff_v(i,:,:), vcres)
    overlay(plot(i), plotv(i))

    ; cnres@cnLineColor = "red"
    ; plotc1(i) = gsn_csm_contour(wks, zave_hig(i,:,:), cnres)
    ; overlay(plot(i), plotc1(i))
    ; cnres@cnLineColor = "blue"
    ; plotc2(i) = gsn_csm_contour(wks, zave_low(i,:,:), cnres)
    ; overlay(plot(i), plotc2(i))

    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/2,3/), pres)
end