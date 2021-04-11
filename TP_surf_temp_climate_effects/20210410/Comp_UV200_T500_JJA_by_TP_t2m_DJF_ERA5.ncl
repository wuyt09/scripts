
begin

yrStrt = 1979
yrLast = 2018
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = (/"JJA"/)

fx = "time_series_TP_t2m_DJF_MAM_ERA5.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(:,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))

; hig_ind = hig_ind+1
; low_ind = low_ind+1

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
filiu = "u_mon_1x1_1979_2019.nc" ;"uwind.monthly.197901-201812.nc"
filiv = "v_mon_1x1_1979_2019.nc" ;"vwind.monthly.197901-201812.nc"
;filiz = "hgt_mon_1x1_1979_2019.nc" ;"zg.monthly.197901-201812.nc"
filiz = "t_mon_1x1_1979_2019.nc"

diri2 = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
filip = "tp_mon_1x1_1979_2019.nc"

fu = addfile(diri+filiu, "r")
fv = addfile(diri+filiv, "r")
fz = addfile(diri+filiz, "r")
fp = addfile(diri2+filip, "r")

Time = fp->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
prec = short2flt(fp->tp(iYYYY,::-1,:))
prec = prec*1000
prec@units = "mm/day"
printVarSummary(prec)

u = short2flt(fu->u(iYYYY,{200},::-1,:))
v = short2flt(fv->v(iYYYY,{200},::-1,:))

z = short2flt(fz->t(iYYYY,{500},::-1,:))
; z = z/9.8
; z@units = "gpm"
printVarSummary(u)

u_sea = month_to_season(u, season)
v_sea = month_to_season(v, season)
z_sea = month_to_season(z, season)
prec_sea = month_to_season(prec, season)
printVarSummary(u_sea)

u_sea_hig = u_sea(hig_ind,:,:)
u_sea_low = u_sea(low_ind,:,:)
v_sea_hig = v_sea(hig_ind,:,:)
v_sea_low = v_sea(low_ind,:,:)
z_sea_hig = z_sea(hig_ind,:,:)
z_sea_low = z_sea(low_ind,:,:)
prec_sea_hig = prec_sea(hig_ind,:,:)
prec_sea_low = prec_sea(low_ind,:,:)


;=================================ttest============================================
usea_hig_ave = dim_avg_n_Wrap(u_sea_hig, 0)
usea_low_ave = dim_avg_n_Wrap(u_sea_low, 0)
vsea_hig_ave = dim_avg_n_Wrap(v_sea_hig, 0)
vsea_low_ave = dim_avg_n_Wrap(v_sea_low, 0)
zsea_hig_ave = dim_avg_n_Wrap(z_sea_hig, 0)
zsea_low_ave = dim_avg_n_Wrap(z_sea_low, 0)
precsea_hig_ave = dim_avg_n_Wrap(prec_sea_hig, 0)
precsea_low_ave = dim_avg_n_Wrap(prec_sea_low, 0)

usea_hig_var = dim_variance_n_Wrap(u_sea_hig, 0)
usea_low_var = dim_variance_n_Wrap(u_sea_low, 0)
vsea_hig_var = dim_variance_n_Wrap(v_sea_hig, 0)
vsea_low_var = dim_variance_n_Wrap(v_sea_low, 0)
zsea_hig_var = dim_variance_n_Wrap(z_sea_hig, 0)
zsea_low_var = dim_variance_n_Wrap(z_sea_low, 0)
precsea_hig_var = dim_variance_n_Wrap(prec_sea_hig, 0)
precsea_low_var = dim_variance_n_Wrap(prec_sea_low, 0)

prob_u = usea_hig_ave
prob_v = vsea_hig_ave
prob_z = zsea_hig_ave
prob_prec = precsea_hig_ave

prob_u = ttest(usea_hig_ave, usea_hig_var, dimsizes(hig_ind), usea_low_ave, usea_low_var, dimsizes(low_ind), False, 0)
prob_v = ttest(vsea_hig_ave, vsea_hig_var, dimsizes(hig_ind), vsea_low_ave, vsea_low_var, dimsizes(low_ind), False, 0)
prob_z = ttest(zsea_hig_ave, zsea_hig_var, dimsizes(hig_ind), zsea_low_ave, zsea_low_var, dimsizes(low_ind), False, 0)
prob_prec = ttest(precsea_hig_ave, precsea_hig_var, dimsizes(hig_ind), precsea_low_ave, precsea_low_var, dimsizes(low_ind), False, 0)
printVarSummary(prob_u)

sig = 0.1
usea_diff = usea_hig_ave
usea_diff = usea_hig_ave - usea_low_ave
vsea_diff = vsea_hig_ave
vsea_diff = vsea_hig_ave - vsea_low_ave
zsea_diff = zsea_hig_ave
zsea_diff = zsea_hig_ave - zsea_low_ave
prec_diff = precsea_hig_ave
prec_diff = precsea_hig_ave - precsea_low_ave

; usea_diff = mask(usea_diff, prob_u.le.sig .and. prob_v.le.sig, True)
; vsea_diff = mask(vsea_diff, prob_u.le.sig .and. prob_v.le.sig, True)
; zsea_diff = mask(zsea_diff, prob_z.le.sig, True)

;===================================================================
; plot
;===================================================================
lonstr = 40
lonend = 180
latstr = 10
latend = 70

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210410/"
figname = "Comp_UV200_T500_JJA_by_TP_t2m_DJF_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)
plotv = plot
plotc = plot
plot_topo = plot
; plot_prob = new(1, graphic)

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey40"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""



res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 4
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -0.7          ; move ref vector up

  vcres@vcGlyphStyle            = "LineArrow"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.06            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
  vcres@gsnLeftString = ""
  vcres@gsnRightString = ""
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
    ; cnres@cnLineLabelInterval = 1
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""

;=====================================================

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 2.
res@cnMaxLevelValF = 270
res@cnMinLevelValF = 250
res@cnFillPalette = "BlueWhiteOrangeRed"
res@pmLabelBarOrthogonalPosF = 0.12

vcres@vcRefMagnitudeF         = 20.0 
; res@gsnRightString = "Comp"

cnres@cnLevelSelectionMode = "ExplicitLevels"
cnres@cnLevels = (/5883/);12520,12560


  res@gsnLeftString = "(a) High DJF TP t2m"
  plot(0) = gsn_csm_contour_map_ce(wks, zsea_hig_ave, res)
  plotv(0) = gsn_csm_vector(wks, usea_hig_ave, vsea_hig_ave, vcres)
  plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(0), plotv(0))
  overlay(plot(0), plot_topo(0))

  res@gsnLeftString = "(b) Low DJF TP t2m"
  plot(1) = gsn_csm_contour_map_ce(wks, zsea_low_ave, res)
  plotv(1) = gsn_csm_vector(wks, usea_low_ave, vsea_low_ave, vcres)
  plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(1), plotv(1))
  overlay(plot(1), plot_topo(1))


    res@cnLevelSelectionMode = "ManualLevels"
    res@cnLevelSpacingF = 0.2
    res@cnMaxLevelValF = 1
    res@cnMinLevelValF = -1
    res@cnFillPalette = "BlueWhiteOrangeRed"

    vcres@vcRefMagnitudeF         = 3.0 

  res@gsnLeftString = "(c) a - b"
  plot(2) = gsn_csm_contour_map_ce(wks, zsea_diff, res)
  plotv(2) = gsn_csm_vector(wks, usea_diff, vsea_diff, vcres)
  plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(2), plotv(2))
  overlay(plot(2), plot_topo(2))



;---------------------画到同一张图------------------------------------------------
    resP                     = True         ; modify the panel plot
    resP@gsnMaximize         = True         ; large format
    ; resP@gsnPanelLabelBar    = True         ; add common colorbar
    resP@lbLabelAutoStride   = True         ; auto stride on labels
    resP@lbLabelFontHeightF = 0.012
    resP@pmLabelBarOrthogonalPosF     = -0.005                  ; 调色标和图的间距
;--------------------------------------------------------------------------------  

resP@gsnPanelMainString = season+" mean UV200_T500" ;_WPSH
gsn_panel(wks, plot, (/1,3/), resP)

end