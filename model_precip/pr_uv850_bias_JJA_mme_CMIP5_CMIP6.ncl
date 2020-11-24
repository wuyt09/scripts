
begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

season = "JJA"
yrStrt = 1979
yrLast = 2014
yearo = ispan(yrStrt, yrLast, 1)
nyro = dimsizes(yearo)
LEV = 850

diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
TIME = fobs->time
YYYYo = cd_calendar(TIME, -1)/100
iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)
PR_OBS = fobs->precip(iYYYYo,:,:)
pr_obs = month_to_season(PR_OBS, season)
 ; = dim_avg_n_Wrap(PR_OBSx, 0)
printVarSummary(pr_obs)

plat = pr_obs&lat
plon = pr_obs&lon
pr_obs_ave = dim_avg_n_Wrap(pr_obs, 0)
pr_obs_var = dim_variance_n_Wrap(pr_obs, 0)


diriobs = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
ufili = "uwnd.mon.mean.nc"
vfili = "vwnd.mon.mean.nc"
fuobs = addfile(diriobs+ufili, "r")
fvobs = addfile(diriobs+vfili, "r")
TIME := fuobs->time
YYYYo := cd_calendar(TIME, -1)/100
iYYYYo := ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)
U_OBS = short2flt(fuobs->uwnd(iYYYYo,{LEV},::-1,:))
V_OBS = short2flt(fvobs->vwnd(iYYYYo,{LEV},::-1,:))
u_obs = month_to_season(U_OBS, season)
v_obs = month_to_season(V_OBS, season)
printVarSummary(u_obs)
ulat = u_obs&lat
ulon = u_obs&lon
u_obs_ave = dim_avg_n_Wrap(u_obs, 0)
u_obs_var = dim_variance_n_Wrap(u_obs, 0)
v_obs_ave = dim_avg_n_Wrap(v_obs, 0)
v_obs_var = dim_variance_n_Wrap(v_obs, 0)
;===================================================================
; GCM 19 data
;===================================================================
yrStrt = 1951
yrLast = 2005
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)
sig = 0.05

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP5/"

; models = (/"ACCESS1-0","bcc-csm1-1","CanCM4"/)

models = (/"ACCESS1-0","bcc-csm1-1","CanCM4","CanESM2","CCSM4","CESM1-BGC",\ ;
          "CESM1-CAM5","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5",\
          "FIO-ESM","GFDL-CM2p1","GFDL-ESM2M","GISS-E2-H","GISS-E2-R-CC",\
          "GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR","IPSL-CM5B-LR",\
          "MIROC4h","MIROC5","MIROC-ESM-CHEM","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR",\
          "MPI-ESM-P","MRI-CGCM3","NorESM1-M"/)

nm = dimsizes(models)
diff_pr5 = new((/nm,72,144/), float)
diff_pr5!0 = "model"
diff_pr5&model = models
; prob_pr5 = diff_pr5
diff_pr6 = diff_pr5

diff_u5 = new((/nm,73,144/), float)
diff_u5!0 = "model"
diff_u5&model = models
diff_v5 = diff_u5
; prob_uv = diff_u
diff_u6 = diff_u5
diff_v6 = diff_v5

do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, plon, plat, 0)
    pr_gcm = month_to_season(PR_GCMx, season)
    pr_gcm_ave = dim_avg_n_Wrap(pr_gcm, 0)
    ; pr_gcm_var = dim_variance_n_Wrap(pr_gcm, 0)

    ; prob_pr5(i,:,:) = ttest(pr_obs_ave, pr_obs_var, nyro, pr_gcm_ave, pr_gcm_var, nyr, False, False)
    diff_pr5(i,:,:) = pr_gcm_ave - pr_obs_ave
    copy_VarCoords(pr_obs_ave, diff_pr5(i,:,:))
    ; diff_pr(i,:,:) = mask(diff_pr(i,:,:), prob_pr(i,:,:).le.0.05, True)

    delete([/PR_GCMs,PR_GCMx,pr_gcm,pr_gcm_ave/])

    filiu = systemfunc("ls "+dirigcm+"ua_Amon_"+models(i)+"_historical_r1i1p1_*.nc")
    fu = addfile(filiu, "r")
    filiv = systemfunc("ls "+dirigcm+"va_Amon_"+models(i)+"_historical_r1i1p1_*.nc")
    fv = addfile(filiv, "r")
    TIME := fu->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    U_GCMs := fu->ua(iYYYY,{85000},:,:)
    U_GCMx = linint2_Wrap(U_GCMs&lon, U_GCMs&lat, U_GCMs, True, ulon, ulat, 0)
    V_GCMs := fv->va(iYYYY,{85000},:,:)
    V_GCMx = linint2_Wrap(V_GCMs&lon, V_GCMs&lat, V_GCMs, True, ulon, ulat, 0)

    u_gcm = month_to_season(U_GCMx, season)
    u_gcm_ave = dim_avg_n_Wrap(u_gcm, 0)
    ; u_gcm_var = dim_variance_n_Wrap(u_gcm, 0)
    ; prob_u = ttest(u_obs_ave, u_obs_var, nyro, u_gcm_ave, u_gcm_var, nyr, False, False)

    v_gcm = month_to_season(V_GCMx, season)
    v_gcm_ave = dim_avg_n_Wrap(v_gcm, 0)
    ; v_gcm_var = dim_variance_n_Wrap(v_gcm, 0)
    ; prob_v = ttest(v_obs_ave, v_obs_var, nyro, v_gcm_ave, v_gcm_var, nyr, False, False)

    diff_u5(i,:,:) = u_gcm_ave - u_obs_ave
    copy_VarCoords(u_obs_ave, diff_u5(i,:,:))
    diff_v5(i,:,:) = v_gcm_ave - v_obs_ave
    copy_VarCoords(v_obs_ave, diff_v5(i,:,:))
    ; prob_uv(i,:,:) = mask(prob_u, (prob_u.le.sig) .and. (prob_v.le.sig), True)

    ; diff_u(i,:,:) = mask(diff_u(i,:,:), (prob_u.le.sig) .and. (prob_v.le.sig), True)
    ; diff_v(i,:,:) = mask(diff_v(i,:,:), (prob_u.le.sig) .and. (prob_v.le.sig), True)

    delete([/U_GCMs,U_GCMx,u_gcm,u_gcm_ave/])
    delete([/V_GCMs,V_GCMx,v_gcm,v_gcm_ave/])
    print("end of "+models(i))
    print("============================================")
end do
printVarSummary(diff_u5)
;===================================================================
yrStrt = 1951
yrLast = 2014
year6 = ispan(yrStrt, yrLast, 1)
nyr6 = dimsizes(year6)

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"

models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
           "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
           "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
; models := (/"CAMS-CSM1-0","CanESM5"/)

nm := dimsizes(models)
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1f*_g*_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, plon, plat, 0)
    pr_gcm = month_to_season(PR_GCMx, season) ;PR_GCMx
    pr_gcm_ave = dim_avg_n_Wrap(pr_gcm, 0)
    ; pr_gcm_var = dim_variance_n_Wrap(pr_gcm, 0)

    ; prob_pr(i,:,:) = ttest(pr_obs_ave, pr_obs_var, nyro, pr_gcm_ave, pr_gcm_var, nyr, False, False)
    diff_pr6(i,:,:) = pr_gcm_ave - pr_obs_ave
    copy_VarCoords(pr_obs_ave, diff_pr6(i,:,:))
    ; diff_pr(i,:,:) = mask(diff_pr(i,:,:), prob_pr(i,:,:).le.0.05, True)

    delete([/PR_GCMs,PR_GCMx,pr_gcm,pr_gcm_ave/])

    filiu = systemfunc("ls "+dirigcm+"ua_Amon_"+models(i)+"_historical_r1i1p1f*_*.nc")
    fu = addfile(filiu, "r")
    filiv = systemfunc("ls "+dirigcm+"va_Amon_"+models(i)+"_historical_r1i1p1f*_*.nc")
    fv = addfile(filiv, "r")
    TIME := fu->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    U_GCMs := fu->ua(iYYYY,{85000},:,:)
    U_GCMx = linint2_Wrap(U_GCMs&lon, U_GCMs&lat, U_GCMs, True, ulon, ulat, 0)
    V_GCMs := fv->va(iYYYY,{85000},:,:)
    V_GCMx = linint2_Wrap(V_GCMs&lon, V_GCMs&lat, V_GCMs, True, ulon, ulat, 0)

    u_gcm = month_to_season(U_GCMx, season)
    u_gcm_ave = dim_avg_n_Wrap(u_gcm, 0)
    ; u_gcm_var = dim_variance_n_Wrap(u_gcm, 0)
    ; prob_u = ttest(u_obs_ave, u_obs_var, nyro, u_gcm_ave, u_gcm_var, nyr, False, False)

    v_gcm = month_to_season(V_GCMx, season)
    v_gcm_ave = dim_avg_n_Wrap(v_gcm, 0)
    ; v_gcm_var = dim_variance_n_Wrap(v_gcm, 0)
    ; prob_v = ttest(v_obs_ave, v_obs_var, nyro, v_gcm_ave, v_gcm_var, nyr, False, False)

    diff_u6(i,:,:) = u_gcm_ave - u_obs_ave
    copy_VarCoords(u_obs_ave, diff_u6(i,:,:))
    diff_v6(i,:,:) = v_gcm_ave - v_obs_ave
    copy_VarCoords(v_obs_ave, diff_v6(i,:,:))
    ; prob_uv(i,:,:) = mask(prob_u, (prob_u.le.sig) .and. (prob_v.le.sig), True)

    ; diff_u(i,:,:) = mask(diff_u(i,:,:), (prob_u.le.sig) .and. (prob_v.le.sig), True)
    ; diff_v(i,:,:) = mask(diff_v(i,:,:), (prob_u.le.sig) .and. (prob_v.le.sig), True)

    delete([/U_GCMs,U_GCMx,u_gcm,u_gcm_ave/])
    delete([/V_GCMs,V_GCMx,v_gcm,v_gcm_ave/])
    print("end of "+models(i))
    print("============================================")
end do
printVarSummary(diff_u6)

diff_pr5_mme = dim_avg_n_Wrap(diff_pr5, 0)
diff_u5_mme = dim_avg_n_Wrap(diff_u5, 0)
diff_v5_mme = dim_avg_n_Wrap(diff_v5, 0)

diff_pr6_mme = dim_avg_n_Wrap(diff_pr6, 0)
diff_u6_mme = dim_avg_n_Wrap(diff_u6, 0)
diff_v6_mme = dim_avg_n_Wrap(diff_v6, 0)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

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

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "pr_uv850_bias_JJA_mme_CMIP5_CMIP6"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(2, graphic)
plot_topo = plot
plotv = plot
; plot_prob = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red"
; res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10
res@cnFillPalette = "MPL_BrBG"
res@gsnRightString = "mm/day"

;===========================================================
 vecres  = True
 vecres@gsnDraw  = False
 vecres@gsnFrame = False

 vecres@vcGlyphStyle     = "LineArrow"   ; curly vectors"Fillarrow" 
 vecres@vcRefLengthF     = 0.05            ; define length of vec ref
 vecres@vcMinDistanceF   = 0.03
 vecres@vcLineArrowThicknessF = 3.0
 vecres@tiXAxisString    = " "             ; turn off axis label
 vecres@vcRefAnnoOn      = True            ; draw the reference vector annotation.
 vecres@vcRefAnnoString2On =  False          ;display a string below or to the right
 vecres@vcRefAnnoOrientation = "horizontal"
 vecres@vcRefAnnoString2  = " "
 vecres@vcRefAnnoOrthogonalPosF = -1.0  ; move ref vector into plot
 vecres@vcFillArrowsOn       = True
 vecres@vcLineArrowHeadMinSizeF = 0.008
 vecres@vcLineArrowHeadMaxSizeF = 0.3
 vecres@gsnLeftString  = ""
 ;vecres@gsnRightString  = "g*s/Kg"
 vecres@gsnRightString  = ""
 vecres@gsnRightStringFontHeightF = 0.04
 
 ; vecres@vcLineArrowColor = "gray52"
 vecres@vcRefAnnoFontHeightF = 0.02
 vecres@vcRefMagnitudeF  = 2.             ; define vector ref mag
 ; vecres@vcRefAnnoString1  = "3 Kg/m/s"


gres = True
gres@gsLineColor = "red"
gres@gsLineThicknessF = 3.0
glat1 = (/40,40,28,28,40/)
glon1 = (/85,105,105,85,85/)
glat2 = (/28,28,18,18,28/)
glon2 = (/72,85,85,72,72/)

res@gsnLeftString = "CMIP5 MME"
plot(0) = gsn_csm_contour_map(wks, diff_pr5_mme, res)
plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
plotv(0) = gsn_csm_vector(wks, diff_u5_mme, diff_v5_mme, vecres)
overlay(plot(0), plotv(0))
dum1 = gsn_add_polyline(wks,plot(0),glon1,glat1,gres)
dum2 = gsn_add_polyline(wks,plot(0),glon2,glat2,gres)

res@gsnLeftString = "CMIP6 MME"
plot(1) = gsn_csm_contour_map(wks, diff_pr6_mme, res)
plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))
plotv(1) = gsn_csm_vector(wks, diff_u6_mme, diff_v6_mme, vecres)
overlay(plot(1), plotv(1))
dum3 = gsn_add_polyline(wks,plot(1),glon1,glat1,gres)
dum4 = gsn_add_polyline(wks,plot(1),glon2,glat2,gres)

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/1, 2/), pres)


end