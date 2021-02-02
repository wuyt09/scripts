
begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS t2m data
;===================================================================
yrStrt = 1979
yrLast = 2005
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "MAM"

;========================================================
; obs prec
;========================================================
diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
TIME = fobs->time
YYYYo = cd_calendar(TIME, -1)/100
iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)
PR_OBS = fobs->precip(iYYYYo,:,:)
timeo = fobs->time(iYYYYo)
PR_OBS = dtrend_msg_n(timeo, PR_OBS, False, False, 0)
PR_OBSx = month_to_season(PR_OBS, season)
pr_obs = dim_avg_n_Wrap(PR_OBSx, 0)
printVarSummary(pr_obs)

lat = pr_obs&lat
lon = pr_obs&lon
;===================================================================
; GCM 19 data
;===================================================================


models1 = (/"ACCESS1-0","bcc-csm1-1","CanCM4","CanESM2","CCSM4","CESM1-BGC",\
          "CESM1-CAM5","CESM1-FASTCHEM","CESM1-WACCM","CMCC-CM","CNRM-CM5",\
          "EC-EARTH","FIO-ESM","GFDL-CM2p1","GFDL-ESM2M","GISS-E2-H","GISS-E2-R-CC",\
          "GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR",\
          "MIROC4h","MIROC5","MIROC-ESM-CHEM","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR",\
          "MRI-CGCM3","NorESM1-M"/)

models2 = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

models = venn2_intersection(models1,models2)

; print(models)
; print(dimsizes(models))

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP5/"
nm = dimsizes(models)
diff = new((/nm,72,144/), float)
diff!0 = "model"
diff&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    times = fgcm->time(iYYYY)
    PR_GCMs = dtrend_msg_n(times, PR_GCMs, False, False, 0)
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, lon, lat, 0)
    PR_GCMs_jja = month_to_season(PR_GCMx, season)
    ; printVarSummary(PR_GCMs_jja)
    PR_GCM6 = dim_avg_n_Wrap(PR_GCMs_jja, 0)
    diff(i,:,:) = PR_GCM6 - pr_obs
    copy_VarCoords(PR_GCM6, diff(i,:,:))
    delete([/PR_GCMs/])
    print("end of prec "+models(i))
end do
printVarSummary(diff)


;==================================================================
; sort by bias
;==================================================================

sst_bias_corr = asciiread("./calc_pattern_corr_SST_bias_GCM_with_MME_CMIP5.txt", nm, "float")
printVarSummary(sst_bias_corr)
print(sst_bias_corr)

ip_tbias = dim_pqsort(sst_bias_corr,-1)
print(models(ip_tbias)+"  "+sprintf("%4.2f",sst_bias_corr(ip_tbias)))

;===================================================================
; topo data
;===================================================================
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
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_FAMIP/"
filiname = "check_Ann_SST_bias_GCM_with_MME_CMIP5"

wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(nm, graphic)
plot_topo = plot

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
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 5
res@cnMinLevelValF = -5
res@cnFillPalette = "CBR_drywet"
res@gsnRightString = "mm/day"

do i = 0, nm-1
    res@gsnLeftString = models(ip_tbias(i))
    res@gsnCenterString = sprintf("%4.2f", sst_bias_corr(ip_tbias(i)))
    plot(i) = gsn_csm_contour_map(wks, diff(ip_tbias(i),:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/4, 5/), pres)
end