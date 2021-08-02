begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

season = "JJA"
yrStrt = 1851
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

;=============================================================
diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
; TIME = fobs->time
; YYYYo = cd_calendar(TIME, -1)/100
; iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)
pr_obs = fobs->precip(0,:,:)
printVarSummary(pr_obs)

lat = pr_obs&lat
lon = pr_obs&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)

;=============================================================

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"

models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
           "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
           "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
; models = (/"CAMS-CSM1-0","CanESM5"/)

nm := dimsizes(models)
corr = new((/nm,nlat,nlon/), float)
pr_gcm = new((/nm,nyear,nlat,nlon/), float)

do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1f*_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, lon, lat, 0)
    pr_jul = PR_GCMx(6::12,:,:)
    pr_aug = PR_GCMx(7::12,:,:)
    ; pr_gcm = month_to_season(PR_GCMs, season)
    pr_gcm_med = (pr_aug + pr_jul)/2
    copy_VarCoords(pr_aug, pr_gcm_med)
    pr_gcm(i,:,:,:) = pr_gcm_med
    pr_inx = wgt_areaave_Wrap(pr_gcm_med(:,{15:30},{70:80}), 1., 1., 0)
    corr(i,:,:) = escorc_n(pr_inx, pr_gcm_med, 0, 0)
    copy_VarCoords(PR_GCMx(0,:,:), corr(i,:,:))

    delete([/PR_GCMs,PR_GCMx,pr_jul,pr_aug,pr_gcm_med,pr_inx/])
    print("end of "+models(i))
end do
printVarSummary(corr)

pr_mme = dim_avg_n_Wrap(pr_gcm, 0)
pr_inx_mme = wgt_areaave_Wrap(pr_mme(:,{15:30},{70:80}), 1., 1., 0)
corr_mme = escorc_n(pr_inx_mme, pr_mme, 0, 0)
copy_VarCoords(pr_mme(0,:,:), corr_mme)
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
filiname = "Corr_IS_pr_CMIP6_1980"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(nm+1, graphic)
plot_topo = plot
dum = plot

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
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.8
res@cnMinLevelValF = -0.8
; res@cnLevelSpacingF = 1.
; res@cnMaxLevelValF = 10
; res@cnMinLevelValF = -10
res@cnFillPalette = "MPL_BrBG"
res@gsnRightString = "1"

;---------------------------------------
gres = True
gres@gsLineColor = "red"
gres@gsLineThicknessF = 5.0
glat1 = (/30,30,15,15,30/)
glon1 = (/70,80,80,70,70/)
;---------------------------------------

do i = 0, nm-1
    res@gsnLeftString = models(i)
    plot(i) = gsn_csm_contour_map(wks, corr(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
    dum(i) = gsn_add_polyline(wks,plot(i),glon1,glat1,gres)
end do

res@gsnLeftString = "MME"
plot(nm) = gsn_csm_contour_map(wks, corr_mme, res)
plot_topo(nm) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(nm), plot_topo(nm))
dum(nm) = gsn_add_polyline(wks,plot(nm),glon1,glat1,gres)


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/4,5/), pres)

end
