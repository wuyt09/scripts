
begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

season = "JJA"
yrStrt = 2015
yrLast = 2099
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)

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

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/ssp585/"

; models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
;            "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
;            "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1",\
           "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G",\
           "HadGEM3-GC31-LL","MPI-ESM1-2-HR","NESM3"/)


nm := dimsizes(models)
rc_gcm = new((/nm,nlat,nlon/), float)
rc_gcm!0 = "model"
rc_gcm&model = models
; pr_TP = new((/nm,nyr/), float)
; pr_TP!0 = "model"
; pr_TP&model = models
; pr_IS = pr_TP

do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_ssp585_*_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, lon, lat, 0)
    pr_gcm = month_to_season(PR_GCMx, season)
    rc_gcm(i,:,:) = regCoef_n(year, pr_gcm, 0, 0)
    copy_VarCoords(pr_gcm(0,:,:), rc_gcm(i,:,:))

    delete([/PR_GCMs,pr_gcm/]) ;
    print("end of "+models(i))
end do

rc_gcm = rc_gcm*10
rc_mme = dim_avg_n_Wrap(rc_gcm, 0)
printVarSummary(rc_mme)

;----------------------------------------


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
filiname = "trend_pr_JJA_CMIP6_ssp585_MME"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(nm+1, graphic)
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
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 1
res@cnMinLevelValF = -1
; res@cnLevelSpacingF = 1.
; res@cnMaxLevelValF = 10
; res@cnMinLevelValF = -10
res@cnFillPalette = "MPL_BrBG"
res@gsnRightString = "mm/day/yr"

do i = 0, nm-1
    res@gsnLeftString = models(i)
    plot(i) = gsn_csm_contour_map(wks, rc_gcm(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

res@gsnLeftString = "MME"
plot(nm) = gsn_csm_contour_map(wks, rc_mme, res)
plot_topo(nm) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(nm), plot_topo(nm))


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/4,4/), pres)

end