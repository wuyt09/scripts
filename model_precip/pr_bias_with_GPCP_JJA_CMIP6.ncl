
begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

season = "JJA"
yrStrt = 1979
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)

diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
TIME = fobs->time
YYYYo = cd_calendar(TIME, -1)/100
iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)
PR_OBS = fobs->precip(iYYYYo,:,:)
PR_OBSx = month_to_season(PR_OBS, season)
pr_obs = dim_avg_n_Wrap(PR_OBSx, 0)
printVarSummary(pr_obs)

;===================================================================
; GCM 19 data
;===================================================================
yrStrt = 1951
yrLast = 2014
dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"

models = (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
           "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
           "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)

nm = dimsizes(models)
diff = new((/nm,72,144/), float)
diff!0 = "model"
diff&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1f*_g*_185001-201*12.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    PR_GCMs_jja = month_to_season(PR_GCMs, season)
    printVarSummary(PR_GCMs_jja)
    PR_GCM6 = dim_avg_n_Wrap(PR_GCMs_jja, 0)
    diff(i,:,:) = PR_GCM6 - pr_obs
    copy_VarCoords(PR_GCM6, diff(i,:,:))
    print("end of "+models(i))
end do
printVarSummary(PR_GCM6)


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
filiname = "pr_bias_with_GPCP_JJA_CMIP6"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(20, graphic)
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
res@cnLevelSpacingF = 1
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10
res@cnFillPalette = "CBR_drywet"
res@gsnRightString = "mm/day"

do i = 0, nm-1
    res@gsnLeftString = models(i)
    plot(i) = gsn_csm_contour_map(wks, diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do


diff_mme = dim_avg_n_Wrap(diff, 0)
res@gsnLeftString = "MME"
plot(19) = gsn_csm_contour_map(wks, diff_mme, res)
plot_topo(19) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(19), plot_topo(19))

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/4, 5/), pres)

end