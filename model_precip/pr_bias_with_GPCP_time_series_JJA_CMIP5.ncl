
begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

season = "JJA"
yrStrt = 1979
yrLast = 2005
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)

diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
TIME = fobs->time
YYYYo = cd_calendar(TIME, -1)/100
iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)
PR_OBS = fobs->precip(iYYYYo,:,:)
pr_obs = month_to_season(PR_OBS, season)
; pr_obs = dim_avg_n_Wrap(PR_OBSx, 0)
printVarSummary(pr_obs)

lat = pr_obs&lat
lon = pr_obs&lon
nlat = dimsizes(lat)
nlon = dimsizes(lon)
;===================================================================
; GCM 19 data
;===================================================================
dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP5/"

models = (/"ACCESS1-0","bcc-csm1-1","CanCM4","CanESM2","CCSM4","CESM1-BGC",\
          "CESM1-CAM5","CESM1-FASTCHEM","CESM1-WACCM","CNRM-CM5",\
          "FIO-ESM","GFDL-CM2p1","GFDL-ESM2M","GISS-E2-H","GISS-E2-R-CC",\
          "GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR","IPSL-CM5B-LR",\
          "MIROC4h","MIROC5","MIROC-ESM-CHEM","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR",\
          "MPI-ESM-P","MRI-CGCM3","NorESM1-M"/)
; models = (/"ACCESS1-0","bcc-csm1-1"/)

nm = dimsizes(models)
diff_pr = new((/nm,nyr,nlat,nlon/), float)
diff_pr!0 = "model"
diff_pr&model = models
prob_pr = diff_pr

do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, lon, lat, 0)
    pr_gcm = month_to_season(PR_GCMx, season)
    diff_pr(i,:,:,:) = pr_gcm - pr_obs
    copy_VarCoords(pr_obs, diff_pr(i,:,:,:))
    ; diff_pr(i,:,:) = mask(diff_pr(i,:,:), prob_pr(i,:,:).le.0.05, True)

    delete([/PR_GCMs,PR_GCMx,pr_gcm/])
    print("end of "+models(i))
end do
printVarSummary(diff_pr)

diff_mme = dim_avg_n_Wrap(diff_pr, 0)
printVarSummary(diff_mme)
diff_tp = wgt_areaave_Wrap(diff_mme(:,{28:40},{85:105}), 1., 1., 0)
diff_is = wgt_areaave_Wrap(diff_mme(:,{18:28},{72:85}), 1., 1., 0)
printVarSummary(diff_tp)

diff_std_tp = dim_standardize_Wrap(diff_tp, 1)
diff_std_is = dim_standardize_Wrap(diff_is, 1)
corr = escorc(diff_tp, diff_is)
n = nyr
df = nyr-2
t = corr*sqrt((n-2)/(1-corr^2))
p = student_t(t, df)
print(sprintf("%4.2f", p))

pp = new((/2,nyr/), float)
pp(0,:) = diff_std_tp
pp(1,:) = diff_std_is
copy_VarCoords(diff_tp, pp(0,:))
printVarSummary(pp)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)


;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "pr_bias_with_GPCP_time_series_JJA_CMIP5"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = year
; res@tmXBLabels = year
res@tmXTOn = False
; res@tiXAxisString = "Month"

res@gsnYRefLine = 0
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineThicknessF = 10
res@trYMinF = -3
res@trYMaxF = 3
; res@tmYROn = False
res@tiYAxisString = "Precipitation diff (mm/day)"
res@gsnLeftString = "Corr: "+sprintf("%4.2f", corr)
res@gsnRightString = "Sig: "+sprintf("%4.2f", (1-p))

Colors = (/"black","red"/)
Patterns = (/0,0/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, year, pp, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"CE-TP","North-India"/)

;============================================================
    gx = (/1980,1982/)
    text_pos = 1983
    interval = 0.5
    gy1 = (/2.5,2.5/)
    gy2 = gy1 - interval

    res_lines@gsLineThicknessF = 15
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)

draw(plot(0))

end