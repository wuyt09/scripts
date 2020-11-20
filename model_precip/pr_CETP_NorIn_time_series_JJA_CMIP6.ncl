
begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

season = "JJA"
yrStrt = 1861
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)

; diriobs = "/home/yangsong3/data-observation/GPCP/"
; filiobs = "precip.mon.mean.nc"
; fobs = addfile(diriobs+filiobs, "r")
; TIME = fobs->time
; YYYYo = cd_calendar(TIME, -1)/100
; iYYYYo = ind(YYYYo.ge.yrStrt .and. YYYYo.le.yrLast)
; PR_OBS = fobs->precip(iYYYYo,:,:)
; pr_obs = month_to_season(PR_OBS, season)
; ; pr_obs = dim_avg_n_Wrap(PR_OBSx, 0)
; printVarSummary(pr_obs)

; lat = pr_obs&lat
; lon = pr_obs&lon
; nlat = dimsizes(lat)
; nlon = dimsizes(lon)
;=============================================================

dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"

models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
           "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
           "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
; models = (/"CAMS-CSM1-0","CanESM5"/)

nm := dimsizes(models)
pr_TP = new((/nm,nyr/), float)
pr_TP!0 = "model"
pr_TP&model = models
pr_IS = pr_TP

do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"pr_Amon_"+models(i)+"_historical_r1i1p1f*_*.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PR_GCMs = fgcm->pr(iYYYY,:,:)
    PR_GCMs = PR_GCMs * 86400
    ; PR_GCMx = linint2_Wrap(PR_GCMs&lon, PR_GCMs&lat, PR_GCMs, True, lon, lat, 0)
    pr_gcm = month_to_season(PR_GCMs, season)
    pr_TP(i,:) = wgt_areaave_Wrap(pr_gcm(:,{28:40},{85:105}), 1., 1., 0)
    pr_IS(i,:) = wgt_areaave_Wrap(pr_gcm(:,{18:28},{72:85}), 1., 1., 0)
    ; diff_pr6(i,:,:,:) = pr_gcm - pr_obs
    ; copy_VarCoords(pr_obs, diff_pr6(i,:,:,:))
    ; diff_pr(i,:,:) = mask(diff_pr(i,:,:), prob_pr(i,:,:).le.0.05, True)

    delete([/PR_GCMs,pr_gcm/]) ;
    print("end of "+models(i))
end do
pr_tp = dim_avg_n_Wrap(pr_TP, 0)
pr_is = dim_avg_n_Wrap(pr_IS, 0)
printVarSummary(pr_tp)

std_tp = dim_standardize_Wrap(pr_tp, 1)
std_is = dim_standardize_Wrap(pr_is, 1)
corr6 = escorc(std_tp, std_is)
printVarSummary(std_tp)

n = nyr
df = n-2
t6 = corr6*sqrt((n-2)/(1-corr6^2))
p6 = student_t(t6, df)
print(sprintf("%4.2f", p6))
;----------------------------------------

pp = new((/2,nyr/), float)
; pp(0,:) = diff5_std_tp
; pp(1,:) = diff5_std_is
pp(0,:) = std_tp
pp(1,:) = std_is
copy_VarCoords(std_tp, pp(0,:))
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
filiname = "pr_CETP_NorIn_time_series_JJA_CMIP5"
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
res@gsnYRefLineThicknessF = 8
res@trYMinF = -3
res@trYMaxF = 3
; res@tmYROn = False
res@tiYAxisString = "Precip std(mm/day)"
res@gsnLeftString = "Corr: "+sprintf("%4.2f", corr6)
res@gsnRightString = "Sig: "+sprintf("%4.0f", 100*(1-p6))+"%"
res@gsnStringFontHeightF = 0.015

Colors = (/"black","red","black","red"/)
Patterns = (/0,0,16,16/)
res@xyLineThicknessF = 5
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
    gx = (/1930,1940/)
    text_pos = 1942
    interval = 0.5
    gy1 = (/2.7,2.7/)
    gy2 = gy1 - interval

    gx2 = (/1990,1992/)
    text_pos2 = 1993
    gy3 = gy1
    gy4 = gy2

    res_lines@gsLineThicknessF = 15
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    ; res_lines@gsLineColor = Colors(2)
    ; res_lines@gsLineDashPattern = Patterns(2)
    ; dum31 = gsn_add_polyline(wks, plot(0), gx2, gy3, res_lines)
    ; dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos2, gy3, res_text)
    ; res_lines@gsLineColor = Colors(3)
    ; res_lines@gsLineDashPattern = Patterns(3)
    ; dum41 = gsn_add_polyline(wks, plot(0), gx2, gy4, res_lines)
    ; dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos2, gy4, res_text)

draw(plot(0))

end