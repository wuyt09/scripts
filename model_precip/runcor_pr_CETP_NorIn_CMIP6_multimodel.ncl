load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"

begin

lonstr = 60
lonend = 135
latstr = 5
latend = 50

season = "JJA"
yrStrt = 1851
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
year@units = "hours since 1-1-1 00:00:0.0"
nyear = dimsizes(year)

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
; trend_tp = new((/nm,nyr/), float)
; trend_tp!0 = "model"
; trend_tp&model = models
; trend_is = trend_tp
; line_tp = trend_tp
; line_is = trend_tp

; corr6 = new(nm, float)
; p6 = corr6
yr = ispan(yrStrt+10, yrLast-10, 1)
nyr = dimsizes(yr)
cor = new((/nm,nyr/), float)

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
    ; pr_TP(i,:) = wgt_areaave_Wrap(pr_gcm(:,{28:40},{85:105}), 1., 1., 0)
    ; ; pr_IS(i,:) = wgt_areaave_Wrap(pr_gcm(:,{18:28},{72:85}), 1., 1., 0)
    ; pr_IS(i,:) = wgt_areaave_Wrap(pr_gcm(:,{15:30},{70:80}), 1., 1., 0)
    pr_tp = wgt_areaave_Wrap(pr_gcm(:,{28:40},{85:105}), 1., 1., 0)
    pr_is = wgt_areaave_Wrap(pr_gcm(:,{15:30},{70:80}), 1., 1., 0)

    cor(i,:) = run_cor(pr_tp, pr_is, year, 21)
    ; rc_tp = regline(year, pr_tp)
    ; rc_is = regline(year, pr_is)
    ; ; printVarSummary(rc_tp)

    ; trend_tp(i,:) = year*rc_tp+rc_tp@yintercept
    ; trend_is(i,:) = year*rc_is+rc_is@yintercept

    ; line_tp(i,:) = runave_n_Wrap(pr_tp, 10, 0, 0)
    ; line_is(i,:) = runave_n_Wrap(pr_is, 10, 0, 0)

    ; corr6_med = escorc(pr_tp, pr_is)
    ; corr6(i) = corr6_med
    ; n = nyr
    ; df = n-2
    ; t6 = corr6_med*sqrt((n-2)/(1-corr6_med^2))
    ; p6(i) = student_t(t6, df)

    ; delete([/PR_GCMs,pr_gcm,pr_tp,pr_is,rc_tp,rc_is,corr6_med,n,df,t6/]) ;
    print("end of "+models(i))
end do

;----------------------------------------

; pp = new((/4,nm,nyr/), float)
; ; pp(0,:) = diff5_std_tp
; ; pp(1,:) = diff5_std_is
; pp(0,:,:) = line_tp
; pp(1,:,:) = line_is
; pp(2,:,:) = trend_tp
; pp(3,:,:) = trend_is
; copy_VarCoords(line_tp, pp(0,:,:))
; printVarSummary(pp)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)


;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "runcor_pr_CETP_NorIn_multimodel" ;_1951
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(nm+1, graphic)

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
res@trYMinF = -0.8
res@trYMaxF = 0.8
; res@tmYROn = False
; res@tiYAxisString = "Precip std(mm/day)"
res@gsnStringFontHeightF = 0.025

Colors = (/"black","red","black","red"/)
Patterns = (/0,0,16,16/)
res@xyLineThicknessF = 8
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

do i = 0, nm-1
    res@gsnLeftString = models(i)
    ; res@gsnCenterString = "Corr: "+sprintf("%4.2f", corr6(i))
    ; res@gsnRightString = "Sig: "+sprintf("%4.0f", 100*(1-p6(i)))+"%"
    ; plot(i) = gsn_csm_xy(wks, year, pp(:,i,:), res)
    plot(i) = gsn_csm_xy(wks, yr, cor(i,:), res)
end do

; ;==================================================================================
;     res_lines                   = True                  ; polyline mods desired
;     res_lines@gsLineThicknessF  = 5                     ; line thicker
    

;     res_text                    = True                  ; text mods desired
;     res_text@txFontHeightF      = 0.018                 ; change text size
;     res_text@txJust             = "CenterLeft"          ; text justification

;     Texts  = (/"CE-TP ("+sprintf("%4.2f", rc_tp*10)+" mm/day/decade)","North-India ("+sprintf("%4.2f", rc_is*10)+" mm/day/decade)"/)

; ;============================================================
;     gx = (/1920,1930/)
;     text_pos = 1932
;     interval = 0.1
;     gy1 = (/4.75,4.75/)
;     gy2 = gy1 - interval

;     gx2 = (/1990,1992/)
;     text_pos2 = 1993
;     gy3 = gy1
;     gy4 = gy2

;     res_lines@gsLineThicknessF = 15
;     res_lines@gsLineColor = Colors(0)
;     res_lines@gsLineDashPattern = Patterns(0)
;     dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
;     dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
;     res_lines@gsLineColor = Colors(1)
;     res_lines@gsLineDashPattern = Patterns(1)
;     dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
;     dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)

pres = True
pres@gsnMaximize = True
do i = 0, 4
    gsn_panel(wks, plot(i*4:i*4+3), (/2,2/), pres)
end do

end