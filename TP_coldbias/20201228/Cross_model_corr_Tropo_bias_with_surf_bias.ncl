
begin

lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

season = (/"DJF","JJA"/)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
lato = topo&lat
lono = topo&lon
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1980
yrLast = 2005
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

ts_obs = new((/nyear*nm,361,576/), float)

diriobs = "/home/yangsong3/data-observation/MERRA2/monthly/surfT/"

do i = 0, nyear-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        filiobs = systemfunc("ls "+diriobs+"MERRA2_*.instM_2d_lfo_Nx."+yr+mn+".SUB.nc4")
        fobs = addfile(filiobs, "r")
        ts_obs(i*12+j,:,:) = fobs->TLML(0,:,:)
    end do
end do
printVarSummary(ts_obs)

time = ispan(1, nyear*12, 1)

ts_obs!0 = "time"
; ts_obs!1 = "lat"
; ts_obs!2 = "lon"
ts_obs&time = time
; ts_obs&lat = lat
; ts_obs&lon = lon
; ts_obs&lat@units = "degrees_north"
; ts_obs&lon@units = "degrees_east"
printVarSummary(ts_obs)

ts_obs_seax = month_to_seasonN(ts_obs, season)
ts_obs_sea = dim_avg_n_Wrap(ts_obs_seax, 1)
printVarSummary(ts_obs_sea)

ts_obs_sea_TP = linint2_Wrap(ts_obs&lon, ts_obs&lat, ts_obs_sea, False, lono, lato, 0)
; print(topox2({35},{90}))

ts_obs_sea_TP = mask(ts_obs_sea_TP, topo.ge.1500, True)
TP_obs_sea = wgt_areaave_Wrap(ts_obs_sea_TP(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs_sea)
print(TP_obs_sea)

;===================================================================
; GCM 20 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/"

models1 = (/"ACCESS1-0","bcc-csm1-1","BNU-ESM","CanCM4","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CESM","CNRM-CM5","GFDL-ESM2M","GISS-E2-H","CMCC-CM","CMCC-CMS","GFDL-CM3","GFDL-ESM2G",\
        "IPSL-CM5A-MR","NorESM1-ME","GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR",\
        "IPSL-CM5B-LR","MIROC5","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
;"CESM1-WACCM","FIO-ESM","MPI-ESM-P","CSIRO-Mk3-6-0",
;

models2 = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

models = venn2_intersection(models1,models2)

nm = dimsizes(models)

Tropo_T = asciiread("./write_Tropo_T_bias.txt", (/nm,2/), "float")
print(Tropo_T(nm-1,:))


TP_GCM5 = new((/nm,2/), float)
TP_GCM5!0 = "model"
TP_GCM5&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+models(i)+"/tas_Amon_"+models(i)+"_historical_r1i1p1_18*01-20*12.nc")
    ; print(filigcm)
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.2005)
    TS_GCMs := fgcm->tas(iYYYY,:,:)
    TS_GCMs!1 = "lat"
    TS_GCMs!2 = "lon"
    ; TS_GCMs = TS_GCMs - 273.15
    TS_GCMx := month_to_seasonN(TS_GCMs, season)
    TS_GCM := dim_avg_n_Wrap(TS_GCMx, 1)

    TOPO := linint2_Wrap(topo&lon, topo&lat, topo, True, TS_GCM&lon, TS_GCM&lat, 0)
    TS_GCM = mask(TS_GCM, TOPO.ge.1500, True)
    TP_GCM5(i,:) = wgt_areaave_Wrap(TS_GCM(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
    print("end of "+models(i))
end do
printVarSummary(TP_GCM5)
;================================================================================

lines = new((/4,nm/), float)
lines(0,:) = TP_GCM5(:,0) - TP_obs_sea(0)
lines(1,:) = TP_GCM5(:,1) - TP_obs_sea(1)
lines(2,:) = Tropo_T(:,0)
lines(3,:) = Tropo_T(:,1)

; lines!0 = "season"
; lines&season = season
lines!1 = "model"
lines&model = models
printVarSummary(lines)

Corr_DJF = escorc(lines(0,:), lines(2,:))
Corr_JJA = escorc(lines(1,:), lines(3,:))
;=========================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/20201228/"
figname = "Corss_model_corr_Tropo_bias_with_surf_bias"
; pwks = "png"
; pwks@wkHeight              = 3048
; pwks@wkWidth               = 6048
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

; res@trXMinF = 1
; res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = ispan(1, nm, 1)
res@tmXBLabels = models
res@tmXTOn = False
res@tmXBLabelAngleF = 45
res@tmXBLabelFontHeightF = 0.007
; res@tiXAxisString = "Month"

res@trYMinF = -8
res@trYMaxF = 4
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-5, 4, 1)
res@tmYLPrecision = 2
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYROn = False
res@tiYAxisString = ""
; res@gsnLeftString = season
res@gsnRightString = "DJF Corr: "+sprintf("%4.2f", Corr_DJF)+", JJA Corr: "+sprintf("%4.2f", Corr_JJA)
res@gsnStringFontHeightF = 0.02

res@gsnYRefLine = 0
res@gsnYRefLineDashPattern = 2
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"blue","red","dodgerblue","firebrick3"/) ;
Patterns = (/0,0,16,16/) ;
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, ispan(1, nm, 1), lines, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/season(0)+" Tas bias",season(1)+" Tas bias",season(0)+" Tropo T bias",season(1)+" Tropo T bias"/)

;============================================================
    gx = (/1.5,2.5/)
    text_pos = 2.8
    gy1 = (/3,3/)
    gy2 = (/2,2/)

    gx2 = (/8.5,9.5/)
    text_pos2 = 9.8


    res_lines@gsLineThicknessF = 10
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx2, gy1, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos2, gy1, res_text)
    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    dum41 = gsn_add_polyline(wks, plot(0), gx2, gy2, res_lines)
    dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos2, gy2, res_text)

draw(plot(0))
end
