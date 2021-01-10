
begin

yrStrt = 1980
yrLast = 2005
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/) ;ispan(mnStrt, mnLast, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
season = "DJF"

lat = ispan(-90, 90, 1)
lat@units = "degree_north"
lon = ispan(0, 359, 1)
lon@units = "degree_east"
; print(lat)
; print(lon)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
LEV = (/1000,925,850,700,600,500,400,300,250,200,150,100,50,10/)
LEV@units = "hPa"
LEV2 = LEV*100
LEV2@units = "Pa"
nlev = dimsizes(LEV)
LAT = 35

lonstr = 60 ; 72 ;75
lonend = 120 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================
t = new((/nyear*nmonth,nlev,nlat,nlon/), float)

diri = "/home/yangsong3/data-observation/MERRA2/monthly/t/"

do i = 0, nyear-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        fili = systemfunc("ls "+diri+"MERRA2_*.instM_3d_asm_Np."+yr+mn+".SUB.nc4")
        ; print(fili)
        f = addfile(fili, "r")
        t_med = f->T(0,{LEV},:,:)
        t_med = lonFlip(t_med)
        t(i*12+j,:,:,:) = linint2_Wrap(t_med&lon, t_med&lat, t_med, True, lon, lat, 0)
    end do
    print("end of year "+year(i))
end do
t!0 = "time"
t&time = ispan(1, nyear*nmonth, 1)
t!1 = "lev"
t!2 = "lat"
t!3 = "lon"
printVarSummary(t)

TS_obsx = month_to_season(t, season) ;clmMonTLL(t)
TS_obs = dim_avg_n_Wrap(TS_obsx, 0)
TS_obs_profile = TS_obs(:,{LAT},:)

; topox = linint2_Wrap(topo&lon, topo&lat, topo, False, t&lon, t&lat, 0)
; print(topox({35},{90}))

; TS_obs = mask(TS_obs, topox.ge.1500, True)
; TP_obs = wgt_areaave_Wrap(TS_obs(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
; printVarSummary(TP_obs)

;===================================================================
; GCM 20 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/ta/"
dirigcm2 = "/home/yangsong3/data-model/wuyt/CMIP5/"


models1 = (/"ACCESS1-0","bcc-csm1-1","BNU-ESM","CanCM4","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CESM","CNRM-CM5","GFDL-ESM2M","GISS-E2-H","GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR",\
        "IPSL-CM5B-LR","MIROC5","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
;"CESM1-WACCM","FIO-ESM","MPI-ESM-P","CSIRO-Mk3-6-0",
;"CMCC-CM","CMCC-CMS","GFDL-CM3","GFDL-ESM2G","IPSL-CM5A-MR",,"NorESM1-ME"

models2 = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

models = venn2_intersection(models1,models2)

nm = dimsizes(models)
print(models)

TP_GCM5 = new((/nm,nlev,nlon/), float)
TP_GCM5!0 = "model"
TP_GCM5&model = models
TP_usea = TP_GCM5

do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"/ta_Amon_"+models(i)+"_historical_r1i1p1_1*-20*12.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    T_GCMs := fgcm->ta(iYYYY,{LEV2},:,:)
    ; printVarSummary(T_GCMs)
    T_GCMs!1 = "lev"
    T_GCMs!2 = "lat"
    T_GCMs!3 = "lon"

    T_GCM := linint2_Wrap(T_GCMs&lon, T_GCMs&lat, T_GCMs, True, t&lon, t&lat, 0)
    printVarSummary(T_GCM)

    t_gcm_clmx = month_to_season(T_GCM, season)
    t_gcm_clm = dim_avg_n_Wrap(t_gcm_clmx, 0)
    TP_GCM5(i,:,:) = t_gcm_clm(:,{LAT},:)
    print("end of "+models(i))
    delete([/T_GCMs,T_GCM,t_gcm_clmx,t_gcm_clm/])


    filiu = systemfunc("ls "+dirigcm2+models(i)+"/ua_Amon_"+models(i)+"_historical_r1i1p1_*.nc")
    print(filiu)
    fu = addfile(filiu, "r")
    TIMEu = fu->time
    YYYYu = cd_calendar(TIMEu, -1)/100
    iYYYYu = ind(YYYYu.ge.yrStrt .and. YYYYu.le.yrLast)
    Us = fu->ua(iYYYYu,{LEV2},:,:)
    U = linint2_Wrap(Us&lon, Us&lat, Us, True, lon, lat, 0)
    u_seax = month_to_season(U, season)
    u_sea = dim_avg_n_Wrap(u_seax, 0)
    TP_usea(i,:,:) = u_sea(:,{LAT},:)
    print("end of "+models(i))
    print(u_sea({20000},{LAT},{lonstr:lonend}))

    delete([/TIMEu,YYYYu,iYYYYu,Us,u_seax,u_sea/])

end do
printVarSummary(TP_GCM5)
printVarSummary(TP_usea)

diff = TP_GCM5
do i = 0, nm-1
    diff(i,:,:) = TP_GCM5(i,:,:) - TS_obs_profile
end do
printVarSummary(diff)

;=========================================================
; outline of tp
;=========================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
topo_z = topo({LAT},:)
topo_zonal = 1013.25*(1-topo_z*0.0065/288.15)^5.25145
copy_VarCoords(topo({LAT},:), topo_zonal)
;========================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/20201228/"
figname = "T_profile_bias_"+season
; pwks = "png"
; pwks@wkHeight              = 3048
; pwks@wkWidth               = 6048
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(nm, graphic)
plot_topo = plot
plotc = plot

;=====================================================
    res_topo                   = True

    res_topo@gsnDraw  = False           ;先别直接画出来
    res_topo@gsnFrame = False           ;也不要翻页

    res_topo@gsnLeftString  = " "              ;不要标示出变量名和单位等信息
    res_topo@gsnRightString = " "

    ;res_topo@gsnMaximize        = True
    ;res_topo@trYMaxF            = 1000.0
    res_topo@gsnYRefLine        = 1000 ;(/100,1000/)          ; create a reference line 
    res_topo@gsnBelowYRefLineColor = "black" ;(/"black","black"/) ; above ref line fill grey
    ;res_topo@gsnXYBelowFillColors = "goldenrod2"
    res_topo@gsnAboveYRefLineColor = "blue"
    res_topo@xyLineColors = "black"
    res_topo@xyLineThicknesses = 4.
;===================================================================
;===================================================================

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@vpHeightF = 0.4
    res@vpWidthF  = 0.7

    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@gsnStringFontHeightF = 0.03
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02

    res@pmLabelBarOrthogonalPosF = -0.06
    res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@tiYAxisString = ""

    res@tmYLMode      = "Explicit"
    res@tmYLValues    = (/"1000","850","750","500","300","250","150","100"/) ;
    res@tmYLLabels    = (/"1000","850","750","500","300","250","150","100"/) ;

    res@trYMinF   = 100
    res@trYMaxF   = 1000
    res@gsnPresHgtHeightLabelOn = False
    res@tmYROn = False

;=====================================================

;----------------------------------------------
; draw prob
;----------------------------------------------
 res3  = True
 res3@gsnDraw   = False
 res3@gsnFrame  = False

 ;res3@cnFillOn  = False 
 res3@cnLinesOn = False
 res3@cnLineLabelsOn = False
 res3@cnInfoLabelOn  = False
 ;res3@lbLabelBarOn   = False

 res3@gsnRightString = ""

 res3@cnLevelSelectionMode = "ManualLevels"
 res3@cnMinLevelValF = 0.00
 res3@cnMaxLevelValF = 0.10
 res3@cnLevelSpacingF = 0.01

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------

 ; ====================================================================================
    cnres = True
    cnres@gsnDraw  = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = False
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "purple"
    cnres@cnLineThicknessF = 4
    cnres@cnInfoLabelOn = True
    cnres@cnInfoLabelOrthogonalPosF = -0.25
    cnres@cnInfoLabelFontHeightF = 0.02
    cnres@cnLabelDrawOrder = "PostDraw"
    cnres@gsnContourNegLineDashPattern = 2
    cnres@gsnContourPosLineDashPattern = 0
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""
    cnres@tmYLMode      = "Explicit"
    cnres@tmYLValues    = (/"1000","850","750","500","300","250","150","100"/) ;
    cnres@tmYLLabels    = (/"1000","850","750","500","300","250","150","100"/) ;

    cnres@trYMinF   = 100
    cnres@trYMaxF   = 1000
    cnres@gsnPresHgtHeightLabelOn = False
    cnres@tmYROn = False
;====================================================================================
 ;----------------------------------------------

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 10.
res@cnMinLevelValF = -10.
res@cnFillPalette = "temp_diff_18lev" ; "MPL_RdBu"

cnres@cnLevelSelectionMode = "ManualLevels"
cnres@cnMinLevelValF = 20
cnres@cnMaxLevelValF = 50
cnres@cnLevelSpacingF = 10

do i = 0, nm-1
    res@gsnLeftString = models(i)
    plot(i) = gsn_csm_pres_hgt(wks, diff(i,:,{lonstr:lonend}), res)
    plot_topo(i) = gsn_csm_xy(wks, topo&lon, topo_zonal, res_topo)
    overlay(plot(i), plot_topo(i))
    plotc(i) = gsn_csm_pres_hgt(wks, TP_usea(i,:,{lonstr:lonend}), cnres)
    overlay(plot(i), plotc(i))
    ; plot_prob(i) = gsn_csm_contour(wks,prob(i,:,{lonstr:lonend}),res3)
    ; plot_prob(i) = ShadeLtContour(plot_prob(i),0.055,17)
    ; overlay(plot(i), plot_prob(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar  = True
pres@gsnPanelMainString = season+" mean: T profiles over "+LAT+"~S~o~N~N"
gsn_panel(wks, plot, (/6,5/), pres)

end
