
begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; GCM 24 data
;===================================================================
dirigcm = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/CMIP5/post-pro/"

models = (/"ACCESS1-0","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\  ;"ACCESS1-3",
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\ ;"HadGEM2-CC","HadGEM2-ES",
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)

season = "DJF"
lon = fspan(0, 359, 360)
lat = fspan(-90, 90, 181)
lev = (/1000, 925, 850, 700, 600, 500, 400, 300/)
lev@units = "hPa"

TS_GCM = new((/nm,8,181,360/), float)
TS_GCM!0 = "model"
TS_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"hus_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TSX = fgcm->hus
    TSX_sea = month_to_season(TSX, season)
    TS_GCM(i,:,:,:) = dim_avg_n_Wrap(TSX_sea, 0)
    print("end of "+models(i))
end do
printVarSummary(TS_GCM)
TS_GCMA = dim_avg_n_Wrap(TS_GCM, 0)
printVarSummary(TS_GCMA)
TS_GCMA!0 = "lev"
TS_GCMA&lev = lev
TS_GCMA!1 = "lat"
TS_GCMA&lat = lat
TS_GCMA!2 = "lon"
TS_GCMA&lon = lon
printVarSummary(TS_GCMA)

;===================================================================
; NCEP data
;===================================================================

yrStrt = 1951
yrLast = 2000
fili = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/shum.mon.mean.nc"
f = addfile(fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

QY = f->shum(iYYYY,:,::-1,:)
QY_sea = month_to_season(QY, season)
Q_ncep = dim_avg_n_Wrap(QY_sea, 0)
printVarSummary(Q_ncep)

TS_GCMa = linint2_Wrap(TS_GCMA&lon, TS_GCMA&lat, TS_GCMA, True, Q_ncep&lon, Q_ncep&lat, 0)
printVarSummary(TS_GCMa)
Q_diff = Q_ncep
q_diff = TS_GCMa - Q_ncep ;conform_dims(dimsizes(TS_GCMa), Q_ncep, (/1,2,3/))
copy_VarCoords(Q_ncep, q_diff)
printVarSummary(q_diff)

; Q_X = new((/17,360/), float)
; Q_Y = new((/17,181/), float)
; q_diff = TS_GCMA
Q_X = new((/8,144/), float)
Q_Y = new((/8,73/), float)
Q_X = dim_avg_n_Wrap(q_diff(:,{latstr:latend},:), 1)
Q_Y = dim_avg_n_Wrap(q_diff(:,:,{lonstr:lonend}), 2)
printVarSummary(Q_X)


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
figname = "CMIP5_"+season+"_wv_airT_profile"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(2, graphic)

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
    res@tmXBLabelFontHeightF = 0.025
    res@tmYLLabelFontHeightF = 0.025

    ;res@pmLabelBarOrthogonalPosF = 0.12
    ;res@lbLabelBarOn = False
    ;res@lbLabelFontHeightF = 0.02
    res@cnFillPalette = "BlueWhiteOrangeRed"
    res@tiYAxisString = ""
    res@gsnPresHgtHeightLabelOn = False

 ;------------------------------------------
     res@trYMinF   = 300
     res@trYMaxF   = 1000
      
     res@tmYLMode      = "Explicit"
     res@tmYLValues    = (/"1000","850","750","500","300","200","100"/)
     res@tmYLLabels    = (/"1000","850","750","500","300","200","100"/)

;=====================================================

    cnres = True
    cnres@gsnDraw = False
    cnres@gsnFrame = False

    cnres@cnFillOn = False
    cnres@cnLinesOn = True
    cnres@cnLineLabelsOn = True
    cnres@cnLineLabelFontColor = "red"
    cnres@cnLineColor = "grey50"
    cnres@cnLineThicknessF = 5
    cnres@cnInfoLabelOn = False
    ;cnres@cnLineDrawOrder = "PostDraw"
    cnres@cnLabelDrawOrder = "PostDraw"
    ;cnres@cnLevelFlag = "LineAndLabel"
    cnres@cnLineLabelInterval = 1
    cnres@gsnLeftString = ""
    cnres@gsnRightString = ""
;-----------------------------------------------------

cnres@gsnContourPosLineDashPattern = 0
cnres@gsnContourNegLineDashPattern = 2
;-----------------------------------------------------

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
 ;----------------------------------------------

;=============================================================================
res@gsnLeftString = "Zonal mean ("+latstr+"-"+latend+"~S~o~N~N)"; (34~S~o~N~N)
res@gsnRightString = "Diff"
res@cnFillPalette = "CBR_drywet"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1
res@cnMaxLevelValF = 10
res@cnMinLevelValF = -10

plot(0) = gsn_csm_pres_hgt(wks, Q_X(:,{60:120}), res)



res@gsnLeftString = "Meridional mean ("+lonstr+"-"+lonend+"~S~o~N~E)";
res@gsnRightString = "Diff"
; res@cnFillPalette = "CBR_drywet"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 0.0001
; res@cnMaxLevelValF = 0.001
; res@cnMinLevelValF = -0.001

plot(1) = gsn_csm_pres_hgt(wks, Q_Y(:,{-20:60}), res)

 

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "CMIP5 ensmble mean Q diff .Vs. NCEP(Kg/Kg, shading)"
gsn_panel(wks, plot, (/2,1/), pres)

end
