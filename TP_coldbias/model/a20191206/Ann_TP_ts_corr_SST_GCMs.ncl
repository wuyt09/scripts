
begin

yrStrt = 1979
yrLast = 2017
nyear = yrLast - yrStrt + 1

lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
;===================================================================
; GCM 24 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/tas/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

; models = (/"bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
;         "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
;         "GISS-E2-H","GISS-E2-R","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
;         "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TS_GCM = new((/nm,50,181,360/), float)
TS_GCM!0 = "model"
TS_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"tas_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TS_GCMs = fgcm->tas
    TS_GCMs = dtrend_msg_n(TS_GCMs&time, TS_GCMs, False, False, 0)
    TS_GCM(i,:,:,:) = month_to_annual(TS_GCMs, 1)
end do
printVarSummary(TS_GCM)
;TP_GCMX = TS_GCM
TS_GCM = mask(TS_GCM, topo.ge.1500, True)
TP_GCM = wgt_areaave_Wrap(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_GCM)

;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/ts/"
SST_GCM = new((/nm,50,181,360/), float)
SST_GCM!0 = "model"
SST_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"ts_mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    SST_GCMs = fgcm->ts
    SST_GCMs = dtrend_msg_n(SST_GCMs&time, SST_GCMs, False, False, 0)
    SST_GCM(i,:,:,:) = month_to_annual(SST_GCMs, 1)
end do
printVarSummary(SST_GCM)

;===================================================================
; land sea mask
;===================================================================
a      = addfile ("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata, SST_GCM&lat, SST_GCM&lon)
;lmask = where(lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, 0, 1)
SST_GCM = mask(SST_GCM, lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, False)
;===================================================================
; regression data
;===================================================================
cc = new((/nm,181,360/), float)
prob = cc
do i = 0, nm-1
    cc(i,:,:) = escorc_n(TP_GCM(i,:), SST_GCM(i,:,:,:), 0, 0)
    prob(i,:,:) = rtest(cc(i,:,:), nm, 0)
end do
copy_VarCoords(SST_GCM(:,0,:,:), cc)
copy_VarCoords(SST_GCM(:,0,:,:), prob)
prob = mask(prob, cc.eq.0, False)
printVarSummary(cc)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20191206/"
filiname = "Ann_TP_ts_corr_SST_GCMs"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(nm, graphic)
plot_prob = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.8
res@cnMinLevelValF = -0.8
res@cnFillPalette = "CBR_coldhot"
;res@gsnRightString = "K"

;------------------------
; draw prob
;------------------------
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

;res3@cnFillPattern = 17 ;(/17/)
res3@cnFillColor  =  "black" ;"green4"
res3@gsnLeftString = ""
res3@cnFillDotSizeF = 0.003
res3@cnFillScaleF  = 1 ;1.2

sig = 0.95
res@gsnRightString = sig*100+"%"

do i = 0, nm-1
    res@gsnLeftString = models(i)
    plot(i) = gsn_csm_contour_map(wks, cc(i,:,:), res)
    plot_prob(i) = gsn_csm_contour(wks, prob(i,:,:), res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),(1-sig),6)
    overlay(plot(i), plot_prob(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Correlation coefficient of annual mean: model Ts & model SST "
gsn_panel(wks, plot, (/7, 4/), pres)

cc_mme = dim_avg_n_Wrap(cc, 0)
prob_mme = dim_avg_n_Wrap(prob, 0)
res@gsnLeftString = "MME"
res@lbLabelBarOn = True
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 0.3
res@cnMinLevelValF = -0.3
plot(0) = gsn_csm_contour_map(wks, cc_mme, res)
plot_prob(0) = gsn_csm_contour(wks, prob_mme, res3)
plot_prob(0) = ShadeLtContour(plot_prob(0),(1-sig),6)
overlay(plot(0), plot_prob(0))
draw(plot(0))
frame(wks)
end