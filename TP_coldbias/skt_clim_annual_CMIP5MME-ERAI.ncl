
begin

;===================================================================
; GCM 24 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/ts/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TS_GCM = new((/nm,12,181,360/), float)
TS_GCM!0 = "model"
TS_GCM&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"ts_12mon_"+models(i)+"_historical_r1i1p1_195101-200012.nc")
    fgcm = addfile(filigcm, "r")
    TS_GCM(i,:,:,:) = fgcm->ts
end do
printVarSummary(TS_GCM)
ts_MME = dim_avg_n_Wrap(TS_GCM, (/0,1/))
printVarSummary(ts_MME)

;===================================================================
; Reanalysis ERAI data
;===================================================================
diriera = "/home/yangsong3/data-observation/wuyt/ERAI/monthly/1x1/"
;filir = "t2m.mon.mean.7901.1804.nc"
filiera = "skt.mon.mean.7917.nc"
fera = addfile(diriera+filiera, "r")
; timer = fera->time
; TIME = cd_calendar(timer, -1)
; print(TIME)
;TS_era = short2flt(fr->ts(0:467,:,:))
TS_era = short2flt(fera->skt(:,::-1,:))
printVarSummary(TS_era)

ts_era = dim_avg_n_Wrap(TS_era, 0)
printVarSummary(ts_era)

;===================================================================
; Minus
;===================================================================
eramme = ts_MME - ts_era
copy_VarCoords(ts_era, eramme)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "skt_clim_annual_CMIP5MME-ERAI"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"

res@mpCenterLonF = 180
res@mpMinLatF = -30
res@mpMaxLatF = 50

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

;res@lbLabelBarOn = False
res@pmLabelBarOrthogonalPosF = 0.25
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 4.
res@cnMinLevelValF = -4.
res@cnFillPalette = "CBR_coldhot"
res@gsnRightString = "K"

res@tiMainString = "Clim Ts Bias"; (MME - ERAI)
res@tiMainFontHeightF = 0.02
plot(0) = gsn_csm_contour_map(wks, eramme, res)

draw(plot)
frame(wks)
end