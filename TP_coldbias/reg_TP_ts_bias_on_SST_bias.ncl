begin

lonstr = 75
lonend = 105
latstr = 25
latend = 40

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
TP_GCM = wgt_areaave(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
TP_GCMA = dim_avg_n_Wrap(TP_GCM, 0)
printVarSummary(TP_GCMA)

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
;print(TS_era&latitude)
ts_era = clmMonTLL(TS_era)
printVarSummary(ts_era)
TP_era = wgt_areaave(ts_era(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_era)

ts_erabias = TP_era - TP_GCMA
printVarSummary(ts_erabias)

;===================================================================
; Reanalysis NCEP_DOE data
;===================================================================
dirin2 = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
;filir = "t2m.mon.mean.7901.1804.nc"
filin2 = "skt.sfc.mon.mean.nc"
fn2 = addfile(dirin2+filin2, "r")
timen2 = fn2->time
YYYY = cd_calendar(timen2, -1)/100
iYYYY = ind(YYYY.ge.1979 .and. YYYY.le.2018)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_n2 = short2flt(fn2->skt(iYYYY,::-1,:))
printVarSummary(TS_n2)
ts_n2 = clmMonTLL(TS_n2)
printVarSummary(ts_n2)
TP_n2 = wgt_areaave(ts_n2(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_n2)

ts_n2bias = TP_n2 - TP_GCMA
printVarSummary(ts_n2bias)

;===================================================================
; Reanalysis NCEP_NCAR data
;===================================================================
diriN1 = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
filiN1 = "skt.sfc.mon.mean.nc"
fN1 = addfile(diriN1+filiN1, "r")
timeN1 = fN1->time
YYYYN1 = cd_calendar(timeN1, -1)/100
iYYYYN1 = ind(YYYYN1.ge.1948 .and. YYYYN1.le.2018)
;TS_era = short2flt(fr->ts(iYYYY,:,:))
TS_N1 = short2flt(fN1->skt(iYYYYN1,::-1,:))
ts_N1 = clmMonTLL(TS_N1)
ts_N1 = ts_N1 + 273.15
printVarSummary(ts_N1)
TP_N1 = wgt_areaave(ts_N1(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0) 
printVarSummary(TP_N1)

ts_n1bias = TP_N1 - TP_GCMA
printVarSummary(ts_n1bias)




ts_bias = new((/3,12/), float)
ts_bias(0,:) = ts_erabias
ts_bias(1,:) = ts_n1bias
ts_bias(2,:) = ts_n2bias


;===================================================================
; SST historical bias data
;===================================================================

dirisst = "/home/yangsong3/huxm/ENSO_decadal/"
; filicor = "sst_HadOIBl_bc_1.9x2.5_1979_2010_clim.nc"
; fcor = addfile(dirisst+filicor, "r")
; latsst = fcor->lat
; lonsst = fcor->lon
filisst = "cmip5_hist_sst_mn_clm_bias_mme.dat"
setfileoption("bin", "ReadByteOrder", "Native")

sst_ano = new((/12,91,180/), float)
lato = fspan(-90, 90, 91)
lono = fspan(0, 358, 180)
lato@units = "degrees_north"
lono@units = "degrees_east"
do i = 0, 11
    sst_ano(i,:,:) = fbindirread(dirisst+filisst, i, (/91,180/), "float")
end do
sst_ano!0 = "month"
sst_ano!1 = "lat"
sst_ano!2 = "lon"
sst_ano&lat = lato
sst_ano&lon = lono
sst_ano@_FillValue = 1e+20
printVarSummary(sst_ano)

;===================================================================
; reg TP ts bias on historical SST bias
;===================================================================

rc = regCoef_n(ts_bias, sst_ano, 1, 0)
copy_VarCoords(sst_ano(0,:,:), rc(0,:,:))
printVarSummary(rc)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "reg_TP_ts_bias_on_SST_bias"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(3, graphic)

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

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.2
res@cnMaxLevelValF = 2.
res@cnMinLevelValF = -2.
res@cnFillPalette = "CBR_coldhot"
res@gsnRightString = "K"

res@gsnLeftString = "ERAI ts bias"
plot(0) = gsn_csm_contour_map(wks, rc(0,:,:), res)

res@gsnLeftString = "NCEP_NCAR ts bias"
plot(1) = gsn_csm_contour_map(wks, rc(1,:,:), res)

res@gsnLeftString = "NCEP_DOE ts bias"
plot(2) = gsn_csm_contour_map(wks, rc(2,:,:), res)

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,1/), pres)

end