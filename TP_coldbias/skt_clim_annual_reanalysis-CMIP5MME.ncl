
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
TS_N1 = TS_N1 + 273.15
ts_N1 = dim_avg_n_Wrap(TS_N1, 0)
printVarSummary(ts_N1)
ts_n1 = linint2_Wrap(ts_N1&lon, ts_N1&lat, ts_N1, True, TS_GCM&lon, TS_GCM&lat, 0)
printVarSummary(ts_n1)

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
ts_N2 = dim_avg_n_Wrap(TS_n2, 0)
printVarSummary(ts_N2)
ts_n2 = linint2_Wrap(ts_N2&lon, ts_N2&lat, ts_N2, True, TS_GCM&lon, TS_GCM&lat, 0)
printVarSummary(ts_n2)

;===================================================================
; Minus
;===================================================================
eramme = ts_era - ts_MME
n1mme = ts_n1 - ts_MME
n2mme = ts_n2 - ts_MME
copy_VarCoords(ts_era, eramme)
copy_VarCoords(ts_n1, n1mme)
copy_VarCoords(ts_n2, n2mme)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "skt_clim_annual_reanalysis-CMIP5MME"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(3, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False

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
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 5.
res@cnMinLevelValF = -5.
res@cnFillPalette = "CBR_coldhot"
res@gsnRightString = "K"

res@gsnLeftString = "ERAI(1979-2017) - MME(1951-2000)"
plot(0) = gsn_csm_contour_map(wks, eramme, res)

res@gsnLeftString = "NCEP_NCAR(1948-2018) - MME(1951-2000)"
plot(1) = gsn_csm_contour_map(wks, n1mme, res)

res@gsnLeftString = "NCEP_DOE(1979-2018) - MME(1951-2000)"
plot(2) = gsn_csm_contour_map(wks, n2mme, res)


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,1/), pres)
end