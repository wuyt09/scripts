
begin

yrStrt = 1979
yrLast = 2017
nyear = yrLast - yrStrt + 1

lonstr = 75
lonend = 105
latstr = 25
latend = 40

;===================================================================
; GCM 24 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/tas/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/);,"NorESM1-ME"

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
TP_GCM = wgt_areaave_Wrap(TS_GCM(:,:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_GCM)

;===================================================================
; regression data
;===================================================================
cc = new((/nm,181,360/), float)
prob = cc
Nr = dimsizes(TS_GCM(0,:,0,0))
do i = 0, nm-1
    cc(i,:,:) = escorc_n(TP_GCM(i,:), TS_GCM(i,:,:,:), 0, 0)
    prob(i,:,:) = rtest(cc(i,:,:), Nr, 0)
end do
prob = mask(prob, cc.eq.0, False)
copy_VarCoords(TS_GCM(:,0,:,:), cc)
copy_VarCoords(TS_GCM(:,0,:,:), prob)
printVarSummary(cc)

;===================================================================
; land sea mask
;===================================================================
a      = addfile ("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm = landsea_mask(lsdata, cc&lat, cc&lon)
;lmask = where(lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, 0, 1)
cc = mask(cc, lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, False)

;===================================================================
; EOFs
;===================================================================
neval = 2
cc_eof = eofunc(cc(lat|:,lon|:,model|:), neval, False)
copy_VarCoords(cc(0,:,:), cc_eof(0,:,:))
printVarSummary(cc_eof)

;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/"
filiname = "TP_tas_corr_SST_GCMs_EOF"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(2, graphic)

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

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

;res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.002
res@cnMaxLevelValF = 0.02
res@cnMinLevelValF = -0.02
res@cnFillPalette = "CBR_coldhot"
;res@gsnRightString = "K"

res@gsnLeftString = "EOF 1"
res@gsnRightString = cc_eof@pcvar(0)+"%"
plot(0) = gsn_csm_contour_map(wks, cc_eof(0,:,:), res)

res@cnLevelSpacingF = 0.002
res@cnMaxLevelValF = 0.02
res@cnMinLevelValF = -0.02
res@gsnLeftString = "EOF 2"
res@gsnRightString = cc_eof@pcvar(1)+"%"
plot(1) = gsn_csm_contour_map(wks, cancor(  , y, option)c_eof(1,:,:), res)

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "EOFs: correlation of Annual mean model Tas & model SST "
gsn_panel(wks, plot, (/2, 1/), pres)


end