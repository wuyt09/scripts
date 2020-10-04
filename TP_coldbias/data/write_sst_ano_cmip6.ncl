
begin
;=========================================================
; Note: the changes has been add as annual cycle, 
;       each year has the same anomalies.
;       but I still confused by how the anomaly data derived.
;                            2019-9-27 Wu Yuting
;=========================================================

fili_ctrl = "sst_HadOIBl_bc_1.9x2.5_1850_2012_c130411.nc"
f_ctrl = addfile(fili_ctrl, "r")

sst = f_ctrl->SST_cpl
printVarSummary(sst)

lat = sst&lat
lon = sst&lon

;==============================================================
fili_landsea = "/home/yangsong3/soft/ncl650/lib/ncarg/data/cdf/landsea.nc"
flandsea = addfile(fili_landsea, "r")
LANDSEA = flandsea->LSMASK
landsea = linint2_Wrap(LANDSEA&lon, LANDSEA&lat, LANDSEA, True, lon, lat, 0)
printVarSummary(landsea)
; sst = mask(sst, landsea.eq.0, True)

;=======================================================
;read hist bias
;=======================================================
dirisst = "/home/yangsong3/wuyt/sysu/scripts/TP_coldbias/model/a20200910/"
filisst = "cmip6_hist_sst_bias_1950_1999_mn_mme.dat"
; setfileoption("bin", "ReadByteOrder", "Native")

sst_A = new((/600,89,180/), float)
lato = fspan(-88, 88, 89)
lono = fspan(0, 358, 180)
lato@units = "degrees_north"
lono@units = "degrees_east"
do i = 0, 599
    sst_A(i,:,:) = fbindirread(dirisst+filisst, i, (/89,180/), "float")
end do
sst_A!0 = "month"
sst_A!1 = "lat"
sst_A!2 = "lon"
sst_A&lat = lato
sst_A&lon = lono
sst_A@_FillValue = 9.96921e+36
printVarSummary(sst_A)

if(any(ismissing(sst_A))) then
    print("Your data contains some missing values. Beware.")
end if
sst_ano = clmMonTLL(sst_A)
;=======================================================

;=======================================================
; add bias to original sst
;=======================================================
sstA = linint2_Wrap(lono, lato, sst_ano, True, lon, lat, 0)
;print(sstA(11,{-30:30},{0}))
sstA = where(ismissing(sstA), 0, sstA) ;! careful, original sst has no FillValue
;sstA = where(sstA.eq.sstA@_FillValue, 0, sstA)  ;! this not works
;print(sstA(11,{-30:30},{0}))
printVarSummary(sstA)

; X_ano = wgt_areaave_Wrap(sst_ano(:,{-30:30},:), 1, 1, 0) ;sst_ano
; print(X_ano)
; x_ano = avg(X_ano)
; print(x_ano)  ; -0.1923212  ;-0.1434768 ;-0.1434765
; exit

hist_sst = sst
nyear = 2012 - 1850 + 1
do i = 0, nyear-1
    do j = 0, 11
        hist_sst(12*i+j,{-30:30},:) = sst(12*i+j,{-30:30},:)+sstA(j,{-30:30},:)
    end do
end do
printVarSummary(hist_sst)

;=======================================================
; check if it works
;=======================================================
test_o = sst(240:251,:,:) ;clmMonTLL(sst)
test_o = mask(test_o, landsea.eq.0, True)
test_h = hist_sst(240:251,:,:) ;clmMonTLL(hist_sst)
test_h = mask(test_h, landsea.eq.0, True)
test_o30 = wgt_areaave_Wrap(test_o(:,{-30:30},:), 1, 1, 0)
test_h30 = wgt_areaave_Wrap(test_h(:,{-30:30},:), 1, 1, 0)
diff_30 = test_h30 - test_o30
; print("ctrl "+test_o30)
; print("hist "+test_h30)
sstA = mask(sstA, landsea.eq.0, True)
ssta30 = wgt_areaave_Wrap(sstA(:,{-30:30},:), 1, 1, 0)
print("diff "+diff_30+" sstA "+ssta30)
print(avg(ssta30))

;=======================================================
; write new sst files
;=======================================================
date = systemfunc("date +%Y%m%d")
Hist_fili = "Hist_bias_sst_cmip6_bc_1.9x2.5_c"+date+".nc"
system("rm "+Hist_fili)
system("cp "+fili_ctrl+" "+Hist_fili)
f_hist = addfile(Hist_fili, "w")
f_hist->SST_cpl = hist_sst

;=======================================================
; plot 
;=======================================================
wks = gsn_open_wks("pdf", "write_sst_ano_cmip6")
plot = new(3, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"

res@mpCenterLonF = 180
res@mpMinLatF = -30
res@mpMaxLatF = 30

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.015
res@tmXBLabelFontHeightF = 0.015
res@tmYLLabelFontHeightF = 0.015

res@lbLabelBarOn = True
res@pmLabelBarOrthogonalPosF = 0.25
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 30
res@cnMinLevelValF = 20
res@cnFillPalette = "CBR_coldhot"
res@gsnRightString = "K"
res@gsnStringFontHeightF = 0.025

res@gsnLeftString = "Ctrl run"
SSTX = month_to_annual(sst, 1)
SST = dim_avg_n_Wrap(SSTX, 0)
plot(0) = gsn_csm_contour_map(wks, SST, res)

res@gsnLeftString = "Hist run"
hist_SSTX = month_to_annual(hist_sst, 1)
hist_SST = dim_avg_n_Wrap(hist_SSTX, 0)
plot(1) = gsn_csm_contour_map(wks, hist_SST, res)

res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 1
res@cnMinLevelValF = -1
res@gsnLeftString = "Hist_bias"
sst_ANOX = month_to_annual(sst_ano, 1)
sst_ANO = dim_avg_n_Wrap(sst_ANOX, 0)
plot(2) = gsn_csm_contour_map(wks, sst_ANO, res)

pres = True
pres@gsnMaximize = True
pres@gsnDraw = True
pres@gsnFrame = True

gsn_panel(wks, plot, (/3,1/), pres)
end