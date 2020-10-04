
begin
;=========================================================
; Note: the changes has been add as annual cycle, 
;       each year has the same anomalies.
;       but I still confused by how the anomaly data derived.
;                            2019-9-27 Wu Yuting
;=========================================================


yrStrt = 1850
yrLast = 2012
year = ispan(yrStrt, yrLast, 1)


;fili_ctrl = "sst_HadOIBl_bc_1.9x2.5_1850_2013_c140701.nc"
fili_ctrl = "sst_HadOIBl_bc_1.9x2.5_1850_2012_c130411.nc"
; fili_hist = "Hist_bias_sst_cmip5_bc_1.9x2.5_c20191204.nc"
fili_hist = "Hist_bias_sst_cmip6_bc_1.9x2.5_c20201004.nc"

f_ctrl = addfile(fili_ctrl, "r")
f_hist = addfile(fili_hist, "r")

sst1 = f_ctrl->SST_cpl
sst2 = f_hist->SST_cpl
printVarSummary(sst1)
printVarSummary(sst2)

; ;==============================================================
; fili_landsea = "/home/yangsong3/soft/ncl650/lib/ncarg/data/cdf/landsea.nc"
; flandsea = addfile(fili_landsea, "r")
; LANDSEA = flandsea->LSMASK
; landsea = linint2_Wrap(LANDSEA&lon, LANDSEA&lat, LANDSEA, True, sst1&lon, sst1&lat, 0)
; printVarSummary(landsea)
; sst1 = mask(sst1, landsea.eq.0, True)
; sst2 = mask(sst2, landsea.eq.0, True)

SST1 = month_to_annual(sst1, 1)
SST2 = month_to_annual(sst2, 1)
printVarSummary(SST1)

x1 = wgt_areaave_Wrap(SST1(:,{-30:30},:), 1, 1, 0)
x2 = wgt_areaave_Wrap(SST2(:,{-30:30},:), 1, 1, 0)
printVarSummary(x1)

pp = new((/2,163/), float)
pp(0,:) = x1
pp(1,:) = x2
pp!1 = "year"
pp&year = year
copy_VarCoords(x1, pp(0,:))
printVarSummary(pp)
test = pp(1,:) - pp(0,:)
print(test) ;-0.1434765
exit
;=======================================================
; read hist bias
;=======================================================
dirisst = "/home/yangsong3/huxm/ENSO_decadal/"
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
; print(sst_ano(:,{-30:30},:))
; exit

X_ano = wgt_areaave_Wrap(sst_ano(:,{-30:30},:), 1, 1, 0)
print(X_ano)
x_ano = avg(X_ano)
print(x_ano)  ; -0.1923212

;=======================================================
; plot 
;=======================================================
wks = gsn_open_wks("pdf", "check_sst_diff_files")
plot = new(4, graphic)

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
; res@tmXBLabels = Month
; res@tmXTOn = False
; res@tiXAxisString = "Month"

res@xyLineThicknessF = 5
res@xyLineColors = (/"black", "blue", "red", "green4"/)
;(/"CTRL", "HIST", "FUTR", "FUBIAS"/)
res@xyDashPatterns = (/0,2/)

plot(0) = gsn_csm_xy(wks, year, pp, res)
draw(plot(0))
frame(wks)


end