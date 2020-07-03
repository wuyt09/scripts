
begin
;=========================================================
; Note: the changes has been add as annual cycle, 
;       each year has the same anomalies.
;       but I still confused by how the anomaly data derived.
;                            2019-9-27 Wu Yuting
;=========================================================


yrStrt = 1850
yrLast = 2013
year = ispan(yrStrt, yrLast, 1)


fili_ctrl = "sst_HadOIBl_bc_1.9x2.5_1850_2013_c140701.nc"
;fili_ctrl = "../sst_HadOIBl_bc_1.9x2.5_1850_2012_c130411.nc"
fili_hist = "sst_cmip5_bc_1.9x2.5_1950_1999.nc"
fili_futr = "sst_cmip5_bc_1.9x2.5_2050_2099.nc"
fili_fubias = "sst_cmip5_bc_1.9x2.5_2050_2099_clm.nc"

f_ctrl = addfile(fili_ctrl, "r")
f_hist = addfile(fili_hist, "r")
f_futr = addfile(fili_futr, "r")
f_fubias = addfile(fili_fubias, "r")

sst1 = f_ctrl->SST_cpl
sst2 = f_hist->SST_cpl
sst3 = f_futr->SST_cpl
sst4 = f_fubias->SST_cpl
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
SST3 = month_to_annual(sst3, 1)
SST4 = month_to_annual(sst4, 1)
printVarSummary(SST1)

x1 = wgt_areaave_Wrap(SST1(:,{-30:30},:), 1, 1, 0)
x2 = wgt_areaave_Wrap(SST2(:,{-30:30},:), 1, 1, 0)
x3 = wgt_areaave_Wrap(SST3(:,{-30:30},:), 1, 1, 0)
x4 = wgt_areaave_Wrap(SST4(:,{-30:30},:), 1, 1, 0)
printVarSummary(x1)

pp = new((/4,164/), float)
pp(0,:) = x1
pp(1,:) = x2
pp(2,:) = x3
pp(3,:) = x4
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
exit


; ;===================================================================
; ; SST historical bias data
; ;===================================================================
; yrStrtmodel1 = 1950
; yrLastmodel1 = 1999
; nm = 28

; dirh="ls /home/yangsong3/jyt/ensocro/sst/ssth/"
; ;read the name of all of models
; dirhname = asciiread("/home/yangsong3/jyt/name/dirssth28.txt",(/28/),"string")

; ssth = new((/600,91,180/), float)
; ssthdt = new((/nm,600,91,180/), float)
; ssthclm = new((/nm,12,91,180/), float)
; time = ispan(1, 600, 1)

; do i = 0, nm-1
;     ;==============read SST==============================
;     filh = systemfunc(dirh+dirhname(i))
;     ssth = fbindirread(filh, 0, (/600,91,180/), "float")
;     ssth!0 = "time"
;     ssth!1 = "lat"
;     ssth!2 = "lon"
;     ssth&time = time
;     ssth&lat  = fspan(-90, 90, 91)
;     ssth&lon  = fspan(0, 358, 180)
;     ssth&lat@units = "degrees_north"
;     ssth&lon@units = "degrees_east"

;     if (any(isnan_ieee(ssth))) then
;         value = 1.e+20
;         replace_ieeenan(ssth, value, 0)
;         ssth@_FillValue = value
;     end if

;     ;------detrend---------
;     ssthdt(i,:,:,:) =  dtrend_msg_n(time,ssth,False,False,0);去趋势
;     copy_VarAtts(ssth,ssthdt(i,:,:,:))
;     copy_VarCoords(ssth,ssthdt(i,:,:,:))

;     ;------anomaly---------
;     ssthclm(i,:,:,:) = clmMonTLL(ssthdt(i,:,:,:))
; end do
; printVarSummary(ssthclm)
; ;===============================================
; ;bias of models
; ;===============================================
; ;----------------read reanalysis----------------
;   dir1 = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
;   fName = addfile(dir1+"sst.mnmean.v5.nc","r")

;   TIME      = fName->time
;   YYYY      = cd_calendar(TIME,-1)/100
;   iYYYY     = ind(YYYY.ge.yrStrtmodel1 .and. YYYY.le.yrLastmodel1)
;   Time      = fName->time(iYYYY)
;   time      = ut_calendar(Time,-1)
;   lon       = fName->lon
;   lat       = fName->lat(::-1)

;   sstobs    = fName->sst(iYYYY,::-1,:)
;   sstobs!0  = "time"
;   sstobs!1  = "lat"
;   sstobs!2  = "lon"
;   sstobs&time = time

;  ;------interpolate sstobs to same resolution to SST from CMIP5---------
;   LON = fspan(0,358,180)
;   LAT = fspan(-90,90,91)
;   copy_VarAtts(lon, LON)
;   copy_VarAtts(lat, LAT)
;   SSTOBS = linint2_Wrap(lon,lat,sstobs,True,LON,LAT,0)

;   ;------detrend---------
;   sstodt  =  dtrend_msg_n(time,SSTOBS,False,False,0);去趋势
;   copy_VarAtts(SSTOBS,sstodt)
;   copy_VarCoords(SSTOBS,sstodt)

;   ;------anomaly---------
;   sstoclm = clmMonTLL(sstodt)

;   ; ;-------model bias------
;   ; sstOBSclm = conform(ssthclm, sstoclm, (/1,2,3/))
;   ; biasclm = ssthclm-sstOBSclm-273.15
;   ; copy_VarMeta(ssthclm, biasclm)
;   ; biasclm!0 = "time"
;   ; printVarSummary(biasclm)

;   ;-------model bias------
;   ssthdt_mme = dim_avg_n_Wrap(ssthdt, 0)
;   printVarSummary(ssthdt_mme)
;   printVarSummary(sstodt)
;   biasclm = ssthdt_mme - sstodt - 273.15
;   copy_VarCoords(sstodt, biasclm)
;   BIASclm = month_to_annual(biasclm, 1)
;   BIASclm&year = ispan(yrStrtmodel1, yrLastmodel1, 1)
;   ;BIASCLM = dim_avg_n_Wrap(BIASclm({1979:2013},:,:), 0)
;   BIASCLM = dim_avg_n_Wrap(BIASclm, 0)
;   test1 = wgt_areaave_Wrap(BIASCLM({-30:30},:), 1, 1, 0)
;   print(test1) ;-0.2008476 -0.1923213
;   exit


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
res@xyDashPattern = 0

plot(0) = gsn_csm_xy(wks, year, pp, res)
draw(plot(0))
frame(wks)


end