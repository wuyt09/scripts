
begin

lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

season = (/"DJF","MAM"/)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo1 = ftopo->topo({14.75:55.25},{69.75:140.25})
printVarSummary(topo1)
topo2 = ftopo->topo(::-1,:)
topo = topo2
;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

x = new((/4,nyear/), float)
;===================================================================
; CN051
;===================================================================
diri = "/home/yangsong3/data-observation/CN05.1/"
fili = "CN05.1_Tm_1961_2018_month_025x025.nc"
f1 = addfile(diri+fili, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
time = f1->time(iYYYY)
ts_obs = short2flt(f1->tm(iYYYY,:,:))
printVarSummary(ts_obs)

nlon = dimsizes(ts_obs&lon)
nlat = dimsizes(ts_obs&lat)

; ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
TOPO = linint2_Wrap(topo1&lon, topo1&lat, topo1, False, ts_obs&lon, ts_obs&lat, 0)
TS_obs = mask(ts_obs, TOPO.ge.1500, True)
copy_VarCoords(ts_obs, TS_obs)
TS_obs_sea = month_to_seasonN(ts_obs, season)
printVarSummary(TS_obs_sea)

TP_cn051 = dim_avg_n_Wrap(TS_obs_sea, 1)
printVarSummary(TP_cn051)

; ts_mask = linint2_Wrap(ts_obs&lon, ts_obs&lat, ts_obs, False, topo2&lon, topo2&lat, 0)
delete([/diri,fili,f1,Time,YYYY,iYYYY,time,ts_obs,nlon,nlat,ts_obs,TOPO,TS_obs,TS_obs_sea/])
;===================================================================
; ERA5
;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/surface/"
fili = "t2m_mon_1x1_1979_2019.nc"
f1 = addfile(diri+fili, "r")

Time = f1->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
time = f1->time(iYYYY)
ts_obs = short2flt(f1->t2m(iYYYY,::-1,:))
printVarSummary(ts_obs)
nlon = dimsizes(ts_obs&longitude)
nlat = dimsizes(ts_obs&latitude)

; ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)
; ts_obs = mask(ts_obs, .not.ismissing(ts_mask), True)
TS_obs = mask(ts_obs, topo2.ge.1500, True)
copy_VarCoords(ts_obs, TS_obs)
TS_obs_sea = month_to_seasonN(ts_obs, season)
printVarSummary(TS_obs_sea)

TP_era5 = dim_avg_n_Wrap(TS_obs_sea, 1)
TP_era5 = TP_era5 - 273.15
printVarSummary(TP_era5)

;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210410/"
figname = "check_spatial_distribution_TP_t2m_DJF_MAM_CN051_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(4, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey40"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;---------------------------------------------------------
res = True
res@gsnDraw = False
res@gsnFrame = False
; res@gsnMaximize = True
res@gsnAddCyclic = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillPalette = "BlueWhiteOrangeRed"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 2
res@cnMaxLevelValF = 20
res@cnMinLevelValF = -20
res@lbLabelBarOn = False

do i = 0, 1
    res@gsnLeftString = "CN05.1"
    res@gsnRightString = season(i)
    plot(i) = gsn_csm_contour_map(wks, TP_cn051(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    res@gsnLeftString = "ERA5"
    plot(i+2) = gsn_csm_contour_map(wks, TP_era5(i,:,:), res)
    plot_topo(i+2) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i+2), plot_topo(i+2))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/2,2/), pres)
end