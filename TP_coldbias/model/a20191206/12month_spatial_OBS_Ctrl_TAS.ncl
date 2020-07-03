
begin
lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================
yrStrt = 1961
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyr = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nm = dimsizes(month)

ts_obsx = new((/nyr*nm,72,128/), float,-99.) ;

diriobs = "/home/yangsong3/data-observation/stationChina2000/Grid/t2m/"

do i = 0, nyr-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        filiobs = "SURF_CLI_CHN_TEM_MON_GRID_0.5-MEAN-"+yr+mn+".txt"
        x = asciiread(diriobs+filiobs, -1, "float")
        ts_obsx(i*12+j,:,:) = reshape(x(6:), (/72,128/))
    end do
end do
printVarSummary(ts_obsx)


opt = True
opt@fout = "./test.txt"
write_matrix(ts_obsx(0,:,:), "128f7.2", opt)

ts_obsx = ts_obsx(:,::-1,:)
time = ispan(1, nyr*12, 1)
ts_obs_dtrend = dtrend_msg_n(time, ts_obsx, False, False, 0)
ts_obs = reshape(ts_obs_dtrend, (/nyr,nm,72,128/))

lat = fspan(18.25, 53.75, 72)
lon = fspan(72.25, 135.75, 128)
print(lat)
print(lon)

ts_obs!0 = "year"
ts_obs!1 = "month"
ts_obs!2 = "lat"
ts_obs!3 = "lon"
ts_obs&year = year
ts_obs&month = month
ts_obs&lat = lat
ts_obs&lon = lon
ts_obs&lat@units = "degrees_north"
ts_obs&lon@units = "degrees_east"
printVarSummary(ts_obs)

TS_OBS = dim_avg_n_Wrap(ts_obs, 0)
printVarSummary(TS_OBS)

;===================================================================
; Ctrl 24 data
;===================================================================
yrStrt = 198102
yrLast = 200601
diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_Ctrl_FAMIP/a20191206/"
fili = "Ctrl_FAMIP.cam.h0.TREFHT.1979-2005.nc"
f = addfile(diri+fili, "r")
TIME = f->time
YYYY = cd_calendar(TIME, -1)
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
TS = f->TREFHT(iYYYY,{18:54},{72:136})
TS_Ctrl = clmMonTLL(TS)
TS_Ctrl = TS_Ctrl - 273.15
printVarSummary(TS_Ctrl)

TS_obs = linint2_Wrap(lon, lat, TS_OBS, False, TS_Ctrl&lon, TS_Ctrl&lat, 0)

ts_diff = TS_Ctrl
ts_diff = TS_Ctrl - TS_obs
printVarSummary(ts_diff)


; ;===================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20191206/"
figname = "12month_spatial_OBS_Ctrl_TAS"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(12, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "green4" ;grey80
 res_topo@cnLineDashPattern      = 12
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

; cmap = (/                     \
;             (/ 255./255, 255./255, 255./255 /),    \  ; 0 - White background.
;             (/ 0./255  , 0./255  , 0./255   /),         \  ; 1 - Black foreground.
;             (/ 255./255, 0./255  , 0./255   /),        \  ; 2 - Red.
;             (/ 0./255  , 0./255  , 255./255 /),        \  ; 3 - Blue.
;             (/ 164./255, 244./255, 131./255 /),      \  ; 4 - Ocean Blue.
;             (/ 0./255  , 0./255  , 255./255   /),         \  ; 5 - Bar 1
;          (/ 0./255  , 153./255, 255./255 /),      \  ; 6 - Bar 2
;             (/ 0./255, 153./255, 153./255 /),        \  ; 7 - Bar 3 
;             (/ 0./255  , 255./255, 0./255  /),       \  ; 8 - Bar 4   
;             (/ 255./255, 255./255  , 102./255 /),        \  ; 9 - Bar 5 
;             (/ 255./255, 153./255  , 102./255   /),        \  ; 10 - Bar 6
;             (/ 255./255, 0./255  , 255./255   /)        \  ; 11 - Bar 7   
;          /)
;  gsn_define_colormap(wks,cmap)

  res = True
  res@gsnAddCyclic = False
  res@gsnDraw  = False
  res@gsnFrame = False
  res@vpHeightF = 0.4
  res@vpWidthF = 0.7

  res@mpOutlineOn = True
  res@mpFillOn = True
  res@mpGeophysicalLineThicknessF = 2
  res@mpGeophysicalLineColor = "grey45" ;"grey25"
  res@mpFillDrawOrder = "PreDraw"

   res@mpFillOn       = True
   res@mpDataSetName = "Earth..4"
   res@mpDataBaseVersion = "MediumRes"
   res@mpAreaMaskingOn  = True
  ; re0@mpCenterLonF = 104
   res@mpMinLonF    = 71 ;97 ;72
   res@mpMaxLonF    = 138
   res@mpMinLatF    = 17
   res@mpMaxLatF    = 55
   ;res@mpMaskAreaSpecifiers = (/"China"/)
   res@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
   res@mpOutlineBoundarySets   = "NoBoundaries"  ;if national then the coutries around will be included
   res@mpFillBoundarySets         = "NoBoundaries"
   res@mpOceanFillColor = 0 
   res@mpLandFillColor    = 0
   res@mpInlandWaterFillColor = 0 
   res@mpNationalLineColor        = "black"
   res@mpProvincialLineColor      = "black"
   res@mpGeophysicalLineColor     = "black"
   res@mpNationalLineThicknessF   = 2
   res@mpProvincialLineThicknessF = 1

  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False
  res@lbLabelBarOn = False

  res@mpMinLatF = latstr
  res@mpMaxLatF = latend
  res@mpMinLonF = lonstr
  res@mpMaxLonF = lonend
  res@mpCenterLonF = (lonstr+lonend)/2

  res@gsnStringFontHeightF = 0.035
  res@tmXBLabelFontHeightF = 0.025
  res@tmYLLabelFontHeightF = 0.025
  res@pmLabelBarOrthogonalPosF = 0.2

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = 1
  res@cnMaxLevelValF = 10.
  res@cnMinLevelValF = -10.
  
  do i = 0, 11
    res@gsnLeftString = month(i)
    plot(i) = gsn_csm_contour_map(wks, ts_diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
  end do


  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelLabelBar = True
  ;pres@pmLabelBarOrthogonalPosF = -0.01
  pres@lbLabelFontHeightF = 0.015
  pres@gsnPanelMainString = "Ctrl - OBS"
  gsn_panel(wks, plot, (/3,4/), pres)
; ;end do

end