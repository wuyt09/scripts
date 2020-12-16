
begin

lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

yrStrt = 1979
yrLast = 2017
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = (/"DJF","MAM"/)

fx = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/multi_standard_onset_BOB_SCS_ISM.txt"
onset_time = asciiread(fx, (/nyear,4/), "integer")
inx = onset_time(:,2)


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================

ts_obs = new((/nyear*nmonth,72,128/), float,-99.)

diriobs = "/home/yangsong3/data-observation/stationChina2000/Grid/t2m/"

do i = 0, nyear-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        filiobs = "SURF_CLI_CHN_TEM_MON_GRID_0.5-MEAN-"+yr+mn+".txt"
        x = asciiread(diriobs+filiobs, -1, "float")
        ts_obs(i*12+j,:,:) = reshape(x(6:), (/72,128/))
    end do
end do
printVarSummary(ts_obs)

ts_obs = ts_obs(:,::-1,:)
time = ispan(1, nyear*12, 1)
ts_obs = dtrend_msg_n(time, ts_obs, False, False, 0)


lat = fspan(18.25, 53.75, 72)
lon = fspan(72.25, 135.75, 128)
; print(lat)
; print(lon)

ts_obs!0 = "time"
ts_obs!1 = "lat"
ts_obs!2 = "lon"
ts_obs&time = ispan(1, nyear*nmonth, 1)
ts_obs&lat = lat
ts_obs&lon = lon
ts_obs&lat@units = "degrees_north"
ts_obs&lon@units = "degrees_east"
printVarSummary(ts_obs)

topox1 = ftopo->topo({18:54},{72:136})
lonx = topox1&lon
latx = topox1&lat
printVarSummary(topox1)
topox2 = linint2_Wrap(lonx, latx, topox1, False, lon, lat, 0)
print(topox2({35},{90}))

; TS_obs = mask(ts_obs, topox2.ge.1500, True)
; copy_VarCoords(ts_obs, TS_obs)
t2m_sea = month_to_seasonN(ts_obs, season)
printVarSummary(t2m_sea)

rc = regCoef_n(inx, t2m_sea, 0, 1)
copy_VarCoords(t2m_sea(:,0,:,:), rc)
printVarSummary(rc)

df = onedtond(rc@nptxy, dimsizes(rc))-2
tval = onedtond(rc@tval, dimsizes(rc))
b = tval
b = 0.5
prob = betainc(df/(df+tval^2), df/2.0, b)
copy_VarCoords(rc, prob)
printVarSummary(prob)


;==============================================================
;plot
;==============================================================
lonstr = 60
lonend = 120
latstr = 10
latend = 50


figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Reg_TP_t2m_DJF_MAM_on_SCS_onset"

wks = gsn_open_wks("pdf", figpath+figname) ;
plot = new(2, graphic)
plot_topo = plot
plot_prob = plot

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

res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False
res@gsnAddCyclic          = False


res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@gsnStringFontHeightF = 0.025
; res@gsnRightString       = "K"

res@pmLabelBarOrthogonalPosF = 0.2
res@pmLabelBarWidthF = 0.7
res@lbLabelFontHeightF = 0.015
res@cnFillPalette = "CBR_coldhot"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 0.6
res@cnMinLevelValF = -0.6
res@cnLevelSpacingF = 0.1
res@lbLabelBarOn = False

;----------------------------------------------
; draw prob
;----------------------------------------------
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

 ; res3@cnFillPattern = 17 ;(/17/)
 res3@cnFillColor  =  "black"
 res3@gsnLeftString = ""
 res3@cnFillDotSizeF = 0.003
 res3@cnFillScaleF  = 1.2
 ;----------------------------------------------
 ;----------------------------------------------

do i = 0, 1
    res@gsnLeftString = "t2m: "+season(i)
    plot(i) = gsn_csm_contour_map(wks, rc(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))

    plot_prob(i) = gsn_csm_contour(wks,prob(i,:,:),res3)
    plot_prob(i) = ShadeLtContour(plot_prob(i),0.05,6)
    overlay(plot(i), plot_prob(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarWidthF = 0.7
pres@gsnPanelMainString = "SCS onset"

gsn_panel(wks, plot, (/1,2/), pres)
end