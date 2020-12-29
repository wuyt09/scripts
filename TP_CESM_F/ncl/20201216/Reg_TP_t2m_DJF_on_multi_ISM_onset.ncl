
begin

yrStrt = 1979
yrLast = 2018
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = "DJF"

nlat = 181
nlon = 360
nlev = 19


diri = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/ISM/"
fili1 = "MTG_onset_pentad_1979_2018.txt"
fili2 = "LBQ_onset_pentad_1979_2017.txt"
fili3 = "IMD_onset_pentad_1970_2020.txt"
fili4 = "TTGI_onset_pentad_1979_2018.txt"
fili5 = "WFI_onset_pentad_1979_2018.txt"
fili6 = "HOWI_onset_pentad_1979_2018.txt"

x1 = asciiread(diri+fili1, (/nyear,2/), "integer")
x2 = asciiread(diri+fili2, (/nyear-1,2/), "integer")
x3 = asciiread(diri+fili3, (/51,3/), "integer")
x4 = asciiread(diri+fili4, (/nyear,2/), "integer")
x5 = asciiread(diri+fili5, (/nyear,2/), "integer")
x6 = asciiread(diri+fili6, (/nyear,2/), "integer")
; print(x3(:,2))

x = new((/6,nyear/), integer)
x(0,:) = x1(:,1)
x(1,0:nyear-2) = x2(:,1)
x(2,:) = x3(9:48,2)
x(3,:) = x4(:,1)
x(4,:) = x5(:,1)
x(5,:) = x6(:,1)
printVarSummary(x)

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
        var = asciiread(diriobs+filiobs, -1, "float")
        ts_obs(i*12+j,:,:) = reshape(var(6:), (/72,128/))
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
t2m_sea = month_to_season(ts_obs, season)
printVarSummary(t2m_sea)


x_djf = x ;(:,1:nyear-1)
t2m_djf = t2m_sea ;(0:nyear-2,:,:)

rc = new((/6,72,128/), float)
prob = rc
do i = 0, 5
    rcx = regCoef_n(x_djf(i,:), t2m_djf, 0, 0)
    rc(i,:,:) = rcx
    df = onedtond(rcx@nptxy, dimsizes(rcx))-2
    tval = onedtond(rcx@tval, dimsizes(rcx))
    b = tval
    b = 0.5
    prob(i,:,:) = betainc(df/(df+tval^2), df/2.0, b)
end do
copy_VarCoords(t2m_sea(0,:,:), rc(0,:,:))
copy_VarCoords(rc, prob)
printVarSummary(rc)
printVarSummary(prob)


;==============================================================
;plot
;==============================================================
lonstr = 60
lonend = 120
latstr = 10
latend = 50

Texts  = (/"MTG","LBQ","IMD","TTGI","WFI","HOWI"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Reg_TP_t2m_DJF_on_multi_ISM_onset"

wks = gsn_open_wks("pdf", figpath+figname) ;
plot = new(6, graphic)
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

do i = 0, 5
    res@gsnLeftString = Texts(i)
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
pres@gsnPanelMainString = "DJF surf temp"

gsn_panel(wks, plot, (/2,3/), pres)
end