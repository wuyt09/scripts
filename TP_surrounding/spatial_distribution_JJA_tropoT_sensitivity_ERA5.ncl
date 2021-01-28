
begin

lonstr = 0
lonend = 130
latstr = 10
latend = 60

yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

; LEVs = 500
; LEVe = 200
; LEVs = 850
; LEVe = 200
; LEVs = 850
; LEVe = 300
; LEVs = 1000
; LEVe = 300
LEVs = (/500,850,850,1000/)
LEVe = (/200,200,300,300/)

latstr_east = 25
latend_east = 45
lonstr_east = 105
lonend_east = 120

latstr_west = 30
latend_west = 55
lonstr_west = 30
lonend_west = 70

latstr_south = 10
latend_south = 25
lonstr_south = 75
lonend_south = 110

latstr_north = 45
latend_north = 60
lonstr_north = 70
lonend_north = 110

latstr_tp = 25 ;28
latend_tp = 40 ;38
lonstr_tp = 75 ;80
lonend_tp = 105

region = (/"West","East","South","North","TP"/)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)


;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili = "temp.monthly.197901-201812.nc"
f = addfile(diri+fili, "r")

Time = f->time
YYYY = cd_calendar(Time, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

t = short2flt(f->t(iYYYY,::-1,::-1,:))
printVarSummary(t)

lat = t&latitude
lon = t&longitude
nlat = dimsizes(lat)
nlon = dimsizes(lon)

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))
printVarSummary(ps)


plev = tofloat(t&level)
plev = plev*100
plev!0 = "lev"
plev&lev = t&level
plev@units = "Pa"

t_tropo = new((/dimsizes(iYYYY),4,nlat,nlon/), float)
do i = 0, 3
    dp := dpres_plevel_Wrap(plev, ps, LEVs(i)*100, 0)
    ; dP = dpres_plevel_Wrap(plev, max(plev), min(plev), 0)
    ; dp = conform_dims(dimsizes(t), dP, 1)
    ; copy_VarCoords(t, dp)
    printVarSummary(dp)

    tdp := t(:,{LEVs(i):LEVe(i)},:,:)*dp(:,{LEVs(i):LEVe(i)},:,:)
    copy_VarCoords(t(:,{LEVs(i):LEVe(i)},:,:), tdp)
    printVarSummary(tdp)

    t_tropo(:,i,:,:) = dim_sum_n_Wrap(tdp, 1)/dim_sum_n_Wrap(dp(:,{LEVs(i):LEVe(i)},:,:), 1)
    copy_VarCoords(t(:,0,:,:), t_tropo(:,i,:,:))
end do

; t_tropo_2 = dim_sum_n_Wrap(tdp(:,{LEVs(1):LEVe(1)},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs(1):LEVe(1)},:,:), 1)
; copy_VarCoords(t(:,0,:,:), t_tropo_2)

; t_tropo_3 = dim_sum_n_Wrap(tdp(:,{LEVs(2):LEVe(2)},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs(2):LEVe(2)},:,:), 1)
; copy_VarCoords(t(:,0,:,:), t_tropo_3)

; t_tropo_4 = dim_sum_n_Wrap(tdp(:,{LEVs(3):LEVe(3)},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs(3):LEVe(3)},:,:), 1)
; copy_VarCoords(t(:,0,:,:), t_tropo_4)
; t_tropo = dim_sum_n_Wrap(tdp, 1)/dim_sum_n_Wrap(dp, 1)

t2m_clmx = month_to_season(t_tropo, season)
t2m_clm = dim_avg_n_Wrap(t2m_clmx, 0)
printVarSummary(t2m_clm)


; lines = new((/4,nlat,nlon/), float)
; lines(0,:,:) = t2m_clm_1
; lines(1,:,:) = t2m_clm_2
; lines(2,:,:) = t2m_clm_3
; lines(3,:,:) = t2m_clm_4
; copy_VarCoords(t2m_clm_1, lines(0,:,:))
; lines!0 = "type"
; lines&type = region
; printVarSummary(lines)
;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/"
figname = "spatial_distribution_"+season+"_tropoT_sensitivity_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(4, graphic)

plot_topo = plot
plot_prob = new(1, graphic)

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 5
 res_topo@cnLineColor            = "forestgreen"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


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

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 5
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
; res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnFillPalette = "CBR_coldhot"
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 250
res@cnMaxLevelValF = 270
res@cnLevelSpacingF = 2
; res@lbLabelBarOn = False

; res@gsnLeftString = season+" trend ("+yrStrt+"-"+yrLast+")"
res@pmLabelBarOrthogonalPosF = 0.12
res@lbLabelFontHeightF = 0.015
; res@tiMainString = "surface temperature"

do i = 0, 3
    res@gsnLeftString = LEVs(i)+"-"+LEVe(i)+" hPa Tropo T"
    plot(i) = gsn_csm_contour_map(wks, t2m_clm(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True

gsn_panel(wks, plot, (/2,2/), pres)

end