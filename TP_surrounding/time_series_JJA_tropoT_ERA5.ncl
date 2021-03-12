
begin

lonstr = 0
lonend = 130 ;360
latstr = 10 ;10
latend = 60

yrStrt = 1979
yrLast = 2019
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

LEVs = 1000
LEVe = 300

latstr_east = 35
latend_east = 55
lonstr_east = 84 ;80
lonend_east = 115 ;120

latstr_west = 40
latend_west = 55
lonstr_west = 10
lonend_west = 50

latstr_south = 15
latend_south = 25
lonstr_south = 60
lonend_south = 110

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)


;===================================================================
diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"
fili = "t_mon_1x1_1979_2019.nc" ;"temp.monthly.197901-201812.nc"
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
; print(plev)

dp = dpres_plevel_Wrap(plev, ps, min(plev), 0) ;LEVs(i)*100
; dP = dpres_plevel_Wrap(plev, max(plev), min(plev), 0)
; dp = conform_dims(dimsizes(t), dP, 1)
; copy_VarCoords(t, dp)
printVarSummary(dp)
; print(dp(90,:,{35},{90}))
; exit

tdp = t*dp
copy_VarCoords(t, tdp)
printVarSummary(tdp)

t_tropo = dim_sum_n_Wrap(tdp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
copy_VarCoords(t(:,0,:,:), t_tropo)

t2m_clm = month_to_season(t_tropo, season)
printVarSummary(t2m_clm)

fa      = addfile ("/home/yangsong3/wuyt/data/lsm_1279l4_0.1x0.1.grb_v4_unpack.nc","r")
; lsdata = doubletofloat(a->lsm(0,:,:))
lsdata = fa->lsm(0,::-1,:)
printVarSummary(lsdata)
lsm = linint2_Wrap(lsdata&longitude, lsdata&latitude, lsdata, True, t&longitude, t&latitude, 0)
printVarSummary(lsm)


lsm1 = conform_dims(dimsizes(t2m_clm), lsm, (/1,2/))
t2m_clm = where(lsm1.eq.1, t2m_clm, t@_FillValue) ;



lines = new((/3,nyear/), float)
lines(0,:) = wgt_areaave_Wrap(t2m_clm(:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
lines(1,:) = wgt_areaave_Wrap(t2m_clm(:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
lines(2,:) = wgt_areaave_Wrap(t2m_clm(:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)

lines_base = dim_avg_n_Wrap(lines(:,0:9), 1)
lines = lines - conform_dims(dimsizes(lines), lines_base, 0)
; lines = dim_standardize_n_Wrap(lines, 1, 1)
lines_runave = runave_n_Wrap(lines, 10, 0, 1)
printVarSummary(lines_runave)

lines_trend = lines
rc = lines(:,0)
do i = 0, 2
    x_med = regline(year, lines(i,:))
    rc(i) = x_med
    lines_trend(i,:) = year*x_med+x_med@yintercept
end do
printVarSummary(lines_trend)

rc = rc*10
rc@units = "K/decade"

pp = new((/9,nyear/), float)
pp(0:2,:) = lines
; pp(3:5,:) = lines_runave
pp(6:8,:) = lines_trend

; rc = regCoef_n(year, t2m_clm, 0, 0)
; copy_VarCoords(t2m_clm(0,:,:), rc)

; tval = onedtond(rc@tval, dimsizes(rc))
; df   = onedtond(rc@nptxy, dimsizes(rc))
; b    = tval
; b    = 0.5
; prob = betainc(df/(df+tval^2), df/2.0, b)
; copy_VarCoords(rc, prob)

; rc@long_name = "regression coefficient"
; prob@long_name = "probability"

; trend = rc*10
; copy_VarCoords(rc, trend)
; trend@units = "K/decade"
; printVarSummary(trend)

; sig = 0.05


;==============================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/"
figname = "time_series_JJA_tropoT_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(1, graphic)


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.8
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
res@trYMaxF = 3 ;1.5
res@trYMinF = -1 ;-0.5
; res@tmXBMode = "Explicit"
; res@tmXBValues = year
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@tmXTOn = False
; res@tiXAxisString = "Month"
res@tmYROn = False
res@tiYAxisString = ""
res@gsnStringFontHeightF = 0.03

res@gsnYRefLine = 0
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineDashPattern = 16
res@gsnYRefLineThicknessF = 6

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"red","blue","green4","red","blue","green4","red","blue","green4"/)
Patterns = (/0,0,0,2,2,2,16,16,16/)
res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns


plot(0) = gsn_csm_xy(wks, year, pp, res)




;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"west region "+sprintf("%8.3f", rc(0))+"K/decade","east region "+sprintf("%8.3f", rc(1))+"K/decade","south region "+sprintf("%8.3f", rc(2))+"K/decade"/)

;============================================================
    gx = (/1981.5,1984.2/)
    text_pos = 1984.5

    interval = 0.3
    gy1 = (/2.7,2.7/)
    gy2 = gy1 - interval
    gy3 = gy2 - interval

    res_lines@gsLineThicknessF = 5
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)


draw(plot(0))

end