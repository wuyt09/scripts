
begin

lonstr = 0
lonend = 130
latstr = 10
latend = 60

yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = (/"MAM","JJA","SON","DJF"/)

LEVs = 500
LEVe = 200

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

filips = "/home/yangsong3/data-observation/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2019.nc"
fps = addfile(filips, "r")
ps = short2flt(fps->sp(iYYYY,::-1,:))
printVarSummary(ps)


plev = t&level
plev = plev*100
plev!0 = "lev"
plev&lev = t&level
plev@units = "Pa"

dp = dpres_plevel_Wrap(plev, ps, min(plev), 0)
printVarSummary(dp)

tdp = t*dp
copy_VarCoords(t, tdp)
printVarSummary(tdp)

t_tropo = dim_sum_n_Wrap(tdp(:,{LEVs:LEVe},:,:), 1)/dim_sum_n_Wrap(dp(:,{LEVs:LEVe},:,:), 1)
; t_tropo = dim_sum_n_Wrap(tdp, 1)/dim_sum_n_Wrap(dp, 1)
copy_VarCoords(t(:,0,:,:), t_tropo)

t2m_clm = month_to_seasonN(t_tropo, season)
printVarSummary(t2m_clm)

lines_west = wgt_areaave_Wrap(t2m_clm(:,:,{latstr_west:latend_west},{lonstr_west:lonend_west}), 1., 1., 0)
lines_east = wgt_areaave_Wrap(t2m_clm(:,:,{latstr_east:latend_east},{lonstr_east:lonend_east}), 1., 1., 0)
lines_south = wgt_areaave_Wrap(t2m_clm(:,:,{latstr_south:latend_south},{lonstr_south:lonend_south}), 1., 1., 0)
lines_north = wgt_areaave_Wrap(t2m_clm(:,:,{latstr_north:latend_north},{lonstr_north:lonend_north}), 1., 1., 0)

x_mask = mask(t2m_clm, topo.ge.1500, True)
lines_tp = wgt_areaave_Wrap(t2m_clm(:,:,{latstr_tp:latend_tp},{lonstr_tp:lonend_tp}), 1., 1., 0)
; lines_tp = wgt_areaave_Wrap(t2m_clm(:,:,{20:latstr_tp},{lonstr_tp:lonstr_tp}), 1., 1., 0)


lines = new((/5,4,nyear/), float)
lines(0,:,:) = lines_west
lines(1,:,:) = lines_east
lines(2,:,:) = lines_south
lines(3,:,:) = lines_north
lines(4,:,:) = lines_tp
copy_VarCoords(t2m_clm(:,:,0,0), lines(0,:,:))
lines!0 = "region"
lines&region = region
printVarSummary(lines)

lines_runave = runave_n_Wrap(lines, 5, 0, 2)
printVarSummary(lines_runave)

; rc = regCoef_n(year, lines, 0, 2)
; copy_VarCoords(lines(:,:,0), rc)
; printVarSummary(rc)
; exit
lines_trend = lines
rc = lines(:,:,0)
do i = 0, 4
    do j = 0, 3
        x_med = regline(year, lines(i,j,:))
        rc(i,j) = x_med
        lines_trend(i,j,:) = year*x_med+x_med@yintercept
    end do
end do
printVarSummary(lines_trend)

; tval = onedtond(rc@tval, dimsizes(rc))
; df   = onedtond(rc@nptxy, dimsizes(rc))
; b    = tval
; b    = 0.5
; prob = betainc(df/(df+tval^2), df/2.0, b)
; copy_VarCoords(rc, prob)

; rc@long_name = "regression coefficient"
; prob@long_name = "probability"


trend = rc*10
copy_VarCoords(rc, trend)
trend@units = "K/decade"
printVarSummary(trend)

pp = new((/3,5,4,nyear/), float)
pp(0,:,:,:) = lines
pp(1,:,:,:) = lines_runave
pp(2,:,:,:) = lines_trend
copy_VarCoords(lines, pp(0,:,:,:))
;===============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/"
figname = "time_series_sea_tropofT_ERA5"
wks = gsn_open_wks("pdf", figpath+figname)

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
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@tmXTOn = False
; res@tiXAxisString = "Month"


; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-20, 20, 5)
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYLPrecision = 4
res@tmYROn = False
res@tiYAxisString = ""
res@gsnStringFontHeightF = 0.03

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"forestgreen","red","orange2","blue"/)
Patterns = (/0,16,16,0,0/)
res@xyLineThicknessF = 10
res@xyLineColors = (/"black","red","blue"/)
res@xyDashPatterns = Patterns


; res@trYMinF = 280
; res@trYMaxF = 305
; res@gsnLeftString = SEA(0)
; plot(0) = gsn_csm_xy(wks, year, lines(:,0,:), res)

; res@trYMinF = 296
; res@trYMaxF = 302
; res@gsnLeftString = SEA(1)
; plot(1) = gsn_csm_xy(wks, year, lines(:,1,:), res)

; res@trYMinF = 280
; res@trYMaxF = 302
; res@gsnLeftString = SEA(2)
; plot(2) = gsn_csm_xy(wks, year, lines(:,2,:), res)

; res@trYMinF = 260
; res@trYMaxF = 302
; res@gsnLeftString = SEA(3)
; plot(3) = gsn_csm_xy(wks, year, lines(:,3,:), res)

pres = True
pres@gsnMaximize = True

do i = 0, 4
    ; res@xyLineColors = Colors(i)
    do j = 0, 3
        res@gsnLeftString = season(j)
        res@gsnRightString = "trend: "+sprintf("%4.2f", trend(i,j))+" K/decade"
        plot(j) = gsn_csm_xy(wks, year, pp(:,i,j,:), res)
    end do

    pres@gsnPanelMainString = region(i)+" ("+yrStrt+"-"+yrLast+")"
    gsn_panel(wks, plot, (/2,2/), pres)
end do



; ;==================================================================================
;     res_lines                   = True                  ; polyline mods desired
;     res_lines@gsLineThicknessF  = 5                     ; line thicker
    

;     res_text                    = True                  ; text mods desired
;     res_text@txFontHeightF      = 0.018                 ; change text size
;     res_text@txJust             = "CenterLeft"          ; text justification

;     Texts  = region

; ;============================================================
;     gx = (/1980,1983/)
;     text_pos = 1984
;     gy1 = (/302,302/)
;     interval = 2
;     gy2 = gy1 - interval
;     gy3 = gy2 - interval
;     gy4 = gy3 - interval
;     gy5 = gy4 - interval
;     ; gy6 = (/8,8/)

;     res_lines@gsLineThicknessF = 10
;     res_lines@gsLineColor = Colors(0)
;     res_lines@gsLineDashPattern = Patterns(0)
;     dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
;     dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
;     res_lines@gsLineColor = Colors(1)
;     res_lines@gsLineDashPattern = Patterns(1)
;     dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
;     dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
;     res_lines@gsLineColor = Colors(2)
;     res_lines@gsLineDashPattern = Patterns(2)
;     dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
;     dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

;     res_lines@gsLineColor = Colors(3)
;     res_lines@gsLineDashPattern = Patterns(3)
;     dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
;     dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
;     ; res_lines@gsLineColor = Colors(4)
;     ; res_lines@gsLineDashPattern = Patterns(4)
;     ; dum51 = gsn_add_polyline(wks, plot(0), gx, gy5, res_lines)
;     ; dum52 = gsn_add_text(wks, plot(0), Texts(4), text_pos, gy5, res_text)


; pres = True
; pres@gsnMaximize = True

; pres@gsnPanelMainString = region(0)
; gsn_panel(wks, plot(0:3), (/2,2/), pres)
; pres@gsnPanelMainString = region(1)
; gsn_panel(wks, plot(4:7), (/2,2/), pres)
; pres@gsnPanelMainString = region(2)
; gsn_panel(wks, plot(8:11), (/2,2/), pres)
; pres@gsnPanelMainString = region(3)
; gsn_panel(wks, plot(12:15), (/2,2/), pres)

end