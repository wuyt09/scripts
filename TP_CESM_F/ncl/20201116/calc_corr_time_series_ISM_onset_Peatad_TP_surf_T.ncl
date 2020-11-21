
begin
yrStrt = 1979
yrLast = 2017
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

month = ispan(1, 12, 1)
nmonth = dimsizes(month)

inx_name = (/"TP_t2m","TTGI","WYI","WFI"/)

fili = "./ISM_onset_pentad.txt"
x = asciiread(fili, (/nyear,4/), "float")
; print(x)
x(38,2) = x@Fill_Value
x(18,2) = 27


x(:,0) = (/33,31,30,30,33,30,30,31,31,29,31,27,32,33,30,30,33,31,33,32,28,30,32,32,32,27,32,29,31,30,31,29,30,32,30,32,32,31,31/)



x!0 = "year"
x!1 = "var"
printVarSummary(x)
; exit
;====read t2m daily data================================
; yrStrt = 1979
; yrLast = 1981 ;2017
; year := ispan(yrStrt, yrLast, 1)
; nyear := dimsizes(year)

nlat = 181
nlon = 360
sdiri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/"
T2M = new((/nyear*nmonth,nlat,nlon/), float)
do j = 0, nyear-1
    YEAR = sprinti("%0.4i", year(j))
    ; print(YEAR)
    do i = 0, nmonth-1
        MONTH = sprinti("%0.2i", month(i))
        ; print(MONTH)
        filist = systemfunc("ls "+sdiri+"t2m/t2m."+YEAR+"-"+MONTH+"*.daily.nc")
        print(filist)
        fst = addfile(filist, "r")
        T2MX = short2flt(fst->t2m(:,::-1,:))
        T2M(j*12+i,:,:) = dim_avg_n_Wrap(T2MX, 0)
        copy_VarCoords(T2MX(0,:,:), T2M(j*12+i,:,:))
        delete([/T2MX/])
    end do
end do
T2M!0 = "time"
printVarSummary(T2M)


ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
T2M = mask(T2M, topo.ge.1500, True)
t2m = wgt_areaave(T2M(:,{25:40},{75:105}), 1, 1, 0)
copy_VarCoords(T2M(:,0,0), t2m)

t2m_spr = month_to_season(t2m, "MAM")
t2m_win = month_to_season(t2m, "DJF")
printVarSummary(t2m_spr)

;====================================================

corr_spr2 = new(4, float)
corr_spr = new(4, float)
corr_win = new(4, float)
do i = 0, 3
    corr_spr(i) = escorc(t2m_spr, x(:,i))    
    corr_win(i) = escorc(t2m_win, x(:,i))
end do

corr_t2m = escorc(t2m_spr, t2m_win)

print("spr: (month2season): "+sprintf("%4.2f", corr_spr))
print("win "+sprintf("%4.2f", corr_win))
print("t2m "+sprintf("%4.2f", corr_t2m))



y = new((/2,nyear/), float)
y!0 = "sea"
y!1 = "year"
y&year = year
y(0,:) = (/t2m_spr/)
y(1,:) = (/t2m_win/)
; copy_VarCoords(t2m_spr, y(0,:))
printVarSummary(y)
;=============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201116/"
figname = "calc_corr_time_series_ISM_onset_Peatad_TP_surf_T"



wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = Month
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; res@tmXTOn = False
; res@tiXAxisString = "Month"

res@trYMinF = 22
res@trYMaxF = 40
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(22, 40, 1)
res@tmYLLabels = "P"+res@tmYLValues
; res@tmYROn = False
; res@tiYAxisString = "Temperature (~S~o~N~C)"

;Colors = (/"black","deeppink2", "orange2", "green4", "black" /);
Colors = (/"red","blue","forestgreen","purple"/)
Patterns = (/0,0,0,0/)
res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

; plot(0) = gsn_csm_xy(wks, year, ptime, res)

res2 = True
res2@gsnDraw  = False
res2@gsnFrame = False
res2@gsnMaximize = True

res2@vpWidthF = 0.7
res2@vpHeightF = 0.4

res2@trXMinF = yrStrt
res2@trXMaxF = yrLast

res2@xyLineThicknessF = 10
res2@xyLineColor = "black"
res2@xyDashPatterns = (/0,16/)


plot(0) = gsn_csm_xy2(wks, year, x(var|:,year|:), y(sea|:,year|:), res, res2)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"LiuBoqi","TTGI","WYI","WFI"/)

;============================================================
    gx = (/1983,1984.5/)
    text_pos = 1985

    interval = 1
    gy1 = (/39,39/)
    gy2 = gy1 - interval
    gy3 = gy2 - interval
    gy6 = gy3 - interval

    gx2 = (/1991.5,1993/)
    text_pos2 = 1993.5
    gy4 = gy1
    gy5 = gy2

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

    res_lines@gsLineColor = Colors(3)
    res_lines@gsLineDashPattern = Patterns(3)
    dum61 = gsn_add_polyline(wks, plot(0), gx, gy6, res_lines)
    dum62 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy6, res_text)



    res_lines@gsLineColor = "black"
    res_lines@gsLineDashPattern = 0
    dum41 = gsn_add_polyline(wks, plot(0), gx2, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(0), "TP_t2m_spring", text_pos2, gy4, res_text)
    res_lines@gsLineColor = "black"
    res_lines@gsLineDashPattern = 2
    dum51 = gsn_add_polyline(wks, plot(0), gx2, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), "TP_t2m_winter", text_pos2, gy5, res_text)


    draw(plot(0))

end