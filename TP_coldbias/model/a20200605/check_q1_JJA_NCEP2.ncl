
begin


lonstr = 75
lonend = 105
latstr = 25
latend = 40

;global constants
G_CPD   =   1004.07     ;Cpd=1004.07 J*kg^-1*K^-1
G_RD    =   287.05      ;Rd=287.05 J*kg^-1*K^-1
G_G     =   9.8         ;g=9.8 m*s^-2


nlev = 17
nlat = 73
nlon = 144

season = "DJF"

f = addfile("/home/yangsong3/wangzq/DATA/jra/JRA55.Q1I.year.mon.1958.2014.nc", "r")
Q1X = f->Q1I

Q1_mon = dim_avg_n_Wrap(Q1X, 0)
q1 = (Q1_mon(0,:,:)+Q1_mon(1,:,:)+Q1_mon(11,:,:))/3
copy_VarCoords(Q1X(0,0,:,:), q1)
; Q1_mon = reshape(Q1X, (/57*12,145,288/))
; copy_VarCoords(Q1X(0,0,:,:), Q1_mon)
; Q1 = month_to_season(Q1_mon, season)
; q1 = dim_avg_n_Wrap(Q1, 0)
; printVarSummary(q1)

; diri = "/home/yangsong3/data-observation/NCEP_DOE/monthly/q1/"

; filis = systemfunc("ls "+diri+"q1.*.mon.mean.nc") 
; print(filis)

; f = addfiles(filis, "r")
; x = f[:]->q1
; printVarSummary(x)


; fp = addfile("/home/yangsong3/data-observation/NCEP_DOE/monthly/pres.mon.mean.nc", "r")
; TIME = fp->time
; YYYY = cd_calendar(TIME, -1)/100
; iYYYY = ind(YYYY.ge.1979 .and. YYYY.le.2013)
; PS = fp->pres(iYYYY,:,:)

; x = x*G_CPD/G_G
; ; x = x * 86400
; ; x9 = month_to_season(x, season)
; ; x8 = dim_avg_n_Wrap(x9, 0)
; ; printMinMax(x, False)

; lev = x&lev
; plev = lev
; plev = plev*100
; plev!0 = "lev"
; plev&lev = lev
; plev@units = "Pa"
; levs = 600
; levt = 300

; ptop = min(plev)
; psfc = PS*100
; dp = dpres_plevel_Wrap(plev, psfc, ptop, 0)
; printVarSummary(dp)

; xdp = x*dp
; copy_VarCoords(x, xdp)
; printMinMax(xdp, False)
; print(xdp&lev)
; printVarSummary(xdp)

; x0 = dim_sum_n_Wrap(xdp(:,{levs:levt},:,:), 1) ;
; copy_VarCoords(xdp(:,0,:,:), x0)
; printMinMax(x0, False)
; printVarSummary(x0)

; vt0 = month_to_season(x0, season)
; vt = dim_avg_n_Wrap(vt0, 0)

; printVarSummary(vt)
; printMinMax(vt, False)
;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = "check_q1_JJA_NCEP2"

lonstr = 60
lonend = 150
latstr = 10
latend = 50
;leftsti = (/"dTdt","Q1_v1","dTdx","dTdp","Q1_v2"/)

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "grey"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn        = False
res@cnInfoLabelOn         = False
res@cnFillPalette         = "CBR_coldhot"
;res@lbLabelBarOn          = False

res@mpMinLatF            = latstr      ; range to zoom in on
res@mpMaxLatF            = latend
res@mpMinLonF            = lonstr
res@mpMaxLonF            = lonend
; res@mpCenterLonF         = 180.
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@gsnStringFontHeightF = 0.035

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -200
res@cnMaxLevelValF = 200
res@cnLevelSpacingF = 20
; res@cnMinLevelValF = -5
; res@cnMaxLevelValF = 5
; res@cnLevelSpacingF = 0.5

; do i = 0, 1
    plot(0) = gsn_csm_contour_map(wks, q1, res) ; x8({300},:,:)
    plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(0), plot_topo(0))
; end do
draw(plot)
; pres = True
; pres@gsnMaximize = True
; ;pres@gsnPanelLabelBar = True
; ; pres@gsnPanelMainString = season+" thermal equation (600hPa to 300hPa)" ;
; ; pres@pmLabelBarOrthogonalPosF = -0.02
; ; pres@lbLabelFontHeightF = 0.015

; gsn_panel(wks, plot, (/2,1/), pres)

end