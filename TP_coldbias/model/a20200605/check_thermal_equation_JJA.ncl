
begin


lonstr = 75
lonend = 105
latstr = 25
latend = 40
yrStrt = 198102
yrLast = 200601

;global constants
G_CPD   =   1004.07     ;Cpd=1004.07 J*kg^-1*K^-1
G_RD    =   287.05      ;Rd=287.05 J*kg^-1*K^-1
G_G     =   9.8         ;g=9.8 m*s^-2

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
ncase = dimsizes(casename)
nlev = 17
nlat = 96
nlon = 144

season = "DJF"

vt = new((/8,96,144/), float)
diri1 = "/home/yangsong3/data-model/wuyt/TPbias/TEMP/"
diri2 = "/home/yangsong3/data-model/wuyt/TPbias/TEMPv2/"

i = 1
; do i = 0, ncase-1
    diri3 = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".dTdt.1981-2005.mon.mean.nc"
    fili2 = casename(i)+".Q1.1981-2005.mon.mean.nc"
    fili3 = casename(i)+".VgradT.1981-2005.mon.mean.nc"
    fili4 = casename(i)+".Vertical.1981-2005.mon.mean.nc"

    fili5 = casename(i)+".dTdt.1981_2005.monmean.nc"
    fili6 = casename(i)+".q1.1981_2005.monmean.nc"
    fili7 = casename(i)+".VgradT.1981_2005.monmean.nc"
    fili8 = casename(i)+".Vertical.1981_2005.monmean.nc"

    fili9 = casename(i)+".cam.h0.PS.1979-2005.nc"

    f1 = addfile(diri1+fili1, "r")
    f2 = addfile(diri1+fili2, "r")
    f3 = addfile(diri1+fili3, "r")
    f4 = addfile(diri1+fili4, "r")

    f5 = addfile(diri2+fili5, "r")
    f6 = addfile(diri2+fili6, "r")
    f7 = addfile(diri2+fili7, "r")
    f8 = addfile(diri2+fili8, "r")

    f9 = addfile(diri3+fili9, "r")

    dTdt = f1->P1
    Q1_v1 = f2->Q1
    dTdx = f3->P2
    dTdx = dTdx*6.371e6
    dTdp = f4->P3
    Q1_v1 = dTdx+dTdt+dTdp
    printVarSummary(Q1_v1)

    dTdt_v2 = f5->dsdt
    Q1_v2 = f6->q1
    dTdx_v2 = f7->adv_s
    dTdp_v2 = f8->dsdp
    dTdp_v2 = dTdp_v2*(-1.)
    Q1_v2 = dTdx_v2+dTdt_v2+dTdp_v2
    printVarSummary(Q1_v2)

    TIME = f9->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    PS = f9->PS(iYYYY,:,:)

    x = new((/300,8,nlev,nlat,nlon/), float)
    x(:,0,:,:,:) = (/Q1_v1/)
    x(:,1,:,:,:) = (/Q1_v2/)
    x(:,2,:,:,:) = (/dTdt/)
    x(:,3,:,:,:) = (/dTdt_v2/)
    x(:,4,:,:,:) = (/dTdx/)
    x(:,5,:,:,:) = (/dTdx_v2/)
    x(:,6,:,:,:) = (/dTdp/)
    x(:,7,:,:,:) = (/dTdp_v2/)
    copy_VarCoords(Q1_v2, x(:,0,:,:,:))
    x = x*G_CPD/G_G

    lev = Q1_v2&lev
    plev = lev
    plev = plev*100
    plev!0 = "lev"
    plev&lev = lev
    plev@units = "Pa"
    levs = 600
    levt = 300

    ptop = min(plev)
    psfc = PS
    dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
    dp = conform_dims(dimsizes(x), dP, (/0,2,3,4/))
    copy_VarCoords(x, dp)
    xdp = x*dp
    copy_VarCoords(x, xdp)
    x0 = dim_sum_n_Wrap(xdp(:,:,{levs:levt},:,:), 2) ;/dim_sum_n_Wrap(dp(:,:,{levs:levt},:,:), 2)
    ; x0 = dim_sum_n_Wrap(xdp, 2)/dim_sum_n_Wrap(dp, 2)
    copy_VarCoords(xdp(:,:,0,:,:), x0)
    vt0 = month_to_season(x0, season)
    vt = dim_avg_n_Wrap(vt0, 0)
    printVarSummary(vt)
; end do

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = "check_thermal_equation_JJA"

lonstr = 60
lonend = 150
latstr = 10
latend = 50
leftsti = (/"Q1","dTdt","dTdx","dTdp"/)
rghtsti = (/"Version1","Version2"/)

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(8, graphic)
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

do i = 0, 7
    res@gsnLeftString = leftsti(i/2)
    res@gsnRightString = rghtsti(mod(i, 2))
    plot(i) = gsn_csm_contour_map(wks, vt(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = season+" thermal equation (600hPa to 300hPa)" ;
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/4,2/), pres)

end