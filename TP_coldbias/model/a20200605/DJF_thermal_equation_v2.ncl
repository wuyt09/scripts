
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

vt = new((/ncase,4,96,144/), float)
diri1 = "/home/yangsong3/data-model/wuyt/TPbias/TEMPv2/"

do i = 0, ncase-1
    diri2 = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".dTdt.1981_2005.monmean.v2.nc"
    fili2 = casename(i)+".q1.1981_2005.monmean.v2.nc"
    fili3 = casename(i)+".VgradT.1981_2005.monmean.v2.nc"
    fili4 = casename(i)+".Vertical.1981_2005.monmean.v2.nc"
    fili5 = casename(i)+".cam.h0.Q1.1979-2005.nc"
    fili6 = casename(i)+".cam.h0.PS.1979-2005.nc"
    f1 = addfile(diri1+fili1, "r")
    f2 = addfile(diri1+fili2, "r")
    f3 = addfile(diri1+fili3, "r")
    f4 = addfile(diri1+fili4, "r")
    f5 = addfile(diri2+fili5, "r")
    f6 = addfile(diri2+fili6, "r")

    dTdt = f1->dsdt
    Q1_v1 = f2->q1
    dTdx = f3->adv_s
    dTdp = f4->dsdp
    ;dTdp = dTdp*(-1.)
    Q1_v1 = dTdx+dTdt+dTdp
    printVarSummary(Q1_v1)

    TIME = f6->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    Q1_v2 = f5->Q1(iYYYY,:,:,:)
    PS = f6->PS(iYYYY,:,:)
    printVarSummary(Q1_v2)

    x = new((/300,4,nlev,nlat,nlon/), float)
    ; x(:,0,:,:,:) = (/dTdt/)
    x(:,2,:,:,:) = (/Q1_v1/)
    x(:,0,:,:,:) = (/-1*dTdx/)
    x(:,1,:,:,:) = (/-1*dTdp/)
    x(:,3,:,:,:) = x(:,1,:,:,:) + x(:,2,:,:,:) +x(:,0,:,:,:)
    ; x(:,3,:,:,:) = (/Q1_v2/)
    copy_VarCoords(Q1_v2, x(:,0,:,:,:))

    x = x*G_CPD/G_G  ; convert units


    lev = Q1_v2&level
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
    vt(i,:,:,:) = dim_avg_n_Wrap(vt0, 0)
end do
vt@units = "W/m~S~2~C~"
printVarSummary(vt)

vt_diff = vt(0,:,:,:)
vt_diff = vt(1,:,:,:) - vt(0,:,:,:)

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)

;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200605/"
figname = season+"_thermal_equation_v2"

; lonstr = 60
; lonend = 120
; latstr = 10
; latend = 50
leftsti = (/"(a) -dTdxdy","(b) -dTdp","(c) Q1","(d) a+b+c"/) ;,"Q1_v2"

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
res@cnMinLevelValF = -20
res@cnMaxLevelValF = 20
res@cnLevelSpacingF = 2.

do i = 0, 2
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_contour_map(wks, vt_diff(i,:,:), res)
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.1
res@gsnLeftString = leftsti(3)
plot(3) = gsn_csm_contour_map(wks, vt_diff(3,:,:), res)
plot_topo(3) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(3), plot_topo(3))


pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = season+" thermal equation (600hPa to 300hPa)" ;
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,2/), pres)

end