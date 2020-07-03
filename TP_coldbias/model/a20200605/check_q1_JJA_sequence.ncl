
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

vt = new((/ncase,2,96,144/), float)

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(i)+"/a20191206/"
    fili1 = casename(i)+".cam.h0.QRL.1979-2005.nc"
    fili2 = casename(i)+".cam.h0.QRS.1979-2005.nc"
    fili3 = casename(i)+".cam.h0.DTV.1979-2005.nc"
    fili4 = casename(i)+".cam.h0.DTCOND.1979-2005.nc"
    fili5 = casename(i)+".cam.h0.Q1.1979-2005.nc"
    fili6 = casename(i)+".cam.h0.PS.1979-2005.nc"

    f1 = addfile(diri+fili1, "r")
    f2 = addfile(diri+fili2, "r")
    f3 = addfile(diri+fili3, "r")
    f4 = addfile(diri+fili4, "r")
    f5 = addfile(diri+fili5, "r")
    f6 = addfile(diri+fili6, "r")
    
    TIME = f6->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    x1 = f1->QRL(iYYYY,:,:,:)
    x2 = f2->QRS(iYYYY,:,:,:)
    x3 = f3->DTV(iYYYY,:,:,:)
    x4 = f4->DTCOND(iYYYY,:,:,:)
    x5 = f5->Q1(iYYYY,:,:,:)

    x = new((/300,2,nlev,nlat,nlon/), float)
    x(:,0,:,:,:) = x1+x2+x3+x4
    x(:,1,:,:,:) = x5
    copy_VarCoords(x3, x(:,0,:,:,:))
    printVarSummary(x)
    
    PS = f6->PS(iYYYY,:,:)

    x = x*G_CPD/G_G

    lev = x1&level
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
    xdp = x*dp
    copy_VarCoords(x, xdp)
    x0 = dim_sum_n_Wrap(xdp(:,:,{levs:levt},:,:), 2) ;/dim_sum_n_Wrap(dp(:,:,{levs:levt},:,:), 2)
    ; x0 = dim_sum_n_Wrap(xdp, 2)/dim_sum_n_Wrap(dp, 2)
    copy_VarCoords(xdp(:,:,0,:,:), x0)
    vt0 = month_to_season(x0, season)
    vt(i,:,:,:) = dim_avg_n_Wrap(vt0, 0)
end do
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
figname = "check_q1_JJA_seqence"

lonstr = 60
lonend = 150
latstr = 10
latend = 50
;leftsti = (/"dTdt","Q1_v1","dTdx","dTdp","Q1_v2"/)

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(2, graphic)
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
; res@cnMinLevelValF = -50
; res@cnMaxLevelValF = 50
; res@cnLevelSpacingF = 5

do i = 0, 1
    plot(i) = gsn_csm_contour_map(wks, vt(0,i,:,:), res) ;vt_diff
    plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    overlay(plot(i), plot_topo(i))
end do

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
; pres@gsnPanelMainString = season+" thermal equation (600hPa to 300hPa)" ;
; pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbLabelFontHeightF = 0.015

gsn_panel(wks, plot, (/2,1/), pres)

end