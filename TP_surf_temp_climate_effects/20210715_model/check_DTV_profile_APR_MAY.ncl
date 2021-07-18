begin

nlat = 96
nlon = 144
nlev = 37
season = "JJA"

lonstr = 30
lonend = 150
latstr = 0
latend = 60

;===================================================================
; topo data
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo ;(::-1,:)
printVarSummary(topo)
tpmask = topo
;===================================================================

casename = (/"NOHEAT_CTRL","NOHEAT_ALL","NOHEAT_APR","NOHEAT_MAY"/) ;
ncase = dimsizes(casename)
season = "JJA"

STRT_ESM=1
END_ESM=30
ensemble = ispan(STRT_ESM, END_ESM, 1)
member = dimsizes(ensemble)

x = new((/ncase,6,member,nlev,nlat,nlon/), float)
x!0 = "case"
x&case = casename
x!1 = "time"
title = (/"March","April","May","June","July","JJA"/) ;,"August""Jan","Feb",
x&time = title
x!2 = "member"
x&member = ensemble
printVarSummary(x)

do i = 0, ncase-1
  do j = 0, member-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/Ensemble_"+ensemble(j)+"/"
    fname = casename(i)+".DTV.plevel.mon.51.nc"
    f = addfile(diri+fname, "r")
    X = f->DTV(:,:,::-1,:)
    X = X*86400 ;K/day

    x(i,0:4,j,:,:,:) = X(2:6,:,:,:)
    x_med = month_to_season(X, season)
    x(i,5,j,:,:,:) = (/x_med(0,:,:,:)/)
    copy_VarCoords(X(0,:,:,:), x(i,0,j,:,:,:))
  end do
end do
printVarSummary(x)

x_clm = dim_avg_n_Wrap(x, 2)
TOPO = linint2_Wrap(topo&lon, topo&lat, topo, True, x&lon, x&lat, 0)
x_clm = mask(x_clm, TOPO.ge.1500, True)
TP_clm = wgt_areaave_Wrap(x_clm(:,:,:,{25:40},{75:105}), 1., 1., 0)
printVarSummary(TP_clm)
lev = f->lev


;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surf_temp_climate_effects/20210715_model/"
figname = "check_DTV_profile_APR_MAY"
Colors = (/"black","red","blue","green4"/)

wks = gsn_open_wks("pdf", figpath+figname)

plot = new(6, graphic)
dum = new(4, graphic)
dum2 = dum


 res = True
 res@gsnDraw  = False
 res@gsnFrame = False

 res@gsnXRefLine = 0
 res@gsnXRefLineDashPattern = 2
 res@gsnXRefLineColor = "grey"
 res@gsnXRefLineThicknessF = 3
 res@trXMinF = -6.
 res@trXMaxF = 6
 res@trYMinF   = 100
 res@trYMaxF   = 1000
 res@trYLog = True
 
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)
 res@xyLineColors = Colors
 res@xyLineThicknessF = 5
 res@xyDashPattern = 0

 res@trYReverse = True
 
 res@gsnRightString = "K/day"
 res@gsnStringFontHeightF = 0.025

 do i = 0,5
    res@gsnLeftString = title(i)
    plot(i) = gsn_csm_xy(wks, TP_clm(:,i,{1000:50}), lev({1000:50}), res)
 end do

;===========================================================
; legend
;===========================================================
xpos = (/-5.5,-4.2/)
ypos = (/120,150,200,300/)
gsres = True
gsres@gsLineThicknessF = 5
gsres@gsLineDashPattern = 0

do i = 0, 3
    gsres@gsLineColor = Colors(i)
    dum(i) = gsn_add_polyline(wks, plot(0), xpos, (/ypos(i),ypos(i)/), gsres)
end do

Text = casename
txres = True
txres@txFontHeightF = 0.02
txres@txJust             = "CenterLeft"
tx_xpos = -4.0
tx_ypos = ypos
do i = 0, 3
    dum2(i) = gsn_add_text(wks, plot(0), Text(i), tx_xpos, tx_ypos(i), txres)
end do

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,3/), pres)

end