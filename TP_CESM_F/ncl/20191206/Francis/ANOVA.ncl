
begin
t = 25
I = 2 ;dynamical
J = 2 ;heating
LEV = 500

casename = (/"CTRL","NOTP","NOQ1","NOTP_NOQ1"/)
ncase = dimsizes(casename)
Z = new((/ncase,25,96,144/), float)

season = "JJA" ;(/"MAM","JJA","SON","DJF"/)
do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    filiz = "F2000_"+casename(i)+".plevel.Z3.mon.0601-3012.nc"
    ;filiz = "F2000_"+casename(i)+".TS.mon.0601-3012.nc"
    fz = addfile(diri+filiz, "r")
    ZX = fz->Z3(:,{LEV},:,:) ;240:299
    ;ZX = fz->TS
    printVarSummary(ZX)
    Z(i,:,:,:) = month_to_season(ZX, season)
end do
lat = ZX&lat
lon = ZX&lon
;Z = Z/9.8
z = new((/I,J,t,96,144/), float)
z(0,0,:,:,:) = Z(0,:,:,:)
z(1,0,:,:,:) = Z(1,:,:,:)
z(0,1,:,:,:) = Z(2,:,:,:)
z(1,1,:,:,:) = Z(3,:,:,:)
printVarSummary(z)

;====================================================================================
; SST = SSA + SSB + SSI + SSE
;====================================================================================
Zijo = dim_avg_n_Wrap(z, 2)

Zooo = dim_avg_n_Wrap(z, (/0,1,2/))

Zioo = dim_avg_n_Wrap(z, (/1,2/))

Zojox = dim_avg_n_Wrap(z, 2)
Zojo = dim_avg_n_Wrap(Zojox, 0)

Zoot = dim_avg_n_Wrap(z, (/0,1/))
printVarSummary(Zioo)
printVarSummary(Zojo)
printVarSummary(Zoot)

zooo = conform_dims(dimsizes(z), Zooo, (/3,4/))
SST = dim_sum_n_Wrap((z - zooo)^2, (/0,1,2/)) 

zioo = conform_dims(dimsizes(Zioo), Zooo, (/1,2/))
SSA = dim_sum_n_Wrap((Zioo - zioo)^2, 0)
SSA = J*t*SSA

zojo = conform_dims(dimsizes(Zojo), Zooo, (/1,2/))
SSB = dim_sum_n_Wrap((Zojo - zojo)^2, 0)
SSB = I*t*SSB

zooox = conform_dims(dimsizes(Zijo), Zooo, (/2,3/))
zioox = conform_dims(dimsizes(Zijo), Zioo, (/0,2,3/))
zojox = conform_dims(dimsizes(Zijo), Zojo, (/1,2,3/))
SSI = dim_sum_n_Wrap((Zijo - zioox - zojox + zooox)^2, (/0,1/))
SSI = t*SSI

zijo = conform_dims(dimsizes(z), Zijo, (/0,1,3,4/))
SSE = dim_sum_n_Wrap((z-zijo)^2, (/0,1,2/))

pp = z(0,0,0:3,:,:)
SST = where(SST.eq.0, SST@_FillValue, SST)
pp(0,:,:) = SSA/SST
pp(1,:,:) = SSB/SST
pp(2,:,:) = SSI/SST
pp(3,:,:) = SSE/SST
copy_VarCoords(z(0,0,0:2,:,:), pp)
printVarSummary(pp)

;====================================================================================
; F value
;====================================================================================
; sig = 0.05
; prob_dyn = 1 - ftest(SSA, I-1, SSE, I*J*(t-1), 0)
; prob_het = 1 - ftest(SSB, J-1, SSE, I*J*(t-1), 0)
; prob_int = 1 - ftest(SSI, (I-1)*(J-1), SSE, I*J*(t-1), 0)
; F_sig = 3.94 ; a = 1, b = 96
; F_dyn = SSA
; F_dyn = SSA*(I-1)/SSE*I*J*(t-1)
; F_het = SSB
; F_het = SSB*(J-1)/SSE*I*J*(t-1)
; F_int = SSI
; F_int = SSI*(I-1)*(J-1)/SSE*I*J*(t-1)
; printVarSummary(F_dyn)
; pp = z(0,0,0:2,:,:)
; pp(0,:,:) = F_dyn
; pp(1,:,:) = F_het
; pp(2,:,:) = F_int
; ;pp(3,:,:) = dim_avg_n_Wrap(Z(0,:,:,:), 0)
; copy_VarCoords(z(0,0,0:2,:,:), pp)
; pp = where(pp.ge.F_sig, pp, pp@_FillValue)
; printVarSummary(pp)
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20191206/Francis/"
figname = "ANOVA_Z500_"+season

lonstr = 0 ;20
lonend = 360 ;220
latstr = -30
latend = 80

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4 ;4*ncase
plot = new(nplot, graphic)
plot_topo = plot

 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "blue"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/0.9/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey45" ;"grey25"
res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend
res@mpCenterLonF = (lonstr+lonend)/2

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 200
; res@cnMaxLevelValF = 2000
; res@cnMinLevelValF = -2000
res@cnFillPalette = "WhiteYellowOrangeRed"

;leftsti = (/"SST","SSA","SSB","SSI","SSE"/)
leftsti = (/"Dynamical effects", "Heating effects", "Interaction effects", "Error"/)
do i = 0, 3
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_contour_map_ce(wks, pp(i,:,:), res)
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelMainString = "ANOVA results "+season+" mean: H"+LEV
pres@gsnPanelLabelBar = True
;pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/3,2/), pres)

end
