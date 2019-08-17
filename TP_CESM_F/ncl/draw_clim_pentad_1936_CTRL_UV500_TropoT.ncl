
begin

LEV = 500
nlat = 96
nlon = 144

lonstr = 50
lonend = 180
latstr = -10
latend = 40

diri = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/daily/"

tname = systemfunc("ls "+diri+"F2000_CTRL.plevel.T.*.daily.nc")
print(tname)
ft = addfiles(tname, "r")
TT = ft[:]->T
printVarSummary(TT)

TTx = reshape(TT, (/25,365,37,nlat,nlon/))
copy_VarCoords(TT(0:364,:,:,:), TTx(0,:,:,:,:))
T = dim_avg_n_Wrap(TTx, 0)
printVarSummary(T)


plev = T&lev*100
plev!0 = "lev"
plev&lev = T&lev
psfc = max(plev)
ptop = min(plev)
dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
dp = conform_dims(dimsizes(T), dP, 1)
Tdp = T
Tdp = T*dp
t = dim_sum_n_Wrap(Tdp(:,{500:200},:,:), 1)/sum(dP({500:200}))
copy_VarCoords(T(:,0,:,:), t)
printVarSummary(t)

pentad_t = t(0:72,:,:)

do p = 0,72  ; a year has 73 penta
    pentad_t(p,:,:) = dim_sum_n_Wrap(t(5*p:5*p+4,:,:), 0)/5.0
end do
;pentad_test = pentad_lp(:,0,:,:)
;pentad_test = (cp(20,0,{20:40},{130:140})+cp(20,1,{20:40},{130:140})+cp(20,2,{20:40},{130:140})+cp(20,3,{20:40},{130:140})+cp(20,4,{20:40},{130:140}))/5.0
;print(pentad_cp(20,0,{20:40},{130:140}))
;print(pentad_lp(20,0,{20:40},{130:140}))
;print(pentad_test)
;exit

clim_pentad_t = pentad_t ;dim_avg_n_Wrap(pentad_t, 0)
printVarSummary(clim_pentad_t)
delete([/TT,T,Tdp,t,pentad_t/])
;--------------------------------------------------------------------

uname = systemfunc("ls "+diri+"F2000_CTRL.plevel.U.*.daily.nc")
vname = systemfunc("ls "+diri+"F2000_CTRL.plevel.V.*.daily.nc")
print(uname)

fu = addfiles(uname, "r")
fv = addfiles(vname, "r")

U = fu[:]->U(:,{LEV},:,:)
V = fv[:]->V(:,{LEV},:,:)

u = reshape(U, (/25,365,nlat,nlon/))
v = reshape(V, (/25,365,nlat,nlon/))
copy_VarCoords(U(0:364,:,:), u(0,:,:,:))
copy_VarCoords(V(0:364,:,:), v(0,:,:,:))

pentad_u = u(:,0:72,:,:)
pentad_v = v(:,0:72,:,:)

do p = 0,72  ; a year has 73 penta
    pentad_u(:,p,:,:) = dim_sum_n_Wrap(u(:,5*p:5*p+4,:,:), 1)/5.0
    pentad_v(:,p,:,:) = dim_sum_n_Wrap(v(:,5*p:5*p+4,:,:), 1)/5.0
end do

clim_pentad_u = dim_avg_n_Wrap(pentad_u, 0)
clim_pentad_v = dim_avg_n_Wrap(pentad_v, 0)
printVarSummary(clim_pentad_u)
delete([/U,V,u,v,pentad_u,pentad_v/])
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_pentad_1936_CTRL_UV500_TropoT"

wks = gsn_open_wks("ps", figpath+figname)

nplot = 18
plot = new(nplot, graphic)
plotv = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
res@mpFillDrawOrder = "PreDraw"

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.035
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025
;------------------------------------------

  vcres                         = True             ; plot mods desired
  vcres@gsnFrame                = False            ; so we can draw time stamp
  vcres@gsnDraw                 = False
  vcres@vcRefAnnoOrthogonalPosF = -1.0  ;  -0.35           ; move ref vector up

  vcres@vcGlyphStyle            = "CurlyVector"  ;"LineArrow"    ; turn on curly vectors
  ;vcres@vcMinDistanceF          = 0.017
  ;vcres@vcGlyphStyle            = "LineArrow"
  vcres@vcMinDistanceF          = 0.03
  vcres@vcRefLengthF            = 0.025            ; define length of vec ref
  
  vcres@vcLineArrowThicknessF   = 3.0
  ;vecres@vcLineArrowHeadMinSizeF = 0.008
  ;vecres@vcLineArrowHeadMaxSizeF = 0.3

  vcres@vcRefAnnoString2On = False
  vcres@vcRefAnnoFontHeightF = 0.02
;------------------------------------------
vcres@vcRefMagnitudeF         = 5.0  ; define vector ref mag
res@lbLabelBarOn = False
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1.
res@cnMaxLevelValF = 250.
res@cnMinLevelValF = 240.
;res@cnFillPalette = "BlueWhiteOrangeRed"
res@gsnRightString = "K"


do i = 0,17

    res@gsnLeftString = "Pentad "+(i+19)

    plot(i) = gsn_csm_contour_map(wks, clim_pentad_t(i+18,:,:), res)
    plotv(i) = gsn_csm_vector(wks, clim_pentad_u(i+18,:,:), clim_pentad_v(i+18,:,:), vcres)
    overlay(plot(i), plotv(i))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "CTRL UV500_TropoT(500-200hPa)_clim"
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/5,4/), pres)


end