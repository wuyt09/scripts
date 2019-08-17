
begin

LEV = 850
nlat = 96
nlon = 144

lonstr = 50
lonend = 180
latstr = -10
latend = 40

diri = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/daily/"

cname = "F2000_NOTP_NOQ1.PRECC.daily.06-30.nc"
lname = "F2000_NOTP_NOQ1.PRECL.daily.06-30.nc"

fc = addfile(diri+cname, "r")
fl = addfile(diri+lname, "r")

cp = fc->PRECC
lp = fl->PRECL
printVarSummary(cp)

pentad_cp = cp(:,0:72,:,:)
;new((/25,73,nlat,nlon/), float)
pentad_lp = pentad_cp

do p = 0,72  ; a year has 73 penta
    pentad_cp(:,p,:,:) = dim_sum_n_Wrap(cp(:,5*p:5*p+4,:,:), 1)/5.0
    pentad_lp(:,p,:,:) = dim_sum_n_Wrap(lp(:,5*p:5*p+4,:,:), 1)/5.0
end do
;pentad_test = pentad_lp(:,0,:,:)
;pentad_test = (cp(20,0,{20:40},{130:140})+cp(20,1,{20:40},{130:140})+cp(20,2,{20:40},{130:140})+cp(20,3,{20:40},{130:140})+cp(20,4,{20:40},{130:140}))/5.0
;print(pentad_cp(20,0,{20:40},{130:140}))
;print(pentad_lp(20,0,{20:40},{130:140}))
;print(pentad_test)
;exit

pentad_p = pentad_lp
pentad_p = (pentad_cp+pentad_lp) * 86400 * 1000. ; convert units from m/s to mm/day
clim_pentad_p = dim_avg_n_Wrap(pentad_p, 0)
printVarSummary(clim_pentad_p)

;--------------------------------------------------------------------

uname = systemfunc("ls "+diri+"F2000_NOTP_NOQ1.plevel.U.*.daily.nc")
vname = systemfunc("ls "+diri+"F2000_NOTP_NOQ1.plevel.V.*.daily.nc")
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
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_pentad_1936_NOTP_NOQ1_UV850_precip"

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
  vcres@vcRefAnnoOrthogonalPosF =  -0.35 ;-1.0             ; move ref vector up

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
res@cnLevelSpacingF = 2.
res@cnMaxLevelValF = 20.
res@cnMinLevelValF = 1.
res@cnFillPalette = "CBR_wet"
res@gsnRightString = "mm/day"


do i = 0,17

    res@gsnLeftString = "Pentad "+(i+19)

    plot(i) = gsn_csm_contour_map(wks, clim_pentad_p(i+18,:,:), res)
    plotv(i) = gsn_csm_vector(wks, clim_pentad_u(i+18,:,:), clim_pentad_v(i+18,:,:), vcres)
    overlay(plot(i), plotv(i))

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "NOTP_NOQ1 UV850_Precip_clim"
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/5,4/), pres)


end