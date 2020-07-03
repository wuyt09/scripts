
begin

LEV = 850
nlat = 96
nlon = 144

lonstr = 50
lonend = 180
latstr = -10
latend = 40

diri = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/daily/"

cname = "F2000_CTRL.TS.daily.06-30.nc"
fc = addfile(diri+cname, "r")
cp = fc->TS
printVarSummary(cp)

pentad_cp = cp(:,0:72,:,:)
do p = 0,72  ; a year has 73 penta
    pentad_cp(:,p,:,:) = dim_sum_n_Wrap(cp(:,5*p:5*p+4,:,:), 1)/5.0
end do

clim_pentad_p = dim_avg_n_Wrap(pentad_cp, 0)
printVarSummary(clim_pentad_p)
;print(clim_pentad_p(25,{20:30},{120:130}))

;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/"
figname = "draw_clim_pentad_1936_CTRL_TS"

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
res@cnMaxLevelValF = 310.
res@cnMinLevelValF = 290.
;res@cnFillPalette = "CBR_wet"
;res@gsnRightString = "mm/day"


do i = 0,17

    res@gsnLeftString = "Pentad "+(i+19)

    plot(i) = gsn_csm_contour_map(wks, clim_pentad_p(i+18,:,:), res)

end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "CTRL TS_clim"
pres@pmLabelBarOrthogonalPosF = -0.01

gsn_panel(wks, plot, (/5,4/), pres)


end