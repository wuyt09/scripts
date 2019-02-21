
begin

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
fili1 = "B2000_f19g16_CP_CTRL.cam.h0.PS.200-349.nc"

diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
fili2 = "TP_noTOPO_hybrid.cam.h0.0250-01.nc"

f1 = addfile(diri1+fili1, "r")
f2 = addfile(diri2+fili2, "r")

ps1 = f1->PS
ps2 = f2->PS

;ps1&lat@units = "degress_north"
;ps1&lon@units = "degrees_east"
;ps2&lat@units = "degress_north"
;ps2&lon@units = "degrees_east"
ps1 = ps1/100
ps2 = ps2/100
printVarSummary(ps1)
diff = ps2(0,:,:) - ps1(50,:,:)
copy_VarCoords(ps1(0,:,:), diff)

;======================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_B/check_heating/"
figname = "check_PS_CTRL_noTOPO"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(3, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@cnFillOn = True
res@cnLineLabelsOn = False

;res@lbLabelBarOn = False
res@gsnRightString = "hPa"
res@gsnStringFontHeightF = 0.035
res@mpGeophysicalLineThicknessF = 3
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

res@cnFillPalette = "cmp_b2r"
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/500,700,850,900,950,1000,1050/)
res@cnFillColors = (/4,9,13,18,38,44,49,56/)

res@gsnLeftString = "CTRL"
plot(0) = gsn_csm_contour_map(wks, ps1(50,:,:), res)
res@gsnLeftString = "noTOPO"
plot(1) = gsn_csm_contour_map(wks, ps2(0,:,:), res)

delete(res@cnFillColors)
res@cnLevelSelectionMode = "ManualLevels"
res@cnMaxLevelValF = 20
res@cnMinLevelValF = -20
res@cnLevelSpacingF = 5
res@gsnLeftString = "diff"
plot(2) = gsn_csm_contour_map(wks, diff, res)

pres = True
pres@gsnMaximize = True
;pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/3,1/), pres)

end