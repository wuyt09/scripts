
begin

lonstr = 50
lonend = 120
latstr = -5
latend = 50

diri = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
fili1 = "TPmask_1000m_C20190410.nc"
fili2 = "India_mask_C20190710.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")

hgt1 = f1->tp
hgt2 = f2->india
printVarSummary(hgt1)

;=================================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/model_need_input_data/"
figname = "check_mask_India_TP"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new(2, graphic)

;---------------------------------------

res                 =  True
res@gsnDraw         = False
res@gsnFrame        = False
res@gsnAddCyclic    = True
;res@gsnSpreadColors = True

res@mpOutlineOn  = True
res@mpFillOn     = True
;res@mpCenterLonF = 0.
res@mpMinLatF    = latstr
res@mpMaxLatF    = latend
res@mpMinLonF    = lonstr
res@mpMaxLonF    = lonend
;res@mpCenterLonF = 180
res@mpGeophysicalLineThicknessF = 1

res@cnFillOn     = False
res@cnLinesOn    = True
res@cnLineLabelsOn         = False
res@cnInfoLabelOn          = False
;res@cnLabelDrawOrder = "PostDraw"
res@cnLineLabelFontHeightF = 0.025
res@cnLineLabelFontColor = "red"

res@cnLineThicknessF      = 3
;res@cnLineColor            = "black"
res@cnLineDashPatterns = (/0,2,4/)
res@cnLineColor = "green4"
;res@cnFillMode = "CellFill"
res@cnFillPalette = "BlueWhiteOrangeRed"

res@gsnStringFontHeightF = 0.035
res@gsnRightString = "m"

res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnLevelSelectionMode   ="ExplicitLevels"
res@cnLevels               = (/0.2,0.4,0.6,0.8,0.99/)
res@gsnLeftString = "Mask data"
plot(0) = gsn_csm_contour_map(wks, hgt1, res)

delete(res@gsnLeftString)
delete(res@gsnRightString)
delete(res@mpOutlineOn)
delete(res@mpFillOn)
delete(res@mpGeophysicalLineThicknessF)
delete(res@mpMaxLatF)
delete(res@mpMinLatF)
delete(res@mpMaxLonF)
delete(res@mpMinLonF)
res@lbLabelBarOn = False
plot(1) = gsn_csm_contour(wks, hgt2, res)

overlay(plot(0), plot(1))
draw(plot(0))
frame(wks)


end


