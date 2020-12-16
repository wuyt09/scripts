
begin

yrStrt = 1979
yrLast = 2017
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = (/"DJF","MAM"/)

fx = "time_series_TP_t2m_DJF_MAM.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(:,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))


;=======================================================================
diri = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
fili1 = "3_onset_pentad_for_grid.nc"
fili2 = "multi_standard_onset_BOB_SCS_ISM.txt"

f1 = addfile(diri+fili1, "r")
peak = f1->peak
printVarSummary(peak)

peak_hig = dim_avg_n_Wrap(peak(hig_ind,:,:), 0)
peak_low = dim_avg_n_Wrap(peak(low_ind,:,:), 0)
peak_diff = peak_hig
peak_diff = peak_hig - peak_low

ASMX = asciiread(diri+fili2, (/nyear,4/), "integer")
; print(ASMX(:,0))
ASM = ASMX(:,1:3)

hig_onset = dim_avg_n_Wrap(ASM(hig_ind,:), 0)
low_onset = dim_avg_n_Wrap(ASM(low_ind,:), 0)
print(hig_onset)
print(low_onset)

;===================================================================
; plot
;===================================================================
lonstr = 50
lonend = 150
latstr = 0
latend = 30

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201216/"
figname = "Comp_TP_t2m_DJF_ASM_onset"
wks = gsn_open_wks("pdf", figpath+figname)
gsn_define_colormap(wks,"so4_23")  ; choose colormap

plot = new(3, graphic)


res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False
res@cnFillMode = "CellFill"

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 20
res@cnMaxLevelValF = 40
res@cnLevelSpacingF = 2
res@lbLabelBarOn = False

res@gsnLeftString = "(a) High DJF TP t2m"
plot(0) = gsn_csm_contour_map(wks, peak_hig, res)

res@gsnLeftString = "(b) Low DJF TP t2m"
plot(1) = gsn_csm_contour_map(wks, peak_low, res)


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF = 1
res@lbLabelBarOn = True
res@cnFillPalette = "NCV_jaisnd"
res@gsnLeftString = "(a) - (b)"
plot(2) = gsn_csm_contour_map(wks, peak_diff, res)


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarOrthogonalPosF = -0.015
gsn_panel(wks, plot, (/3,1/), pres)
end