
begin

bob_lat1 = 5
bob_lat2 = 15
bob_lon1 = 90
bob_lon2 = 100

scs_lat1 = 10
scs_lat2 = 20
scs_lon1 = 110
scs_lon2 = 120

ism_lat1 = 10
ism_lat2 = 20
ism_lon1 = 60
ism_lon2 = 85


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

nlat = 181
nlon = 360
nlev = 19

diri = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
fili1 = "2_U850_onset_pentad_1979_2017.nc"
fili2 = "2_MTG_onset_pentad_1979_2017.nc"
fili3 = "2_precip_onset_pentad_GPCP_1979_2017.nc"

f1 = addfile(diri+fili1, "r")
f2 = addfile(diri+fili2, "r")
f3 = addfile(diri+fili3, "r")
x1 = f1->u850_peak
x2 = f2->MTG_peak
x3 = f3->prec_peak
printVarSummary(x1)

x = new((/3,nyear,nlat,nlon/), float)
x(0,:,:,:) = x1
x(1,:,:,:) = x2
x(2,:,:,:) = x3

x_max = dim_max_n_Wrap(x, 0)
; print("BOB "+x_max(:,0))
; print("SCS "+x_max(:,1))
; print("ISM "+x_max(:,2))
; x_max(20,0) = 19
; x_max(21,0) = 20

x_max@units = "penta"
x_max@long_name = "the onset pentad define by multi standard"
printVarSummary(x_max)

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "3_onset_pentad_for_grid.nc"
system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->peak = x_max
;=====================================================================
; plot
;=====================================================================
xmax_clm = dim_avg_n_Wrap(x_max, 0)

lonstr = 50
lonend = 150
latstr = 0
latend = 30

figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "3_onset_pentad_for_grid"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

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

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 20
res@cnMaxLevelValF = 40
res@cnLevelSpacingF = 2

plot(0) = gsn_csm_contour_map(wks, xmax_clm, res)
draw(plot)

end