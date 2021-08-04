begin

yrStrt = 198212
yrLast = 198301

fin = addfile("/home/yangsong3/huxm/antarctic_melt/199012_201001_ta_4time_daily.nc", "r")
lat = fin->latitude(::-1)
lon = fin->longitude
printVarSummary(lat)

filiname = (/"rsds","rss"/)
varname = (/"ssrd","ssr"/)
nvar = dimsizes(varname)

outpath = "./"
outname = (/"ssrd_warm.dat","ssru_warm.dat"/)

diri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/"

nlat = 181
nlon = 360
x_med = new((/nvar,62,nlat,nlon/), float)
do i = 0, nvar-1
    fili1 = diri+filiname(i)+"/"+filiname(i)+".1982-12.daily.nc"
    fili2 = diri+filiname(i)+"/"+filiname(i)+".1983-01.daily.nc"
    f1 = addfile(fili1, "r")
    f2 = addfile(fili2, "r")
    x1 = short2flt(f1->$varname(i)$(:,::-1,:))
    x2 = short2flt(f2->$varname(i)$(:,::-1,:))
    printVarSummary(x1)
    printVarSummary(x2)
    x_med(i,0:30,:,:) = x1
    x_med(i,31:61,:,:) = x2
    copy_VarCoords(x1(0,:,:), x_med(0,0,:,:))

end do
x_med = x_med/3600.
x0 = linint2_Wrap(x_med&longitude, x_med&latitude, x_med, True, lon, lat, 0)
x = x0(:,:,::-1,:)
printVarSummary(x)
ssrd = x(0,:,:,:)
ssru = x(0,:,:,:) - x(1,:,:,:)

system("rm -f "+outpath+outname(0))
fbindirwrite(outpath+outname(0),ssrd)
system("rm -f "+outpath+outname(1))
fbindirwrite(outpath+outname(1),ssru)

;===========================================
wks = gsn_open_wks("pdf", "ssrd_ssru")

plot = new(2, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnAddCyclic = False
res@gsnMaximize = True

; res@mpMinLatF = latstr
; res@mpMaxLatF = latend
; res@mpMinLonF = lonstr
; res@mpMaxLonF = lonend

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red"
; res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 0
res@cnMaxLevelValF = 320
res@cnLevelSpacingF = 20
res@cnFillPalette = "BlueWhiteOrangeRed"

res@gsnLeftString = "ssrd"
plot(0) = gsn_csm_contour_map(wks, ssrd(22,:,:), res)
res@gsnLeftString = "ssru"
plot(1) = gsn_csm_contour_map(wks, ssru(22,:,:), res)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,1/), pres)
end