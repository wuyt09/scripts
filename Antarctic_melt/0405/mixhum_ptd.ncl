
begin

yrStrt = 200412
yrLast = 200501

fin = addfile("/home/yangsong3/huxm/antarctic_melt/199012_201001_ta_4time_daily.nc", "r")
lat = fin->latitude
lon = fin->longitude
printVarSummary(lat)

filiname = (/"dwt2m","ps"/)
varname = (/"d2m","sp"/)
nvar = dimsizes(varname)

outpath = "./"
outname = (/"huss_warm.dat","ps_warm.dat"/)

diri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/"

nlat = 181
nlon = 360
x_med = new((/2,62,nlat,nlon/), float)
do i = 0, nvar-1
    fili1 = diri+filiname(i)+"/"+filiname(i)+".2004-12.daily.nc"
    fili2 = diri+filiname(i)+"/"+filiname(i)+".2005-01.daily.nc"
    f1 = addfile(fili1, "r")
    f2 = addfile(fili2, "r")
    x1 = short2flt(f1->$varname(i)$)
    x2 = short2flt(f2->$varname(i)$)
    printVarSummary(x1)
    printVarSummary(x2)
    x_med(i,0:30,:,:) = x1
    x_med(i,31:61,:,:) = x2
    copy_VarCoords(x1(0,:,:), x_med(0,0,:,:))

    ; x = linint2_Wrap(x_med&longitude, x_med&latitude, x_med, True, lon, lat, 0)
    x = f2fsh_Wrap(x_med, (/dimsizes(lat),dimsizes(lon)/))
end do

printVarSummary(x)

d2m = x(0,:,:,:)
ps  = x(1,:,:,:)
huss = mixhum_ptd(ps, d2m, 2)
printVarSummary(ps)
system("rm -f "+outpath+outname(0))
fbindirwrite(outpath+outname(0),huss)
system("rm -f "+outpath+outname(1))
fbindirwrite(outpath+outname(1),ps)

;===========================================
wks = gsn_open_wks("pdf", "mixhum_ptd")

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

res@gsnLeftString = "huss"
plot(0) = gsn_csm_contour_map(wks, huss(22,:,:), res)
res@gsnLeftString = "ps"
plot(1) = gsn_csm_contour_map(wks, ps(22,:,:), res)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,1/), pres)

end