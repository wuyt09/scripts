begin

yrStrt = 198212
yrLast = 198301

fin = addfile("/home/yangsong3/huxm/antarctic_melt/199012_201001_ta_4time_daily.nc", "r")
lat = fin->latitude(::-1)
lon = fin->longitude
printVarSummary(lat)

filiname = (/"rsdt"/) ;,"t2m"
varname = (/"tisr"/) ;,"t2m"
nvar = dimsizes(varname)

outpath = "./"
outname = (/"tisr_warm.dat"/) ;,"t2m_warm.dat"

diri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/"

nlat = 181
nlon = 360
x_med = new((/62,nlat,nlon/), float)
do i = 0, nvar-1
    fili1 = diri+filiname(i)+"/"+filiname(i)+".1982-12.daily.nc"
    fili2 = diri+filiname(i)+"/"+filiname(i)+".1983-01.daily.nc"
    f1 = addfile(fili1, "r")
    f2 = addfile(fili2, "r")
    x1 = short2flt(f1->$varname(i)$(:,::-1,:))
    x2 = short2flt(f2->$varname(i)$(:,::-1,:))
    printVarSummary(x1)
    printVarSummary(x2)
    x_med(0:30,:,:) = x1
    x_med(31:61,:,:) = x2
    copy_VarCoords(x1(0,:,:), x_med(0,:,:))

    x_med = x_med/3600.
    x0 = linint2_Wrap(x_med&longitude, x_med&latitude, x_med, True, lon, lat, 0)
    x = x0(:,::-1,:)
    printVarSummary(x)

    system("rm -f "+outpath+outname(i))
    fbindirwrite(outpath+outname(i),x)

;===========================================
wks = gsn_open_wks("pdf", varname(i))

plot = new(1, graphic)

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

; res@gsnLeftString = "APHRO trend ("+yrStrt+"-"+yrLast+")"
plot(0) = gsn_csm_contour_map(wks, x(22,:,:), res)
draw(plot)
end do
end
