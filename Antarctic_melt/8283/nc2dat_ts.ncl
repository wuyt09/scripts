begin

yrStrt = 198212
yrLast = 198301

; fin = addfile("/home/yangsong3/huxm/antarctic_melt/199012_201001_ta_4time_daily.nc", "r")
; lat = fin->latitude
; lon = fin->longitude
; printVarSummary(lat)

filiname = (/"skt"/)
varname = (/"skt"/)
nvar = dimsizes(varname)

outpath = "./"
outname = (/"ts_warm.dat"/)

diri = "/home/yangsong3/wuyt/sysu/scripts/Antarctic_melt/data/"

nlat = 121
nlon = 240
x = new((/62,nlat,nlon/), float)
do i = 0, nvar-1
    fili1 = diri+filiname(i)+".1982_1991_2004_2015_12.daily.nc"
    fili2 = diri+filiname(i)+".1983_1992_2005_2016_01.daily.nc"
    f1 = addfile(fili1, "r")
    f2 = addfile(fili2, "r")
    TIME1 := f1->time
    YYYY1 := cd_calendar(TIME1, -1)
    iYYYY1 := ind(YYYY1.eq.yrStrt)
    printVarSummary(iYYYY1)
    TIME2 := f2->time
    YYYY2 := cd_calendar(TIME2, -1)
    iYYYY2 := ind(YYYY2.eq.yrLast)
    printVarSummary(iYYYY2)
    x1 = short2flt(f1->$varname(i)$(iYYYY1,:,:))
    x2 = short2flt(f2->$varname(i)$(iYYYY2,:,:))
    printVarSummary(x1)
    printVarSummary(x2)
    x(0:30,:,:) = x1
    x(31:61,:,:) = x2
    copy_VarCoords(x1(0,:,:), x(0,:,:))

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

plot(0) = gsn_csm_contour_map(wks, x(22,:,:), res)
draw(plot)
end do
end
