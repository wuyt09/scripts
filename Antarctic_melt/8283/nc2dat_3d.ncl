begin

yrStrt = 198212
yrLast = 198301

fin = addfile("/home/yangsong3/huxm/antarctic_melt/199012_201001_ta_4time_daily.nc", "r")
lat = fin->latitude
lon = fin->longitude
printVarSummary(lat)

filiname = (/"cc","cli","clw","t","q","o3"/)
varname = (/"cc","ciwc","clwc","t","q","o3"/)
nvar = dimsizes(varname)

outpath = "./"
outname = (/"cc_warm.dat","ciwc_warm.dat","clwc_warm.dat","ta_warm.dat","hus_warm.dat","o3_warm.dat"/)

diri = "/home/yangsong3/wuyt/sysu/scripts/Antarctic_melt/data/"

nlat = 181
nlon = 360
nlev = 37
x_med = new((/62,nlev,nlat,nlon/), float)
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
    x1 = short2flt(f1->$varname(i)$(iYYYY1,:,:,:))
    x2 = short2flt(f2->$varname(i)$(iYYYY2,:,:,:))
    ; printVarSummary(x1)
    ; printVarSummary(x2)
    x_med(0:30,:,:,:) = x1
    x_med(31:61,:,:,:) = x2
    copy_VarCoords(x1(0,:,:,:), x_med(0,:,:,:))

    x = linint2_Wrap(x_med&longitude, x_med&latitude, x_med, True, lon, lat, 0)
    printVarSummary(x)

    system("rm -f "+outpath+outname(i))
    fbindirwrite(outpath+outname(i),x)
end do
end
