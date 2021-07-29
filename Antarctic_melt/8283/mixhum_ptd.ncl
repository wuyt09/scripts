
begin

yrStrt = 198212
yrLast = 198301

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
    fili1 = diri+filiname(i)+"/"+filiname(i)+".1982-12.daily.nc"
    fili2 = diri+filiname(i)+"/"+filiname(i)+".1983-01.daily.nc"
    f1 = addfile(fili1, "r")
    f2 = addfile(fili2, "r")
    x1 = short2flt(f1->$varname(i)$)
    x2 = short2flt(f2->$varname(i)$)
    printVarSummary(x1)
    printVarSummary(x2)
    x_med(i,0:30,:,:) = x1
    x_med(i,31:61,:,:) = x2
    copy_VarCoords(x1(0,:,:), x_med(0,0,:,:))

    x = linint2_Wrap(x_med&longitude, x_med&latitude, x_med, True, lon, lat, 0)
end do

printVarSummary(x)

d2m = x(0,:,:,:)
ps  = x(1,:,:,:)
huss = mixhum_ptd(ps, d2m, 2)

system("rm -f "+outpath+outname(0))
fbindirwrite(outpath+outname(0),huss)
system("rm -f "+outpath+outname(1))
fbindirwrite(outpath+outname(1),ps)


end