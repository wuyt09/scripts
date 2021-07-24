begin

yrStrt = 198212
yrLast = 198301


filiname = (/"cl","cli","clw","temp","q","ozone"/)
varname = (/"cc","ciwc","clwc","t","q","o3"/)
nvar = dimsizes(varname)

outpath = "./"
outname = (/"cc_warm.dat","ciwc_warm.dat","clwc_warm.dat","ta_warm.dat","hus_warm.dat","o3_warm.dat"/)

diri = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/"

nlat = 181
nlon = 360
nlev = 37
x = new((/62,nlev,nlat,nlon/), float)
do i = 0, nvar-1
    fili1 = diri+filiname(i)+"/"+filiname(i)+".1982-12.daily.nc"
    fili2 = diri+filiname(i)+"/"+filiname(i)+".1983-01.daily.nc"
    f1 = addfile(fili1, "r")
    f2 = addfile(fili2, "r")
    x1 = short2flt(f1->$varname(i)$)
    x2 = short2flt(f2->$varname(i)$)
    printVarSummary(x1)
    printVarSummary(x2)
    x(0:30,:,:,:) = x1
    x(31:61,:,:,:) = x2
    copy_VarCoords(x1(0,:,:,:), x(0,:,:,:))
    printVarSummary(x)

    system("rm -f "+outpath+outname(i))
    fbindirwrite(outpath+outname(i),x)
end do
end