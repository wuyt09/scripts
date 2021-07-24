
begin

yrStrt = 198212
yrLast = 198301


filiname = (/"dwt2m","ps"/)
varname = (/"d2m","sp"/)
nvar = dimsizes(varname)

outpath = "./"
outname = (/"huss_warm.dat","ps_warm.dat"/)

diri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/"

nlat = 181
nlon = 360
x = new((/2,62,nlat,nlon/), float)
do i = 0, nvar-1
    fili1 = diri+filiname(i)+"/"+filiname(i)+".1982-12.daily.nc"
    fili2 = diri+filiname(i)+"/"+filiname(i)+".1983-01.daily.nc"
    f1 = addfile(fili1, "r")
    f2 = addfile(fili2, "r")
    x1 = short2flt(f1->$varname(i)$)
    x2 = short2flt(f2->$varname(i)$)
    printVarSummary(x1)
    printVarSummary(x2)
    x(i,0:30,:,:) = x1
    x(i,31:61,:,:) = x2
    copy_VarCoords(x1(0,:,:), x(0,0,:,:))
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