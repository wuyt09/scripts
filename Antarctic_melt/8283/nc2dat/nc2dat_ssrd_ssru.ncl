begin

yrStrt = 198212
yrLast = 198301


filiname = (/"rsds","rss"/)
varname = (/"ssrd","ssr"/)
nvar = dimsizes(varname)

outpath = "./"
outname = (/"ssrd_warm.dat","ssru_warm.dat"/)

diri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/"

nlat = 181
nlon = 360
x = new((/nvar,62,nlat,nlon/), float)
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
    printVarSummary(x)
end do
ssrd = x(0,:,:,:)
ssru = x(0,:,:,:) - x(1,:,:,:)

system("rm -f "+outpath+outname(0))
fbindirwrite(outpath+outname(0),ssrd)
system("rm -f "+outpath+outname(1))
fbindirwrite(outpath+outname(1),ssru)
end