begin

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

mdiri = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/uwind/"
LEV = 850

u = new((/nyear,npenta,nlat,nlon/), float)

do i = 0, nyear-1
    filiu = systemfunc("ls "+mdiri+"uwind."+year(i)+"*.daily.nc")
    fu = addfiles(filiu, "r")    
    UX = short2flt(fu[:]->u(:,{LEV},::-1,:))

    if (dimsizes(UX&time).eq.366) then
        U = UX(0:364,:,:)
        U(0:58,:,:) = UX(0:58,:,:)
        U(59:364,:,:) = UX(60:365,:,:)
    else
        U = UX
    end if
    ; printVarSummary(U)
    delete([/UX/])

    do p = 0, 72
        u(i,p,:,:) = dim_sum_n_Wrap(U(5*p:5*p+4,:,:), 0)/5.0
    end do
    copy_VarCoords(U(0,:,:), u(0,0,:,:))
end do
u!0 = "year"
u&year = year
u!1 = "pentad"
u&pentad = penta
printVarSummary(u)

u@units = "m/s"
u@long_name = "zonal wind at 850 hPa derived from ERA5 daily data"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "1_U850_pentad_ERA5_1979_2017.nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->u850 = u


end