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

mdiri = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/temp/"

t = new((/npenta,nlev,nlat,nlon/), float)
t_tropo = new((/nyear,npenta,nlat,nlon/), float)

do i = 0, nyear-1
    filit = systemfunc("ls "+mdiri+"temp."+year(i)+"*.daily.nc")
    ft = addfiles(filit, "r")    
    TX = short2flt(ft[:]->t(:,:,::-1,:))

    if (dimsizes(TX&time).eq.366) then
        T = TX(0:364,:,:,:)
        T(0:58,:,:,:) = TX(0:58,:,:,:)
        T(59:364,:,:,:) = TX(60:365,:,:,:)
    else
        T = TX
    end if
    ; printVarSummary(T)
    delete([/TX/])

    do p = 0, 72
        t(p,:,:,:) = dim_sum_n_Wrap(T(5*p:5*p+4,:,:,:), 0)/5.0
    end do
    copy_VarCoords(T(0,:,:,:), t(0,:,:,:))
    t!0 = "pentad"
    t&pentad = penta

    lev = tofloat(t&level)
    plev = lev*100
    plev@units = "Pa"
    psfc = max(plev)
    ptop = min(plev)

    dP = dpres_plevel(plev, psfc, ptop, 0)
    dP!0 = "lev"
    dP&lev = lev
    dp = conform_dims(dimsizes(t), dP, 1)
    tdp = t*dp
    copy_VarCoords(t, tdp)
    t_tropo(i,:,:,:) = dim_sum_n_Wrap(tdp(:,{500:200},:,:), 1)/sum(dP({500:200}))
    copy_VarCoords(tdp(:,0,:,:), t_tropo(i,:,:,:))
end do
t_tropo!0 = "year"
t_tropo&year = year
printVarSummary(t_tropo)
t_tropo@units = "K"
t_tropo@long_name = "troposphere (500hPa-200hPa) temperature derived from ERA5 daily data"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "1_MTG_pentad_ERA5_1979_2017.nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->MTG = t_tropo

end

