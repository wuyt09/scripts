
begin

yrStrt = 1979
yrLast = 2018
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

T = new((/365,nlev,nlat,nlon/), float)
t = new((/npenta,nlev,nlat,nlon/), float)
tropot = new((/nyear,npenta,nlat,nlon/), float)

do i = 0, nyear-1
    irecS = 0
    do j = 0, nmonth-1
        mth = sprinti("%0.2i", month(j))
        filiu = "temp."+year(i)+"-"+mth+".daily.nc"
        ; print(filiu)

        ft = addfile(mdiri+filiu, "r")    
        UX := short2flt(ft->t(:,::-1,::-1,:))

        ndays := dimsizes(UX&time)
        if (ndays.eq.29) then
            irecL = irecS+ndays-2
            T(irecS:irecL,:,:,:) = UX(0:27,:,:,:)
            irecS = irecS+ndays-1
        else
            irecL = irecS+ndays-1
            T(irecS:irecL,:,:,:) = UX 
            irecS = irecS+ndays
        end if
        ; print(irecL)
    end do
    ; print(irecL)
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
    tropot(i,:,:,:) = dim_sum_n_Wrap(tdp(:,{600:200},:,:), 1)/sum(dP({600:200}))
    copy_VarCoords(tdp(:,0,:,:), tropot(i,:,:,:))
    print("end of year "+year(i))
end do
tropot!0 = "year"
tropot&year = year
tropot@units = "K"
tropot@long_name = "troposphere (600hPa-200hPa) temperature derived from ERA5 daily data"
printVarSummary(tropot)

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/ISM/"
outname = "1_TTGI_pentad_"+yrStrt+"_"+yrLast+".nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->tropoT = tropot

end