
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

diriu = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/uwind/"
diriq = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/q/"

U = new((/365,nlev,nlat,nlon/), float)
Q = U

qu = new((/npenta,nlev,nlat,nlon/), float)
column_qu = new((/nyear,npenta,nlat,nlon/), float)

do i = 0, nyear-1
    irecS = 0
    do j = 0, nmonth-1
        mth = sprinti("%0.2i", month(j))
        filiu = "uwind."+year(i)+"-"+mth+".daily.nc"
        filiq = "q."+year(i)+"-"+mth+".daily.nc"
        ; print(filiu)

        fu = addfile(diriu+filiu, "r")    
        UX := short2flt(fu->u(:,::-1,::-1,:))

        fq = addfile(diriq+filiq, "r")    
        QX := short2flt(fq->q(:,::-1,::-1,:))

        ndays := dimsizes(UX&time)
        if (ndays.eq.29) then
            irecL = irecS+ndays-2
            U(irecS:irecL,:,:,:) = UX(0:27,:,:,:)
            Q(irecS:irecL,:,:,:) = QX(0:27,:,:,:)
            irecS = irecS+ndays-1
        else
            irecL = irecS+ndays-1
            U(irecS:irecL,:,:,:) = UX 
            Q(irecS:irecL,:,:,:) = QX 
            irecS = irecS+ndays
        end if
        ; print(irecL)
    end do
    QU = Q
    QU = Q*U
    do p = 0, 72
        qu(p,:,:,:) = dim_sum_n_Wrap(QU(5*p:5*p+4,:,:,:), 0)/5.0
    end do
    copy_VarCoords(U(0,:,:,:), qu(0,:,:,:))
    qu!0 = "pentad"
    qu&pentad = penta


    lev = tofloat(qu&level)
    plev = lev*100
    plev@units = "Pa"
    psfc = max(plev)
    ptop = min(plev)

    dP = dpres_plevel(plev, psfc, ptop, 0)
    dP!0 = "lev"
    dP&lev = lev
    dp = conform_dims(dimsizes(qu), dP, 1)
    qudp = qu*dp
    copy_VarCoords(qu, qudp)
    column_qu(i,:,:,:) = dim_sum_n_Wrap(qudp(:,{1000:300},:,:), 1)
    copy_VarCoords(qudp(:,0,:,:), column_qu(i,:,:,:))
    print("end of year "+year(i))
end do
column_qu!0 = "year"
column_qu&year = year
column_qu@units = "K"
column_qu@long_name = "column integrated q*u derived from ERA5 daily data"
printVarSummary(column_qu)

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/ISM/"
outname = "1_column_qu_pentad_"+yrStrt+"_"+yrLast+".nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->column_qu = column_qu

end

