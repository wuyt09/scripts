
begin

yrStrt = 6
yrLast = 30
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)

nlat = 96
nlon = 144
nlev = 37



casename = (/"CTRL","TPQ1_double"/) ;"NOTP"
ncase = dimsizes(casename)

outpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/ncl/20210102/ISM_onset/"



do ll = 0,ncase-1
    t = new((/npenta,nlev,nlat,nlon/), float)
    tropot = new((/nyear,npenta,nlat,nlon/), float)

    ps = new((/nyear,npenta,nlat,nlon/), float)

    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(ll)+"_F2000_F19G16_CAM4/pro/daily/"

    filips = "F2000_"+casename(ll)+".PS.daily.06-30.nc"
    fps = addfile(diri+filips, "r")
    PS = fps->PS
    do p = 0, 72
        ps(:,p,:,:) = dim_sum_n_Wrap(PS(:,5*p:5*p+4,:,:), 1)/5.0
    end do
    copy_VarCoords(PS(:,0,:,:), ps(:,0,:,:))
    ps!1 = "pentad"
    ps&pentad = penta

    do i = 0, nyear-1
        yr = sprinti("%0.2i", year(i))
        filit = "F2000_"+casename(ll)+".plevel.T."+yr+".daily.nc"
        ft = addfile(diri+filit, "r")
        T = ft->T

        do p = 0, 72
            t(p,:,:,:) = dim_sum_n_Wrap(T(5*p:5*p+4,:,:,:), 0)/5.0
        end do
        copy_VarCoords(T(0,:,:,:), t(0,:,:,:))
        t!0 = "pentad"
        t&pentad = penta


        lev = tofloat(t&lev)
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
    tropot@long_name = "troposphere (600hPa-200hPa) temperature derived from "+casename(ll)+" daily data"
    printVarSummary(tropot)

    outname = "1_TTGI_pentad_"+yrStrt+"_"+yrLast+"_"+casename(ll)+".nc"

    system("rm -rf "+outpath+outname)

    fout = addfile(outpath+outname, "c")
    setfileoption(fout, "DefineMode", True)

    fout->tropoT = tropot
end do
end