
begin

yrStrt = 1979
yrLast = 2005
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
nlev = 17


casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_PO","Hist_FAMIP_CMIP6"/) ;"NOTP"
ncase = dimsizes(casename)

outpath = "/home/yangsong3/wuyt/sysu/scripts/TP_FAMIP/"


do ll = 0, 0 ;ncase-1
    ps = new((/nyear,npenta,nlat,nlon/), float)
    column_qu = new((/nyear,npenta,nlat,nlon/), float)

    diri = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(ll)+"/a20191206/"

    filips = casename(ll)+".cam.h1.PS."+yrStrt+"-"+yrLast+".nc"
    fps = addfile(diri+filips, "r")
    PS = fps->PS
    do p = 0, 72
        ps(:,p,:,:) = dim_sum_n_Wrap(PS(:,5*p:5*p+4,:,:), 1)/5.0
    end do
    copy_VarCoords(PS(:,0,:,:), ps(:,0,:,:))
    ps!1 = "pentad"
    ps&pentad = penta

    do i = 0, nyear-1
        qu = new((/npenta,nlev,nlat,nlon/), float)
        yr = sprinti("%0.2i", year(i))
        filiu = casename(ll)+".cam.h1.U."+yr+".nc"
        filiq = casename(ll)+".cam.h1.Q."+yr+".nc"
        fu = addfile(diri+filiu, "r")
        fq = addfile(diri+filiq, "r")

        Q = fq->Q
        U = fu->U
        QU = Q*U
        copy_VarCoords(Q, QU)

        do p = 0, 72
            qu(p,:,:,:) = dim_sum_n_Wrap(QU(5*p:5*p+4,:,:,:), 0)/5.0
        end do
        copy_VarCoords(U(0,:,:,:), qu(0,:,:,:))
        qu!0 = "pentad"
        qu&pentad = penta

        lev = tofloat(qu&level)
        plev = lev*100
        plev@units = "Pa"
        psfc = ps(i,:,:,:)
        ptop = min(plev)
        dp = dpres_plevel_Wrap(plev, psfc, ptop, 0)

        qudp = qu*dp
        copy_VarCoords(qu, qudp)
        column_qu(ll,:,:,:) = dim_sum_n_Wrap(qudp(:,{1000:300},:,:), 1)
        copy_VarCoords(qudp(:,0,:,:), column_qu(ll,:,:,:))
        print("end of year "+year(i))
        delete([/Q,U,QU,qu,dp,qudp/])
    end do
    print("end of "+casename(ll))
    delete([/PS,ps/])


    column_qu!0 = "year"
    column_qu&year = year
    column_qu@units = "K"
    column_qu@long_name = "column integrated q*u derived from "+casename(ll)+" daily data"
    printVarSummary(column_qu)

    outname = "1_column_qu_pentad_"+yrStrt+"_"+yrLast+"_"+casename(ll)+".nc"

    system("rm -rf "+outpath+outname)

    fout = addfile(outpath+outname, "c")
    setfileoption(fout, "DefineMode", True)

    fout->column_qu = column_qu

    delete(column_qu)
end do

end

