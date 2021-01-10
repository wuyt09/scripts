
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


pens = 18
diri = "/home/yangsong3/wuyt/sysu/scripts/TP_FAMIP/"

outpath = diri

casename = (/"Ctrl_FAMIP","Hist_FAMIP","Hist_PO"/) ;,"Hist_FAMIP_CMIP6"
ncase = dimsizes(casename)

do i = 0, ncase-1
    fili = "1_MTG_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".nc"
    f = addfile(diri+fili, "r")

    MTG = f->MTG(2:nyear-1,pens:npenta-1,:,:)
    printVarSummary(MTG)

    mtg = dim_avg_n_Wrap(MTG, 0) ;wgt_areaave_Wrap(qu(:,:,{5:20},{45:80}), 1, 1, 0)
    printVarSummary(mtg)


    ; howi = new((/npenta-pens/), float)
    diff = new((/npenta-pens/), float)
    last = diff
    onset_time = new((/nlat,nlon/), integer)
    copy_VarCoords(mtg(0,:,:), onset_time)

    do ilat = 0,nlat-1
        do ilon = 0,nlon-1
            var = mtg(:,ilat,ilon)
            do n = 1, npenta-pens-2
                diff(n) = var(n-1)*var(n)
                last(n) = var(n+1)*var(n)
            end do
            a= ind((diff.lt.0) .and. (var.gt.0) .and. (last.gt.0))
            onset_time(ilat,ilon) =(/a(0)/) + pens +1
            delete(a)

        end do
    end do
    ; print(onset_time)

    onset_time@units = "penta"
    onset_time@long_name = "the onset pentad define by MTG reverse"

    outname = "2_MTG_onset_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".nc"

    system("rm -rf "+outpath+outname)

    fout = addfile(outpath+outname, "c")
    setfileoption(fout, "DefineMode", True)

    fout->MTG = onset_time
end do
end
