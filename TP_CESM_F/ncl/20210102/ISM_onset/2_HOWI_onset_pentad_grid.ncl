
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


pens = 18
diri = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/ncl/20210102/ISM_onset/"

outpath = diri

casename = (/"CTRL","TPQ1_double"/) ;"NOTP"
ncase = dimsizes(casename)

do i = 0,ncase-1
    fili = "1_column_qu_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".nc"
    f = addfile(diri+fili, "r")

    qu = f->column_qu(:,pens:npenta-1,:,:)
    printVarSummary(qu)

    qdp_wgt = dim_avg_n_Wrap(qu, 0) ;wgt_areaave_Wrap(qu(:,:,{5:20},{45:80}), 1, 1, 0)
    printVarSummary(qdp_wgt)

    ; howi = new((/npenta-pens/), float)
    diff = new((/npenta-pens/), float)
    last = diff
    onset_time = new((/nlat,nlon/), integer)
    copy_VarCoords(qdp_wgt(0,:,:), onset_time)

    do ilat = 0,nlat-1
        do ilon = 0,nlon-1
            var = qdp_wgt(:,ilat,ilon)
            qdpmax = max(var)
            qdpmin = min(var)
            ; printVarSummary(qdpmin)

            howi = 2*((var - qdpmin)/(qdpmax - qdpmin))-1
            copy_VarCoords(var, howi)
            ; printVarSummary(howi)


            do n = 1, npenta-pens-2
                diff(n) = howi(n-1)*howi(n)
                last(n) = howi(n+1)*howi(n)
            end do
            ; printVarSummary(diff)


            a= ind((diff.lt.0) .and. (howi.gt.0) .and. (last.gt.0))
            onset_time(ilat,ilon) =(/a(0)/) + pens +1
            delete(a)

        end do
    end do
    ; print(onset_time)

    onset_time@units = "penta"
    onset_time@long_name = "the onset pentad define by HOWI reverse"

    outname = "2_HOWI_onset_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".nc"

    system("rm -rf "+outpath+outname)

    fout = addfile(outpath+outname, "c")
    setfileoption(fout, "DefineMode", True)

    fout->HOWI = onset_time
end do
end
