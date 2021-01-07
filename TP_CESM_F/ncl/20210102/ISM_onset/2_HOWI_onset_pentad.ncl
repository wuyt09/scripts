
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

filipath = "./"

casename = (/"CTRL","TPQ1_double"/) ;"NOTP"
ncase = dimsizes(casename)

do i = 0,ncase-1
    fili = "1_column_qu_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".nc"
    f = addfile(diri+fili, "r")

    qu = f->column_qu(:,pens:npenta-1,:,:)
    printVarSummary(qu)

    qdp_wgt = wgt_areaave_Wrap(qu(:,:,{5:20},{45:80}), 1, 1, 0)

    qdpmax = max(qdp_wgt)
    qdpmin = min(qdp_wgt)

    howi = 2*((qdp_wgt - qdpmin)/(qdpmax - qdpmin))-1
    copy_VarCoords(qdp_wgt, howi)
    printVarSummary(howi)

    diff = howi
    last = diff
    onset_time = new(nyear, integer)

    do n = 1, npenta-pens-2
        diff(:,n) = howi(:,n-1)*howi(:,n)
        last(:,n) = howi(:,n+1)*howi(:,n)
    end do
    printVarSummary(diff)

    do yy = 0, nyear-1
        a= ind((diff(yy,:).lt.0) .and. (howi(yy,:).gt.0) .and. (last(yy,:).gt.0))
        onset_time(yy) =(/a(0)/) + pens +1
        delete(a)
    end do
    print(onset_time)

    filiname = "HOWI_onset_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".txt"
    write_table(filipath+filiname, "w", [/"YEAR","Pentad"/], "%4s %4s")
    write_table(filipath+filiname, "a", [/year, onset_time/], "%0.2i %0.2i")
end do
end
