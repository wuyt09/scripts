
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

nlat = 181
nlon = 360
nlev = 19


pens = 18
diri = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/ncl/20210102/ISM_onset/"

filipath = "./"

casename = (/"CTRL","TPQ1_double"/) ;"NOTP"
ncase = dimsizes(casename)

do i = 0,ncase-1
    fili = "1_TTGI_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".nc"
    f = addfile(diri+fili, "r")

    t = f->tropoT(:,pens:npenta-1,:,:)
    printVarSummary(t)

    tN = wgt_areaave_Wrap(t(:,:,{5:35},{40:100}), 1, 1, 0)
    tS = wgt_areaave_Wrap(t(:,:,{-15:5},{40:100}), 1, 1, 0)
    ttgi = tN - tS
    copy_VarCoords(tN, ttgi)

    diff = ttgi
    last = diff
    onset_time = new(nyear, integer)

    do n = 1, npenta-pens-2
        diff(:,n) = ttgi(:,n-1)*ttgi(:,n)
        last(:,n) = ttgi(:,n+1)*ttgi(:,n)
    end do
    printVarSummary(diff)

    do yy = 0, nyear-1
        a= ind((diff(yy,:).lt.0) .and. (ttgi(yy,:).gt.0) .and. (last(yy,:).gt.0))
        onset_time(yy) =(/a(0)/) + pens +1
        delete(a)
    end do
    print(onset_time)

    filiname = "TTGI_onset_pentad_"+yrStrt+"_"+yrLast+"_"+casename(i)+".txt"
    write_table(filipath+filiname, "w", [/"YEAR","Pentad"/], "%4s %4s")
    write_table(filipath+filiname, "a", [/year, onset_time/], "%0.2i %0.2i")
end do



end