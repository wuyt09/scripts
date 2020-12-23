
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


pens = 18
diri = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
fili = "1_U850_pentad_ERA5_1979_2018.nc"
f = addfile(diri+fili, "r")

u850 = f->u850(:,pens:npenta-1,:,:)
printVarSummary(u850)

uS = wgt_areaave_Wrap(u850(:,:,{5:15},{40:80}), 1, 1, 0)
uN = wgt_areaave_Wrap(u850(:,:,{20:30},{70:90}), 1, 1, 0)
wfi = uS - uN
copy_VarCoords(uS, wfi)
printVarSummary(wfi)

diff = wfi
last = diff
onset_time = new(nyear, integer)

do n = 1, npenta-pens-2
    diff(:,n) = wfi(:,n-1)*wfi(:,n)
    last(:,n) = wfi(:,n+1)*wfi(:,n)
end do
printVarSummary(diff)

do yy = 0, nyear-1
    a= ind((diff(yy,:).lt.0) .and. (wfi(yy,:).gt.0) .and. (last(yy,:).gt.0))
    onset_time(yy) =(/a(0)/) + pens +1
    delete(a)
end do
print(onset_time)

filipath = "./"
filiname = "WFI_onset_pentad_"+yrStrt+"_"+yrLast+".txt"
write_table(filipath+filiname, "w", [/"YEAR","Pentad"/], "%4s %4s")
write_table(filipath+filiname, "a", [/year, onset_time/], "%0.2i %0.2i")

end