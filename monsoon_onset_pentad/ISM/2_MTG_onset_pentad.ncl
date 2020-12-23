
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

ism_lat1 = 10
ism_lat2 = 20
ism_lon1 = 60
ism_lon2 = 85


pens = 18
f = addfile("../1_MTG_pentad_ERA5_"+yrStrt+"_"+yrLast+".nc", "r")
MTGx = f->MTG
MTG = MTGx(:,pens:npenta-1,:,:)
printVarSummary(MTG)
delete(MTGx)

mtg = wgt_areaave_Wrap(MTG(:,:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)



diff = mtg
last = diff
onset_time = new(nyear, integer)

do n = 1, npenta-pens-2
    diff(:,n) = mtg(:,n-1)*mtg(:,n)
    last(:,n) = mtg(:,n+1)*mtg(:,n)
end do
printVarSummary(diff)

do yy = 0, nyear-1
    a= ind((diff(yy,:).lt.0) .and. (mtg(yy,:).gt.0) .and. (last(yy,:).gt.0))
    onset_time(yy) =(/a(0)/) + pens +1
    delete(a)
end do
print(onset_time)

filipath = "./"
filiname = "MTG_onset_pentad_"+yrStrt+"_"+yrLast+".txt"
write_table(filipath+filiname, "w", [/"YEAR","Pentad"/], "%4s %4s")
write_table(filipath+filiname, "a", [/year, onset_time/], "%0.2i %0.2i")

end