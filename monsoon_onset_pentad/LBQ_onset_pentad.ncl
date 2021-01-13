
begin
;===========================================
; this was derived from JRA55 daily data,
; 3 standards
;===========================================
yrStrt = 1979
yrLast = 2017
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

bob_mtg_time = (/26,26,26,23,27,21,21,25,25,25,19,25,24,26,26,25,26,23,28,27,13,18,24,26,27,25,26,23,24,24,21,28,24,23,27,26,27,28,27/)
scs_mtg_time = (/27,28,30,26,29,28,29,27,29,29,28,28,27,28,29,25,28,26,26,29,24,25,26,27,28,28,27,26,28,26,22,29,26,25,27,27,28,29,29/)
ism_mtg_time = (/32,28,28,30,31,30,30,31,31,28,28,27,31,33,30,30,27,30,32,32,28,27,27,26,26,26,31,28,31,30,28,29,30,29,27,32,32,29,29/)

bob_pr_time = (/21,18,19,20,24,21,22,20,22,20,24,24,23,25,27,23,24,23,27,26,18,21,24,24,24,24,25,21,23,23,19,27,22,20,24,25,24,28,18/)
scs_pr_time = (/27,28,25,30,28,24,22,27,29,27,28,27,31,27,27,25,26,26,28,28,18,27,26,28,27,26,29,27,26,24,22,29,27,27,25,31,33,32,26/)
ism_pr_time = (/33,31,30,30,33,30,28,31,31,29,31,26,31,32,30,30,33,31,33,32,28,30,32,32,32,25,32,29,29,29,31,28,30,32,30,32,31,31,31/)

bob_wind_time = (/27,23,21,23,25,22,20,26,27,25,22,21,24,26,24,25,26,22,25,27,20,22,24,26,27,26,25,23,25,25,21,26,24,22,27,26,25,28,21/)
scs_wind_time = (/29,28,30,30,27,28,30,27,32,29,32,29,32,28,28,30,27,34,28,28,28,27,28,27,28,34,25,30,28,29,29,29,26,27,27,27,28,29,28/)
ism_wind_time = (/32,28,27,30,29,28,29,26,31,27,27,27,32,28,28,28,27,30,28,28,24,26,26,26,28,27,30,27,26,25,27,29,29,29,26,30,28,29,28/)

time = new((/3,3,nyear/), integer)
time(0,0,:) = bob_mtg_time
time(0,1,:) = bob_pr_time
time(0,2,:) = bob_wind_time
time(1,0,:) = scs_mtg_time
time(1,1,:) = scs_pr_time
time(1,2,:) = scs_wind_time
time(2,0,:) = ism_mtg_time
time(2,1,:) = ism_pr_time
time(2,2,:) = ism_wind_time

onset_time = dim_max_n(time, 1)
printVarSummary(onset_time)

bob_time = (/27,26,26,23,27,22,22,26,27,25,24,25,24,26,27,25,26,23,28,27,20,22,24,26,27,26,26,23,25,25,21,28,24,23,27,26,27,28,27/)
scs_time = (/29,28,30,30,29,28,30,27,32,29,32,29,32,28,29,30,28,34,28,29,28,27,28,28,28,34,29,30,28,29,29,29,27,27,27,31,33,32,29/)
ism_time = (/33,31,30,30,33,30,30,31,31,29,31,27,32,33,30,30,33,31,33,32,28,30,32,32,32,27,32,29,31,30,31,29,30,32,30,32,32,31,31/)

bob = onset_time(0,:) - bob_time
scs = onset_time(1,:) - scs_time
ism = onset_time(2,:) - ism_time
print(year+"    "+bob+"    "+scs+"    "+ism)

filipath = "./"
filiname = "LBQ_onset_pentad_"+yrStrt+"_"+yrLast+".txt"
write_table(filipath+filiname, "w", [/"YEAR","BOB","SCS","ISM"/], "%4s %4s %4s %4s")
write_table(filipath+filiname, "a", [/year, bob_time, scs_time, ism_time/], "%0.2i %0.2i %0.2i %0.2i")


filiname = "LBQ_onset_pentad_"+yrStrt+"_"+yrLast+"_MTG.txt"
write_table(filipath+filiname, "w", [/"YEAR","BOB","SCS","ISM"/], "%4s %4s %4s %4s")
write_table(filipath+filiname, "a", [/year, time(0,0,:), time(1,0,:), time(2,0,:)/], "%0.2i %0.2i %0.2i %0.2i")

filiname = "LBQ_onset_pentad_"+yrStrt+"_"+yrLast+"_prec.txt"
write_table(filipath+filiname, "w", [/"YEAR","BOB","SCS","ISM"/], "%4s %4s %4s %4s")
write_table(filipath+filiname, "a", [/year, time(0,1,:), time(1,1,:), time(2,1,:)/], "%0.2i %0.2i %0.2i %0.2i")

filiname = "LBQ_onset_pentad_"+yrStrt+"_"+yrLast+"_wind.txt"
write_table(filipath+filiname, "w", [/"YEAR","BOB","SCS","ISM"/], "%4s %4s %4s %4s")
write_table(filipath+filiname, "a", [/year, time(0,2,:), time(1,2,:), time(2,2,:)/], "%0.2i %0.2i %0.2i %0.2i")

end