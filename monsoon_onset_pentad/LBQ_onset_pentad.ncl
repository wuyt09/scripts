
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

bob_time = (/27,26,26,23,27,22,22,26,27,25,24,25,24,26,27,25,26,23,28,27,20,22,24,26,27,26,26,23,25,25,21,28,24,23,27,26,27,28,27/)
scs_time = (/29,28,30,30,29,28,30,27,32,29,32,29,32,28,29,30,28,34,28,29,28,27,28,28,28,34,29,30,28,29,29,29,27,27,27,31,33,32,29/)
ism_time = (/33,31,30,30,33,30,30,31,31,29,31,27,32,33,30,30,33,31,33,32,28,30,32,32,32,27,32,29,31,30,31,29,30,32,30,32,32,31,31/)
print(year+"    "+scs_time)

filipath = "./"
filiname = "LBQ_onset_pentad_"+yrStrt+"_"+yrLast+".txt"
write_table(filipath+filiname, "w", [/"YEAR","BOB","SCS","ISM"/], "%4s %4s %4s %4s")
write_table(filipath+filiname, "a", [/year, bob_time, scs_time, ism_time/], "%0.2i %0.2i %0.2i %0.2i")

end