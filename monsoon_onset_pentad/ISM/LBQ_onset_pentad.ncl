
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

onset_time = (/33,31,30,30,33,30,30,31,31,29,31,27,32,33,30,30,33,31,33,32,28,30,32,32,32,27,32,29,31,30,31,29,30,32,30,32,32,31,31/)
print(year+"    "+onset_time)

filipath = "./"
filiname = "LBQ_onset_pentad_"+yrStrt+"_"+yrLast+".txt"
write_table(filipath+filiname, "w", [/"YEAR","Pentad"/], "%4s %4s")
write_table(filipath+filiname, "a", [/year, onset_time/], "%0.2i %0.2i")

end