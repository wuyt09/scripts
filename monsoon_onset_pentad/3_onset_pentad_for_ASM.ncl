
begin

bob_lat1 = 5
bob_lat2 = 15
bob_lon1 = 85
bob_lon2 = 100

scs_lat1 = 5
scs_lat2 = 20
scs_lon1 = 110
scs_lon2 = 120

ism_lat1 = 10
ism_lat2 = 20
ism_lon1 = 60
ism_lon2 = 75


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

diri = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
filiu = "2_U850_onset_pentad.nc"
filim = "2_MTG_onset_pentad.nc"
filip = "2_precip_onset_pentad_GPCP.nc"

fu = addfile(diri+filiu, "r")
fm = addfile(diri+filim, "r")
fp = addfile(diri+filip, "r")

u = fu->u850_peak
mtg = fm->MTG_peak
precx = fp->prec_peak
printVarSummary(u)
printVarSummary(mtg)
prec = linint2_Wrap(precx&lon, precx&lat, precx, True, u&longitude, u&latitude, 0)
printVarSummary(prec)

x = new((/3,nyear,nlat,nlon/), float)
x(0,:,:,:) = mtg
x(1,:,:,:) = u
x(2,:,:,:) = prec

onset_pentad = dim_max_n_Wrap(x, 0)
printVarSummary(onset_pentad)
onset_pentad_clm = dim_avg_n_Wrap(onset_pentad, 0)

;========================================================
; individual monsoon component
;========================================================

bob = wgt_areaave_Wrap(onset_pentad(:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
scs = wgt_areaave_Wrap(onset_pentad(:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
ism = wgt_areaave_Wrap(onset_pentad(:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)
print(bob)

; write_table(filipath+filiname, "w", [/"Year","BOB","SCS","ISM"/], "%10s %10s %10s %10s")
; write_table(filipath+filiname, "a", [/year,bob,scs,ism/], "%2.4i %10.2f %10.2f %10.2f %10.2f")

end