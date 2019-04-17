
begin

yrStat = 260
yrLast = 304 

diri2 = "/home/yangsong3/data-model/wuyt/TP_NOSH_B2000_F19G16_CAM4/pro/"

fili2 = "TP_NOSH.T.mon.26001-30712.nc"

f2 = addfile(diri2+fili2, "r")

time2 = f2->time
TIME2 = cd_calendar(time2, -2)/100
;print(TIME1(0:40))
;print(TIME2)

iYYYY2 = ind(TIME2.ge.26002 .and. TIME2.le.30501)

x2 = f2->T(iYYYY2,:,:,:)
printVarSummary(x2)


;=========================================================================
; hyam hybm PS
;----------------------------------

filia2 = "TP_NOSH.hyam.mon.30001.nc"
filib2 = "TP_NOSH.hybm.mon.30001.nc"
filip2 = "TP_NOSH.PS.mon.26001-30712.nc"

fa2 = addfile(diri2+filia2, "r")
fb2 = addfile(diri2+filib2, "r")
fp2 = addfile(diri2+filip2, "r")

hyam2 = fa2->hyam
hybm2 = fb2->hybm
ps2   = fp2->PS(iYYYY2,:,:)
printVarSummary(ps2)


p0  = 1000.
p0@units = "hPa"

;plevo = (/1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700,\
;         650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175,\
;         150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1/)
plevo = (/1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200,\
         225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,\
         775, 800, 825, 850, 875, 900, 925, 950, 975, 1000/)
plevo@units = "hPa"
;==========================================================================

x2p = vinth2p(x2, hyam2, hybm2, plevo, ps2, 1, p0, 1, False)
x2p!1 = "lev"
x2p&lev = plevo
printVarSummary(x2p)

outpath2 = diri2
outname2 = "TP_NOSH.plevel.T.mon.260-304.nc"

system("rm "+outpath2+outname2)
fout2 = addfile(outpath2+outname2, "c")
fout2->T=x2p

end