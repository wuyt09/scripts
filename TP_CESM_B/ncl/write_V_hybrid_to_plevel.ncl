
begin

yrStat = 260
yrLast = 304 

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_B2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_heat_B2000_F19G16_CAM4/pro/"

fili1 = "B2000_f19g16_CP_CTRL.cam.h0.V.200-349.nc"
fili2 = "TP_NOTOPO.V.mon.26001-30508.nc"
fili3 = "TP_HEAT.V.mon.26001-34912.nc"

f1 = addfile(diri1+fili1, "r")
f2 = addfile(diri2+fili2, "r")
f3 = addfile(diri3+fili3, "r")

time1 = f1->time
time2 = f2->time
time3 = f3->time
TIME1 = cd_calendar(time1, -2)/100
TIME2 = cd_calendar(time2, -2)/100
TIME3 = cd_calendar(time3, -2)/100
;print(TIME1(0:40))
;print(TIME2)

iYYYY1 = ind(TIME1.ge.26002 .and. TIME1.le.30501)
iYYYY2 = ind(TIME2.ge.26002 .and. TIME2.le.30501)
iYYYY3 = ind(TIME3.ge.26002 .and. TIME3.le.30501)

x1 = f1->V(iYYYY1,:,:,:) ;
x2 = f2->V(iYYYY2,:,:,:)
x3 = f3->V(iYYYY3,:,:,:)
printVarSummary(x1)
printVarSummary(x2)
printVarSummary(x3)

;=========================================================================
; hyam hybm PS
;----------------------------------

filia1 = "B2000_f19g16_CP_CTRL.cam.h0.hyam.30001.nc"
filib1 = "B2000_f19g16_CP_CTRL.cam.h0.hybm.30001.nc"
filip1 = "B2000_f19g16_CP_CTRL.cam.h0.PS.200-349.nc"

filia2 = "TP_NOTOPO.hyam.mon.30001.nc"
filib2 = "TP_NOTOPO.hybm.mon.30001.nc"
filip2 = "TP_NOTOPO.PS.mon.26001-30508.nc"

filia3 = "TP_HEAT.hyam.mon.30001.nc"
filib3 = "TP_HEAT.hybm.mon.30001.nc"
filip3 = "TP_HEAT.PS.mon.26001-34912.nc"

fa1 = addfile(diri1+filia1, "r")
fb1 = addfile(diri1+filib1, "r")
fp1 = addfile(diri1+filip1, "r")

fa2 = addfile(diri2+filia2, "r")
fb2 = addfile(diri2+filib2, "r")
fp2 = addfile(diri2+filip2, "r")

fa3 = addfile(diri3+filia3, "r")
fb3 = addfile(diri3+filib3, "r")
fp3 = addfile(diri3+filip3, "r")

hyam1 = fa1->hyam
hybm1 = fb1->hybm
ps1   = fp1->PS(iYYYY1,:,:)

hyam2 = fa2->hyam
hybm2 = fb2->hybm
ps2   = fp2->PS(iYYYY2,:,:)

hyam3 = fa3->hyam
hybm3 = fb3->hybm
ps3   = fp3->PS(iYYYY3,:,:)
printVarSummary(ps1)
printVarSummary(ps2)
printVarSummary(ps3)

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

x1p = vinth2p(x1, hyam1, hybm1, plevo, ps1, 1, p0, 1, False)
x2p = vinth2p(x2, hyam2, hybm2, plevo, ps2, 1, p0, 1, False)
x3p = vinth2p(x3, hyam3, hybm3, plevo, ps3, 1, p0, 1, False)
x1p!1 = "lev"
x1p&lev = plevo
copy_VarCoords(x1p, x2p)
copy_VarCoords(x1p, x3p)
printVarSummary(x1p)
printVarSummary(x2p)
printVarSummary(x3p)

outpath1 = diri1
outpath2 = diri2
outpath3 = diri3

outname1 = "B2000_f19g16_CP_CTRL.plevel.V.mon.260-304.nc"
outname2 = "TP_NOTOPO.plevel.V.mon.260-304.nc"
outname3 = "TP_HEAT.plevel.V.mon.260-304.nc"

system("rm "+outpath1+outname1)
fout1 = addfile(outpath1+outname1, "c")
fout1->V=x1p

system("rm "+outpath2+outname2)
fout2 = addfile(outpath2+outname2, "c")
fout2->V=x2p

system("rm "+outpath3+outname3)
fout3 = addfile(outpath3+outname3, "c")
fout3->V=x3p

end