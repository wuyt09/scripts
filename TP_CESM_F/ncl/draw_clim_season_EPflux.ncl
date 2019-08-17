
begin

LEV = 200

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
diri2 = "/home/yangsong3/data-model/wuyt/TP_NOQ1_F2000_F19G16_CAM4/pro/"
diri3 = "/home/yangsong3/data-model/wuyt/TP_NOTP_NOQ1_F2000_F19G16_CAM4/pro/"

filiu1 = "F2000_CTRL.plevel.U.mon.0601-3012.nc"
filiu2 = "F2000_NOQ1.plevel.U.mon.0601-3012.nc"
filiu3 = "F2000_NOTP_NOQ1.plevel.U.mon.0601-3012.nc"

filiv1 = "F2000_CTRL.plevel.V.mon.0601-3012.nc"
filiv2 = "F2000_NOQ1.plevel.V.mon.0601-3012.nc"
filiv3 = "F2000_NOTP_NOQ1.plevel.V.mon.0601-3012.nc"

filic1 = "F2000_CTRL.plevel.T.mon.0601-3012.nc"
filic2 = "F2000_NOQ1.plevel.T.mon.0601-3012.nc"
filic3 = "F2000_NOTP_NOQ1.plevel.T.mon.0601-3012.nc"

fu1 = addfile(diri1+filiu1, "r")
fu2 = addfile(diri2+filiu2, "r")
fu3 = addfile(diri3+filiu3, "r")
fv1 = addfile(diri1+filiv1, "r")
fv2 = addfile(diri2+filiv2, "r")
fv3 = addfile(diri3+filiv3, "r")
u1 = fu1->U;(:,{LEV},:,:)
u2 = fu2->U;(:,{LEV},:,:)
u3 = fu3->U;(:,{LEV},:,:)
v1 = fv1->V;(:,{LEV},:,:)
v2 = fv2->V;(:,{LEV},:,:)
v3 = fv3->V;(:,{LEV},:,:)
printVarSummary(u1)

fc1 = addfile(diri1+filic1, "r")
fc2 = addfile(diri2+filic2, "r")
fc3 = addfile(diri3+filic3, "r")
t1 = fc1->T;(:,{LEV},:,:)
t2 = fc2->T;(:,{LEV},:,:)
t3 = fc3->T;(:,{LEV},:,:)
printVarSummary(t1)

lat = fu1->lat
lon = fu1->lon
lev = fu1->lev
;=========================================
plvl = lev
ep1 = epflux(u1, v1, t1, plvl, lat, False)
Fphi1 = ep1[0]
Fp1 = ep1[1]
EPdiv1 = ep1[2]
dudt1 = ep1[3]
printVarSummary(EPdiv1)

end