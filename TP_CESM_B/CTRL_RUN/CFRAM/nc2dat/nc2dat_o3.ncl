
begin

yrStrt = 200
yrLast = 349
year = ispan(yrStrt, yrLast, 1)

yrbs = 200
yrbe = 249
yrws = 300
yrwe = 349

basename = "o3_base.dat"
warmname = "o3_warm.dat"

;===========================================================
diri = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"

filia = "B2000_f19g16_CP_CTRL.cam.h0.hyam.200-349.nc"
filib = "B2000_f19g16_CP_CTRL.cam.h0.hybm.200-349.nc"
fa = addfile(diri+filia, "r")
fb = addfile(diri+filib, "r")
hyam = fa->hyam
hybm = fb->hybm
printVarSummary(hyam)

filip = "B2000_f19g16_CP_CTRL.cam.h0.PS.200-349.nc"
fp = addfile(diri+filip, "r")
ps = fp->PS
printVarSummary(ps)
ps0 = clmMonTLL(ps)
printVarSummary(ps0)

fili = "ozone_1.9x2.5_L26_2000clim_c091112.nc"
f = addfile(diri+fili, "r")
Cc = f->O3
printVarSummary(Cc)
printVarSummary(Cc&lev)

p0 = 1000.
plevo = (/1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150,\
         175, 200, 225,250, 300, 350, 400, 450, 500, 550, 600,\
         650, 700, 750, 775, 800, 825,850, 875, 900, 925, 950, 975, 1000/)
CC = vinth2p(Cc, hyam, hybm, plevo, ps0, 1, p0, 1, False)
copy_VarMeta(Cc(:,0,:,:), CC(:,0,:,:))
CC!1 = "lev"
CC&lev = plevo
CC&lev@units = "hPa"
printVarSummary(CC)

cc = dim_avg_n_Wrap(CC, 0)
cc_base = cc
cc_warm = cc

;=========write data===========================================
outpath = "/home/yangsong3/wuyt/sysu/data/TP_CESM/CTRL_RUN/nc2dat/"
system("rm -f "+outpath+basename)
fbindirwrite(outpath+basename, cc_base)

system("rm -f "+outpath+warmname)
fbindirwrite(outpath+warmname, cc_warm)


;============data checking====================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM/CTRL_RUN/nc2dat/"
figname = "O3"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(2, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False

res@cnFillOn = True
res@cnLineLabelsOn = False
res@cnInfoLabelOn  = False
res@cnFillPalette  = "precip2_17lev"

res@lbLabelAutoStride = True

res@tmYLMode      = "Explicit"
res@tmYLValues    = (/"1000","850","750","500","300","100","70","50"/)
res@tmYLLabels    = (/"1000","850","750","500","300","100","70","50"/)

res@gsnRightString = "lat :35N"

res@gsnLeftString = "O3: 200-249"
plot(0) = gsn_csm_pres_hgt(wks, cc_base({1000:50},{35},:), res)

res@gsnLeftString = "O3: 300-349"
plot(1) = gsn_csm_pres_hgt(wks, cc_warm({1000:50},{35},:), res)


pres = True
pres@gsnMaximize = True

gsn_panel(wks, plot, (/2,1/), pres)
end
