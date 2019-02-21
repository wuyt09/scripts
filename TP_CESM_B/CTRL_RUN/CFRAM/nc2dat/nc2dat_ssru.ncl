
begin

yrStrt = 200
yrLast = 349
year = ispan(yrStrt, yrLast, 1)

yrbs = 200
yrbe = 249
yrws = 300
yrwe = 349

basename = "ssru_base.dat"
warmname = "ssru_warm.dat"


diri = "/home/yangsong3/data-model/wuyt/TP_CTRL_B2000_F19G16_CAM4/pro/"
fil1 = "B2000_f19g16_CP_CTRL.cam.h0.FSDS.200-349.nc"
fil2 = "B2000_f19g16_CP_CTRL.cam.h0.FSNS.200-349.nc"
f1 = addfile(diri+fil1, "r")
f2 = addfile(diri+fil2, "r")

C1 = f1->FSDS
C2 = f2->FSNS
CC = C1 - C2
copy_VarMeta(C1, CC)
CC@long_name = "surface upward shortwave flux"
printVarSummary(CC)

cc = month_to_annual(CC, 1)
cc&year = year

cc_base = dim_avg_n_Wrap(cc({yrbs:yrbe},:,:), 0)
cc_warm = dim_avg_n_Wrap(cc({yrws:yrwe},:,:), 0)

;=========write data===========================================
outpath = "/home/yangsong3/wuyt/sysu/data/TP_CESM/CTRL_RUN/nc2dat/"
system("rm -f "+outpath+basename)
fbindirwrite(outpath+basename, cc_base)

system("rm -f "+outpath+warmname)
fbindirwrite(outpath+warmname, cc_warm)

;==============check data=============================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM/CTRL_RUN/nc2dat/"
figname = "ssru"
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(2, graphic)

 res = True
 res@gsnDraw          = False
 res@gsnFrame         = False
 ;res@gsnSpreadColors  = True
 
 res@cnFillOn         = True
 res@cnLineLabelsOn   = False
 res@cnInfoLabelOn    = False
 res@cnFillPalette    = "CBR_wet"
 
 res@lbLabelAutoStride = True


 res@mpFillOn  = True
 res@mpOutlineOn = True

 res@gsnLeftString = "FSDS-FSNS: 200-249"
 plot(0) = gsn_csm_contour_map(wks, cc_base, res)

 res@gsnLeftString = "FSDS-FSNS: 300-349"
 plot(1) = gsn_csm_contour_map(wks, cc_warm, res)

 pres = True
 pres@gsnMaximize = True

 gsn_panel(wks, plot, (/2,1/), pres)
end
