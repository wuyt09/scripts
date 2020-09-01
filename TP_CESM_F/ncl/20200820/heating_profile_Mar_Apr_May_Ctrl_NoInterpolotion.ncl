
begin

lonstr = 40
lonend = 150
latstr = -15
latend = 50

casename = (/"CTRL","NOTPQ1"/) ;,"NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"
ncase = dimsizes(casename)
nlev = 26
nlat = 96 ;34
nlon = 144 ;45

dtc = new((/ncase,12,nlev,nlat,nlon/), float)
dtv = dtc

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/"
    ; filis = "F2000_"+casename(i)+".plevel.QRS.mon.0601-3012.nc"
    ; filir = "F2000_"+casename(i)+".plevel.QRL.mon.0601-3012.nc"
    filic = "F2000_"+casename(i)+".DTCOND.mon.0601-3012.nc"
    filiv = "F2000_"+casename(i)+".DTV.mon.0601-3012.nc"

    ; fs = addfile(diri+filis, "r")
    ; fr = addfile(diri+filir, "r")
    fc = addfile(diri+filic, "r")
    fv = addfile(diri+filiv, "r")
    ; s1 = fs->QRS(:,:,{latstr:latend},{lonstr:lonend})
    ; r1 = fr->QRL(:,:,{latstr:latend},{lonstr:lonend})
    c1 = fc->DTCOND ;(:,:,{latstr:latend},{lonstr:lonend})
    v1 = fv->DTV ;(:,:,{latstr:latend},{lonstr:lonend})
    c1 = c1*86400
    v1 = v1*86400
    printVarSummary(c1)

    dtc(i,:,:,:,:) = clmMonTLLL(c1)
    dtv(i,:,:,:,:) = clmMonTLLL(v1)
end do
dtc!0 = "case"
dtc&case = casename
copy_VarCoords(dtc, dtv)
printVarSummary(dtc)

lev = dtc&lev
lev!0 = "level"
lev&level = lev

print(dtc(0,2,20,{25:30},{100:110}))
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
; tpmaskname = "TPmask_1000m_C20190410.nc"
tpmaskname = "TPmask_1500m_C20200713.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
printVarSummary(tpmask)
;-------------------------------------------------------------------

dtc = mask(dtc, tpmask.ge.0.5, True)
dtv = mask(dtv, tpmask.ge.0.5, True)

profile_dtc = wgt_areaave_Wrap(dtc, 1, 1, 0) 
profile_dtv = wgt_areaave_Wrap(dtv, 1, 1, 0)
printVarSummary(profile_dtv)

print(profile_dtc(0,2,:))
; exit
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200820/"
figname = "heating_profile_Mar_Apr_May_Ctrl_NoInterpolotion"

 Colors = (/"red","green4"/) ;"black","goldenrod",
 Texts = (/"condensation heating","vertical diffusion"/)  ;"total","QRL+QRS",
 leftsti = (/"March","April","May"/)

wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True
res@vpHeightF = 0.5
res@vpWidthF = 0.8

 res@gsnXRefLine = 0
 res@gsnXRefLineDashPattern = 2
 res@gsnXRefLineColor = "grey"
 res@gsnXRefLineThicknessF = 3
 res@trXMinF = -6.
 res@trXMaxF = 6
 res@trYMinF   = 100
 res@trYMaxF   = 1000
 res@trYLog = True
 
 res@tmYLMode      = "Explicit"
 res@tmYLValues    = (/"1000","850","750","500","300","100"/)
 res@tmYLLabels    = (/"1000","850","750","500","300","100"/)
 res@xyLineColors = Colors
 res@xyLineThicknessF = 5
 res@xyDashPattern = 0

 res@trYReverse = True
 res@gsnRightString = "K"
 res@gsnStringFontHeightF = 0.025



levs = 1000
levt = 100
do i = 0, 2
    res@gsnLeftString = leftsti(i)
    plot(i) = gsn_csm_xy(wks, profile_dtc(:,i+2,{levs:levt}), lev({levs:levt}), res)
end do

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/1,3/), pres)



end