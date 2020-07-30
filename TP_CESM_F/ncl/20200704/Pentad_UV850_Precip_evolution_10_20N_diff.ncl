
begin

nlat = 96
nlon = 144

lonstr = 60
lonend = 150
latstr = 0 ;-10
latend = 50

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)

LEV = 850

;-------------------------------------------------------------------
;-------------------------------------------------------------------
latS = 10
latN = 20
yy = ispan(6, 30, 1)

U = new((/nyear,365,nlon/), float)
V = U
PREC = U
pentad_u = new((/ncase,73,nlon/), float)
pentad_v = pentad_u
pentad_prec = pentad_u
do i = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
  
  do j = 0, dimsizes(yy)-1
    yr = sprinti("%0.2i", yy(j))
    cname = "F2000_"+casename(i)+".plevel.U."+yr+".daily.nc"
    fc = addfile(diri+cname, "r")
    UX = fc->U(:,{LEV},:,:)
    U(j,:,:) = dim_avg_n_Wrap(UX(:,{latS:latN},:), 1)
    dname = "F2000_"+casename(i)+".plevel.V."+yr+".daily.nc"
    fd = addfile(diri+dname, "r")
    VX = fd->V(:,{LEV},:,:)
    V(j,:,:) = dim_avg_n_Wrap(VX(:,{latS:latN},:), 1)
  end do


    aname = "F2000_"+casename(i)+".PRECC.daily.06-30.nc"
    fa = addfile(diri+aname, "r")
    AX = fa->PRECC
    bname = "F2000_"+casename(i)+".PRECL.daily.06-30.nc"
    fb = addfile(diri+bname, "r")
    BX = fb->PRECL
    PRECX = AX
    PRECX = AX + BX
    PRECX = PRECX * 1000. * 86400
    PRECX@units = "mm/day"
    PREC = dim_avg_n_Wrap(PRECX(:,:,{latS:latN},:), 2)


  print("end of reading")
  u = dim_avg_n_Wrap(U, 0) ; 365,96,144
  printVarSummary(u)
  pentad_ux = u(0:72,:)

  v = dim_avg_n_Wrap(V, 0) ; 365,96,144
  printVarSummary(v)
  pentad_vx = v(0:72,:)

  prec = dim_avg_n_Wrap(PREC, 0) ; 365,96,144
  printVarSummary(prec)
  pentad_precx = prec(0:72,:)

  do p = 0,72  ; a year has 73 penta
    pentad_ux(p,:) = dim_sum_n_Wrap(u(5*p:5*p+4,:), 0)/5.0
    pentad_vx(p,:) = dim_sum_n_Wrap(v(5*p:5*p+4,:), 0)/5.0
    pentad_precx(p,:) = dim_sum_n_Wrap(prec(5*p:5*p+4,:), 0)/5.0
  end do
  pentad_u(i,:,:) = pentad_ux
  pentad_v(i,:,:) = pentad_vx  
  pentad_prec(i,:,:) = pentad_precx
end do
pentad_u!0 = "case"
pentad_u&case= casename
pentad_u!1 = "time"
pentad_u&time = ispan(1, 73, 1)
copy_VarCoords(pentad_u, pentad_v)
copy_VarCoords(pentad_u, pentad_prec)
printVarSummary(pentad_u)

u_diff = pentad_u(1:ncase-1,:,:)
v_diff = pentad_v(1:ncase-1,:,:)
prec_diff = pentad_prec(1:ncase-1,:,:)
do i = 1, ncase-1
  u_diff(i-1,:,:) = pentad_u(i,:,:) - pentad_u(0,:,:)
  v_diff(i-1,:,:) = pentad_v(i,:,:) - pentad_v(0,:,:)
  prec_diff(i-1,:,:) = pentad_prec(i,:,:) - pentad_prec(0,:,:)
end do
;----------------------------------------------------------------
tpmaskpath = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_F/model_need_input_data/"
tpmaskname = "TPmask_1000m_C20190410.nc"
ftp = addfile(tpmaskpath+tpmaskname, "r")

tpmask = ftp->tp
tpmask = smth9_Wrap(tpmask, 0.5, -0.25, True)
printVarSummary(tpmask)
;-------------------------------------------------------------------
;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
title = casename ;(/"CTRL","NOHEAT"/) ;,"NOTP"

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/"
figname = "Pentad_UV"+LEV+"_Precip_evolution_10_20N_diff"

wks = gsn_open_wks("pdf", figpath+figname)

xaxis = ispan(1, 73, 2)
pens = 13
pene = 48

nplot = ncase-1
plot = new(nplot, graphic)
plotd = plot


res          = True                ; plot mods desired
res@gsnDraw  = False
res@gsnFrame = False
res@vpHeightF = 0.6
res@vpWidthF  = 0.8


res@gsnScalarContour = True
res@cnFillOn         = True                ; color on
res@cnFillPalette    = "precip_diff_12lev"    ; set color map
res@lbLabelBarOn = False
res@cnInfoLabelOn = False
res@cnLinesOn  = False
res@cnLineLabelsOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -5
res@cnMaxLevelValF = 5
res@cnLevelSpacingF  = 0.5                  ; contour spacing

; res@trYMinF = 1
; res@trYMaxF = 73
; res@tmYLMode = "Explicit"
; res@tmYLLabels = yaxis
; res@tmYLValues = yaxis
; res@tmYLLabelStride = 3
; res@tmYLMinorOn = True

res@tmXBMinorOn = False
res@tmXBMode = "Explicit"
res@tmXBValues = ispan(lonstr, lonend, 10)
res@tmXBLabels = (/"60E","70E","80E","90E","100E","110E","120E","130E","140E","150E"/)  ;(/"60E","80E","100E","120E","140E"/)
res@tmXBLabelFontHeightF = 0.02


res@tmYLPrecision = 2
res@tiYAxisString    = "Pentad"              ; y axis title
res@gsnStringFontHeightF = 0.03

vcres = True
vcres@gsnDraw  = False
vcres@gsnFrame = False
vcres@vcRefLengthF            = 0.04            ;the length of the reference magnitude
vcres@vcRefMagnitudeF         = 2.0
vcres@vcMinDistanceF          = 0.018            ;the minimum distance of neighboring vectors
vcres@vcGlyphStyle            = "CurlyVector"   ;the style of vector
vcres@vcLineArrowThicknessF   = 3.0

vcres@vcRefAnnoOrthogonalPosF = -0.4          ;move ref vector down

vcres@vcRefAnnoString2On = False
vcres@vcRefAnnoFontHeightF = 0.02
vcres@vcLineArrowColor = "black" ;grey50
; res@vcVectorDrawOrder = "PostDraw"

do i = 0, ncase-2
    res@gsnLeftString = title(i+1)+" - Ctrl"
    ; plot(i) = gsn_csm_vector_scalar(wks, pentad_u(i,pens-1:pene-1,{lonstr:lonend}), pentad_v(i,pens-1:pene-1,{lonstr:lonend}), pentad_prec(i,pens-1:pene-1,{lonstr:lonend}), res)
    plot(i) = gsn_csm_contour(wks,  prec_diff(i,pens-1:pene-1,{lonstr:lonend}), res)
    plotd(i) = gsn_csm_vector(wks, u_diff(i,pens-1:pene-1,{lonstr:lonend}), v_diff(i,pens-1:pene-1,{lonstr:lonend}), vcres)
    overlay(plot(i), plotd(i))
end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "UV"+LEV+" & Precip ("+latS+"-"+latN+"N): Diff with Ctrl"

gsn_panel(wks, plot, (/2,2/), pres)


end