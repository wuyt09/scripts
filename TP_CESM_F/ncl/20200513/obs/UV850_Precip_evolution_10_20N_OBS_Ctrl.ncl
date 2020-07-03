
begin

nlat = 96
nlon = 144

lonstr = 60
lonend = 150
latstr = 0 ;-10
latend = 50

casename = (/"CTRL"/) ;"NOTP","NOQ1"
ncase = dimsizes(casename)
year = ispan(6, 30, 1)
nyear = dimsizes(year)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
nmonth = dimsizes(month)

LEV = 850
penA = 19
penB = 60

;-------------------------------------------------------------------
; Ctrl data
;-------------------------------------------------------------------
lonW = 110
lonE = 120
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

; ud = pentad_u(case|:,lat|:,time|:)
; vd = pentad_v(case|:,lat|:,time|:)
; precd = pentad_prec(case|:,lat|:,time|:)
; printVarSummary(ud)
; ==============================================================
; OBS 
; ==============================================================
    yrStrt = 1980 ;1979
    yrLast = 1989 ;2015
    year := ispan(yrStrt, yrLast, 1)
    nyr  = dimsizes(year)
    npnt = 48 ;pentad 19 to pentad 60

   ; filu  = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/pressure/uwnd/uwnd.interim.19*.nc /home/yangsong3/data-observation/ERA-interim/pressure/uwnd/uwnd.interim.200*.nc  /home/yangsong3/data-observation/ERA-interim/pressure/uwnd/uwnd.interim.201[0-7].nc")
   ; filv  = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/pressure/vwnd/vwnd.interim.19*.nc /home/yangsong3/data-observation/ERA-interim/pressure/vwnd/vwnd.interim.200*.nc  /home/yangsong3/data-observation/ERA-interim/pressure/vwnd/vwnd.interim.201[0-7].nc")
   filu  = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/pressure/uwnd/uwnd.interim.198*.nc")
   filv  = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/pressure/vwnd/vwnd.interim.198*.nc")
   print(filu)

   fu  = addfiles(filu,"r")
   fv  = addfiles(filv,"r")
   u0 = short2flt(fu[:]->uwnd(:,{LEV},:,:))
   v0 = short2flt(fv[:]->vwnd(:,{LEV},:,:))
   printVarSummary(u0)

   yr   := new(nyr,integer,-999)
   mn   = new((/nyr,12/),integer,-999)
   dy   = new((/nyr,31/),integer,-999)

   ; yr   = ispan( 1979,2016,1 )
   yr = year
   do n=0,nyr-1
     mn(n,:)   = ispan(1,12,1)
     dy(n,:)   = ispan(1,31,1)
   end do
   ddd4   = day_of_year(yr,mn(:,3),dy(:,0))
   ddd6   = day_of_year(yr,mn(:,10),dy(:,25))
   ddd12  = day_of_year(yr,mn(:,11),dy(:,30))
   ; print(ddd4)
   ; print(ddd6)
   ; print(ddd12)
   ; exit
;************************************************

   u_TLL  = u0(:,::-1,:);13149
   v_TLL  = v0(:,::-1,:)
   printVarSummary(u_TLL)
   nlat = dimsizes(u_TLL(1,:,1))
   nlon = dimsizes(u_TLL(1,1,:))

   u_AMJJASO = new((/nyr,npnt*5,nlat,nlon/),float)
   v_AMJJASO = new((/nyr,npnt*5,nlat,nlon/),float)
;   t_Pnt     = t_TLLL(0:1079,:,:,:)

   days = 0
   do n=0,nyr-1
    ddd6(n) = ddd6(n)-1+days
    ddd4(n) = ddd4(n)-1+days
    u_AMJJASO(n,:,:,:) = u_TLL(ddd4(n):ddd6(n),:,:)
    v_AMJJASO(n,:,:,:) = v_TLL(ddd4(n):ddd6(n),:,:)
    days    = days+ddd12(n)
   end do
   printVarSummary(u_AMJJASO)

   u_Pnt     = u_AMJJASO(:,0:npnt-1,:,:)
   u_Pnt     = 0
   v_Pnt     = v_AMJJASO(:,0:npnt-1,:,:)
   v_Pnt     = 0
   do n=0,nyr-1
    do k=0,npnt-1
      do i=0,4
        u_Pnt(n,k,:,:)   = u_Pnt(n,k,:,:) + u_AMJJASO(n,k*5+i,:,:)/5.0
        v_Pnt(n,k,:,:)   = v_Pnt(n,k,:,:) + v_AMJJASO(n,k*5+i,:,:)/5.0
      end do
    end do
   end do
u_Pnt&time = ispan(19, 66, 1)
v_Pnt&time = ispan(19, 66, 1)
uclim_pentad = dim_avg_n_Wrap(u_Pnt, 0)
vclim_pentad = dim_avg_n_Wrap(v_Pnt, 0)
printVarSummary(uclim_pentad)
uclim = dim_avg_n_Wrap(uclim_pentad(:,{latS:latN},:), 1)
vclim = dim_avg_n_Wrap(vclim_pentad(:,{latS:latN},:), 1)
; uclim = uxclim(lat|:,time|:)
; vclim = vxclim(lat|:,time|:)
printVarSummary(uclim)
printMinMax(uclim({10:50},0:26), False)
printMinMax(vclim({10:50},0:26), False)
;---------------------------------------------
dirip = "/home/yangsong3/data-observation/GPCP/"
filip = "precip.pentad.clim.mean.nc"
fp = addfile(dirip+filip, "r")
pclim_pentad = fp->precip(:,::-1,:)
pclim_pentad&pentad = ispan(1, 73, 1)
pclim_pentad@units = "mm/day"
printVarSummary(pclim_pentad)
pclim = dim_avg_n_Wrap(pclim_pentad(:,{latS:latN},:), 1)
; pclim = pxclim(lat|:,pentad|:)
;-------------------------------------------------------------------
;-------------------------------------------------------------------
Month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; title = (/"CTRL","NOHEAT"/) ;,"NOTP"
title = (/"OBS","CTRL"/)

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/obs/"
figname = "UV"+LEV+"_Precip_evolution_10_20N_OBS_Ctrl"

wks = gsn_open_wks("pdf", figpath+figname)

xaxis = ispan(1, 73, 2)
pens = 19
pene = 65

nplot = 2
plot = new(nplot, graphic)
plotd = plot


res          = True                ; plot mods desired
res@gsnDraw  = False
res@gsnFrame = False
res@vpHeightF = 0.6
res@vpWidthF  = 0.8

res@gsnScalarContour = True
res@cnFillOn         = True                ; color on
res@cnFillPalette    = "precip2_17lev"    ; set color map
res@lbLabelBarOn = False
res@cnInfoLabelOn = False
res@cnLinesOn  = False
res@cnLineLabelsOn = False

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 2
res@cnMaxLevelValF = 20
res@cnLevelSpacingF  = 2                  ; contour spacing
; res@cnFillDrawOrder = "PreDraw"
; res@trYMinF = 1
; res@trYMaxF = 73
; res@tmYLMode = "Explicit"
; res@tmYLLabels = yaxis
; res@tmYLValues = yaxis
; res@tmYLLabelStride = 3
; res@tmYLMinorOn = True

res@tmXBMinorOn = False
res@tmXBMode = "Explicit"
; res@tmXBValues =(/50,60,70,80,90,100,110,120,130,140/)  ;(/60,80,100,120,140/)
res@tmXBValues = ispan(lonstr, lonend, 10)
res@tmXBLabels = (/"60E","70E","80E","90E","100E","110E","120E","130E","140E","150E"/)  ;(/"60E","80E","100E","120E","140E"/)
res@tmXBLabelFontHeightF = 0.02


res@tmXBPrecision = 2
res@tiYAxisString    = "Pentad"              ; y axis title
res@gsnStringFontHeightF = 0.03

; res@lbLabelStride    = 2                   ; every other label
; res@cnLinesOn        = False               ; turn off contour lines

; vcres = True
; vcres@gsnDraw  = False
; vcres@gsnFrame = False

res@vcRefLengthF            = 0.04            ;the length of the reference magnitude
res@vcRefMagnitudeF         = 5.0
res@vcMinDistanceF          = 0.018            ;the minimum distance of neighboring vectors
res@vcGlyphStyle            = "CurlyVector"   ;the style of vector
res@vcLineArrowThicknessF   = 3.0

res@vcRefAnnoOrthogonalPosF = -0.4          ;move ref vector down

res@vcRefAnnoString2On = False
res@vcRefAnnoFontHeightF = 0.02
res@vcLineArrowColor = "black" ;grey50
; res@vcVectorDrawOrder = "PostDraw"

res@gsnLeftString = title(0)
; plot(0) = gsn_csm_contour(wks, pclim({10:50},pens-1:pene-1), res)
; ; plotd(0) = gsn_csm_vector(wks, uclim({10:50},pens-19:pene-19), vclim({10:50},pens-19:pene-19), vcres)
; plotd(0) = gsn_csm_vector(wks, uclim({10:50},0:26), vclim({10:50},0:26), vcres)
; overlay(plot(0), plotd(0))
plot(0) = gsn_csm_vector_scalar(wks, uclim(pens-19:pene-19,{lonstr:lonend}), vclim(pens-19:pene-19,{lonstr:lonend}), pclim(pens-1:pene-1,{lonstr:lonend}), res)

res@gsnLeftString = title(1)
; plot(1) = gsn_csm_contour(wks, precd(0,{10:50},pens-1:pene-1), res)
; plotd(1) = gsn_csm_vector(wks, ud(0,{10:50},pens-1:pene-1), vd(0,{10:50},pens-1:pene-1), vcres)
; overlay(plot(1), plotd(1))
plot(1) = gsn_csm_vector_scalar(wks, pentad_u(0,pens-1:pene-1,{lonstr:lonend}), pentad_v(0,pens-1:pene-1,{lonstr:lonend}), pentad_prec(0,pens-1:pene-1,{lonstr:lonend}), res)


; res@cnFillPalette = "BlueWhiteOrangeRed"
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnMinLevelValF = -20
; res@cnMaxLevelValF = 20
; res@cnLevelSpacingF  = 4  
; res@gsnLeftString = title(2)
; plot(2) = gsn_csm_hov(wks, x(2,18:35,{90:160}), res)


  ; delete(res@gsnLeftString)
  ; delete(res@gsnRightString)
  ; res@cnLinesOn = True
  ; res@cnFillOn  = False
  ; res@cnLineLabelsOn = False
  ; res@cnLevelSelectionMode = "ExplicitLevels"
  ; res@cnLevels = (/0/)
  ; res@cnMonoLineThickness = True
  ; res@cnLineThicknessF = 5
  ; res@gsnLeftString = ""
  ; res@gsnRightString = ""
  ; ;res@gsnContourZeroLineThicknessF = 5
  ; do i = 0, ncase-1
  ;   plots(i) = gsn_csm_hov(wks, prec_tlat(i,18:35,{50:140}), res)
  ;   overlay(plot(i), plots(i))
  ; end do

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "UV"+LEV+" & Precip ("+latS+"-"+latN+"N)"

gsn_panel(wks, plot, (/1,2/), pres)


end