
begin

lonstr = 50
lonend = 150
latstr = 0
latend = 30

  latS1  = 10
  latN1  = 20
  lonL1  = 60
  lonR1  = 85.0

  latS2  = 5
  latN2  = 15
  lonL2  = 90
  lonR2  = 100.0

  latS3  = 10
  latN3  = 20
  lonL3  = 110
  lonR3  = 120.0

casename = (/"CTRL","NOTPQ1","NOTPQ1_MAR","NOTPQ1_APR","NOTPQ1_MAY"/) ;
ncase = dimsizes(casename)
pens = 13 ;19
pene = 42 ;36
npen = pene - pens + 1
nlat = 96
nlon = 144
nlev = 37

yrStrt = 6
yrLast = 30
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

clim_pentad_ps = new((/nyear,npen,nlat,nlon/), float)
clim_pentad_t = new((/nyear,npen,nlev,nlat,nlon/), float)

do ll = 0, ncase-1
  diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(ll)+"_F2000_F19G16_CAM4/pro/daily/"

  pname = "F2000_"+casename(ll)+".PS.daily.06-30.nc"
  fp = addfile(diri+pname, "r")

  ps = fp->PS
  printVarSummary(ps)

  pentad_ps = ps(:,0:72,:,:)

  do p = 0,72  ; a year has 73 penta
    pentad_ps(:,p,:,:) = dim_sum_n_Wrap(ps(:,5*p:5*p+4,:,:), 1)/5.0
  end do

  ; clim_pentad_PS = dim_avg_n_Wrap(pentad_ps, 0)
  clim_pentad_ps = pentad_ps(:,pens-1:pene-1,:,:)
  delete([/ps,pentad_ps/])

  ;--------------------------------------------------------------------
  tname = systemfunc("ls "+diri+"F2000_"+casename(ll)+".plevel.T.*.daily.nc")
  print(tname)

  ft = addfiles(tname, "r")

  T = ft[:]->T(:,::-1,:,:)
  t = reshape(T, (/25,365,37,nlat,nlon/))
  copy_VarCoords(T(0:364,:,:,:), t(0,:,:,:,:))

  pentad_t = t(:,0:72,:,:,:)

  do p = 0,72  ; a year has 73 penta
    pentad_t(:,p,:,:,:) = dim_sum_n_Wrap(t(:,5*p:5*p+4,:,:,:), 1)/5.0
  end do

  ; clim_pentad_T = dim_avg_n_Wrap(pentad_t, 0)
  clim_pentad_t = pentad_t(:,pens-1:pene-1,:,:,:)
  delete([/T,t,pentad_t/])

printVarSummary(clim_pentad_ps)
printVarSummary(clim_pentad_t)
;-------------------------------------------------------------------
; partial T / partial y
rearth = 6371229.0 ; meters

lon = fp->lon
dlon = (lon(2) - lon(1))*0.0174533
dY = tofloat(rearth*dlon) 
;dX = cos(lat*0.0174533)*dlon
printVarSummary(dY)

MTG = clim_pentad_t
MTG = center_finite_diff_n(clim_pentad_t, dY, True, 0, 3)
MTG!0 = "year"
MTG&year = year
printVarSummary(MTG)

; vertical integration

lev = MTG&lev
plev = lev*100.
plev@units = "Pa"
plev!0 = "lev"
plev&lev = lev
;psfc = max(plev)
;ptop = min(plev)
;dP = dpres_plevel_Wrap(plev, psfc, ptop, 0)
;dp = conform_dims(dimsizes(MTG), dP, 1)

;tdp = MTG
;tdp = MTG*dp
;printVarSummary(tdp)

;mtg = dim_sum_n_Wrap(tdp(:,{500:200},:,:), 1)/sum(dP({500:200}))
;copy_VarCoords(MTG(:,0,:,:), mtg)

Cmp_MTG = MTG(:,:,0,:,:)
pbot = 50000
ptop = 20000
do k = 0,nyear-1
  do p = 0, npen-1
    xxx  = MTG(year|k,time|p,lat|:,lon|:,lev|:)
    Cmp_MTG(k,p,:,:) = vibeta(plev, xxx, 1, clim_pentad_ps(k,p,:,:), pbot, ptop)
  end do
end do

printVarSummary(Cmp_MTG)


;-------------------------------------------------------------------
; Reversal Pentad

  bob = wgt_areaave_Wrap(Cmp_MTG(:,:,{latS2:latN2},{lonL2:lonR2}),1.0,1.0,1)
  scs = wgt_areaave_Wrap(Cmp_MTG(:,:,{latS3:latN3},{lonL3:lonR3}),1.0,1.0,1)
  ism = wgt_areaave_Wrap(Cmp_MTG(:,:,{latS1:latN1},{lonL1:lonR1}),1.0,1.0,1)
  printVarSummary(bob)

  diff  = Cmp_MTG
  last  = Cmp_MTG
  diff1 = bob
  last1 = bob
  diff2 = bob
  last2 = bob
  diff3 = bob
  last3 = bob
  do n=1,npen-1
    diff(:,n,:,:) = Cmp_MTG(:,n,:,:)*Cmp_MTG(:,n-1,:,:)
    diff1(:,n) = bob(:,n)*bob(:,n-1)
    diff2(:,n) = scs(:,n)*scs(:,n-1)
    diff3(:,n) = ism(:,n)*ism(:,n-1)
    end do
  do n=0,npen-3
    last(:,n,:,:) = Cmp_MTG(:,n,:,:)*Cmp_MTG(:,n+1,:,:)
    last1(:,n) = bob(:,n)*bob(:,n+1)
    last2(:,n) = scs(:,n)*scs(:,n+1)
    last3(:,n) = ism(:,n)*ism(:,n+1)
   end do
;************************************************
  u_peak = Cmp_MTG(:,0,:,:)
  BOB = new(nyear,integer)
  SCS = BOB
  ISM = BOB
  
  do k = 0, nyear-1
    do j=0,nlat-1
      do i=0,nlon-1
        a= ind((diff(k,:,j,i).lt.0).and.(Cmp_MTG(k,:,j,i).gt.0).and.(last(k,:,j,i).gt.0))
        u_peak(k,j,i)=(/a(0)/)  +pens-1
        delete(a)
      end do
    end do

    a= ind((diff1(k,:).lt.0).and.(bob(k,:).gt.0).and.(last1(k,:).gt.0))
    BOB(k) =(/a(0)/)  +pens
    delete(a)
    a= ind((diff2(k,:).lt.0).and.(scs(k,:).gt.0).and.(last2(k,:).gt.0))
    SCS(k) =(/a(0)/)  +pens
    delete(a)
    a= ind((diff3(k,:).lt.0).and.(ism(k,:).gt.0).and.(last3(k,:).gt.0))
    ISM(k) =(/a(0)/)  +pens
    delete(a)
end do

  print(u_peak(:,{10},{100}))

  print(BOB)
  print(SCS)
  print(ISM)

write_table("onset_pentad_by_MTG_"+casename(ll)+".txt", "w", [/"year","BOB","SCS","ISM"/], "%2s %2s %2s %2s")
write_table("onset_pentad_by_MTG_"+casename(ll)+".txt", "a", [/year,BOB,SCS,ISM/], "%0.2i %0.2i %0.2i %0.2i")


pp = new((/3,nyear/), integer)
pp(0,:) = BOB
pp(1,:) = SCS 
pp(2,:) = ISM

; ==============================================================
;  plot (a)
; ==============================================================
 figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/" 
 figname = "calc_onset_pentad_by_MTG_"+casename(ll)

 Colors = (/"green4","red","blue"/)

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 1
plot = new(nplot, graphic)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPattern = 0

res@trXMinF = yrStrt
res@trXMaxF = yrLast
res@tmXBPrecision = 4
res@trYMinF = 20
res@trYMaxF = 35
res@tmYLPrecision = 2

res@gsnRightString = ""
res@gsnStringFontHeightF = 0.025

; res@gsnYRefLine = (/0.8,-0.8/)
; res@gsnYRefLine = (/1,-1/)
; res@gsnYRefLineColor = "gray54"
; res@gsnYRefLineThicknesses = 3
; res@gsnYRefLineDashPattern = 2
res@gsnLeftString = "Onset pentad: MTG reverse"


plot(0) = gsn_csm_xy(wks, year, pp, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"BOB","SCS","ISM"/)

; ;============================================================
    interval = 1
    gx = (/2002,2005/)
    text_pos = 2006
    gy1 = (/34,34/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    gy4 = gy3-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    ;res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ;res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ;res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

draw(plot(0))
frame(wks)

end do
end