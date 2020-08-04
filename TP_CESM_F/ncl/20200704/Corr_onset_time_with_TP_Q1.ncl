;*************************************************
; NCL tutorial script: vert_1.ncl
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
; x = 2 .gt. 1
; print(x)
; exit
   latS = 0.
   latN = 30.
   lonL = 50.
   lonR = 150.

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

  yrStrt = 1979
  yrLast = 2015
  year = ispan(yrStrt, yrLast, 1)
  nyr  = dimsizes(year)

  levs = 500
  levt = 200

   nptS  = 19
   npnt = 24
   x      = ispan(1,12,1)
; ==============================================================
; (I) Basic Field ==> VP200
; ==============================================================
   ; f3  = addfile("./data/sp.interim.1979.2017.nc","r")
   f3  = addfile("/home/yangsong3/data-observation/wuyt/ERAI/monthly/1.5x1.5/sp_era_mn_197901_201507.nc","r")
   PS0 = short2flt(f3->sp(0:431,::-1,:))

   fils  = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/pressure/air/air.interim.19*.nc /home/yangsong3/data-observation/ERA-interim/pressure/air/air.interim.200*.nc /home/yangsong3/data-observation/ERA-interim/pressure/air/air.interim.201[0-5].nc")
   ; fils = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/pressure/air/air.interim.200*.nc")
   print(fils)

   f  = addfiles(fils,"r")
   t0 = short2flt(f[:]->air(:,::-1,::-1,:))
   printVarSummary(t0)

   ; Ps0 = linint2_Wrap(PS0&longitude, PS0&latitude, PS0, True, t0&lon, t0&lat, 0)
   ; ps0 = clmMonTLL(Ps0)
   ps0 = clmMonTLL(PS0)
   printVarSummary(ps0)


   yr   = new(nyr,integer,-999)
   mn   = new((/nyr,12/),integer,-999)
   dy   = new((/nyr,31/),integer,-999)

   ; yr   = ispan( 1979,2016,1 )
   yr = year
   do n=0,nyr-1
       mn(n,:)   = ispan(1,12,1)
       dy(n,:)   = ispan(1,31,1)
   end do
   ddd4   = day_of_year(yr,mn(:,3),dy(:,0)) ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;要减一！！！
   ddd6   = day_of_year(yr,mn(:,6),dy(:,28))   ;13-42pentad
   ddd12  = day_of_year(yr,mn(:,11),dy(:,30))
   print(ddd4)
   print(ddd6)
   print(ddd12)
   ; exit
;************************************************

   nlev = dimsizes(t0(1,:,1,1))
   nlat = dimsizes(t0(1,1,:,1))
   nlon = dimsizes(t0(1,1,1,:))

   tYPLLL_AMJ = new((/nyr,npnt*5,nlev,nlat,nlon/),float,-9999)

   days = 0
   do n=0,nyr-1
        ddd6(n) = ddd6(n)-1+days
        ddd4(n) = ddd4(n)-1+days
        tYPLLL_AMJ(n,:,:,:,:) = t0(ddd4(n):ddd6(n),:,:,:)
        days    = days+ddd12(n)
   end do
   printVarSummary(tYPLLL_AMJ)
   
   delete(t0)

   t_Pnt     = tYPLLL_AMJ(:,0:npnt-1,:,:,:)
   t_Pnt    = 0
   do n=0,nyr-1
        do k=0,npnt-1
            do i=0,4
;               t_Pnt(n*npnt+k,:,:,:)   = t_Pnt(n*npnt+k,:,:,:) + tYPLLL_AMJ(n,k*5+i,:,:,:)/5.0
                t_Pnt(n,k,:,:,:)   = t_Pnt(n,k,:,:,:) + tYPLLL_AMJ(n,k*5+i,:,:,:)/5.0
            end do
        end do
   end do
   t_Pnt!0 = "year"
   t_Pnt&year = year
   printVarSummary(t_Pnt)

;***************************************************************************
; calc MTG reverse pentad
;***************************************************************************

   lon  = t_Pnt&lon
   dlon = (lon(2)-lon(1))*0.0174533
   dY   = 6371229.0*dlon
   MTG   = center_finite_diff_n (t_Pnt, dY , True,0,3)
   copy_VarCoords(t_Pnt, MTG)
   printVarSummary(MTG)
; 垂直积分
    linlog =1
    pbot   =50000
    ptop   =20000
    lev = t_Pnt&level
    plev = lev*100
    print(plev)

    MTG_column = MTG(:,:,0,:,:)
    MTG_column = 0
    ppppp = ps0(6,:,:)
    do k= 0, nyr-1
    do n=0,npnt-1
        kkkkk=MTG(year|k,time|n,lat|:,lon|:,level|:)
        MTG_column(k,n,:,:) = vibeta(plev,kkkkk,linlog,ppppp,pbot,ptop)
    end do
    end do
    printVarSummary(MTG_column)

; ==============================================================
; (IV) Reversal Pentad
; ==============================================================

    bob = wgt_areaave_Wrap(MTG_column(:,:,{latS2:latN2},{lonL2:lonR2}),1.0,1.0,1)
    scs = wgt_areaave_Wrap(MTG_column(:,:,{latS3:latN3},{lonL3:lonR3}),1.0,1.0,1)
    ism = wgt_areaave_Wrap(MTG_column(:,:,{latS1:latN1},{lonL1:lonR1}),1.0,1.0,1)
    printVarSummary(bob)

    diff  = MTG_column
    last  = MTG_column
    diff1 = bob
    last1 = bob
    boblast = new(dimsizes(bob), logical)
    ; last11 = bob    
    ; last111 = bob
    diff2 = bob
    last2 = bob
    scslast = boblast
    ; last22 = bob
    ; last222 = bob
    diff3 = bob
    last3 = bob
    ismlast = boblast
    ; last33 = bob
    ; last333 = bob

    do n=1,npnt-1
        diff(:,n,:,:) = MTG_column(:,n,:,:)*MTG_column(:,n-1,:,:)
        diff1(:,n) = bob(:,n)*bob(:,n-1)
        diff2(:,n) = scs(:,n)*scs(:,n-1)
        diff3(:,n) = ism(:,n)*ism(:,n-1)
    end do
    do n=0,npnt-4
        last(:,n,:,:) = MTG_column(:,n,:,:)*MTG_column(:,n+1,:,:)
        last1(:,n) = bob(:,n)*bob(:,n+1)
        boblast(:,n) = bob(:,n+2).gt.0 .and. bob(:,n+3).gt.0
        ; last11(:,n) = bob(:,n)*bob(:,n+1)*bob(:,n+2)
        ; last111(:,n) = bob(:,n)*bob(:,n+1)*bob(:,n+2)*bob(:,n+3)
        last2(:,n) = scs(:,n)*scs(:,n+1)
        scslast(:,n) = scs(:,n+2).gt.0 .and. scs(:,n+3).gt.0
        ; last22(:,n) = scs(:,n)*scs(:,n+1)*scs(:,n+2)
        ; last222(:,n) = scs(:,n)*scs(:,n+1)*scs(:,n+2)*scs(:,n+3)
        last3(:,n) = ism(:,n)*ism(:,n+1)
        ismlast(:,n) = ism(:,n+2).gt.0 .and. ism(:,n+3).gt.0
        ; last33(:,n) = ism(:,n)*ism(:,n+1)*ism(:,n+2)
        ; last333(:,n) = ism(:,n)*ism(:,n+1)*ism(:,n+2)*ism(:,n+3)
   end do
;************************************************
    u_peak = MTG_column(:,0,:,:)
    BOB = new( nyr , integer )
    SCS = BOB
    ISM = BOB

    do k=0,nyr-1
        do j=0,nlat-1
        do i=0,nlon-1
    ;       a= ind(MTG(:,j,i).gt.0)
            a= ind((diff(k,:,j,i).lt.0).and.(MTG_column(k,:,j,i).gt.0).and.(last(k,:,j,i).gt.0))
            u_peak(k,j,i)=(/a(0)/)  +nptS-1
            delete(a)
        end do
        end do
        a= ind((diff1(k,:).lt.0).and.(bob(k,:).gt.0).and.(last1(k,:).gt.0).and.boblast(k,:))
        BOB(k) =(/a(0)/)    +nptS
        delete(a)
        a= ind((diff2(k,:).lt.0).and.(scs(k,:).gt.0).and.(last2(k,:).gt.0).and.scslast(k,:))
        SCS(k) =(/a(0)/)    +nptS
        delete(a)
        a= ind((diff3(k,:).lt.0).and.(ism(k,:).gt.0).and.(last3(k,:).gt.0).and.ismlast(k,:))
        ISM(k) =(/a(0)/)    +nptS
        delete(a)
    end do
    print(u_peak(:,{10},{50}))


    print(BOB)
    print(SCS)
    print(ISM)

write_table("onset_pentad_ERAI.txt", "w", [/"year","BOB","SCS","ISM"/], "%4s %2s %2s %2s")
write_table("onset_pentad_ERAI.txt", "a", [/year,BOB,SCS,ISM/], "%0.4i %0.2i %0.2i %0.2i")

pp = new((/3,nyr/), integer)
pp(0,:) = BOB
pp(1,:) = SCS 
pp(2,:) = ISM

; ==============================================================
;  plot (a)
; ==============================================================
 figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200704/" 
 figname = "Corr_onset_time_with_TP_Q1"

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
    interval = 2
    gx = (/2002,2005/)
    text_pos = 2006
    gy1 = (/33,33/)
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



end