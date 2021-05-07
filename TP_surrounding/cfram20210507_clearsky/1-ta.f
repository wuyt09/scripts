      program base
      implicit none

      integer,parameter:: xt=360,yt=181,z1=37,z18=38,nn=31
      real             :: plev(1:z1),co2(1:nn),co2ts
      real             :: pres(1:xt,1:yt),tro3(1:xt,1:yt,1:z1)
      real             :: q(1:xt,1:yt,1:z1),tem_a(1:xt,1:yt,1:z1)
      real             :: camt(1:xt,1:yt,1:z1),cice(1:xt,1:yt,1:z1)
      real             :: solar(1:xt,1:yt),cliq(1:xt,1:yt,1:z1)
      real             :: swdn_surf(1:xt,1:yt),swup_surf(1:xt,1:yt)
      real             :: t_surf(1:xt,1:yt),hus_s(1:xt,1:yt)
      integer          :: irec,i,j,k,n1,n2,n3,n4,nnn

      plev = (/1.,2.,3.,5.,7.,10.,20.,30.,50.,70.,100.,125.,150.,175.,
     &       200.,225.,250.,300.,350.,400.,450.,500.,550.,
     &       600.,650.,700.,750.,775.,800.,825.,850.,875.,900.,
     &       925.,950.,975.,1000./)

***********************data input*************************************
      open ( unit = 11, file = './data/solarin_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 12, file = './data/ssrd_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 13, file = './data/ssru_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 14, file = './data/t2m_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 15, file = './data/huss_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 16, file = './data/sp_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 17, file = './data/o3_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 18, file = './data/cc_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 19, file = './data/clwc_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 110, file = './data/ciwc_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 111, file = './data/hus_base.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 112, file = './data/t_warm.dat',
     & form='unformatted', access='direct',recl = xt*yt )
      open ( unit = 113, file = './data/co2_base.dat',
     & form='unformatted', access='direct',recl = nn )

      irec = 1
      read(113,rec=irec)co2
      print*,co2
      close(113)

      irec = 1
      read(11,rec=irec)((solar(i,j),i=1,xt),j=1,yt)
      read(12,rec=irec)((swdn_surf(i,j),i=1,xt),j=1,yt)
      read(13,rec=irec)((swup_surf(i,j),i=1,xt),j=1,yt)
      read(14,rec=irec)((t_surf(i,j),i=1,xt),j=1,yt)
      read(15,rec=irec)((hus_s(i,j),i=1,xt),j=1,yt)
      read(16,rec=irec)((pres(i,j),i=1,xt),j=1,yt)

      irec=1
      do k = 1,z1,1
        read(17,rec=irec)((tro3(i,j,k),i=1,xt),j=1,yt)
        irec=irec+1
      enddo

      irec=1
      do k = 1,z1,1
        read(18,rec=irec)((camt(i,j,k),i=1,xt),j=1,yt)
        irec=irec+1
      enddo

      irec=1
      do k = 1,z1,1
        read(19,rec=irec)((cliq(i,j,k),i=1,xt),j=1,yt)
        irec=irec+1
      enddo

      irec=1
      do k = 1,z1,1
        read(110,rec=irec)((cice(i,j,k),i=1,xt),j=1,yt)
        irec=irec+1
      enddo

      irec = 1
      do k = 1,z1,1
        read(111,rec=irec)((q(i,j,k),i=1,xt),j=1,yt)
        irec=irec+1
      enddo

      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(17)
      close(18)
      close(19)
      close(110)

      do nnn = 1,nn

      irec= (nnn-1)*37+1
      do k = 1,z1,1
        read(112,rec=irec)((tem_a(i,j,k),i=1,xt),j=1,yt)
        irec=irec+1
      enddo

      print*,"end of input for case ", nnn

      n1=xt
      n2=yt
      n3=z1
      n4=z18
      co2ts = co2(1)

      call baseline(n1,n2,n3,n4,plev,pres,tro3,
     &  tem_a,q,camt,co2ts,cice,cliq,
     &  solar,swdn_surf,swup_surf,t_surf,hus_s,nnn)
      print*,'case ', nnn, 'finished!'
      end do

      close(112)
      end program
*---------------------------------------------------------------------
      include 'cas_fu_radiation.f'

C     NCLFORTSTART
      subroutine baseline(ix,iy,zd1,zd18,plevel,ps,o3,ta,hus,
     & cld_amt,co2mass,ice_wat,liq_wat,solar_in,
     & swdn_sfc,swup_sfc,tsurf,huss,mm)

      INTEGER ix,iy,zd1,zd18,mm
      character*20 mm_ch
      real plevel(zd1),ps(ix,iy),o3(ix,iy,zd1),ta(ix,iy,zd1)
      real hus(ix,iy,zd1)
      real cld_amt(ix,iy,zd1),ice_wat(ix,iy,zd1),liq_wat(ix,iy,zd1)
      real solar_in(ix,iy)
      real swdn_sfc(ix,iy),swup_sfc(ix,iy),tsurf(ix,iy)
      real huss(ix,iy)
C     NCLEND


c        program main
C       calculate 3D temperature response to  radiative forcing derived
C       from AR4 data, Uses huss ctl run and hus co2 run, and surface pressure co2 run
      include 'para.file'

      integer iseed
      real plev(zd18),co2mass
      real as(mbs), albedo(IX,IY),ee(mbir)
      real dp(IX,IY,zd18)
      real X(IX,IY)
c      real tas(IX,IY),huss(IX,IY),rlus(IX,IY)

      real rad_base(zd18),SWD(IX,IY,zd18),rad_conv(ix,iy,zd18)
      real LWU(IX,IY,zd18),SWU(IX,IY,zd18),LWD(IX,IY,zd18)
      real rad_conv_clr(ix,iy,zd18),swd_clr(ix,iy,zd18)
      real swu_clr(ix,iy,zd18)
      real lwu_clr(ix,iy,zd18),lwd_clr(ix,iy,zd18)
      real sw_conv(ix,iy,zd18),lw_conv(ix,iy,zd18)
      real rad_conv_mc(ix,iy,zd18),swd_mc(ix,iy,zd18)
      real swu_mc(ix,iy,zd18)
      real lwu_mc(ix,iy,zd18),lwd_mc(ix,iy,zd18)
      real sw_conv_mc(ix,iy,zd18),lw_conv_mc(ix,iy,zd18)

      real sw_base(zd18),lw_base(zd18)


      real var(100),var2(100),var3(100),var6(100),area_c(100)
      real water_c(100),ice_c(100),var4(100),var5(100)
      integer no_cloud_out(100,100)
      common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
      common /atmosp/ pp(100), pt(100),ptc(100),ph(100),po(100)
      common /clouds/ pre(100), plwc(100), pde(100), piwc(100)
      common /rains/ prwc(100)
      common /graups/ pgwc(100)
      common /umcon/ umco2, umch4, umn2o
      common /radiat/ fds(100), fus(100), dts(100),
     1                  fdir(100), fuir(100), dtir(100),
     1                  fd(100), fu(100), ht(100)
	  common /radnew/ fdsljh(100), fusljh(100),
     1                  fdirljh(100), fuirljh(100),
     1                  fuljh(100),fdljh(100),
     1                  fsljh(100),firljh(100)
        common /radflx/ fd_u(100), cd_u(100)

        common /clearsky/fuir_clr(100),fdir_clr(100),
     1                   fus_clr(100),fds_clr(100)
        common /totalsky/fuir_tot(100),fdir_tot(100),
     1                   fus_tot(100),fds_tot(100)

        data as / mbs * 1/
        data ee / mbir * 1.0 /

      PI=ATAN(1.)*4.
*
      ss = 1360.89
      sbc= 5.67e-8
      write(mm_ch,*)mm
      print*,mm_ch
*-------------------------------------------------------------------------------
*     OUTPUT
*
       open ( unit = 22, file =
     & './ta_radranc_'//Trim(AdjustL(mm_ch))//'.grd',
     & form='unformatted', access='direct',recl= ix*iy )
       open ( unit = 32, file =
     & './ta_radsfc_ranc_'//Trim(AdjustL(mm_ch))//'.grd',
     & form='unformatted', access = 'direct',recl = ix*iy )

       open ( unit = 51, file =
     & './ta_input_'//Trim(AdjustL(mm_ch))//'.dat',
     & form='unformatted', access = 'direct', recl = ix*iy )
       open ( unit = 52, file =
     & './base_no_cloud_out_1.dat',
     & form='unformatted', access = 'direct', recl = 100*100 )

       print*, "begining"

c      do i = 1, ix
c         do j = 1, iy  !ps is in unit of pa (not in hPa) L/Rv=2.5*1000*1000/461.
c            qs=(0.622*6.11/ps(i,j))*
c    &   exp((2.5*(tsair(i,j)-273.16)/(4.61*tsair(i,j)*2.7316/100.)))
c            huss(i,j)=qs*rhs_air(i,j)/100.
c         enddo
c      enddo


       do i = 1, ix
          do j = 1, iy
             if(swdn_sfc(i,j) .eq. 0.)then
                albedo(i,j)=0.0
                print*,i,j,swup_sfc(i,j)
             else
                albedo(i,j)=swup_sfc(i,j)/swdn_sfc(i,j)
             endif
          enddo
       enddo

cai       do i = 1, ix
cai          do j = 1, iy
c             if(solar_in(i,j) .lt. 20)
c     &          print*,"solar ",i,j,solar_in(i,j)
c             if(swdn_sfc(i,j) .lt. 5)
c     &          print*,"swdn_sfc ",i,j,swdn_sfc(i,j)
c             if(swup_sfc(i,j) .gt. 0)
c     &          print*,"swup_sfc ",i,j,swup_sfc(i,j)
c             if (tsurf(i,j) .lt. 200)
c     &          print*,"tsurf ",i,j,tsurf(i,j)
c             if (ps(i,j) .lt. 500)
c     &          print*,"ps ",i,j,ps(i,j)
cai             if(huss(i,j) .gt. 0.1)then
c                print*,i,j,huss(i,j)
cai                huss(i,j)=0.
cai             endif
c             do k = 1, 26
c              if(hus(i,j,k) .gt. 0.1)
c     &           print*,"air q ",i,j,k,hus(i,j,k)
c           enddo
cai          enddo
cai       enddo

       irec=1
       write(51,rec=irec)solar_in
       irec=irec+1
       write(51,rec=irec)swdn_sfc
       irec=irec+1
       write(51,rec=irec)swup_sfc
       irec=irec+1
       write(51,rec=irec)tsurf
       irec=irec+1
       write(51,rec=irec)huss
       irec=irec+1
       write(51,rec=irec)ps

       irec=irec+1
       call out3d(o3,ix,iy,zd1,51,irec)
       irec=irec+zd1
       call out3d(ta,ix,iy,zd1,51,irec)
       irec=irec+zd1
       call out3d(hus,ix,iy,zd1,51,irec)
       irec=irec+zd1
       call out3d(cld_amt,ix,iy,zd1,51,irec)
       irec=irec+zd1
       call out3d(liq_wat,ix,iy,zd1,51,irec)
       irec=irec+zd1
       call out3d(ice_wat,ix,iy,zd1,51,irec)

       irec22 = 1
       irec32 = 1
       irec52 = 1

       do 3000 job = 1,1
       do 2000 ILAT= 1,IY
       do 2000 ILONG= 1,IX

c       do 2000 ilat = 52, 52
c          do 2000 ilong = 123, 124
c          print*,ilat,ilong


          do 20 i = 1, 100
             pre(i) = 0.0
             plwc(i) = 0.0
             pde(i) = 0.0
             piwc(i) = 0.0
             prwc(i) = 0.0
             pgwc(i) = 0.0
 20       continue              !first calculate clear-sky radiation

          u0 = solar_in(ILONG,ILAT)/ss
          ts = tsurf(ilong,ilat)

          level_lowest = 1
          plev(zd18)= ps(ilong,ilat)/100.0
          do l = 1, zd1
             if(plevel(l).lt.plev(zd18))then
                plev(l) = plevel(l)
                level_lowest = l
             else
                exit
             end if
          enddo
          if(level_lowest .gt. zd1)then
          print*,ilong,ilat,level_lowest,zd1
          endif
          plev(level_lowest+1) = ps(ilong,ilat)/100.0

          p_surf     = ps(ilong,ilat)/100.0

          nv1        = level_lowest+1
          nv         = nv1-1
          ndfs       = nv
          mdfs       = nv+1
          ndfs4      = 4 * ndfs
          ndfs2      = ndfs * 2

          umco2 = co2mass
          umch4 = 1.6
          umn2o = 0.28

c          print*,p_surf,ts
          do l = 1, nv
             pt(l)=ta(ilong,ilat,l)
             ph(l)=hus(ilong,ilat,l)
             pp(l)=plev(l)
             po(l)=o3(ILONG,ILAT,l)
             water_c(l)=liq_wat(ilong,ilat,l)
             if(water_c(l) .lt. 0)water_c(l)=0.0
             ice_c(l)=ice_wat(ilong,ilat,l)
             if(ice_c(l) .lt. 0)ice_c(l) = 0.0
             area_c(l)=cld_amt(ilong,ilat,l)
             if(area_c(l) .lt. 0)area_c(l)=0.0
             if(area_c(l) .lt. 1e-5)then
                water_c(l)=0.0
                ice_c(l)=0.0
                area_c(l)=0.0
             endif
c             xxx2=liq_wat2(ilong,ilat,l)
c             xxx3=liq_wat3(ilong,ilat,l)
c             if(area_c(l) .gt. 0)then !convert from area-average to in_cloud amount
c	         water_c(l)=water_c(l)/area_c(l)
c                 ice_c(l)=ice_c(l)/area_c(l)
c              endif
c	     print*,l,pt(l),ph(l),pp(l),po(l)
c            print*,l,water_c(l),xxx2,xxx3,area_c(l),ice_c(l)
c            xxx1=xxx1*(plev(l+1)-plev(l))*100.*1000/9.81
c            xxx4=xxx4*(plev(l+1)-plev(l))*100.*1000/9.81
          enddo

          pt(nv1) = ts
          ptc(nv1) = pt(nv1)
          pts = pt(nv1)
          ph(nv1) = huss(ILONG,ILAT)

          pp(nv1) = p_surf
          po(nv1) = o3(ILONG,ILAT,nv)

          do imbs = 1,mbs
             as(imbs) = albedo(ILONG,ILAT)
          enddo

          do l = 1, nv
             pre(l) = 15  !10-15 micron meters for water clouds (assumed)
             pde(l) = 40  !40-45 micro meters for ice clouds
c             pmean=0.5*(pp(l)+pp(l+1))
             pmean = pp(l+1)
             dry_density = pmean*100./(287*pt(l))  ! convert back to pascal
             plwc(l) = 1000.*water_c(l)*dry_density  !converting to g/m^3 from kg/kg
             piwc(l) = 1000.*ice_c(l)*dry_density
          enddo

 ! cloudy sky calculation

          read(52,rec=irec52)no_cloud_out
          irec52 = irec52+1
          call S_R_cloudy (u0,as,ss,pts,rad_base,area_c,sw_base,
     &            lw_base,water_c,ice_c,iseed,no_cloud_out)


          DO l=1,nv
c             print*,fuir_tot(l),fdir_tot(l),l
             LWU(ILONG,ILAT,l)=fuir_tot(l)
             LWD(ILONG,ILAT,l)=fdir_tot(l)
             SWU(ILONG,ILAT,l)=fus_tot(l)
             SWD(ILONG,ILAT,l)=fds_tot(l)
             rad_conv(ILONG,ILAT,l)=rad_base(l)
             sw_conv(ilong,ilat,l)=sw_base(l)
             lw_conv(ilong,ilat,l)=lw_base(l)
          ENDDO
          IF(nv.LT.zd1) then
             do l=nv1,zd1
                LWU(ILONG,ILAT,l)=-999
                LWD(ILONG,ILAT,l)=-999
                SWU(ILONG,ILAT,l)=-999
                SWD(ILONG,ILAT,l)=-999
                rad_conv(ILONG,ILAT,l)=-999
                sw_conv(ilong,ilat,l)=-999
                lw_conv(ilong,ilat,l)=-999

             enddo
          ENDIF
          LWU(ILONG,ILAT,zd18)=fuir_tot(nv1)
          LWD(ILONG,ILAT,zd18)=fdir_tot(nv1)
          SWU(ILONG,ILAT,zd18)=fus_tot(nv1)
          SWD(ILONG,ILAT,zd18)=fds_tot(nv1)
          rad_conv(ILONG,ILAT,zd18)=rad_base(nv1)
          sw_conv(ilong,ilat,zd18)=sw_base(nv1)
          lw_conv(ilong,ilat,zd18)=lw_base(nv1)
 2000  continue

       write(888,*)"start to write",job
       call out3d(rad_conv,ix,iy,zd18,22,irec22)
       irec22=irec22+zd18
       call out3d(sw_conv,ix,iy,zd18,22,irec22)
       irec22=irec22+zd18
       call out3d(lw_conv,ix,iy,zd18,22,irec22)
       irec22=irec22+zd18
       call out3d(lwu,ix,iy,zd18,22,irec22)
       irec22=irec22+zd18
       call out3d(lwd,ix,iy,zd18,22,irec22)
       irec22=irec22+zd18
       call out3d(swu,ix,iy,zd18,22,irec22)
       irec22=irec22+zd18
       call out3d(swd,ix,iy,zd18,22,irec22)
       irec22=irec22+zd18

       write(888,*)swd(1,128,1), "ranc"

c surface radiation
       call out3_2d(rad_conv,ix,iy,zd18,zd18,32,irec32)
       irec32=irec32+1
       call out3_2d(lwu,ix,iy,zd18,zd18,32,irec32)
       irec32=irec32+1
       call out3_2d(lwd,ix,iy,zd18,zd18,32,irec32)
       irec32=irec32+1
       call out3_2d(swu,ix,iy,zd18,zd18,32,irec32)
       irec32=irec32+1
       call out3_2d(swd,ix,iy,zd18,zd18,32,irec32)
       irec32=irec32+1
!at toa
       call out3_2d(lwu,ix,iy,zd18,1,32,irec32)
       irec32=irec32+1
       call out3_2d(swd-swu-lwu,ix,iy,zd18,1,32,irec32)
       irec32=irec32+1
       call out3_2d(swu,ix,iy,zd18,1,32,irec32)
       irec32=irec32+1
       call out3_2d(swd,ix,iy,zd18,1,32,irec32)
       irec32=irec32+1

 3000 continue

 9999 continue

      return
      end

      subroutine out3d(out,ii,jj,kk,iunit,irec0)
      real out(ii,jj,kk)
      real x(ii,jj)
      integer irec

      irec=irec0
      do k = 1,kk
         do i = 1, ii
            do j = 1, jj
               x(i,j)=out(i,j,k)
            enddo
         enddo
         write(iunit,rec=irec)x
         irec=irec+1
      enddo
      return
      end

      subroutine out3_2d(out,ii,jj,kk,k0,iunit,irec)
      real out(ii,jj,kk)
      real x(ii,jj)

      do i = 1, ii
         do j = 1, jj
            x(i,j)=out(i,j,k0)
         enddo
      enddo
      write(iunit,rec=irec)x

      return
      end

*------------------------------------------------
      SUBROUTINE S_R(u0,as,ss,pts,rad_base)

*    Subroutine for the calculation of S-R

        include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
        common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
        common /clouds/ pre(100),plwc(100),pde(100), piwc(100)
        common /rains/ prwc(100)
        common /graups/ pgwc(100)
        common /umcon/ umco2, umch4, umn2o
        common /radiat/ fds(100), fus(100), dts(100),
     1                  fdir(100), fuir(100), dtir(100),
     1                  fd(100), fu(100), ht(100)
        common /radnew/ fdsljh(100), fusljh(100),
     1                  fdirljh(100), fuirljh(100),
     1                  fuljh(100),fdljh(100),
     1                  fsljh(100),firljh(100)

        dimension as(mbs), ee(mbir)
        real  rad_base(nv1)
        data ee / mbir * 1.0 /

        call rad (as,u0,ss,pts,ee)

c        print*,fds(1)-fus(1)-fuir(1),fus(1),fuir(1),"clear sky"

        do l=1,nv
           rad_base(l)=ht(l)
        enddo
        rad_base(nv1)=fd(nv1)-fu(nv1)

       return
       end


      SUBROUTINE S_R_MC_cloud(u0,as,ss,pts,rad_base,area_c,sw_base,
     &                      lw_base)

*    Subroutine for the calculation of S-R

        include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
        common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
        common /clouds/ pre(100),plwc(100),pde(100), piwc(100)
        common /rains/ prwc(100)
        common /graups/ pgwc(100)
        common /umcon/ umco2, umch4, umn2o
        common /radiat/ fds(100), fus(100), dts(100),
     1                  fdir(100), fuir(100), dtir(100),
     1                  fd(100), fu(100), ht(100)
        common /radnew/ fdsljh(100), fusljh(100),
     1                  fdirljh(100), fuirljh(100),
     1                  fuljh(100),fdljh(100),
     1                  fsljh(100),firljh(100)

      common /clearsky/fuir_clr(100),fdir_clr(100),
     1                 fus_clr(100),fds_clr(100)
      common /totalsky/fuir_tot(100),fdir_tot(100),
     1                 fus_tot(100),fds_tot(100)

      dimension as(mbs), ee(mbir), area_c(100)
      real  rad_base(nv1), fuir_conv(nv1),fdir_conv(nv1)
      real                 fus_conv(nv1),fds_conv(nv1)
      real sw_base(nv1),lw_base(nv1)
      data ee / mbir * 1.0 /

      call rad (as,u0,ss,pts,ee)

      do l = 1, nv
         clr=fuir_clr(l+1)-fuir_clr(l)
         cld=fuir(l+1)-fuir(l)
         fuir_conv(l)=(1.-area_c(l))*clr+area_c(l)*cld

         clr=fdir_clr(l)-fdir_clr(l+1)
         cld=fdir(l)-fdir(l+1)
         fdir_conv(l)=(1.-area_c(l))*clr+area_c(l)*cld

         clr=fus_clr(l+1)-fus_clr(l)
         cld=fus(l+1)-fus(l)
         fus_conv(l)=(1.-area_c(l))*clr+area_c(l)*cld

         clr=fds_clr(l)-fds_clr(l+1)
         cld=fds(l)-fds(l+1)
         fds_conv(l)=(1.-area_c(l))*clr+area_c(l)*cld
      enddo

      fds_tot(1)=fds(1)    !incoming solar radiation which is not influenced by clouds
      fdir_tot(1)=fdir(1)  ! should be zero at all time
      do l = 1, nv
         fds_tot(l+1)=fds_tot(l)-fds_conv(l)
         fdir_tot(l+1)=fdir_tot(l)-fdir_conv(l)
      enddo

      fuir_tot(nv1)=fuir(nv1)   !upward LW wave radiation at surface is not influenced by clouds
      fus_tot(nv1) = as(1)*fds_tot(nv1)    !upward SW radiation at surface = reflected solar rad.

      do l = nv, 1, -1
         fuir_tot(l)=fuir_tot(l+1)-fuir_conv(l)
         fus_tot(l)=fus_tot(l+1)-fus_conv(l)
      enddo

      do l=1,nv
         rad_base(l) = fds_conv(l)+fus_conv(l)+fdir_conv(l)+fuir_conv(l)
         sw_base(l)=fds_conv(l)+fus_conv(l)
         lw_base(l)=fdir_conv(l)+fuir_conv(l)
      enddo
      l = nv1
      rad_base(l)=fds_tot(l)-fus_tot(l)+fdir_tot(l)-fuir_tot(l)
      sw_base(l)=fds_tot(l)-fus_tot(l)
      lw_base(l)=fdir_tot(l)-fuir_tot(l)
      return
      end


      SUBROUTINE S_R_cloudy(u0,as,ss,pts,rad_base,area_c,sw_base,
     &                      lw_base,water_c,ice_c,iseed,no_cloud_out)

*    Subroutine for the calculation of S-R

        include 'para.file'
        real water_c(100),ice_c(100)

        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
        common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
        common /clouds/ pre(100),plwc(100),pde(100), piwc(100)
        common /rains/ prwc(100)
        common /graups/ pgwc(100)
        common /umcon/ umco2, umch4, umn2o
        common /radiat/ fds(100), fus(100), dts(100),
     1                  fdir(100), fuir(100), dtir(100),
     1                  fd(100), fu(100), ht(100)
        common /radnew/ fdsljh(100), fusljh(100),
     1                  fdirljh(100), fuirljh(100),
     1                  fuljh(100),fdljh(100),
     1                  fsljh(100),firljh(100)

      common /clearsky/fuir_clr(100),fdir_clr(100),
     1                 fus_clr(100),fds_clr(100)
      common /totalsky/fuir_tot(100),fdir_tot(100),
     1                 fus_tot(100),fds_tot(100)

      dimension as(mbs), ee(mbir), area_c(100)

      integer no_cloud(100,100),iseed,iarea_c(nv1),mran
      integer yes_c(nv1),no_c(nv1),no_cloud_out(100,100)

      real  rad_base(nv1)
      real sw_base(nv1),lw_base(nv1)
      data ee / mbir * 1.0 /

      do l = 1, nv
         iarea_c(l)=int(area_c(l)*100)
c         print*,l,iarea_c(l),area_c(l)
      enddo

      yes_c(:)=0
      no_c(:)=0
      no_cloud(:,:)=no_cloud_out(:,:)
c      no_cloud_out(:,:)=-999

c      do igrid = 1, 100
c         do l = 1, nv
c            if (iarea_c(l) .eq. 0) then
c               no_c(l)=no_c(l)+1
c               no_cloud(igrid,l)=1
c               go to 10
c            endif
c            mran=int(100*ran3(iseed))
c            if(mran .lt. iarea_c(l))then
c               if (yes_c(l) .lt. iarea_c(l))then
c                  no_cloud(igrid,l)=0
c                  yes_c(l) = yes_c(l) + 1
c               else
c                  no_c(l)=no_c(l)+1
c                  no_cloud(igrid,l)=1  !in this case, regardless of mran, no clouds
c               endif
c            else
c               if (no_c(l) .lt. (100 - iarea_c(l)))then
c                  no_cloud(igrid,l)=1
c                  no_c(l)=no_c(l)+1
c               else
c                  yes_c(l)=yes_c(l)+1    !in this case, yes cloud for the remaining call of ran3
c                  no_cloud(igrid,l)=0
c               endif
c            endif
c 10         continue
c         enddo
c      enddo

c        no_cloud_out(:,:)=no_cloud(:,:)

c       print*,"no_cloud",no_cloud(1,:)
c       print*,"no_cloud_out",no_cloud_out(1,:)
c      do l = 1, nv
c         print*,yes_c(l),yes_c(l)+no_c(l),iarea_c(l),l
c         test=0.0
c         do igrid = 1, 100
c            test=test+no_cloud(igrid,l)
c         enddo
c         print*,test,no_c(l)
c      enddo

c      do l1= 1, 5
c         l2 = (l1-1)*20+1
c         l3 = l2 + 19
c         print*,l2,l3
c         do l = 1, nv
c            write(6,*)"l=",l,(no_cloud(i,l),i=l2, l3)
c         enddo
c      enddo

      rad_base(:)=0.0
      sw_base(:)=0.0
      lw_base(:)=0.0
      fus_tot(:)=0.0
      fds_tot(:)=0.0
      fuir_tot(:)=0.0
      fdir_tot(:)=0.0

      do igrid = 1, 100
         test=0.0
         do l = 1, nv
            test=test+no_cloud(igrid,l)
            if(no_cloud(igrid,l) .eq. 0)then
             pre(l) = 15  !10-15 micron meters for water clouds (assumed)
             pde(l) = 40  !40-45 micro meters for ice clouds
c             pmean=0.5*(pp(l)+pp(l+1))
             pmean = pp(l+1)
             dry_density = pmean*100./(287*pt(l))  ! convert back to pascal
             plwc(l) = 1000.*water_c(l)*dry_density  !converting to g/m^3 from kg/kg
             piwc(l) = 1000.*ice_c(l)*dry_density
             pgwc(l)=0.0
c               print*,l,no_cloud(igrid,l),igrid
            else
               pre(l) = 0.0
               plwc(l) = 0.0
               pde(l) = 0.0
               piwc(l) = 0.0
               prwc(l) = 0.0
               pgwc(l) = 0.0
c               print*,l,no_cloud(igrid,l),"no clouds"
            endif
         enddo

         call rad (as,u0,ss,pts,ee)

c         print*,fds(1)-fus(1)-fuir(1),fus(1),fuir(1),test,igrid
         do l=1,nv
            rad_base(l)=rad_base(l)+ht(l)/100.0
            sw_base(l)=sw_base(l)+fsljh(l)/100.0
            lw_base(l)=lw_base(l)+firljh(l)/100.0
            fus_tot(l)=fus_tot(l)+fus(l)/100.0
            fds_tot(l)=fds_tot(l)+fds(l)/100.0
            fdir_tot(l)=fdir_tot(l)+fdir(l)/100.
            fuir_tot(l)=fuir_tot(l)+fuir(l)/100.
         enddo
         rad_base(nv1)=rad_base(nv1)+
     &                     (fd(nv1)-fu(nv1))/100.0
         sw_base(nv1)=sw_base(nv1)+fsljh(nv1)/100.0
         lw_base(nv1)=lw_base(nv1)+firljh(nv1)/100.0
         l = nv1
         fus_tot(l)=fus_tot(l)+fus(l)/100.0
         fds_tot(l)=fds_tot(l)+fds(l)/100.0
         fdir_tot(l)=fdir_tot(l)+fdir(l)/100.
         fuir_tot(l)=fuir_tot(l)+fuir(l)/100.
      enddo

c      print*,fds_tot(1)-fus_tot(1)-fuir_tot(1),fus_tot(1),fuir_tot(1)
      return
      end

      function ran3(idum)  !uniform random generator between (0 and 1)
      integer idum         !set idum = any negative value to initialize the sequence
      integer mbig,mseed,mz
      real ran3, fac
      parameter(mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
      integer i,iff,ii,inext,inextp,k
      integer mj,mk,ma(55)
      save iff,inext,inextp,ma
      data iff /0/
      if(idum .lt. 0 .or. iff .eq. 0)then !initialization
         iff=1
c        mj = mseed - iabs(int(real(idum)))  !initialize mas(55) using the seed idum and the large number mseed
         mj = mseed - iabs(idum) !initialize mas(55) using the seed idum and the large number mseed
         mj = mod(mj,mbig)
         ma(55)=mj
         mk=1
         do 11 i = 1, 54
            ii=mod(21*i,55) !now intialize the rest of the table, in a slightly random order
            ma(ii)=mk       !with numbers that are not especially random
            mk=mj-mk
            if(mk .lt. mz)mk = mk + mbig
            mj = ma(ii)
 11      continue
         do 13 k = 1, 4  !we randomize by "warming up the generator"
            do 12 i = 1, 55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
 12         continue
 13         continue
         inext=0 !prepare indices for our first generated number
         inextp=31  !31 is special and cannot be changed.
         idum=1
      endif
      inext=inext+1 !here is where we start, except on initialization.  increment inext, wrapping around 56 to 1
      if(inext .eq. 56)inext = 1
      inextp=inextp+1   !ditto for inextp
      if(inextp .eq. 56)inextp = 1
      mj=ma(inext)-ma(inextp) !now generate a new random number subtractively.
      if(mj .lt. mz)mj = mj+mbig  !be sure that it is in range
      ma(inext)=mj                 !store it
      ran3=mj*fac                   !and output the derived uniform deviate.
      return
      end
