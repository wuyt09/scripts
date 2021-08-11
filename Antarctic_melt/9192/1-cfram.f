      program main

      INTEGER,PARAMETER :: zd1=37,IX=240,IY=121,zd18=38,nn=62

c  input data
      real rht_cloud(ix,iy,zd18),rht_base_sw(ix,iy,zd18)
      real rht_base_lw(ix,iy,zd18),rht_base(ix,iy,zd18)
      real rht_warm(ix,iy,zd18),rht_warm_sw(ix,iy,zd18)
      real rht_warm_lw(ix,iy,zd18),rht_o3(ix,iy,zd18)
      real rht_wv(IX,IY,zd18),rht_albedo(IX,IY,zd18)
      real rht_cloud_sw(ix,iy,zd18),rht_cloud_lw(ix,iy,zd18)
      real rht_solar(ix,iy,zd18),rht_co2(ix,iy,zd18)
      real lhflx_base(ix,iy),shflx_base(ix,iy)
      real lhflx_warm(ix,iy),shflx_warm(ix,iy)

!forcing output
      real fc_dyn(ix,iy,zd18),fc_cloud(ix,iy,zd18)
      real fc_wv(ix,iy,zd18),fc_albedo(ix,iy,zd18)
      real fc_o3(ix,iy,zd18),fc_co2(ix,iy,zd18)
      real fc_cloud_sw(ix,iy,zd18),fc_cloud_lw(ix,iy,zd18)
      real fc_solar(ix,iy,zd18)
      real fc_atm_dyn(ix,iy,zd18),fc_sfc_dyn(ix,iy,zd18)
      real fc_lhflx(ix,iy,zd18),fc_shflx(ix,iy,zd18)

!partial temp. output
      real dt_wv(ix,iy,zd18),dt_albedo(ix,iy,zd18)
      real dt_dyn(ix,iy,zd18),dt_shflx(ix,iy,zd18)
      real dt_cloud(ix,iy,zd18),dt_cloud_sw(ix,iy,zd18)
      real dt_cloud_lw(ix,iy,zd18),dt_lhflx(ix,iy,zd18)
      real dt_o3(ix,iy,zd18),dt_ocean(ix,iy,zd18)
      real dt_co2(ix,iy,zd18),dt_solar(ix,iy,zd18)
      real dt_atm_dyn(ix,iy,zd18),dt_sfc_dyn(ix,iy,zd18)

      real drdt(100,100),fc(100),x(ix,iy)
      integer irec,i,j,k,n1,n2,n3,n4,nnn,k1,k2
      character*20 nn_ch


      open ( unit = 21, file =
     & './drdt_ranc_1.dat',
     & form='unformatted', access='direct',recl=100*100)

       do nnn=1,nn

       write(nn_ch,*)nnn
       print*,nn_ch

       open ( unit = 11, file =
     & './baseline_radranc_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access='direct',recl=240*121 )

       open ( unit = 12, file =
     & './albedo_radranc_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access='direct',recl=240*121 )

       open ( unit = 13, file =
     & './wv_radranc_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access='direct',recl=240*121 )

       open ( unit = 14, file =
     & './cloud_radranc_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access='direct',recl=240*121 )

       open ( unit = 16, file =
     & './o3_radranc_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access='direct',recl=240*121 )

       open ( unit = 18, file =
     & './warm_radranc_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access='direct',recl=240*121 )

       open(51,file = './partial_T_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access = 'direct', recl = 240*121)

       open(52,file = './forcing_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access = 'direct', recl = 240*121)

       nt = 1
       xnt= float(nt)

       rht_base(:,:,:)=0.0
       rht_base_sw(:,:,:)=0.0
       rht_base_lw(:,:,:)=0.0
       rht_o3(:,:,:)=0.0
       rht_wv(:,:,:)=0.0
       rht_albedo(:,:,:)=0.0
       rht_cloud(:,:,:)=0.0
       rht_warm(:,:,:)=0.0
       rht_warm_sw(:,:,:)=0.0
       rht_warm_lw(:,:,:)=0.0
       rht_cloud_sw(:,:,:)=0.0
       rht_cloud_lw(:,:,:)=0.0

       do 100 it = 1, nt
          irec=(it-1)*7*zd18
       do l = 1, zd18

          irec=irec+1
          read(11,rec=irec)x !base
          do i = 1, ix
             do j = 1, iy
                rht_base(i,j,l)=rht_base(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(12,rec=irec)x !albedo
          do i = 1, ix
             do j = 1, iy
                rht_albedo(i,j,l)=rht_albedo(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(13,rec=irec)x    !wv
          do i = 1, ix
             do j = 1, iy
                rht_wv(i,j,l)=rht_wv(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(14,rec=irec)x !cloud
          do i = 1, ix
             do j = 1, iy
                rht_cloud(i,j,l)=rht_cloud(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(16,rec=irec)x !o3
          do i = 1, ix
             do j = 1, iy
                rht_o3(i,j,l)=rht_o3(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(18,rec=irec)x !2c02
          do i = 1, ix
             do j = 1, iy
                rht_warm(i,j,l)=rht_warm(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

       enddo

       do l = 1, zd18

          irec=irec+1
          read(11,rec=irec)x !base
          do i = 1, ix
             do j = 1, iy
                rht_base_sw(i,j,l)=rht_base_sw(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(14,rec=irec)x ! clouds
          do i = 1, ix
             do j = 1, iy
                rht_cloud_sw(i,j,l)=rht_cloud_sw(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(18,rec=irec)x !warm
          do i = 1, ix
             do j = 1, iy
                rht_warm_sw(i,j,l)=rht_warm_sw(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

       enddo

       do l = 1, zd18

          irec=irec+1
          read(11,rec=irec)x
          do i = 1, ix
             do j = 1, iy
                rht_base_lw(i,j,l)=rht_base_lw(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(14,rec=irec)x ! clouds
          do i = 1, ix
             do j = 1, iy
                rht_cloud_lw(i,j,l)=rht_cloud_lw(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

          read(18,rec=irec)x !warm
          do i = 1, ix
             do j = 1, iy
                rht_warm_lw(i,j,l)=rht_warm_lw(i,j,l)+x(i,j)!/xnt
             enddo
          enddo

       enddo

 100   continue

       print*, "end of 100 loop"
       print*, "end of input"
       close(11)
       close(12)
       close(13)
       close(14)
       close(16)
       close(18)

       do 2000 j = 1,IY
       do 2000 i = 1,IX

          irec=(j-1)*ix+i
          read(21,rec=irec)drdt
          if(isnan(drdt(1,1)))then
             print*,"drdt in ",i,j,"is NaN"
             irec=irec-1 !drdt_ran has infinite values this location
             read(21,rec=irec)drdt
          endif

          do k = 1,zd18
              n1 = rht_base(i,j,k).eq.(-999.0)
              n2 = rht_warm(i,j,k).eq.(-999.0)
              if(n1.or.n2)then
                  exit
              end if
          end do
          if(k.eq.39)then
              k=k-1
              print*,i,j,"k=39"
          end if
          k1 = k

          do k = 1, zd18+1
            n1 = drdt(k,k).eq.(-999.0)
            if(n1)then
              exit
            end if
          end do
          k2 = k-1

          if(k1.le.k2)then
             nv1=k1
           else
             nv1=k2
          end if

          nv=nv1-1


          rht_base(i,j,nv1)    =rht_base(i,j,zd18)
          rht_base_sw(i,j,nv1) =rht_base_sw(i,j,zd18)
          rht_base_lw(i,j,nv1) =rht_base_lw(i,j,zd18)
          rht_albedo(i,j,nv1)  =rht_albedo(i,j,zd18)
          rht_wv(i,j,nv1)      =rht_wv(i,j,zd18)
          rht_cloud(i,j,nv1)   =rht_cloud(i,j,zd18)
          rht_cloud_sw(i,j,nv1)=rht_cloud_sw(i,j,zd18)
          rht_cloud_lw(i,j,nv1)=rht_cloud_lw(i,j,zd18)
          rht_o3(i,j,nv1)      =rht_o3(i,j,zd18)
          rht_warm(i,j,nv1)    =rht_warm(i,j,zd18)
          rht_warm_sw(i,j,nv1) =rht_warm_sw(i,j,zd18)
          rht_warm_lw(i,j,nv1) =rht_warm_lw(i,j,zd18)

          do k = 1, nv1
             fc(k)=rht_albedo(i,j,k)-rht_base(i,j,k)
             fc_albedo(i,j,k)=fc(k)
          enddo
          call delt_gauss(drdt,fc,nv1)
          do k = 1, nv1
             dt_albedo(i,j,k)=fc(k)
          enddo

c          print*,"after albedo"
          do k = 1, nv1
             fc(k)=rht_wv(i,j,k)-rht_base(i,j,k)
             fc_wv(i,j,k)=fc(k)
          enddo
          call delt_gauss(drdt,fc,nv1)
          do k = 1, nv1
             dt_wv(i,j,k)=fc(k)
          enddo

          do k = 1, nv1
             fc(k)=rht_cloud(i,j,k)-rht_base(i,j,k)
             fc_cloud(i,j,k)=fc(k)
          enddo
          call delt_gauss(drdt,fc,nv1)
          do k = 1, nv1
             dt_cloud(i,j,k)=fc(k)
          enddo

c          print*, "after cloud"
          do k = 1, nv1
             fc(k)=rht_cloud_sw(i,j,k)-rht_base_sw(i,j,k)
             fc_cloud_sw(i,j,k)=fc(k)
          enddo
          call delt_gauss(drdt,fc,nv1)
          do k = 1, nv1
             dt_cloud_sw(i,j,k)=fc(k)
          enddo

          do k = 1, nv1
             fc(k)=rht_cloud_lw(i,j,k)-rht_base_lw(i,j,k)
             fc_cloud_lw(i,j,k)=fc(k)
          enddo
          call delt_gauss(drdt,fc,nv1)
          do k = 1, nv1
             dt_cloud_lw(i,j,k)=fc(k)
          enddo

          do k = 1, nv1
             fc(k)=rht_o3(i,j,k)-rht_base(i,j,k)
             fc_o3(i,j,k)=fc(k)
          enddo
          call delt_gauss(drdt,fc,nv1)
          do k = 1, nv1
             dt_o3(i,j,k)=fc(k)
          enddo

          do k = 1, nv
             fc(k)=-(rht_warm(i,j,k)-rht_base(i,j,k))
             fc_dyn(i,j,k)=fc(k)
          enddo
          fc_dyn(i,j,nv1)=0
          call delt_gauss(drdt,fc,nv1)
          do k = 1, nv1
             dt_dyn(i,j,k)=fc(k)
          enddo

           dt_albedo(i,j,zd18)=dt_albedo(i,j,nv1)
           fc_albedo(i,j,zd18)=fc_albedo(i,j,nv1)
           dt_wv(i,j,zd18)=dt_wv(i,j,nv1)
           fc_wv(i,j,zd18)=fc_wv(i,j,nv1)
           dt_cloud(i,j,zd18)=dt_cloud(i,j,nv1)
           fc_cloud(i,j,zd18)=fc_cloud(i,j,nv1)
           dt_cloud_sw(i,j,zd18)=dt_cloud_sw(i,j,nv1)
           fc_cloud_sw(i,j,zd18)=fc_cloud_sw(i,j,nv1)
           dt_cloud_lw(i,j,zd18)=dt_cloud_lw(i,j,nv1)
           fc_cloud_lw(i,j,zd18)=fc_cloud_lw(i,j,nv1)
           dt_o3(i,j,zd18)=dt_o3(i,j,nv1)
           fc_o3(i,j,zd18)=fc_o3(i,j,nv1)
           dt_dyn(i,j,zd18)=dt_dyn(i,j,nv1)
           fc_dyn(i,j,zd18)=fc_dyn(i,j,nv1)
           dt_shflx(i,j,zd18)=dt_shflx(i,j,nv1)
           fc_shflx(i,j,zd18)=fc_shflx(i,j,nv1)
           dt_lhflx(i,j,zd18)=dt_lhflx(i,j,nv1)
           fc_lhflx(i,j,zd18)=fc_lhflx(i,j,nv1)
           IF(nv.LT.zd1) then
             do l=nv1,zd1
                dt_albedo(i,j,l)=-999
                dt_wv(i,j,l)=-999
                dt_cloud(i,j,l)=-999
                dt_cloud_sw(i,j,l)=-999
                dt_cloud_lw(i,j,l)=-999
                dt_o3(i,j,l)=-999
                dt_dyn(i,j,l)=-999
                fc_albedo(i,j,l)=-999
                fc_wv(i,j,l)=-999
                fc_cloud(i,j,l)=-999
                fc_cloud_sw(i,j,l)=-999
                fc_cloud_lw(i,j,l)=-999
                fc_o3(i,j,l)=-999
                fc_dyn(i,j,l)=-999
             enddo
          ENDIF
 2000  continue

       irec=1
       call out3d(dt_albedo,ix,iy,zd18,51,irec)
       call out3d(fc_albedo,ix,iy,zd18,52,irec)
       irec=irec+zd18

       call out3d(dt_wv,ix,iy,zd18,51,irec)
       call out3d(fc_wv,ix,iy,zd18,52,irec)
       irec=irec+zd18

       call out3d(dt_cloud,ix,iy,zd18,51,irec)
       call out3d(fc_cloud,ix,iy,zd18,52,irec)
       irec=irec+zd18

       call out3d(dt_cloud_sw,ix,iy,zd18,51,irec)
       call out3d(fc_cloud_sw,ix,iy,zd18,52,irec)
       irec=irec+zd18

       call out3d(dt_cloud_lw,ix,iy,zd18,51,irec)
       call out3d(fc_cloud_lw,ix,iy,zd18,52,irec)
       irec=irec+zd18

       call out3d(dt_o3,ix,iy,zd18,51,irec)
       call out3d(fc_o3,ix,iy,zd18,52,irec)
       irec=irec+zd18

       call out3d(dt_dyn,ix,iy,zd18,51,irec)
       call out3d(fc_dyn,ix,iy,zd18,52,irec)

       close(51)
       close(52)
       enddo
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

      subroutine delt_gauss(rad_kern,frc,nv1)
c     calculate delt T (n) by individual forcings by
c     solving the matrix equation
c     frc : input as forcing and out as delt T

      real rad_kern(100,100),frc(100)
      real drdt(nv1,nv1),frc_def(nv1)

      do i=1,nv1
         frc_def(i)=frc(i)
         do j=1,nv1
            drdt(i,j)=rad_kern(i,j)
         enddo
      enddo

      call gaussj(drdt,nv1,nv1,frc_def,1,1)
      do l = 1, nv1
         frc(l)=frc_def(l)
      enddo

      return
      end

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),
     &        indxr(NMAX),ipiv(NMAX)
      REAL  big,dum,pivinv

      do j=1,n
         ipiv(j)=0
      enddo

      do i=1,n
         big=0.
         do j=1,n
            if(ipiv(j).ne.1) then
              do k=1,n
                 if(ipiv(k).eq.0) then
                  if(abs(a(j,k)).ge.big) then
                     big=abs(a(j,k))
                     irow=j
                     icol=k
                  endif
                 endif
              enddo
            endif
         enddo
         ipiv(icol)=ipiv(icol)+1

         if(irow.ne.icol) then
           do l=1,n
              dum=a(irow,l)
              a(irow,l)=a(icol,l)
              a(icol,l)=dum
           enddo
           do l=1,m
              dum=b(irow,l)
              b(irow,l)=b(icol,l)
              b(icol,l)=dum
           enddo
         endif

         indxr(i)=irow
         indxc(i)=icol

         if(a(icol,icol).eq.0.) a(icol,icol)= -3.0517578E-05 ! pause

         pivinv=1./a(icol,icol)
         a(icol,icol)=1.

         do l=1,n
            a(icol,l)=a(icol,l)*pivinv
         enddo

         do l=1,m
            b(icol,l)=b(icol,l)*pivinv
         enddo

         do ll=1,n
            if(ll.ne.icol) then
              dum=a(ll,icol)
              a(ll,icol)=0.
              do l=1,n
               a(ll,l)=a(ll,l)-a(icol,l)*dum
              enddo

              do l=1,m
                b(ll,l)=b(ll,l)-b(icol,l)*dum
              enddo
            endif
          enddo
      enddo

      do l=n,1,-1
         if(indxr(l).ne.indxc(l)) then
           do k=1,n
              dum=a(k,indxr(l))
              a(k,indxr(l))=a(k,indxc(l))
              a(k,indxc(l))=dum
           enddo
         endif
      enddo
      return
      end
