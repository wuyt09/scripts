      program main

      INTEGER,PARAMETER :: zd1=37,IX=240,IY=121,zd18=38,var=7,nn=62

      real ptc(IX,IY,zd18,var),tempb(IX,IY,zd1,nn),tempw(IX,IY,zd1,nn)
	  real force(IX,IY,var),force_err(IX,IY),force_tot(IX,IY)
      real drdt(100,100),fc,dt(100),tsum(100)
      integer nnn,lev_lowest,lev_lowest1,lev_lowest2
      character*20 nn_ch


      open ( unit = 21, file =
     & 'ta_base.dat',
     & form='unformatted', access='direct',recl=240*121)

      open ( unit = 22, file =
     & 'ta_warm.dat',
     & form='unformatted', access='direct',recl=240*121)

      open ( unit = 31, file =
     & 'drdt_ranc_1.dat',
     & form='unformatted', access='direct',recl=100*100)

      irec = 1
      do nnn = 1,nn
      do k = 1,zd1
        read(21,rec=irec)((tempb(i,j,k,nnn),i=1,IX),j=1,IY)
        read(22,rec=irec)((tempw(i,j,k,nnn),i=1,IX),j=1,IY)
        irec = irec+1
      end do
      end do
      tempw = tempw-tempb

      close(21)
      close(22)

      do nnn=1,nn

      write(nn_ch,*)nnn
      print*,nn_ch

       open ( unit = 11, file =
     & 'partial_T_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access='direct',recl=240*121)

       open(51,file =
     & './forcing_by_temp_fb_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access = 'direct', recl=240*121)

       open(52,file =
     & './forcing_by_temp_fb_total_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access = 'direct', recl=240*121)

       open(53,file =
     & './forcing_by_air_err_'//Trim(AdjustL(nn_ch))//'.grd',
     & form='unformatted', access = 'direct', recl=240*121)

       irec = 1
       do v = 1,var
          do k = 1,zd18
	        read(11,rec=irec)((ptc(i,j,k,v),i=1,IX),j=1,IY)
	        irec = irec+1
          end do
       end do
       close(11)

       do 2000 j = 1,IY
       do 2000 i = 1,IX

          irec=(j-1)*IX+i
          read(31,rec=irec)drdt
          if(isnan(drdt(1,1)))then
             irec=irec-1                       !drdt_ran has infinite values at this location
             read(31,rec=irec)drdt
          endif

          do k = 1,zd18+1
	        if(drdt(1,k).eq.-999.0)then
            lev_lowest1 = k
            exit
          end if
          end do

          do k = 1,zd1
            lev_lowest2 = zd1
            if(ptc(i,j,k,1).eq.-999.0)then
              lev_lowest2 = k
              exit
            end if
          end do

          if(lev_lowest1.le.lev_lowest2)then
            lev_lowest = lev_lowest1
          else
            lev_lowest = lev_lowest2
          end if

        print*,"lev_lowest",lev_lowest

        nv1 = lev_lowest-1
        nv  = nv1-1
c	    print*,"nv=",nv

        do v = 1,var

        dt(1:nv) = ptc(i,j,1:nv,v)
        dt(nv1)  = 0 !ptc(i,j,zd18,v)

	      fc = 0
        do k = 1,nv
	       fc = fc+drdt(nv1,k)*dt(k)
	      enddo

        force(i,j,v) = -fc

        end do

        do k = 1,nv
          tsum(k) = ptc(i,j,k,1)+ptc(i,j,k,2)+ptc(i,j,k,3)+ptc(i,j,k,6)+
     &           ptc(i,j,k,7)
        end do

        dt(1:nv) = tempw(i,j,1:nv,nnn)-tsum(1:nv)
        dt(nv1)  = 0

        fc = 0
        do k = 1,nv
           fc = fc+drdt(nv1,k)*dt(k)
        end do

        force_err(i,j) = -fc

        dt(1:nv) = tempw(i,j,1:nv,nnn)
        dt(nv1)  = 0

        fc = 0
        do k = 1,nv
           fc = fc+drdt(nv1,k)*dt(k)
        end do

        force_tot(i,j) = -fc


 2000  continue

       irec = 1
       call out3d(force,IX,IY,var,51,irec)
       irec = 1
       call out2d(force_tot,IX,IY,52,irec)
       irec = 1
       call out2d(force_err,IX,IY,53,irec)

       end do
       end

      subroutine out2d(out,ii,jj,iunit,irec0)
      real out(ii,jj)
      real x(ii,jj)
      integer irec

      irec=irec0
      do i = 1, ii
         do j = 1, jj
            x(i,j)=out(i,j)
         enddo
      enddo
         write(iunit,rec=irec)x
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
