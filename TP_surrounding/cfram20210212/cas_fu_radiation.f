CCCC!!!   SUBROUTINE THICKS  CHANGED BY LU JIANHUA              CCCC
CCCCCCC    PT(1~NV)   IS LAYER TEMPERATURE        NOW           CCCC 
c 6-24-98 (1)	
 	subroutine rad (as,u0,ss,pts,ee)
c	subroutine rad ( as, u0, ss, pts, ee, ur )
c 6-24-98 (1)	
c *********************************************************************
c In this radiation scheme,  six  and  12 bands are selected for solar 
c and thermal IR regions, respectively. The spectral division is below: 
c 0.2 - 0.7 um, 0.7 - 1.3 um, 1.3 - 1.9 um, 1.9 - 2.5 um, 2.5 -3.5 um,
c 3.5 - 4.0 um, and 2200 - 1900 cm**-1, 1900 - 1700 cm**-1, 1700 -1400
c cm**-1,  1400 - 1250 cm**-1,  1250 - 1100 cm**-1, 1100 - 980 cm**-1,
c 980 - 800 cm**-1,  800 - 670 cm**-1,  670 - 540 cm**-1, 540 - 400 cm
c **-1,  400 - 280 cm**-1,  280 - 0 cm**-1,  where  the index  for the
c spectral band ( ib = 1, 2, ..., 18 ) is defined.
c
c                       **********************
c                       *  INPUT PARAMETERS  *
c                       **********************
c              as(mbs)   solar surface albedo, mbs = 6
c              u0        cosine of solar zenith angle
c              ss        solar constant ( W / m ** 2 )
c              pts       surface temperature ( K )
c              ee(mbir)  IR surface emissivity, mbir = 12
c              pp(nv1)   atmospheric pressure ( mb )
c              pt(nv1)   atmospheric temperature ( K )
c              ph(nv1)   water vapor mixing ratio ( kg / kg )
c              po(nv1)   ozone mixing ratio ( kg / kg )
c              pre(nv)   effective radius of water cloud ( um )
c              plwc(nv)  liquid water content ( g / m ** 3 )
c              pde(nv)   effective size of ice cloud ( um )
c              piwc(nv)  ice water content ( g / m ** 3 )
c              prwc(nv)  rain water content ( g / m ** 3 )
c              pgwc(nv)  graupel water content ( g / m ** 3 )
c              umco2     concentration of CO2 (ppmv)
c              umch4     concentration of CH4 (ppmv)
c              umn2o     concentration of N2O (ppmv)
c
c Note:  (1)  as(mbs) and ee(mbir) consider the substantial wavelength
c             dependence of surface albedos and emissivities.
c        (2)  For CO2, CH4 and N2O, uniform mixing is assumed  through
c             the atmosphere with concentrations of 330, 1.6 and  0.28
c             ppmv, respectively.  The  concentrations  can be changed
c             through 'common /umcon/ umco2, umch4, umn2o '.
c        (3)  nv, nv1, ndfs, mdfs, ndfs4, mb, mbs, mbir,  and  nc  are  
c             given through 'para.file'. 
c        (4)  nv1 and 1 are the surface and top levels, respectively.
c
c                       **********************
c                       *  OUTPUT PARAMETERS  *
c                       **********************
c              fds(nv1)   downward solar flux ( W / m ** 2 )
c              fus(nv1)   upward solar flux ( W / m **2 )
c              dts(nv)    solar heating rate ( K / day )
c              fdir(nv1)  downward IR flux ( W / m ** 2 )
c              fuir(nv1)  upward IR flux ( W / m **2 )
c              dtir(nv)   IR heating rate ( K / day )
c              fd(nv1)    downward net flux ( W / m ** 2 )
c              fu(nv1)    upward net flux ( W / m **2 )
c              dt(nv)     net heating rate ( K / day )
c
c Note:  Solar, IR, and net represent 0.2 - 0.4 um, 2200 - 0 cm**-1,
c        and  entire spectral regions, respectively.
c
c 6-24-98 (2)	
c                  ********************************* 
c                  * OUTPUT FOR IR WINDOW RADIANCE *
c                  *********************************
c              fiurw(nv1) upward radiance at ur in IR window (W/m**2/Sr)
c
c Note: The IR window is defined as the spectral interval between
c       800 to 1250 cm**-1.
c 6-24-98 (2)	
c
c Fu 07-08-98
c The improved parameterization of cirrus radiative properties in
c both solar and IR spectra (Fu 1996; Fu et al. 1998) has been 
c incorporated into the radiation model.  Note that the definition
c of the generalized effective size used in the new parameterization
c (Eq. 3.10 or Eq. 2.3 in Fu 1996) is different from the mean
c effective size defined in Eq. 2.1 of Fu and Liou (1993).  Now 
c you can make choice between the two versions of cirrus para-
c meterization through the logical variable "fl93i".  Use appropriate
c effective sizes of ice clouds through input "pde" for the two
c different versions of parameterization.
c Fu 07-08-98
c
c *********************************************************************
	include 'para.file'
c Fu 07-08-98
        logical fl93i
c Fu 07-08-98
c 11/4/95 (begin) 8/22/96
c This program is different from "radiatom.f" by adding the modified
c two-stream module (Fu et al. 1996) through "mquadr".
	logical dfsasl, dtsasl
	logical dfsair, dtsair
	logical edding, quadra, hemisp, mquadr
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common / tsfslog / dfsasl, dtsasl, dfsair, dtsair
	common / gtslog / edding, quadra, hemisp, mquadr
c 11/4/95 (end) 8/22/96
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	common /clouds/ pre(100), plwc(100), pde(100), piwc(100)
	common /rains/ prwc(100)
	common /graups/ pgwc(100)
        common /umcon/ umco2, umch4, umn2o 
	common /radiat/ fds(100), fus(100), dts(100),
     1                  fdir(100), fuir(100), dtir(100),
c 6-24-98 (3)	
     1	                fd(100), fu(100), ht(100) 
c    1                  fiurw(100) 
 
	common /radnew/ fdsljh(100), fusljh(100),
     1                  fdirljh(100), fuirljh(100),
     1                  fuljh(100),fdljh(100),
     1                  fsljh(100),firljh(100)
c 3-20-07-ljh 	
	common /dfsout/ fu1(100), fd1(100)
c 6-24-98 (4)	
c 	common /radiance/ fiur(100)  
c 6-24-98 (4)	
	common /planci/ bf(100), bs
	dimension as(mbs), ee(mbir)
c kg(mb) is the number of intervals to perform the g-quadrature in
c each band to consider the nongray gaseous absorption.  In total,
c we need to perform 121 spectral calculations in  the  scattering
c problem for each atmospheric profile.
	dimension kg(mb)
      SAVE kg
	data kg / 10, 8, 12, 7, 12, 5, 
     1            2, 3, 4, 4, 3, 5, 2, 10, 12, 7, 7, 8 /
c 11/4/95 (begin)
	dfsasl = .false.
	dtsasl = .true.
	dfsair = .false.
	dtsair = .true.
c 11/4/95 (end)
c Fu 07-08-98
 	fl93i = .false.
c Fu 07-08-98
	f0 = 1.0 / 3.14159
	do 10 i = 1, nv1
	   fds(i) = 0.0
	   fus(i) = 0.0
	   fdir(i) = 0.0
	   fuir(i) = 0.0
c 6-24-98 (5)	
c          fiurw(i) = 0.0
c 6-24-98 (5)	
10	continue
	call thicks
	call rayle2
	if ( u0 .le. 1.0e-4 ) then
          mbn = mbs + 1
        else
          mbn = 1
        endif
	do 20 ib = mbn, mb
c Fu 07-08-98
           if ( fl93i ) then
	      call ice ( ib )
	   else
              call icenew ( ib )
           endif
c Fu 07-08-98
	   call water ( ib )
	   call rain ( ib )
	   call graup ( ib )
	   call rayle ( ib, u0 )
	   call gascon ( ib )
           if ( ib .gt. mbs ) then
             call planck ( ib, pts )
           endif
	   do 30 ig = 1, kg(ib)
	      call gases ( ib, ig, hk )
	      call comscp 
c 11/4/95 (begin) 8/22/96
	      if ( ib .le. mbs ) then
	        if ( dfsasl ) then
                   call qfts(ib,as(ib),u0,f0 )
                endif
	        if ( dtsasl ) then
	           quadra = .false.
	           hemisp = .false.
                   edding = .true.
		   mquadr = .false.
                   call qftsts(ib,as(ib),u0,f0 )
                endif
	        do 40 i = 1, nv1
                   fds(i) = fds(i) + fd1(i) * hk
                   fus(i) = fus(i) + fu1(i) * hk
40              continue
              else
	        if ( dfsair ) then
                   call qfti (ib,ee(ib-mbs) )
                endif
	        if ( dtsair ) then
	           quadra = .false.
                   edding = .false.
                   hemisp = .false.
		   mquadr = .true.
c 6-24-98 (6)	
c                  call qftisf ( ib,ee(ib-mbs))
c                  call qftisf ( ib, ee(ib-mbs), ur )
c 6-24-98 (6)	
                   call qftits ( ib, ee(ib-mbs) )
                endif
c 11/4/95 (end) 8/22/96
	        do 50 i = 1, nv1
                   fdir(i) = fdir(i) + fd1(i) * hk
                   fuir(i) = fuir(i) + fu1(i) * hk
c 6-24-98 (7)	
c                  fiurw(i) = fiurw(i) + fiur(i) * hk * 0.1591549
c 6-24-98 (7)	
50              continue
              endif
30	   continue  
20	continue
	fuq1 = ss / 1340.0
c In this model, we used the solar spectral irradiance determined by
c Thekaekara (1973), and 1340.0 W/m**2 is the solar energy contained 
c in the spectral region 0.2 - 4.0 um.
c	fuq2 = bs * 0.03 * 3.14159 * ee(12)
c fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
c hk of 0.03.
        fuq2 = 0.0
	do 60 i = 1, nv1
           fds(i) = fds(i) * fuq1
           fus(i) = fus(i) * fuq1
*******************************************************************
***        for STRAT0                            ******************
*           fds(i) = fds(1) 
*           fus(i) = fus(nv1) 
*******************************************************************
           fuir(i) = fuir(i) + fuq2
	   fd(i) = fds(i) + fdir(i)
	   fu(i) = fus(i) + fuir(i)
60	continue

	do 65 i = 1, nv
           fuirljh(i)=fuir(i)
           fdirljh(i+1)=fdir(i)
c	   fdljh(i) = fds(i) + fdirljh(i)
c	   fuljh(i) = fus(i) + fuirljh(i)
65      continue

	fdirljh(1)=0
	fdirljh(nv1)=fdir(nv)
	fuirljh(nv1)=fuir(nv1)

	do i = 1, nv
	   fdljh(i) = fds(i) + fdirljh(i)
           fuljh(i) = fus(i) + fuirljh(i)
	enddo

	fdljh(nv1) = fds(nv1) + fdirljh(nv1)
	fuljh(nv1) = fus(nv1) + fuirljh(nv1)

*Change 5/30/08
          DO i=1,nv1
           fdir(i)=fdirljh(i)
           fuir(i)=fuirljh(i)
           fd(i)=fdljh(i)
           fu(i)=fuljh(i)
          ENDDO

	do 70 i = 1, nv
	   xx = fds(i) -fus(i) - fds(i+1) + fus(i+1)
           fsljh(i)=xx
c          dts(i) = 8.4392 * xx / ( pp(i+1) - pp(i) )
*          xx = fdir(i) -fuir(i) - fdir(i+1) + fuir(i+1)
           xx = fdirljh(i) -fuirljh(i) - fdirljh(i+1) + fuirljh(i+1)
           firljh(i)=xx
c          dtir(i) = 8.4392 * xx / ( pp(i+1) - pp(i) )
c          ht(i) = dts(i) + dtir(i)
cai:  ht now is energy flux convergence in each layer
           ht(i)=fds(i) -fus(i) - fds(i+1) + fus(i+1)
     &          +fdirljh(i)-fuirljh(i)-fdirljh(i+1)+fuirljh(i+1)
70	continue
           fsljh(nv1)=fds(nv1)-fus(nv1)
           firljh(nv1) =  fdirljh(nv1) - fuirljh(nv1)

	return
	end



	subroutine thicks 
c *********************************************************************
c dz is the thickness of a layer in units of km.
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	common /thick/ dz(100)
	do 100 i = 1, nv
*           dz(i) = 0.0146337 * ( pt(i) + pt(i+1) ) 
            dz(i) = 0.0146337 *  2.* pt(i)  
     1	   * alog( pp(i+1) / pp(i) )
100	continue
	return
	end



	block data ice1
c *********************************************************************
c ap and bp are empirical coefficients of Eqs. (2.9) and (2.10) to
c calculate the extiction coefficient (1/m) and single scattering 
c albedo, cps and dps are empirical coefficients of Eq. (2.13) to
c compute the expansion coefficients of the phase function (1, 2, 
c 3, 4) in the solar bands, cpir is the empirical coefficients of 
c Eq. (2.15) to calculate the asymmetry factor in the IR bands (Fu
c and Liou, 1992). The units of mean effective size and ice water
c content are um and g/m*m*m, respectively, in these equations.
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /ic1/ ap(3,mb), bp(4,mb), cps(4,4,mbs), dps(4,mbs),
     1               cpir(4,mbir)
	data ap / -6.656e-3,          3.686,           0.00,
     1            -6.656e-3,          3.686,           0.00,
     1            -6.656e-3,          3.686,           0.00,
     1            -6.656e-3,          3.686,           0.00,
     1            -6.656e-3,          3.686,           0.00,
     1            -6.656e-3,          3.686,           0.00,
     1            -7.770e-3,          3.734,          11.85,
     1            -8.088e-3,          3.717,          17.17,
     1            -8.441e-3,          3.715,          19.48,
     1            -9.061e-3,          3.741,          26.48,
     1            -9.609e-3,          3.768,          34.11,
     1            -1.153e-2,          4.109,          17.32,
     1            -8.294e-3,          3.925,          1.315,
     1            -1.026e-2,          4.105,          16.36,
     1            -1.151e-2,          4.182,          31.13,
     1            -1.704e-2,          4.830,          16.27,
     1            -1.741e-2,          5.541,         -58.42,
     1            -7.752e-3,          4.624,         -42.01 /
	data bp / .10998E-05, -.26101E-07,  .10896E-08, -.47387E-11,
     1            .20208E-04,  .96483E-05,  .83009E-07, -.32217E-09,
     1            .13590E-03,  .73453E-03,  .28281E-05, -.18272E-07,
     1           -.16598E-02,  .20933E-02, -.13977E-05, -.18703E-07,
     1            .46180E+00,  .24471E-03, -.27839E-05,  .10379E-07,
     1            .42362E-01,  .86425E-02, -.75519E-04,  .24056E-06,
     1            .19960E+00,  .37800E-02, -.14910E-04,  .00000E+00,
     1            .30140E+00,  .26390E-02, -.11160E-04,  .00000E+00,
     1            .39080E+00,  .12720E-02, -.55640E-05,  .00000E+00,
     1            .31050E+00,  .26030E-02, -.11390E-04,  .00000E+00,
     1            .20370E+00,  .42470E-02, -.18100E-04,  .00000E+00,
     1            .23070E+00,  .38300E-02, -.16160E-04,  .00000E+00,
     1            .56310E+00, -.14340E-02,  .62980E-05,  .00000E+00,
     1            .52070E+00, -.97780E-03,  .37250E-05,  .00000E+00,
     1            .32540E+00,  .34340E-02, -.30810E-04,  .91430E-07,
     1            .10280E+00,  .50190E-02, -.20240E-04,  .00000E+00,
     1            .39640E+00, -.31550E-02,  .64170E-04, -.29790E-06,
     1            .80790E+00, -.70040E-02,  .52090E-04, -.14250E-06 /
	data cps / .22110E+01, -.10398E-02,  .65199E-04, -.34498E-06,
     1             .32201E+01,  .94227E-03,  .80947E-04, -.47428E-06,
     1             .41610E+01,  .74396E-03,  .82690E-04, -.45251E-06,
     1             .51379E+01,  .51545E-02,  .11881E-04, -.15556E-06,
     1             .22151E+01, -.77982E-03,  .63750E-04, -.34466E-06,
     1             .31727E+01,  .15597E-02,  .82021E-04, -.49665E-06,
     1             .40672E+01,  .25800E-02,  .71550E-04, -.43051E-06,
     1             .49882E+01,  .86489E-02, -.18318E-04, -.59275E-07,
     1             .22376E+01,  .10293E-02,  .50842E-04, -.30135E-06,
     1             .31549E+01,  .47115E-02,  .70684E-04, -.47622E-06,
     1             .39917E+01,  .82830E-02,  .53927E-04, -.41778E-06,
     1             .48496E+01,  .15998E-01, -.39320E-04, -.43862E-07,
     1             .23012E+01,  .33854E-02,  .23528E-04, -.20068E-06,
     1             .31730E+01,  .93439E-02,  .36367E-04, -.38390E-06,
     1             .39298E+01,  .16424E-01,  .10502E-04, -.35086E-06,
     1             .47226E+01,  .25872E-01, -.77542E-04, -.21999E-07,
     1             .27975E+01,  .29741E-02, -.32344E-04,  .11636E-06,
     1             .43532E+01,  .11234E-01, -.12081E-03,  .43435E-06,
     1             .56835E+01,  .24681E-01, -.26480E-03,  .95314E-06,
     1             .68271E+01,  .42788E-01, -.45615E-03,  .16368E-05,
     1             .19655E+01,  .20094E-01, -.17067E-03,  .50806E-06,
     1             .28803E+01,  .36091E-01, -.28365E-03,  .79656E-06,
     1             .34613E+01,  .58525E-01, -.46455E-03,  .13444E-05,
     1             .39568E+01,  .81480E-01, -.64777E-03,  .19022E-05 /
	data dps / .12495E+00, -.43582E-03,  .14092E-04, -.69565E-07,
     1             .12363E+00, -.44419E-03,  .14038E-04, -.68851E-07,
     1             .12117E+00, -.48474E-03,  .12495E-04, -.62411E-07,
     1             .11581E+00, -.55031E-03,  .98776E-05, -.50193E-07,
     1            -.15968E-03,  .10115E-04, -.12472E-06,  .48667E-09, 
     1             .13830E+00, -.18921E-02,  .12030E-04, -.31698E-07 /
	data cpir / .79550,     2.524e-3,    -1.022e-5,     0.000e+0,
     1              .86010,     1.599e-3,    -6.465e-6,     0.000e+0,
     1              .89150,     1.060e-3,    -4.171e-6,     0.000e+0,
     1              .87650,     1.198e-3,    -4.485e-6,     0.000e+0,
     1              .88150,     9.858e-4,    -3.116e-6,     0.000e+0,
     1              .91670,     5.499e-4,    -1.507e-6,     0.000e+0,
     1              .90920,     9.295e-4,    -3.877e-6,     0.000e+0,
     1              .84540,     1.429e-3,    -5.859e-6,     0.000e+0,
     1              .76780,     2.571e-3,    -1.041e-5,     0.000e+0,
     1              .72900,     2.132e-3,    -5.584e-6,     0.000e+0,
     1              .70240,     4.581e-3,    -3.054e-5,     6.684e-8,
     1              .22920,     1.724e-2,    -1.573e-4,     4.995e-7 /
	end

	subroutine ice ( ib)
c *********************************************************************
c ti, wi, and wwi are the optical depth, single scattering albedo,
c and expansion coefficients of the phase function ( 1, 2, 3, and
c 4) due to the scattering of ice clouds for a given layer.
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /clouds/ pre(100), plwc(100), pde(100), piwc(100)
	common /ic1/ ap(3,mb), bp(4,mb), cps(4,4,mbs), dps(4,mbs),
     1               cpir(4,mbir)
	common /thick/ dz(100)
	common /ic/ ti(100), wi(100), wwi(100,4)
	do 10 i = 1, nv
	   if ( piwc(i) .lt. 1.0e-5 ) then
	     ti(i) = 0.0
	     wi(i) = 0.0
	     wwi(i,1) = 0.0
	     wwi(i,2) = 0.0
	     wwi(i,3) = 0.0
	     wwi(i,4) = 0.0
           else
c The constant 1000.0 below is to consider the units of dz(i) is km.
	     fw1 = pde(i)
	     fw2 = fw1 * pde(i)
	     fw3 = fw2 * pde(i)
	     ti(i) = dz(i) * 1000.0 * piwc(i) * ( ap(1,ib) +
     1	     ap(2,ib) / fw1 + ap(3,ib) / fw2 )
	     wi(i) = 1.0 - ( bp(1,ib) + bp(2,ib) * fw1 +
     1	     bp(3,ib) * fw2 + bp(4,ib) * fw3 )
	     if ( ib .le. mbs ) then
	       fd = dps(1,ib) + dps(2,ib) * fw1 +
     1         dps(3,ib) * fw2 + dps(4,ib) * fw3
	       wf1 = cps(1,1,ib) + cps(2,1,ib) * fw1 +
     1         cps(3,1,ib) * fw2 + cps(4,1,ib) * fw3
               wwi(i,1) = ( 1.0 - fd ) * wf1 + 3.0 * fd
	       wf2 = cps(1,2,ib) + cps(2,2,ib) * fw1 +
     1         cps(3,2,ib) * fw2 + cps(4,2,ib) * fw3
               wwi(i,2) = ( 1.0 - fd ) * wf2 + 5.0 * fd
	       wf3 = cps(1,3,ib) + cps(2,3,ib) * fw1 +
     1         cps(3,3,ib) * fw2 + cps(4,3,ib) * fw3
               wwi(i,3) = ( 1.0 - fd ) * wf3 + 7.0 * fd
	       wf4 = cps(1,4,ib) + cps(2,4,ib) * fw1 +
     1         cps(3,4,ib) * fw2 + cps(4,4,ib) * fw3
               wwi(i,4) = ( 1.0 - fd ) * wf4 + 9.0 * fd
             else
	       ibr = ib - mbs
               gg = cpir(1,ibr) + cpir(2,ibr) * fw1 +
     1         cpir(3,ibr) * fw2 + cpir(4,ibr) * fw3
	       x1 = gg
               x2 = x1 * gg
               x3 = x2 * gg
               x4 = x3 * gg
               wwi(i,1) = 3.0 * x1
	       wwi(i,2) = 5.0 * x2
               wwi(i,3) = 7.0 * x3
               wwi(i,4) = 9.0 * x4
	     endif
           endif
10	continue
	return
	end

c Fu 07-08-98	
	block data ice1new
c *********************************************************************
c Following Fu (1996; J. Climate) and Fu et al. (1998; J. Climate),
c ap is the empirical coefficients of Eq. (3.9a) of Fu (1996) and
c Eq. (3.1) of Fu et al. (1998) to calculate the extiction coefficient
c (1/m).  bps is for the single scattering albedo in the solar bands
c (3.9b in Fu) and bpir is for the absorption coefficient (1/m) in the
c IR bands (3.2 in Fu et al.).  cp is the empirical coefficients of
c Eq. (3.9c) in Fu or Eq. (3.3) in Fu et al. to compute the asymmetry
c factor of the phase function.  dps is the empirical coefficients of
c Eq. (3.9d) of Fu to calculate the forward delta-fraction in the 
c solar bands.  The units of generalized effective size and ice water
c content are um and g/m**3, respectively, in these equations.
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /ic1new/ ap(3,mb), bps(4,mbs), bpir(4, mbir),
     1                  cp(4,mb), dps(4,mbs)
        data ap / 
     1           -2.9172062e-05,  2.5192544e+00,  0.0,
     2           -2.2948980e-05,  2.5212550e+00,  0.0,
     3           -2.9772840e-04,  2.5400320e+00,  0.0,
     4            4.2668223e-04,  2.4933372e+00,  0.0,
     5            4.3226531e-04,  2.4642946e+00,  0.0,
     6            9.5918990e-05,  2.5232218e+00,  0.0,
     7		 -2.308881e-03, 2.814002e+00, 1.072211e+00,
     8		 -2.465236e-03, 2.833187e+00,-4.227573e-01,
     9		 -3.034573e-03, 2.900043e+00,-1.849911e+00,
     x		 -4.936610e-03, 3.087764e+00,-3.884262e+00,
     1		 -8.178608e-03, 3.401245e+00,-8.812820e+00,
     2		 -8.372696e-03, 3.455018e+00,-1.516692e+01,
     3		 -1.691632e-03, 2.765756e+00,-8.331033e+00,
     4		 -4.159424e-03, 3.047325e+00,-5.061568e+00,
     5		 -9.524174e-03, 3.587742e+00,-1.068895e+01,
     6		 -1.334860e-02, 4.043808e+00,-2.171029e+01,
     7		  3.325756e-03, 2.601360e+00,-1.909602e+01,
     8		  4.919685e-03, 2.327741e+00,-1.390858e+01 /
	data bps / 
     1  1.3540265e-07,  9.9282217e-08, -7.3843168e-11,  3.3111862e-13,
     2 -2.1458450e-06,  2.1984010e-05, -4.4225520e-09,  1.0711940e-11,
     3  1.4027890e-04,  1.3919010e-03, -5.1005610e-06,  1.4032930e-08,
     4  5.7801650e-03,  2.4420420e-03, -1.1985030e-05,  3.3878720e-08,
     5  2.7122737e-01,  1.9809794e-03, -1.5071269e-05,  5.0103900e-08,
     6  1.6215025e-01,  6.3734393e-03, -5.7740959e-05,  1.9109300e-07 /
	data bpir /
     7	4.346482e-01, 1.721457e-02,-1.623227e-04, 5.561523e-07,
     8  7.428957e-01, 1.279601e-02,-1.391803e-04, 5.180104e-07,
     9  8.862434e-01, 1.226538e-02,-1.523076e-04, 6.000892e-07,
     x  7.152274e-01, 1.621734e-02,-1.868544e-04, 7.078738e-07,
     1  5.874323e-01, 1.876628e-02,-2.045834e-04, 7.510080e-07,
     2  5.409536e-01, 1.949649e-02,-2.050908e-04, 7.364680e-07,
     3  1.195515e+00, 3.350616e-03,-5.266996e-05, 2.233377e-07,
     4  1.466481e+00,-2.129226e-03,-1.361630e-05, 1.193649e-07,
     5  9.551440e-01, 1.309792e-02,-1.793694e-04, 7.313392e-07,
     6  3.003701e-01, 2.051529e-02,-1.931684e-04, 6.583031e-07,
     7  2.005578e-01, 2.132614e-02,-1.751052e-04, 5.355885e-07,
     8  8.869787e-01, 2.118409e-02,-2.781429e-04, 1.094562e-06 /
	data cp / 
     1  7.4812728e-01,  9.5684492e-04, -1.1151708e-06, -8.1557303e-09, 
     2  7.5212480e-01,  1.1045100e-03, -2.9157100e-06, -1.3429900e-09,
     3  7.5320460e-01,  1.8845180e-03, -9.7571460e-06,  2.2428270e-08,
     4  7.7381780e-01,  2.2260760e-03, -1.4052790e-05,  3.7896870e-08,
     5  8.7020490e-01,  1.6645530e-03, -1.4886030e-05,  4.9867270e-08,
     6  7.4212060e-01,  5.2621900e-03, -5.0877550e-05,  1.7307870e-07,
     7  7.962716e-01, 3.003488e-03,-2.082376e-05, 5.366545e-08,
     8  8.472918e-01, 2.559953e-03,-2.182660e-05, 6.879977e-08,
     9  8.741665e-01, 2.455409e-03,-2.456935e-05, 8.641223e-08,
     x  8.522816e-01, 2.523627e-03,-2.149196e-05, 6.685067e-08,
     1  8.609604e-01, 2.200445e-03,-1.748105e-05, 5.176616e-08,
     2  8.906280e-01, 1.903269e-03,-1.733552e-05, 5.855071e-08,
     3  8.663385e-01, 2.797934e-03,-3.187011e-05, 1.217209e-07,
     4  7.984021e-01, 3.977117e-03,-4.471984e-05, 1.694919e-07,
     5  7.363466e-01, 4.798266e-03,-4.513292e-05, 1.525774e-07,
     6  7.260484e-01, 2.664334e-03,-1.251136e-05, 2.243377e-08,
     7  6.891414e-01, 6.192281e-03,-6.459514e-05, 2.436963e-07,
     8  4.949276e-01, 1.186174e-02,-1.267629e-04, 4.603574e-07 /
	data dps / 
     1  1.1572963e-01,  2.5648064e-04,  1.9131293e-06, -1.2460341e-08,
     2  1.1360752e-01,  2.4156171e-04,  2.0185942e-06, -1.2876106e-08,
     3  1.1241170e-01, -1.7635186e-07,  2.1499248e-06, -1.2949304e-08,
     4  1.0855775e-01, -3.2496217e-04,  3.4207304e-06, -1.6247759e-08,
     5  5.7783360e-02, -4.1158260e-04,  4.2361240e-06, -1.7204950e-08,
     6  1.1367129e-01, -1.9711061e-03,  1.6078010e-05, -5.1736898e-08 /
c **********************************************************************
	end

	subroutine icenew ( ib )
c *********************************************************************
c ti, wi, and wwi are the optical depth, single scattering albedo,
c and expansion coefficients of the phase function ( 1, 2, 3, and
c 4) due to the scattering of ice clouds for a given layer.
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /clouds/ pre(100), plwc(100), pdge(100), piwc(100)
	common /ic1new/ ap(3,mb), bps(4,mbs), bpir(4,mbir), 
     1                  cp(4,mb), dps(4,mbs)
	common /thick/ dz(100)
	common /ic/ ti(100), wi(100), wwi(100,4)
	do 10 i = 1, nv
	   if ( piwc(i) .lt. 1.0e-5 ) then
	     ti(i) = 0.0
	     wi(i) = 0.0
	     wwi(i,1) = 0.0
	     wwi(i,2) = 0.0
	     wwi(i,3) = 0.0
	     wwi(i,4) = 0.0
           else
	     fw1 = pdge(i)
	     fw2 = fw1 * pdge(i)
	     fw3 = fw2 * pdge(i)
	     if ( ib .le. mbs ) then
	        tau = dz(i) * 1000.0 * piwc(i) * ( ap(1,ib) +
     1	             ap(2,ib) / fw1 )
	        omega = 1.0 - ( bps(1,ib) + bps(2,ib) * fw1 +
     1	                bps(3,ib) * fw2 + bps(4,ib) * fw3 )
	        asy = cp(1,ib) + cp(2,ib) * fw1 +
     1                cp(3,ib) * fw2 + cp(4,ib) * fw3
	        fd = dps(1,ib) + dps(2,ib) * fw1 +
     1               dps(3,ib) * fw2 + dps(4,ib) * fw3
	        f = 0.5 / omega + fd
	        fw = f * omega
                ti(i) = ( 1.0 - fw ) * tau
                wi(i) = ( 1.0 - f ) * omega / ( 1.0 - fw )
                gg = ( asy - f ) / ( 1.0 - f )
	        x1 = gg
                x2 = x1 * gg
                x3 = x2 * gg
                x4 = x3 * gg
                wwi(i,1) = 3.0 * x1
	        wwi(i,2) = 5.0 * x2
                wwi(i,3) = 7.0 * x3
                wwi(i,4) = 9.0 * x4
             else
	        ibr = ib - mbs
	        betae = piwc(i) * ( ap(1,ib) +
     1	                ap(2,ib) / fw1 + ap(3,ib) / fw2 )
	        betaa = piwc(i) / fw1 * ( bpir(1,ibr) + bpir(2,ibr) * 
     1                  fw1 + bpir(3,ibr) * fw2 + bpir(4,ibr) * fw3 )
	        asy = cp(1,ib) + cp(2,ib) * fw1 +
     1                cp(3,ib) * fw2 + cp(4,ib) * fw3
                ti(i) = dz(i) * 1000.0 * betae
                wi(i) = 1.0 - betaa / betae
                gg = asy
	        x1 = gg
                x2 = x1 * gg
                x3 = x2 * gg
                x4 = x3 * gg
                wwi(i,1) = 3.0 * x1
	        wwi(i,2) = 5.0 * x2
                wwi(i,3) = 7.0 * x3
                wwi(i,4) = 9.0 * x4
	     endif
           endif
10	continue
	return
	end
c Fu 07-08-98	



	block data water1
c *********************************************************************
c bz, wz and gz are the extinction coefficient(1/km), single scattering
c albedo and asymmetry factor for the water clouds (St II, Sc I, St I,
c As, Ns, Sc II, Cu, and Cb) in different bands.   re is the effective 
c radius and fl is the liquid water content (LWC).  See Tables 4.2-4.4 
c of Fu (1991).
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /wat1/ re(nc), fl(nc), bz(nc,mb), wz(nc,mb), gz(nc,mb)
	data re /  4.18,  5.36,  5.89,  6.16, 
     1             9.27,  9.84, 12.10, 31.23 /
	data fl / 0.05, 0.14, 0.22, 0.28,
     1            0.50, 0.47, 1.00, 2.50 /
	data bz /  15.11,  40.25,  59.81,  72.43,
     1             83.69,  73.99, 128.17, 120.91,
     1             15.74,  41.70,  61.52,  74.47,
     1             85.78,  75.59, 130.46, 121.84,
     1             16.38,  43.52,  64.84,  77.97,
     1             87.31,  77.36, 134.30, 124.06,
     1             17.57,  45.78,  66.44,  80.15,
     1             90.49,  79.90, 137.56, 125.92,
     1             18.19,  46.63,  69.39,  82.20,
     1             91.46,  79.99, 138.21, 126.08,
     1             21.30,  51.88,  77.77,  87.02,
     1             94.91,  83.55, 143.46, 128.45,
     1             22.44,  57.35,  84.41, 103.50,
     1            103.49,  84.17, 152.77, 132.07,
     1             18.32,  52.69,  76.67, 100.31,
     1            105.46,  92.86, 157.82, 133.03,
     1             17.27,  50.44,  74.18,  96.76,
     1            105.32,  95.25, 158.07, 134.48,
     1             13.73,  44.90,  67.70,  90.85,
     1            109.16, 105.48, 163.11, 136.21,
     1             10.30,  36.28,  57.23,  76.43,
     1            106.45, 104.90, 161.73, 136.62,
     1              7.16,  26.40,  43.51,  57.24,
     1             92.55,  90.55, 149.10, 135.13,
     1              6.39,  21.00,  33.81,  43.36,
     1             66.90,  63.58, 113.83, 125.65,
     1             10.33,  30.87,  47.63,  60.33,
     1             79.54,  73.92, 127.46, 128.21,
     1             11.86,  35.64,  54.81,  69.85,
     1             90.39,  84.16, 142.49, 135.25,
     1             10.27,  33.08,  51.81,  67.26,
     1             93.24,  88.60, 148.71, 140.42,
     1              6.72,  24.09,  39.42,  51.68,
     1             83.34,  80.72, 140.14, 143.57,
     1              3.92,  14.76,  25.32,  32.63,
     1             60.85,  58.81, 112.30, 145.62 /
	data wz / .999999, .999999, .999999, .999999,
     1            .999998, .999999, .999998, .999997,
     1            .999753, .999700, .999667, .999646,
     1            .999492, .999470, .999344, .998667,
     1            .995914, .994967, .994379, .993842,
     1            .991385, .990753, .988908, .974831,
     1            .983761, .978981, .976568, .974700,
     1            .963466, .959934, .953865, .897690,
     1            .702949, .683241, .679723, .669045,
     1            .642616, .632996, .629776, .588820,
     1            .947343, .929619, .924806, .914557,
     1            .877169, .867047, .853661, .737426,
     1            .919356, .896274, .885924, .881097,
     1            .812772, .781637, .775418, .637341,
     1            .874717, .861122, .847850, .851677,
     1            .787171, .772952, .753143, .618656,
     1            .764750, .752410, .736529, .743435,
     1            .671272, .659392, .639492, .549941,
     1            .807536, .808700, .795994, .805489,
     1            .750577, .755524, .709472, .571989,
     1            .753346, .772026, .767273, .777079,
     1            .751264, .760973, .712536, .568286,
     1            .632722, .676332, .684631, .693552,
     1            .707986, .717724, .682430, .552867,
     1            .288885, .348489, .371653, .380367,
     1            .454540, .465769, .475409, .493881,
     1            .261827, .306283, .321340, .333051,
     1            .392917, .406876, .417450, .484593,
     1            .295804, .339929, .352494, .365502,
     1            .416229, .430369, .435267, .491356,
     1            .301214, .354746, .369346, .381906,
     1            .433602, .447397, .447406, .486968,
     1            .243714, .318761, .344642, .352770,
     1            .427906, .438979, .445972, .477264,
     1            .109012, .187230, .226849, .224976,
     1            .331382, .335917, .374882, .457067 /
	data gz / .838, .839, .844, .847,
     1            .849, .860, .853, .859,
     1            .809, .810, .819, .823,
     1            .823, .849, .833, .843,
     1            .774, .787, .781, .792,
     1            .812, .836, .815, .833,
     1            .801, .802, .793, .793,
     1            .814, .829, .818, .832,
     1            .877, .873, .879, .880,
     1            .885, .899, .891, .908,
     1            .783, .769, .777, .756,
     1            .764, .776, .770, .797,
     1            .818, .805, .824, .830,
     1            .815, .801, .820, .845,
     1            .810, .802, .826, .840,
     1            .829, .853, .840, .868,
     1            .774, .766, .799, .818,
     1            .815, .869, .834, .869,
     1            .734, .728, .767, .797,
     1            .796, .871, .818, .854,
     1            .693, .688, .736, .772,
     1            .780, .880, .808, .846,
     1            .643, .646, .698, .741,
     1            .759, .882, .793, .839,
     1            .564, .582, .637, .690,
     1            .719, .871, .764, .819,
     1            .466, .494, .546, .609,
     1            .651, .823, .701, .766,
     1            .375, .410, .455, .525,
     1            .583, .773, .637, .710,
     1            .262, .301, .334, .406,
     1            .485, .695, .545, .631,
     1            .144, .181, .200, .256,
     1            .352, .562, .413, .517,
     1            .060, .077, .088, .112,
     1            .181, .310, .222, .327 /
	end

	subroutine water ( ib )
c *********************************************************************
c tw, ww, and www are the optical depth, single scattering albedo,
c and expansion coefficients of the phase function ( 1, 2, 3, and
c 4) due to the Mie scattering of water clouds for a given layer. 
c By using the mean single scattering properties of the eight drop
c size distributions in each spectral band, the single scattering
c properties of a water cloud with the given liquid water content
c and effective radius are obtained by interpolating (Eqs. 4.25 -
c 4.27 of Fu, 1991). 
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /clouds/ pre(100), plwc(100), pde(100), piwc(100)
	common /wat1/ re(nc), fl(nc), bz(nc,mb), wz(nc,mb), gz(nc,mb)
	common /thick/ dz(100)
	common /wat/ tw(100), ww(100), www(100,4)
	do 10 i = 1, nv
	   if ( plwc(i) .lt. 1.0e-5 ) then
             tw(i) = 0.0
             ww(i) = 0.0
             www(i,1) = 0.0
             www(i,2) = 0.0
             www(i,3) = 0.0
             www(i,4) = 0.0
           else
	     if ( pre(i) .lt. re(1) ) then
c A cloud with the effective radius smaller than 4.18 um is assumed
c to have an effective radius of 4.18 um with respect to the single
c scattering properties.  
	       tw(i) = dz(i) * plwc(i) * bz(1,ib) / fl(1)
 	       ww(i) = wz(1,ib)
               x1 = gz(1,ib)
               x2 = x1 * gz(1,ib)
               x3 = x2 * gz(1,ib)
	       x4 = x3 * gz(1,ib)
	       www(i,1) = 3.0 * x1
	       www(i,2) = 5.0 * x2
	       www(i,3) = 7.0 * x3
	       www(i,4) = 9.0 * x4
	     elseif ( pre(i) .gt. re(nc) ) then
c A cloud with the effective radius larger than 31.23 um is assumed
c to have an effective radius of 31.18 um with respect to the single
c scattering properties.  
	       tw(i) = dz(i) * plwc(i) * bz(nc,ib) / fl(nc)
	       ww(i) = wz(nc,ib)
	       x1 = gz(nc,ib)
               x2 = x1 * gz(nc,ib)
	       x3 = x2 * gz(nc,ib)
               x4 = x3 * gz(nc,ib)
	       www(i,1) = 3.0 * x1
	       www(i,2) = 5.0 * x2
	       www(i,3) = 7.0 * x3
	       www(i,4) = 9.0 * x4
	     else
	       j = 1
20	       continue
	       if (pre(i) .ge. re(j) .and. pre(i). le. re(j+1)) goto 30
	       j = j + 1
	       goto 20
30	       continue
	       tw(i) = dz(i) * plwc(i) * ( bz(j,ib) / fl(j) + 
     1	       ( bz(j+1,ib) / fl(j+1) - bz(j,ib) / fl(j) ) / 
     1	       ( 1.0 / re(j+1) - 1.0 / re(j) ) * ( 1.0 / pre(i)
     1	       - 1.0 / re(j) ) )
	       ww(i) = wz(j,ib) + ( wz(j+1,ib) - wz(j,ib) ) /
     1	       ( re(j+1) - re(j) ) * ( pre(i) - re(j) )
	       gg = gz(j,ib) + ( gz(j+1,ib) - gz(j,ib) ) /
     1         ( re(j+1) - re(j) ) * ( pre(i) - re(j) )
               x1 = gg
               x2 = x1 * gg
               x3 = x2 * gg
               x4 = x3 * gg
	       www(i,1) = 3.0 * x1
	       www(i,2) = 5.0 * x2
	       www(i,3) = 7.0 * x3
	       www(i,4) = 9.0 * x4
	     endif
           endif
10	continue
	return
	end



	block data rayle1
c *********************************************************************
c ri is the coefficient in Eq.(4.8) of Fu (1991) to compute the optical
c depth due to Rayleigh scattering in the solar bands.
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /ray1/ ri(mbs)
	data ri / 0.9022e-5, 0.5282e-6, 0.5722e-7,
     1	          0.1433e-7, 0.4526e-8, 0.1529e-8 /
	end

	subroutine rayle2 
c *********************************************************************
c trp is P(mb)/T(K)*DZ(m) and the constant 14.6337=R(287)/g(9.806)/2.
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	common /ray2/ trp(100)
	do 100 i = 1, nv
	   trp(i) = 14.6337 * ( pp(i) + pp(i+1) )
     1     * alog( pp(i+1) / pp(i) ) 
100	continue
	return
	end

	subroutine rayle ( ib,u0 )
c *********************************************************************
c tr, wr, and wwr are the optical depth, single scattering albedo,
c and expansion coefficients of the phase function ( 1, 2, 3, and
c 4 ) due to the Rayleigh scattering for a given layer.
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /ray1/ ri(mbs)
	common /ray2/ trp(100)
	common /ray/ tr(100), wr(100), wwr(100,4)
	if ( ib .le. mbs ) then
	  if ( ib .eq. 1 ) then
            x = -3.902860e-6*u0*u0+6.120070e-6*u0+4.177440e-6
	  else
	    x = ri(ib)
	  endif
          do 100 i = 1, nv
	     tr(i) = trp(i) * x
	     wr(i) = 1.0
	     wwr(i,1) = 0.0
	     wwr(i,2) = 0.5
	     wwr(i,3) = 0.0
	     wwr(i,4) = 0.0
100       continue
	else
	  do 200 i = 1, nv
	     tr(i) = 0.0
	     wr(i) = 0.0
	     wwr(i,1) = 0.0
	     wwr(i,2) = 0.0
	     wwr(i,3) = 0.0
	     wwr(i,4) = 0.0
200       continue
	endif
	return
	end



	block data rain1
c *********************************************************************
c brn,  wrnf and  grn  are  the extinction coefficient (1/km),  single
c scattering  albedo  and  asymmetry  factor  for  the rain.  The size
c distribution of  rain  is  in the form of a truncated constant-slope 
c gamma function (Manton and Cotton, 1977)  where rmin = 60 um, rmax =
c 1800 um,  rc = 162 um,  density of water = 1 g/cm**3, and rain water
c content (rwc) = 0.5 g/m**3.
c                        Jan. 19, 1993
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /rai1/ rwc, brn(mb), wrnf(mb), grn(mb)
	data rwc / 0.5 /
	data brn /  1.5377, 1.5377, 1.5379, 1.5385, 1.5396, 1.5417,
     1              1.5454, 1.5478, 1.5512, 1.5559, 1.5600, 1.5642,
     1              1.5647, 1.5741, 1.5862, 1.5993, 1.6149, 1.6765 /
	data wrnf /.999932, .97096, .74627, .56719, .53023, .53815,
     1              .53233, .52884, .53192, .52969, .52716, .52321,
     1              .51904, .53859, .55169, .55488, .55334, .55218 /
	data grn / .88323, .89067, .92835, .96626, .97553, .96626,
     1             .97226, .97663, .97216, .97467, .97745, .98156,
     1             .98584, .96374, .94218, .93266, .92990, .90729 /
	end

	subroutine rain ( ib )
c *********************************************************************
c trn, wrn, and wwrn are the optical depth, single scattering albedo,
c and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
c due to the Mie scattering of rain for a given layer. 
c                        Jan. 19, 1993
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /rains/ prwc(100)
	common /rai1/ rwc, brn(mb), wrnf(mb), grn(mb)
	common /thick/ dz(100)
	common /rai/ trn(100), wrn(100), wwrn(100,4)
        x1 = grn(ib)
        x2 = x1 * grn(ib)
        x3 = x2 * grn(ib)
	x4 = x3 * grn(ib)
	y1 = 3.0 * x1
	y2 = 5.0 * x2
	y3 = 7.0 * x3
	y4 = 9.0 * x4
	do 10 i = 1, nv
	   if ( prwc(i) .lt. 1.0e-5 ) then
             trn(i) = 0.0
             wrn(i) = 0.0
             wwrn(i,1) = 0.0
             wwrn(i,2) = 0.0
             wwrn(i,3) = 0.0
             wwrn(i,4) = 0.0
           else
	     trn(i) = dz(i) * prwc(i) * brn(ib) / rwc
	     wrn(i) = wrnf(ib)
	     wwrn(i,1) = y1
	     wwrn(i,2) = y2
	     wwrn(i,3) = y3
	     wwrn(i,4) = y4
           endif
10	continue
	return
	end



	block data graup1
c *********************************************************************
c bg,  wgf  and  gg are  the  extinction  coefficient (1/km),   single
c scattering  albedo  and  asymmetry  factor for the graupel. The size
c distribution of graupel is in the form of a truncated constant-slope 
c gamma function (Manton and Cotton, 1977)  where rmin = 60 um, rmax =
c 5000 um, rc = 500 um, density of graupel = 0.6 g/cm**3, and  graupel
c water content (gwc) = 0.5 g/m**3.
c                        Jan. 19, 1993
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /gra1/ gwc, bg(mb), wgf(mb), gg(mb)
	data gwc / 0.5 /
	data bg /  0.83939,0.83940,0.83940,0.83941,0.83946,0.83951,
     1             0.83967,0.83979,0.83995,0.84029,0.84058,0.84097,
     1             0.84143,0.84286,0.84418,0.84825,0.85421,0.87477 /
	data wgf / 0.999911,0.97115,0.56192,0.53156,0.52579,0.53846,
     1              0.53296,0.53017,0.53182,0.53180,0.52959,0.52446,
     1              0.52342,0.54914,0.55258,0.54307,0.53160,0.55474 /
	data gg /  0.89218,0.89940,0.96820,0.97816,0.98141,0.96373,
     1             0.97173,0.97559,0.97330,0.97327,0.97626,0.98274,
     1             0.98396,0.94673,0.94213,0.95539,0.97097,0.93183 /
	end

	subroutine graup ( ib )
c *********************************************************************
c tgr, wgr, and wwgr are the optical depth, single scattering albedo,
c and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
c due to the Mie scattering of graupel for a given layer. 
c                        Jan. 19, 1993
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /graups/ pgwc(100)
	common /gra1/ gwc, bg(mb), wgf(mb), gg(mb)
	common /thick/ dz(100)
	common /gra/ tgr(100), wgr(100), wwgr(100,4)
        x1 = gg(ib)
        x2 = x1 * gg(ib)
        x3 = x2 * gg(ib)
	x4 = x3 * gg(ib)
	y1 = 3.0 * x1
	y2 = 5.0 * x2
	y3 = 7.0 * x3
	y4 = 9.0 * x4
	do 10 i = 1, nv
	   if ( pgwc(i) .lt. 1.0e-5 ) then
             tgr(i) = 0.0
             wgr(i) = 0.0
             wwgr(i,1) = 0.0
             wwgr(i,2) = 0.0
             wwgr(i,3) = 0.0
             wwgr(i,4) = 0.0
           else
	     tgr(i) = dz(i) * pgwc(i) * bg(ib) / gwc
             wgr(i) = wgf(ib)
	     wwgr(i,1) = y1
	     wwgr(i,2) = y2
	     wwgr(i,3) = y3
	     wwgr(i,4) = y4
           endif
10	continue
	return
	end



	subroutine gascon ( ib )
c *********************************************************************
c tgm(nv) are the optical depthes due to water vapor continuum absorp-
c tion in nv layers for a given band ib. We include continuum absorp-
c tion in the 280 to 1250 cm**-1 region. vv(11)-vv(17) are the central
c wavenumbers of each band in this region. 
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /con/ tgm(100)
	dimension vv(18)
	data vv / 10*0.0, 1175.0, 1040.0, 890.0, 735.0, 
     1	          605.0, 470.0, 340.0, 0.0 /
	if ( ib .gt. 10 .and. ib .lt. 18 ) then
	   call qopcon ( vv(ib), tgm )
	else
	   do 10 i = 1, nv
              tgm(i) = 0.0
10	   continue
	endif
	return
	end

	subroutine gases ( ib, ig, hk )
c *********************************************************************
c tg(nv) are the optical depthes due to nongray gaseous absorption, in
c nv layers for a given band ib and cumulative probability ig. 
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
        common /umcon/ umco2, umch4, umn2o 
	common /band1/ hk1(10), fk1o3(10)
	common /band2/ hk2(8), c2h2o(3,11,8)
	common /band3/ hk3(12), c3h2o(3,11,12)
	common /band4/ hk4(7), c4h2o(3,11,7)
	common /band5/ hk5(12), c5h2o(3,11,12)
	common /band6/ hk6(5), c6h2o(3,11,5)
	common /band7/ hk7(2), c7h2o(3,19,2)
	common /band8/ hk8(3), c8h2o(3,19,3)
	common /band9/ hk9(4), c9h2o(3,19,4)
	common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
	common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
	common /band12/ hk12(5), c12o3(3,19,5), c12h2o(3,19)
	common /band13/ hk13(2), c13h2o(3,19,2)
	common /band14/ hk14(10), c14hca(3,19,10), c14hcb(3,19,10)
	common /band15/ hk15(12), c15hca(3,19,12), c15hcb(3,19,12)
	common /band16/ hk16(7), c16h2o(3,19,7)
	common /band17/ hk17(7), c17h2o(3,19,7)
	common /band18/ hk18(8), c18h2o(3,19,8)
	common /gas/ tg(100)
	dimension fkg(nv1), fkga(nv1), fkgb(nv1), pq(nv1)
	dimension tg1(nv), tg2(nv), tg3(nv)
	goto ( 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 ) ib 
	stop
1	fk = fk1o3(ig)
	call qopo3s ( fk, tg )
	hk = 619.618 * hk1(ig)
c In this band ( 50000 - 14500 cm**-1 ), we have considered the nongray
c gaseous absorption of O3.    619.618 is the solar energy contained in
c the band in units of Wm**-2.
	goto 20
2	call qks ( c2h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = 484.295 * hk2(ig)
c In this band ( 14500 - 7700 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.  484.295 is the solar energy contained in
c the band in units of Wm**-2.
	goto 20
3	call qks ( c3h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = 149.845 * hk3(ig)
c In this band ( 7700 - 5250 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O. 149.845 is the solar energy contained in
c the band in units of Wm**-2.
	goto 20
4	call qks ( c4h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = 48.7302 * hk4(ig)
c In this band ( 5250 - 4000 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O. 48.7302 is the solar energy contained in
c the band in units of Wm**-2.
	goto 20
5	call qks ( c5h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = 31.6576 * hk5(ig)
c In this band ( 4000 - 2850 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O. 31.6576 is the solar energy contained in
c the band in units of Wm**-2.
	goto 20
6	call qks ( c6h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = 5.79927 * hk6(ig)
c In this band ( 2850 - 2500 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O. 5.79927 is the solar energy contained in
c the band in units of Wm**-2.
	goto 20
7	call qki ( c7h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = hk7(ig)
c In this band ( 2200 - 1900 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
	goto 20
8	call qki ( c8h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = hk8(ig)
c In this band ( 1900 - 1700 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
	goto 20
9	call qki ( c9h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = hk9(ig)
c In this band ( 1700 - 1400 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
	goto 20
10	call qki ( c10h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg1 )
	call qki ( c10ch4, fkg )
	call qopch4 ( fkg, tg2 )
	call qki ( c10n2o, fkg )
	call qopn2o ( fkg, tg3 )
	do 205 i = 1, nv
	   tg(i) = tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
205	continue
	hk = hk10(ig)
c In this band ( 1400 - 1250 cm**-1 ), we have considered the overlapping
c absorption of H2O, CH4, and N2O by approach one of Fu(1991).
	goto 20
11	call qki ( c11h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg1 )
	call qki ( c11ch4, fkg )
	call qopch4 ( fkg, tg2 )
	call qki ( c11n2o, fkg )
	call qopn2o ( fkg, tg3 )
	do 215 i = 1, nv
	   tg(i) = tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
215	continue
	hk = hk11(ig)
c In this band ( 1250 - 1100 cm**-1 ), we have considered the overlapping
c absorption of H2O, CH4, and N2O by approach one of Fu(1991).
	goto 20
12	call qkio3 ( c12o3(1,1,ig), fkg )
	call qopo3i ( fkg, tg1 )
	call qki ( c12h2o, fkg )
	call qoph2o ( fkg, tg2 )
	do 225 i = 1, nv
	   tg(i) = tg1(i) + tg2(i)
225	continue
	hk = hk12(ig)
c In this band ( 1100 - 980 cm**-1 ), we have considered the overlapping
c absorption of H2O and O3 by approach one of Fu(1991).
	goto 20
13	call qki ( c13h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = hk13(ig)
c In this band ( 980 - 800 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
	goto 20
14      do 333 i = 1, nv1
	   if ( pp(i) .ge. 63.1 ) then
             pq(i) = ph(i)
           else
             pq(i) = 0.0
       	   endif
333	continue
	call qki ( c14hca(1,1,ig), fkga )
	call qki ( c14hcb(1,1,ig), fkgb )
	do 343 i = 1, nv1
	   fkg(i) = fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
343	continue
	call qophc ( fkg, tg )
	hk = hk14(ig)
c In this band ( 800 - 670 cm**-1), we have considered the overlapping
c absorption of H2O and CO2 by approach two of Fu(1991).
	goto 20
15      do 353 i = 1, nv1
	   if ( pp(i) .ge. 63.1 ) then
             pq(i) = ph(i)
           else
             pq(i) = 0.0
       	   endif
353	continue
	call qki ( c15hca(1,1,ig), fkga )
	call qki ( c15hcb(1,1,ig), fkgb )
	do 363 i = 1, nv1
	   fkg(i) = fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
363	continue
	call qophc ( fkg, tg )
	hk = hk15(ig)
c In this band ( 670 - 540 cm**-1), we have considered the overlapping
c absorption of H2O and CO2 by approach two of Fu(1991).
	goto 20
16	call qki ( c16h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = hk16(ig)
c In this band ( 540 - 400 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
	goto 20
17	call qki ( c17h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = hk17(ig)
c In this band ( 400 - 280 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
	goto 20
18	call qki ( c18h2o(1,1,ig), fkg )
	call qoph2o ( fkg, tg )
	hk = hk18(ig)
c In this band ( 280 - 000 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
20	continue
	return
	end

	subroutine qks ( coefks, fkg )
c *********************************************************************
c fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
c **-1 for a given cumulative probability in nv1 layers. coefks(3,11)
c are the coefficients to calculate the absorption coefficient at the
c temperature t for the 11 pressures by
c         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
c and the absorption coefficient at conditions other than those eleven
c pressures is interpolated linearly with pressure (Fu, 1991).
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	dimension coefks(3,11)
	dimension fkg(nv1)
	dimension stanp(11)
      	data stanp / 10.0, 15.8, 25.1, 39.8, 63.1, 100.0,
     1	             158.0, 251.0, 398.0, 631.0, 1000.0 /
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefks(1,1) + coefks(2,1) * ( pt(i) - 245.0 )
     1       + coefks(3,1) * ( pt(i) - 245.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(11) ) then
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefks(1,10) + coefks(2,10) * ( pt(i) - 245.0 )
     1	     + coefks(3,10) * y1 )
    	     x2 = exp ( coefks(1,11) + coefks(2,11) * ( pt(i) - 245.0 )
     1	     + coefks(3,11) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(11) - stanp(10) )
     1	     * ( pp(i) - stanp(10) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefks(1,i1-1) + coefks(2,i1-1) * (pt(i)-245.0)
     1	     + coefks(3,i1-1) * y1 )
	     x2 = exp ( coefks(1,i1) + coefks(2,i1) * ( pt(i) - 245.0 )
     1	     + coefks(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     1	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end

	subroutine qki ( coefki, fkg )
c *********************************************************************
c fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
c **-1 for a given cumulative probability in nv1 layers. coefki(3,19)
c are the coefficients to calculate the absorption coefficient at the
c temperature t for the 19 pressures by
c         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
c and the absorption coefficient at  conditions  other  than  those 19
c pressures is interpolated linearly with pressure (Fu, 1991).
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	dimension coefki(3,19)
	dimension fkg(nv1)
	dimension stanp(19)
	data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 
     1	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,
     1	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefki(1,1) + coefki(2,1) * ( pt(i) - 245.0 )
     1       + coefki(3,1) * ( pt(i) - 245.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(19) ) then
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefki(1,18) + coefki(2,18) * ( pt(i) - 245.0 )
     1	     + coefki(3,18) * y1 )
    	     x2 = exp ( coefki(1,19) + coefki(2,19) * ( pt(i) - 245.0 )
     1	     + coefki(3,19) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )
     1	     * ( pp(i) - stanp(18) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefki(1,i1-1) + coefki(2,i1-1) * (pt(i)-245.0)
     1	     + coefki(3,i1-1) * y1 )
	     x2 = exp ( coefki(1,i1) + coefki(2,i1) * ( pt(i) - 245.0 )
     1	     + coefki(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     1	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end

	subroutine qkio3 ( coefki, fkg )
c *********************************************************************
c fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
c **-1 for a given cumulative probability in nv1 layers. coefki(3,19)
c are the coefficients to calculate the absorption coefficient at the
c temperature t for the 19 pressures by
c         ln k = a + b * ( t - 250 ) + c * ( t - 250 ) ** 2
c and the absorption coefficient at  conditions  other  than  those 19
c pressures is interpolated linearly with pressure (Fu, 1991).
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	dimension coefki(3,19)
	dimension fkg(nv1)
	dimension stanp(19)
	data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 
     1	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,
     1	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefki(1,1) + coefki(2,1) * ( pt(i) - 250.0 )
     1       + coefki(3,1) * ( pt(i) - 250.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(19) ) then
	     y1 = ( pt(i) - 250.0 ) * ( pt(i) - 250.0 )
	     x1 = exp ( coefki(1,18) + coefki(2,18) * ( pt(i) - 250.0 )
     1	     + coefki(3,18) * y1 )
    	     x2 = exp ( coefki(1,19) + coefki(2,19) * ( pt(i) - 250.0 )
     1	     + coefki(3,19) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )
     1	     * ( pp(i) - stanp(18) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( pt(i) - 250.0 ) * ( pt(i) - 250.0 )
	     x1 = exp ( coefki(1,i1-1) + coefki(2,i1-1) * (pt(i)-250.0)
     1	     + coefki(3,i1-1) * y1 )
	     x2 = exp ( coefki(1,i1) + coefki(2,i1) * ( pt(i) - 250.0 )
     1	     + coefki(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     1	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end

	subroutine qopo3s ( fk, tg )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	dimension tg(nv)
	fq = 238.08 * fk
	do 10 i = 1, nv
 	   tg(i) = ( po(i) + po(i+1) ) * ( pp(i+1) - pp(i) ) * fq
10	continue
c	do 20 i = 1, nv
c	   tg(i) = tg(i) * 476.16 * fk
c20	continue
c 476.16 = 2.24e4 / M * 10.0 / 9.8, where M = 48 for O3.
	return
	end

	subroutine qoph2o ( fkg, tg )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100),po(100)
	dimension fkg(nv1)
	dimension tg(nv)
	do 10 i = 1, nv
**ljh     tg(i) = ( fkg(i) * ph(i) + fkg(i+1) * ph(i+1) )
	   tg(i) =  2. * fkg(i) * ph(i) 
     1	   * ( pp(i+1) - pp(i) ) * 634.9205
10	continue
c	do 20 i = 1, nv
c	   tg(i) = tg(i) * 1269.841
c20	continue
c 1269.841 = 2.24e4 / M * 10.0 / 9.8, where M = 18 for H2O.
	return
	end

	subroutine qopch4 ( fkg, tg )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	dimension fkg(nv1)
	dimension tg(nv)
	do 10 i = 1, nv
	   tg(i) = ( fkg(i)+fkg(i+1) ) * ( pp(i+1)-pp(i) ) * 6.3119e-4
10	continue
c	do 20 i = 1, nv
c	   tg(i) = tg(i) * 1.26238e-3
c20	continue
c 1.26238e-3 = 2.24e4 / M * 10.0 / 9.8 * 1.6e-6 * M / 28.97, where 
c M = 16 for CH4.
	return
	end

	subroutine qopn2o ( fkg, tg )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100),po(100)
	dimension fkg(nv1)
	dimension tg(nv)
	do 10 i = 1, nv
	   tg(i) = ( fkg(i)+fkg(i+1) ) * ( pp(i+1)-pp(i) ) * 1.10459e-4
10	continue
c	do 20 i = 1, nv
c	   tg(i) = tg(i) * 2.20918e-4
c20	continue
c 2.20918e-4 = 2.24e4 / M * 10.0 / 9.8 * 0.28e-6 * M / 28.97, where
c M = 44 for N2O.
	return
	end

	subroutine qopo3i ( fkg, tg )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100),po(100)
	dimension fkg(nv1)
	dimension tg(nv)
	do 10 i = 1, nv
	   tg(i) = ( fkg(i) * po(i) + fkg(i+1) * po(i+1) )
     1	   * ( pp(i+1) - pp(i) ) * 238.08
10	continue
c	do 20 i = 1, nv
c	   tg(i) = tg(i) * 476.16
c20	continue
	return
	end

	subroutine qophc ( fkg, tg )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100),po(100)
	dimension fkg(nv1)
	dimension tg(nv)
	do 10 i = 1, nv
	   tg(i) = ( fkg(i) + fkg(i+1) ) * ( pp(i+1) - pp(i) ) * 0.5
10	continue
c See page 86 of Fu (1991).
	return
	end

	subroutine qopcon ( vv, tg )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	dimension ff(nv1), pe(nv1)
	dimension tg(nv)
	x = 4.18
	y = 5577.8
	z = 0.00787
	r = 0.002
	s = ( x + y * exp ( - z * vv ) ) / 1013.25
	do 3 i = 1, nv1
 	   pe(i) = pp(i) * ph(i) / ( 0.622 + 0.378 * ph(i) )
	   w = exp ( 1800.0 / pt(i) - 6.08108 )
	   ff(i) = s * ( pe(i) + r * pp(i) ) * w
3	continue
	do 5 i = 1, nv
**ljh    tg(i) = ( ff(i) * ph(i) + ff(i+1) * ph(i+1) )*
	   tg(i) = 2.* ff(i) * ph(i) *
     1	   ( pp(i+1) - pp(i) ) * 0.5098835
5	continue
c	do 7 i = 1, nv
c	   tg(i) = tg(i) * 10.0 / 9.80616
c7	continue
	return
	end
c	function fk ( v, e, p, t )
c The units of fk is cm**2/g. See Eq. (A.19) of Fu (1991).
c	x = 4.18
c	y = 5577.8
c	z = 0.00787
c	r = 0.002
c	w = exp ( 1800.0 / t - 6.08108 )
c	fk = ( x + y * exp ( -z * v ) ) * ( e + r * p ) * w / 1013.25
c	return
c	end



	subroutine comscp
c *********************************************************************
c This subroutine is used to  COMbine Single-Scattering Properties  due
c to  ice crystals,  water droplets, and  Rayleigh molecules along with
c H2O continuum absorption and nongray gaseous absorption.  See Section
c 3.4 of Fu (1991). wc, wc1, wc2, wc3, and wc4, are total (or combined)
c single - scattering  albedo,  and   expansion   coefficients  of  the
c phase function ( 1, 2, 3, and 4 ) in nv layers. tt(nv) are the normal
c optical depth ( from the top of the atmosphere to a given level ) for
c level 2 - level nv1( surface ). The single-scattering  properties  of
c rain and graupel are also incorporated in ( Jan. 19, 1993 ).
c *********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /ic/ ti(100), wi(100), wwi(100,4)
	common /wat/ tw(100), ww(100), www(100,4)
	common /rai/ trn(100), wrn(100), wwrn(100,4)
	common /gra/ tgr(100), wgr(100), wwgr(100,4)
	common /ray/ tr(100), wr(100), wwr(100,4)
        common /con/ tgm(100)
	common /gas/ tg(100)
	common /dfsin/ wc1(100), wc2(100), wc3(100), wc4(100), 
     1	               wc(100), tt(100)
	dimension tc(nv)
	do 10 i = 1, nv
	   tc(i) = ti(i) + tw(i) + tr(i) + tgm(i) + tg(i) +
     1             trn(i) + tgr(i)
	   tis = ti(i) * wi(i)
	   tws = tw(i) * ww(i) 
           trns = trn(i) * wrn(i)
           tgrs = tgr(i) * wgr(i)
           fw = tis + tws + tr(i) + trns + tgrs
	   wc(i) =  fw / tc(i)
           if ( fw .lt. 1.0e-20 ) then
             wc1(i) = 0.0
             wc2(i) = 0.0
             wc3(i) = 0.0
             wc4(i) = 0.0
	   else
             wc1(i) = ( tis * wwi(i,1) + tws * www(i,1) + 
     1	     tr(i) * wwr(i,1) + trns * wwrn(i,1) + tgrs * wwgr(i,1) )/fw
             wc2(i) = ( tis * wwi(i,2) + tws * www(i,2) +
     1	     tr(i) * wwr(i,2) + trns * wwrn(i,2) + tgrs * wwgr(i,2) )/fw
             wc3(i) = ( tis * wwi(i,3) + tws * www(i,3) +
     1	     tr(i) * wwr(i,3) + trns * wwrn(i,3) + tgrs * wwgr(i,3) )/fw
             wc4(i) = ( tis * wwi(i,4) + tws * www(i,4) +
     1	     tr(i) * wwr(i,4) + trns * wwrn(i,4) + tgrs * wwgr(i,4) )/fw
	   endif
10	continue
	tt(1) = tc(1)
	do 20 i = 2, nv
	   tt(i) = tt(i-1) + tc(i)
20	continue
	return
	end



c	function planck1 ( t, w )
c **********************************************************************
c t is the temperature (K), w is the wavenumber (cm-1), and planck1 is
c the blackbody intensity function (W/m**2/Sr/cm-1).  See Eq. (2.8) of
c Fu (1991).
c **********************************************************************
c	a = 1.19107e-8
c	b = 1.43884
c	planck1 = a * w * w * w / ( exp ( b * w / t ) - 1.0 )
c	return
c	end

c	function bt ( t, ve, nd )
c **********************************************************************
c bt (W/m**2/Sr) is the blackbody intensity function integrated over a
c given band, which has a band width of nd*10 (cm-1) from the ve (cm-1).
c **********************************************************************
c	v1 = ve
c	bt = 0.0
c	do 10 j = 1, nd
c	   v2 = v1 - 10.0
c	   w = ( v1 + v2 ) * 0.5
c	   x = planck1 ( t, w )
c	   bt = bt + x
c	   v1 = v2
c10	continue
c	bt = bt * 10.0
c	return
c	end

	subroutine planck ( ib, pts )
c **********************************************************************
c bf and bs are the blackbody intensity function integrated over the
c band ib at the nv1 levels and at the surface, respectively.    The
c units of bf and bs are W/m**2/Sr. nd*10 is the band width from ve.
c **********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /atmosp/ pp(100), pt(100),ptc(100),ph(100), po(100)
	common /planci/ bf(100), bs
	dimension ve(mbir), nd(mbir), bt(nv1)
	data ve / 2200.0, 1900.0, 1700.0, 1400.0, 1250.0, 1100.0,
     1            980.0, 800.0, 670.0, 540.0, 400.0, 280.001 /
	data nd / 30, 20, 30, 15, 15, 12,
     1            18, 13, 13, 14, 12, 28 /
        nv11 = nv1 + 1
	ibr = ib - mbs
        bts = 0.0
	do 10 i = 1, nv1
           bt(i) = 0.0
10	continue
	v1 = ve(ibr)
	do 20 j = 1, nd(ibr)
	   v2 = v1 - 10.0
	   w = ( v1 + v2 ) * 0.5
	   fq1 = 1.19107e-8 * w * w * w
	   fq2 = 1.43884 * w
	   do 30 i = 1, nv11
              if ( i .eq. nv11 ) then
	        x = fq1 / ( exp ( fq2 / pts ) - 1.0 )
	        bts = bts + x
              else
	        x = fq1 / ( exp ( fq2 / pt(i) ) - 1.0 )
                bt(i) = bt(i) + x
              endif
30	    continue
              v1 = v2
20	continue
	do 40 i = 1, nv1
          bf(i) = bt(i) * 10.0
40	continue
        bs = bts * 10.0
	return
	end



c **********************************************************************
c Double-Gauss quadratures and weights (Sykes, 1951).
c **********************************************************************
	block data
	common /dis/ a(4)
	common /point/ u(4)
	data a / 0.5, 0.5, 0.5, 0.5 /
	data u / -0.7886752, -0.2113247, 0.2113247, 0.7886752 /
	end

c *********************************************************************
c p0, p1, p2 and p3 are Legendre polynomials for l = 1, 2, 3.
c *********************************************************************
c	function p0 ( x )
c	p0 = 1.0
c	return
c	end
c	function p1 ( x )
c	p1 = x
c	return
c	end
c	function p2 ( x )
c	p2 = 1.5 * x * x - 0.5
c	return
c	end
c	function p3 ( x )
c	p3 = ( 2.5 * x * x - 1.5 ) * x
c	return
c	end

c **********************************************************************
c p0d(4), p1d(4), p2d(4), and p3d(4) are Legendre polynomials p0(x), 
c p1(x), p2(x), and p3(x) when x = u(1), u(2), u(3), and u(4).
c **********************************************************************
      block data legend
      common /legen/ p0d(4), p1d(4), p2d(4), p3d(4)
      data p0d /  .100000E+01,  .100000E+01,  .100000E+01, .100000E+01 /
      data p1d / -.788675E+00, -.211325E+00,  .211325E+00, .788675E+00 /
      data p2d /  .433013E+00, -.433013E+00, -.433013E+00, .433013E+00 /
      data p3d / -.433940E-01,  .293394E+00, -.293394E+00, .433940E-01 /
      end

c *********************************************************************
c p11d(4,4), p22d(4,4), and p33d(4,4) are defined as 0.5*p1d(i)*p1d(j),
c 0.5*p2d(i)*p2d(j), and 0.5*p3d(i)*p3d(j), respectively.
c *********************************************************************
      block data legenf
      common /legen1/ p11d(4,4), p22d(4,4), p33d(4,4)
      data p11d / .311004E+00, .833334E-01,-.833334E-01,-.311004E+00,
     1            .833334E-01, .223291E-01,-.223291E-01,-.833334E-01,
     1           -.833334E-01,-.223291E-01, .223291E-01, .833334E-01,
     1           -.311004E+00,-.833334E-01, .833334E-01, .311004E+00 /
      data p22d / .937501E-01,-.937501E-01,-.937501E-01, .937501E-01,
     1           -.937501E-01, .937501E-01, .937501E-01,-.937501E-01,
     1           -.937501E-01, .937501E-01, .937501E-01,-.937501E-01,
     1            .937501E-01,-.937501E-01,-.937501E-01, .937501E-01 /
      data p33d / .941520E-03,-.636577E-02, .636577E-02,-.941520E-03,
     1           -.636577E-02, .430400E-01,-.430400E-01, .636577E-02,
     1            .636577E-02,-.430400E-01, .430400E-01,-.636577E-02,
     1           -.941520E-03, .636577E-02,-.636577E-02, .941520E-03 /
      end

c **********************************************************************
c coefficient calculations for four first-order differential equations.
c **********************************************************************
	subroutine coeff1 
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dis/ a(4)
	common /point/ u(4)
        common /legen/ p0d(4), p1d(4), p2d(4), p3d(4)
        common /legen1/ p11d(4,4), p22d(4,4), p33d(4,4)
	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
	common /coedf1/ b(4,3)
	dimension c(4,5)
	x = 0.5 * w
	w0w = x
	w1w = x * w1
	w2w = x * w2
	w3w = x * w3
	if ( ib .le. mbs ) then
	  fw = u0 * u0
	  q1 = - w1w * u0
	  q2 = w2w * ( 1.5 * fw - 0.5 )
	  q3 = - w3w * ( 2.5 * fw - 1.5 ) * u0
	endif
	fq = 0.5 * w0w
	do 10 i = 3, 4
	   do 20 j = 1, 4
	      c(i,j) = fq + w1w * p11d(i,j) +
     1	      w2w * p22d(i,j) + w3w * p33d(i,j) 
	      if ( i .eq. j ) then 
   	        c(i,j) = ( c(i,j) - 1.0 ) / u(i)
	      else
	        c(i,j) = c(i,j) / u(i)
	      endif
20	   continue
10	continue
	do 30 i = 1, 4
	   if ( ib .le. mbs ) then
	     c(i,5) = w0w + q1 * p1d(i) +
     1	     q2 * p2d(i) + q3 * p3d(i) 
	   else
	     c(i,5) = 1.0
           endif
	   c(i,5) = c(i,5) / u(i)
30	continue
	b(1,1) = c(4,4) - c(4,1)
	b(1,2) = c(4,4) + c(4,1)
	b(2,1) = c(4,3) - c(4,2)
	b(2,2) = c(4,3) + c(4,2)
	b(3,1) = c(3,4) - c(3,1)
	b(3,2) = c(3,4) + c(3,1)
	b(4,1) = c(3,3) - c(3,2)
	b(4,2) = c(3,3) + c(3,2)
	b(1,3) = c(4,5) - c(1,5)
	b(2,3) = c(3,5) - c(2,5)
	b(3,3) = c(3,5) + c(2,5)
	b(4,3) = c(4,5) + c(1,5)
	return
	end

c **********************************************************************
c coefficient calculations for second order differential equations.
c **********************************************************************
	subroutine coeff2
	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
	common /coedf1/ b(4,3)
	common /coedf2/ a(2,2,2), d(4)
	fw1 = b(1,1) * b(1,2)
	fw2 = b(2,1) * b(3,2)
	fw3 = b(3,1) * b(2,2)
	fw4 = b(4,1) * b(4,2)
	a(2,2,1) = fw1 + fw2
	a(2,1,1) = b(1,1) * b(2,2) + b(2,1) * b(4,2)
	a(1,2,1) = b(3,1) * b(1,2) + b(4,1) * b(3,2)
	a(1,1,1) = fw3 + fw4
	a(2,2,2) = fw1 + fw3
	a(2,1,2) = b(1,2) * b(2,1) + b(2,2) * b(4,1)
	a(1,2,2) = b(3,2) * b(1,1) + b(4,2) * b(3,1)
	a(1,1,2) = fw2 + fw4
	d(1) = b(3,2) * b(4,3) + b(4,2) * b(3,3) + b(2,3) / u0
	d(2) = b(1,2) * b(4,3) + b(2,2) * b(3,3) + b(1,3) / u0
	d(3) = b(3,1) * b(1,3) + b(4,1) * b(2,3) + b(3,3) / u0
	d(4) = b(1,1) * b(1,3) + b(2,1) * b(2,3) + b(4,3) / u0
	return
	end

c **********************************************************************
c coefficient calculations for fourth-order differential equations.
c **********************************************************************
	subroutine coeff4  
	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
	common /coedf2/ a(2,2,2), d(4)
	common /coedf4/ b1, c1, z(4)
	x = u0 * u0
	b1 = a(2,2,1) + a(1,1,1)
	c1 = a(2,1,1) * a(1,2,1) - a(1,1,1) * a(2,2,1)
	z(1) = a(2,1,1) * d(3) + d(4) / x - a(1,1,1) * d(4)
	z(2) = a(1,2,1) * d(4) - a(2,2,1) *d(3) + d(3) / x
	z(3) = a(2,1,2) * d(1) + d(2) / x - a(1,1,2) * d(2)
	z(4) = a(1,2,2) * d(2) - a(2,2,2) * d(1) + d(1) / x
	return
	end

c **********************************************************************
c fk1 and fk2 are the eigenvalues.
c **********************************************************************
	subroutine coeffl  
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
	common /coedf1/ b(4,3)
	common /coedf2/ a(2,2,2), d(4)
	common /coedf4/ b1, c1, z(4)
	common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2
	dt = t1 - t0
	x = sqrt ( b1 * b1 + 4.0 * c1 )
	fk1 = sqrt ( ( b1 + x ) * 0.5 )
	fk2 = sqrt ( ( b1 - x ) * 0.5 )
	fw = u0 * u0
	x = 1.0 / ( fw * fw ) - b1 / fw - c1
	fw = 0.5 * f0 / x
	z(1) = fw * z(1) 
	z(2) = fw * z(2) 
	z(3) = fw * z(3) 
	z(4) = fw * z(4) 
	z1(1) = 0.5 * ( z(1) + z(3) )
	z1(2) = 0.5 * ( z(2) + z(4) )
	z1(3) = 0.5 * ( z(2) - z(4) )
	z1(4) = 0.5 * ( z(1) - z(3) )
	a2 = ( fk1 * fk1 - a(2,2,1) ) / a(2,1,1)
	b2 = ( fk2 * fk2 - a(2,2,1) ) / a(2,1,1)
	x = b(1,1) * b(4,1) - b(3,1) * b(2,1)
	fw1 = fk1 / x
	fw2 = fk2 / x
	y = fw2 * ( b2 * b(2,1) - b(4,1) ) 
	zx = fw1 * ( a2 * b(2,1) - b(4,1) )
	a1(1,1) = 0.5 * ( 1 - y )
   	a1(1,2) = 0.5 * ( 1 - zx )
	a1(1,3) = 0.5 * ( 1 + zx )
	a1(1,4) = 0.5 * ( 1 + y )
	y = fw2 * ( b(3,1) - b2 * b(1,1) ) 
	zx = fw1 * ( b(3,1) - a2 * b(1,1) ) 
	a1(2,1) = 0.5 * ( b2 - y )
	a1(2,2) = 0.5 * ( a2 - zx )
	a1(2,3) = 0.5 * ( a2 + zx )
	a1(2,4) = 0.5 * ( b2 + y )
        a1(3,1) = a1(2,4)
        a1(3,2) = a1(2,3)
        a1(3,3) = a1(2,2)
        a1(3,4) = a1(2,1)
        a1(4,1) = a1(1,4)
        a1(4,2) = a1(1,3)
        a1(4,3) = a1(1,2)
        a1(4,4) = a1(1,1)
	if ( ib .le. mbs ) then
	  fq0 = exp ( - t0 / u0 )
          fq1 = exp ( - t1 / u0 )
	else
	  fq0 = 1.0
	  fq1 = exp ( - dt / u0 )
	endif
	x = exp ( - fk1 * dt )
	y = exp ( - fk2 * dt )
	do 40 i = 1, 4
	   zz(i,1) = z1(i) * fq0
	   zz(i,2) = z1(i) * fq1
	   aa(i,1,1) = a1(i,1)
	   aa(i,2,1) = a1(i,2)
	   aa(i,3,1) = a1(i,3) * x
	   aa(i,4,1) = a1(i,4) * y
	   aa(i,3,2) = a1(i,3)
	   aa(i,4,2) = a1(i,4)
	   aa(i,1,2) = a1(i,1) * y
	   aa(i,2,2) = a1(i,2) * x
40	continue
	return
	end

c **********************************************************************
c See the paper by Liou, Fu and Ackerman (1988) for the formulation of
c the delta-four-stream approximation in a homogeneous layer.
c **********************************************************************
	subroutine coefft 
	call coeff1
	call coeff2
	call coeff4
	call coeffl
	return
	end

c **********************************************************************
c In the limits of no scattering ( Fu, 1991 ), fk1 = 1.0 / u(3) and
c fk2 = 1.0 / u(4).
c **********************************************************************
	subroutine coefft0 
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /point/ u(4)
	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
	common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2
	fk1 = 4.7320545
	fk2 = 1.2679491
	y = exp ( - ( t1 - t0 ) / u0 )
	fw = 0.5 * f0
	do 10 i = 1, 4
	   if ( ib .le. mbs ) then
             z1(i) = 0.0
             zz(i,1) = 0.0
             zz(i,2) = 0.0
           else
	     jj = 5 - i
	     z1(i) = fw / ( 1.0 + u(jj) / u0 )
	     zz(i,1) = z1(i) 
	     zz(i,2) = z1(i) * y
	   endif
	do 10 j = 1, 4
	   a1(i,j) = 0.0
	do 10 k = 1, 2
	   aa(i,j,k) = 0.0
10	continue
	do 20 i = 1, 4
	   j = 5 - i
	   a1(i,j) = 1.0
20	continue
	dt = t1 - t0
	x = exp ( - fk1 * dt )
	y = exp ( - fk2 * dt )
	aa(1,4,1) = y
	aa(2,3,1) = x
	aa(3,2,1) = 1.0
	aa(4,1,1) = 1.0
	aa(1,4,2) = 1.0
	aa(2,3,2) = 1.0
	aa(3,2,2) = x
	aa(4,1,2) = y
	return
	end

c **********************************************************************
c In the solar band  asbs is the surface albedo, while in the infrared
c band asbs is  blackbody intensity emitted at the surface temperature
c times surface emissivity.  In this subroutine, the delta-four-stream
c is applied to nonhomogeneous atmospheres. See comments in subroutine
c 'qcfel' for array AB(13,4*n).
c **********************************************************************
	subroutine qccfe ( ib, asbs, ee )  
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dis/ a(4)
	common /point/ u(4)
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0(100), f0(100)
	common /coedfi/ ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n
	common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4),
     1                  fk1t, fk2t
	common /qccfeo/ fk1(100), fk2(100), a4(4,4,100), 
     1	                z4(4,100), g4(4,100)
	common /qcfelc/ ab(13,100), bx(100), xx(100)
	dimension fu(4,4), wu(4)
	n = ndfs
	n4 = ndfs4
	do 333 i = 1, n4
	do 333 j = 1, 13
	   ab(j,i) = 0.0
333	continue
	ibn = ib
	wn = w(1)
	w1n = w1(1)
	w2n = w2(1)
	w3n = w3(1)
	t0n = 0.0
	t1n = t(1)
	u0n = u0(1)
	f0n = f0(1)
	if ( wn .ge. 0.999999 ) then
          wn = 0.999999
        endif
	if ( wn .le. 1.0e-4 ) then
 	  call coefft0 
	  fk1(1) = fk1t
	  fk2(1) = fk2t
        else
       	  call coefft
	  fk1(1) = fk1t
	  fk2(1) = fk2t
	endif
	do 10 i = 1, 4
	   z4(i,1) = z1(i)
	do 10 j = 1, 4
	   a4(i,j,1) = a1(i,j)
10	continue
	do 20 i = 1, 2
	   bx(i) = - zz(i+2,1)
           i8 = i + 8
	do 20 j = 1, 4
	   ab(i8-j,j) = aa(i+2,j,1)
20	continue
	do 30 i = 1, 4
	   wu(i) = zz(i,2)
	do 30 j = 1, 4
	   fu(i,j) = aa(i,j,2)
30	continue
	do 40 k = 2, n
	   wn = w(k)
	   w1n = w1(k)
	   w2n = w2(k)
	   w3n = w3(k)
	   t0n = t(k-1)
	   t1n = t(k)
	   u0n = u0(k)
	   f0n = f0(k)
	   if ( wn .ge. 0.999999 ) then
             wn = 0.999999
           endif
	   if ( wn .le. 1.0e-4 ) then
       	     call coefft0
	     fk1(k) = fk1t
	     fk2(k) = fk2t
	   else
	     call coefft 
	     fk1(k) = fk1t
	     fk2(k) = fk2t
           endif
	   do 50 i = 1, 4
	      z4(i,k) = z1(i)
	   do 50 j = 1, 4
	      a4(i,j,k) = a1(i,j)
50	   continue
	   kf = k + k + k + k
	   i1 = kf - 5
	   i2 = i1 + 3
	   j1 = kf - 7
	   j2 = j1 + 3
	   i3 = 0
	   do 55 i = i1, i2
	      i3 = i3 + 1
	      bx(i) = - wu(i3) + zz(i3,1)
	      j3 = 0
              i8 = i + 8
	      do 60 j = j1, j2
	         j3 = j3 + 1
	         ab(i8-j,j) = fu(i3,j3)
60	      continue
	      j3 = 0
	      do 65 j = j2 + 1, j2 + 4
	         j3 = j3 + 1
	         ab(i8-j,j) = - aa(i3,j3,1)
65	      continue
55	   continue
	   do 70 i = 1, 4
	      wu(i) = zz(i,2)
	   do 70 j = 1, 4
	      fu(i,j) = aa(i,j,2)
70	   continue
40	continue
	if ( ib .le. mbs ) then
	  v1 = 0.2113247 * asbs
	  v2 = 0.7886753 * asbs
	  v3 = asbs * u0(1) * f0(1) * exp ( - t(n) / u0(1) )
	  m1 = n4 - 1
	  m2 = n4
          m18 = m1 + 8
          m28 = m2 + 8
	  fw1 = v1 * wu(3)
	  fw2 = v2 * wu(4)
	  bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
	  bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
	  do 80 j = 1, 4
	     j1 = n4 - 4 + j
	     fw1 = v1 * fu(3,j)
	     fw2 = v2 * fu(4,j)
	     ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
	     ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
80	  continue
        else
	  v1 = 0.2113247 * ( 1.0 - ee )
	  v2 = 0.7886753 * ( 1.0 - ee )
	  v3 = asbs
	  m1 = n4 - 1
	  m2 = n4
          m18 = m1 + 8
          m28 = m2 + 8
	  fw1 = v1 * wu(3)
	  fw2 = v2 * wu(4)
	  bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
	  bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
	  do 85 j = 1, 4
	     j1 = n4 - 4 + j
	     fw1 = v1 * fu(3,j)
	     fw2 = v2 * fu(4,j)
	     ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
	     ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
85	  continue
	endif
	call qcfel
	do 90 k = 1, n
	   j = k + k + k + k - 4
	do 90 i = 1, 4
	   j = j + 1
	   g4(i,k) = xx(j)
90	continue
	return
	end

c **********************************************************************
	subroutine qcfel 
c **********************************************************************
c 1. `qcfel' is the abbreviation of ` qiu constants for each layer'.
c 2. The inhomogeneous atmosphere is divided into n adjacent homogeneous
c    layers where the  single scattering properties are constant in each
c    layer and allowed to vary from one to another. Delta-four-stream is
c    employed for each homogeneous layer. The boundary conditions at the
c    top and bottom of the atmosphere,  together with  continuity condi-
c    tions  at  layer interfaces lead to a system of algebraic equations
c    from which 4*n unknown constants in the problom can be solved.
c 3. This subroutine is used for solving the 4*n unknowns of A *X = B by
c    considering the fact that the coefficient matrix is a sparse matrix
c    with the precise pattern in this special problom.
c 4. The method is not different in principle from the general scheme of
c    Gaussian elimination with backsubstitution, but carefully optimized
c    so as to minimize arithmetic operations.  Partial  pivoting is used
c    to quarantee  method's numerical stability,  which will  not change
c    the basic pattern of sparsity of the matrix.
c 5. Scaling special problems so as to make  its nonzero matrix elements
c    have comparable magnitudes, which will ameliorate the stability.
c 6. a, b and x present A, B and X in A*X=B, respectively. and n4=4*n.
c 7. AB(13,4*n) is the matrix A in band storage, in rows 3 to 13; rows 1
c    and 2 and other unset elements should be set to zero on entry.
c 8. The jth column of A is stored in the jth column of the array AB  as
c    follows:
c            AB(8+i-j,j) = A(i,j) for max(1,j-5) <= i <= min(4*n,j+5).
c    Reversedly, we have
c            A(ii+jj-8,jj) = AB(ii,jj).
c **********************************************************************
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /qcfelc/ ab(13,100), b(100), x(100)
	n = ndfs
	n4 = ndfs4
	do 5 k = 1, n - 1
           k44 = 4 * k - 4
	   do 3 l= 1, 4
	      m1 = k44 + l
	      p = 0.0
	      do 10 i = 8, 14 - l
	         if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
	           p = ab(i,m1)
	           i0 = i
                 endif
10	      continue
              i0m1 = i0 + m1
              m18 = m1 + 8
	      if ( i0 .eq. 8 ) goto 20
	      do 15 j = m1, m1 + 8 - l
                 i0f = i0m1 - j
                 m1f = m18 - j
	         t = ab(i0f,j)
	         ab(i0f,j) = ab(m1f,j)
  	         ab(m1f,j) = t
15            continue
              i0f = i0m1 - 8
	      t = b(i0f)
	      b(i0f) = b(m1)
	      b(m1) = t
20	      continue
	      yy = ab(8,m1)
	      ab(8,m1) = 1.0
	      do 25 j = m1 + 1, m1 + 8 - l
                 m1f = m18 - j
	         ab(m1f,j) = ab(m1f,j) / yy
25	      continue
	      b(m1) = b(m1) / yy
	      do 30 i = 9, 14 - l
	         xx = ab(i,m1)
                 ab(i,m1) = 0.0
                 im1 = i + m1
	         do 35 j = m1 + 1, m1 + 8 - l
                    ifq = im1 - j
                    m1f = m18 - j
	            ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
35	         continue
	         ifq = im1 - 8
	         b(ifq) = b(ifq) - b(m1) * xx
30	      continue
3	   continue
5	continue
	n44 = n4 - 4
	do 40 l = 1, 3
	   m1 = n44 + l
	   p = 0.0
	   do 45 i = 8, 12 - l
	      if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
	        p = ab(i,m1)
	        i0 = i
              endif
45	   continue
           i0m1 = i0 + m1
           m18 = m1 + 8
	   if( i0 .eq. 8 ) goto 55
	   do 50 j = m1, m1 + 4 - l
              i0f = i0m1 - j
              m1f = m18 - j
	      t = ab(i0f,j)
	      ab(i0f,j) = ab(m1f,j)
  	      ab(m1f,j) = t
50	   continue
           i0f = i0m1 - 8
	   t = b(i0f)
	   b(i0f) = b(m1)
	   b(m1) = t
55	   continue
	   yy = ab(8,m1)
           ab(8,m1) = 1.0
	   do 60 j = m1 + 1, m1 + 4 - l
              m1f = m18 - j
	      ab(m1f,j) = ab(m1f,j) / yy
60	   continue
	   b(m1) = b(m1) / yy
	   do 65 i = 9, 12 - l
	      xx = ab(i,m1)
              ab(i,m1) = 0.0
              im1 = i + m1
	      do 70 j = m1 + 1, m1 + 4 - l
                 ifq = im1 - j
                 m1f = m18 - j
	         ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
70	      continue
              ifq = im1 - 8
	      b(ifq) = b(ifq) - b(m1) * xx
65	   continue
40	continue
	yy = ab(8,n4)
	ab(8,n4) = 1.0
	b(n4) = b(n4) / yy
	n3 = n4 - 1
	n2 = n3 - 1
	n1 = n2 - 1
	x(n4) = b(n4)
	x(n3) = b(n3) - ab(7,n4) * x(n4)
	x(n2) = b(n2) - ab(7,n3) * x(n3) - ab(6,n4) * x(n4)
	x(n1) = b(n1) - ab(7,n2) * x(n2) - ab(6,n3) * x(n3) -
     1	ab(5,n4) * x(n4)
	do 80 k = 1, n - 1
	   m4 = 4 * ( n - k )
	   m3 = m4 - 1
	   m2 = m3 - 1
	   m1 = m2 - 1
           m48 = m4 + 8
           m38 = m3 + 8
           m28 = m2 + 8
           m18 = m1 + 8
	   x(m4) = b(m4)
	   do 85 m = m4 + 1, m4 + 4
  	      x(m4) = x(m4) - ab(m48-m,m) * x(m)
85	   continue
	   x(m3) = b(m3)
	   do 90 m = m3 + 1, m3 + 5
  	      x(m3) = x(m3) - ab(m38-m,m) * x(m)
90	   continue
	   x(m2) = b(m2)
	   do 95 m = m2 + 1, m2 + 6
  	      x(m2) = x(m2) - ab(m28-m,m) * x(m)
95	   continue
	   x(m1) = b(m1)
	   do 100 m = m1 + 1, m1 + 7
   	      x(m1) = x(m1) - ab(m18-m,m) * x(m)
100	   continue
80	continue
	return
	end

c **********************************************************************
c In this subroutine, we incorporate a delta-function adjustment to
c account for the  forward  diffraction  peak in the context of the 
c four-stream or two stream approximations ( Liou, Fu and Ackerman,
c 1988 ).  The w1(n), w2(n), w3(n), w(n), and t(n) are the adjusted
c parameters.
c **********************************************************************
	subroutine adjust4 
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dfsin/ ww1(100), ww2(100), ww3(100), ww4(100),
     1	               ww(100), tt(100)
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0a(100), f0a(100)
	dimension dtt(ndfs), dt(ndfs)
	n = ndfs
	tt0 = 0.0
	do i = 1, n
	   f = ww4(i) / 9.0
 	   fw = 1.0 - f * ww(i) 
	   w1(i) = ( ww1(i) - 3.0 * f ) / ( 1.0 - f )
	   w2(i) = ( ww2(i) - 5.0 * f ) / ( 1.0 - f )
	   w3(i) = ( ww3(i) - 7.0 * f ) / ( 1.0 - f )
	   w(i) = ( 1.0 - f ) * ww(i) / fw
	   dtt(i) = tt(i) - tt0
	   tt0 = tt(i)
	   dt(i) = dtt(i) * fw
	enddo
	t(1) = dt(1)
	do i = 2, n
	   t(i) = dt(i) + t(i-1)
	enddo
	return
	end

c **********************************************************************
c The delta-four-stream approximation for nonhomogeneous atmospheres
c in the solar wavelengths (Fu, 1991). The input parameters are ndfs,
c mdfs, and ndfs4 through 'para.file',  ib, as, u0, f0 for solar and
c ib, bf, bs, ee for IR through arguments of  'qfts' and 'qfti', and
c ww1(ndfs), ww2(ndfs), ww3(ndfs), ww4(ndfs), ww(ndfs), and tt(ndfs)
c through common statement 'dfsin'.
c **********************************************************************
	subroutine qfts ( ib, as, u0, f0 )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dis/ a(4)
	common /point/ u(4)
	common /dfsin/ ww1(100), ww2(100), ww3(100), ww4(100),
     1	               ww(100), tt(100)
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0a(100), f0a(100)
	common /qccfeo/ fk1(100), fk2(100), a4(4,4,100), 
     1	                z4(4,100), g4(4,100)
	common /dfsout/ ffu(100), ffd(100)
	dimension x(4), fi(4)
	n = ndfs
	m = mdfs
        ee = 0.0
	asbs = as
	call adjust4
	do 5 i = 1, n
	   u0a(i) = u0
	   f0a(i) = f0
5	continue
	call qccfe ( ib, asbs, ee ) 
	fw1 = 0.6638961
	fw2 = 2.4776962
        fw3 = u0 * 3.14159 * f0 
	do 10 i = 1, m
	   if ( i .eq. 1 ) then
             x(1) = 1.0
	     x(2) = 1.0
	     x(3) = exp ( - fk1(1) * t(1) )
	     x(4) = exp ( - fk2(1) * t(1) )
             k = 1
	     y = 1.0
	   elseif ( i .eq. 2 ) then
             x(1) = exp ( - fk2(1) * t(1) )
	     x(2) = exp ( - fk1(1) * t(1) )
	     x(3) = 1.0
	     x(4) = 1.0
	     k = 1
	     y = exp ( - t(1) / u0 )
	   else
	     k = i - 1
	     y1 = t(k) - t(k-1)
             x(1) = exp ( - fk2(k) * y1 )
             x(2) = exp ( - fk1(k) * y1 )
	     x(3) = 1.0
	     x(4) = 1.0
	     y = exp ( - t(k) / u0 )
	   endif
	   do 37 jj = 1, 4
	      fi(jj) = z4(jj,k) * y
37	   continue
	   do 40 ii = 1, 4
	      fw4 = g4(ii,k) * x(ii)
	      do 45 jj = 1, 4
	         fi(jj) = fi(jj) + a4(jj,ii,k) * fw4
45            continue
40	   continue
	   ffu(i)= fw1 * fi(2) + fw2 * fi(1) 
	   ffd(i)= fw1 * fi(3) + fw2 * fi(4) + fw3 * y
10	continue
	return
	end

c **********************************************************************
c The exponential approximation for the Planck function in optical depth
c is used for the infrared ( Fu, 1991). Since the direct solar radiation
c source has an exponential function form in terms of optical depth, the
c formulation of the delta-four-stream approximation for infrared  wave-
c lengths is the same as that for solar wavelengths. 
c **********************************************************************
	subroutine qfti ( ib, ee )
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dis/ a(4)
	common /point/ u(4)
	common /dfsin/ ww1(100), ww2(100), ww3(100), ww4(100),
     1	               ww(100), tt(100)
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0(100), f0(100)
	common /qccfeo/ fk1(100), fk2(100), a4(4,4,100), 
     1	                z4(4,100), g4(4,100)
	common /dfsout/ ffu(100), ffd(100)
	common /planci/ bf(100), bs

	dimension x(4), fi(4)
	n = ndfs
	m = mdfs
	asbs = bs * ee
	call adjust4 
	t0 = 0.0
	do 3 i = 1, n
	   q1 = alog ( bf(i+1) / bf(i) )
c Paul Stackhouse: July 1, 98
c	   q2 = 1.0 / ( t(i) - t0 )
	   timt0 = t(i) - t0 
	   if (timt0 .lt. 1.0e-10) then
              timt0 = 1.0e-10
	   endif
           q2 = 1.0 / timt0
c Paul Stackhouse: July 1, 98
	   f0(i) = 2.0 * ( 1.0 - w(i) ) * bf(i)
	   if ( abs(q1) .le. 1.0e-10 ) then
	     u0(i) = - 1.0e+10 / q2
           else
	     u0(i) = - 1.0 / ( q1 * q2 )
           endif
	   t0 = t(i)
3	continue
	call qccfe ( ib, asbs, ee ) 
	fw1 = 0.6638958
	fw2 = 2.4776962
	do 10 i = 1, m
	   if ( i .eq. 1 ) then
             x(1) = 1.0
	     x(2) = 1.0
	     x(3) = exp ( - fk1(1) * t(1) )
	     x(4) = exp ( - fk2(1) * t(1) )
             k = 1
             xy = 1.0
	   elseif ( i .eq. 2 ) then
             x(1) = exp ( - fk2(1) * t(1) )
	     x(2) = exp ( - fk1(1) * t(1) )
	     x(3) = 1.0
	     x(4) = 1.0
             k = 1
	     xy =  exp ( - t(1) / u0(1) )
	   else
             k = i - 1
	     y1 = t(k) - t(k-1)
             x(1) = exp ( - fk2(k) * y1 )
             x(2) = exp ( - fk1(k) * y1 )
	     x(3) = 1.0
	     x(4) = 1.0
	     xy =  exp ( - y1 / u0(k) )
	   endif
	   do 37 jj = 1, 4
	      fi(jj) = z4(jj,k) * xy
37	   continue
           do 40 ii = 1, 4
	      fw3 = g4(ii,k) * x(ii)
	      do 45 jj = 1, 4
	         fi(jj) = fi(jj) + a4(jj,ii,k) * fw3
45	      continue
40	   continue
	   ffu(i)= fw1 * fi(2) + fw2 * fi(1)
	   ffd(i)= fw1 * fi(3) + fw2 * fi(4)
10	continue
	return
	end


c 11/4/95 (begin) 8/22/96
c 8/1/95 (begin)
	subroutine cfgts0 ( gamma1, gamma2, gamma3,gamma4,ugts1)
c **********************************************************************
c This subroutine is used to calculate the Coefficients For Generalized
c Two-Stream scheme.  We can make choices between Eddington, quadrature
c and  hemispheric  mean  schemes  through  logical variables 'edding',
c 'quadra', and 'hemisp'.  The  Eddington  and  quadrature  schemes are 
c discussed in detail by Liou (1992).  The  hemispheric  mean scheme is 
c derived by assuming that the phase function is equal to 1 + g in  the 
c forward scattering hemisphere and 1 - g  in  the  backward scattering 
c hemisphere where g is the asymmetry factor.   The hemispheric mean is
c only used for infrared wavelengths (Toon et al. 1989).
c We add the modified discrete-ordinate two-stream scheme for  infrared
c radiative transfer (Fu et al.1996). This scheme can be chosen through
c the logical variable 'mquadr'. 8/19/96.
c **********************************************************************
	include 'para.file'
	logical edding, quadra, hemisp, mquadr
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common / gtslog / edding, quadra, hemisp, mquadr
	common /coedfi / ib, w, w1, w2, w3, t0, t1, u0, f0
	if ( edding ) then
	   x = 0.25 * w1
	   y = w * x
           gamma1 = 1.75 - w - y
	   gamma2 = - 0.25 + w - y
           gamma3 = 0.0
           gamma4 = 0.0
	   if ( ib .le. mbs ) then
              gamma3 = 0.5 - x * u0
              gamma4 = 1.0 - gamma3
           endif
           ugts1 = 0.5
	endif
	if ( quadra ) then
	   x = 0.866 * w
           y = 0.2887 * w1
	   z = y * w
           gamma1 = 1.732 - x - z
	   gamma2 = x - z
           gamma3 = 0.0
           gamma4 = 0.0
	   if ( ib .le. mbs ) then
              gamma3 = 0.5 - y * u0
              gamma4 = 1.0 - gamma3
           endif
           ugts1 = 0.57735
	endif
	if ( hemisp ) then
	   x = w * w1 / 3.0
           gamma1 = 2.0 - w - x
           gamma2 = w - x
           gamma3 = 0.0
           gamma4 = 0.0
           ugts1 = 0.5
        endif
	if ( mquadr ) then
           y = 0.2767 * w
	   x = y + y + y
	   z = y * w1
           gamma1 = 1.66 - x - z
	   gamma2 = x - z
           gamma3 = 0.0
           gamma4 = 0.0
           ugts1 = 0.6024
	endif
	return
	end


	subroutine cfgts ( lamda, gamma, cadd0, cadd1, cmin0, cmin1,
     1                     g1g2, fkb )
c **********************************************************************
c This subroutine is used to calculate the Coefficients For Generalized
c Two-Stream scheme. 
c **********************************************************************
	include 'para.file'
	real lamda
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
	call cfgts0 ( gamma1, gamma2, gamma3, gamma4, ugts1)
	lamda = sqrt ( ( gamma1 + gamma2 ) * ( gamma1 - gamma2 ) )
	gamma = gamma2 / ( gamma1 + lamda )
	g1g2 = gamma1 + gamma2
        fq = 1.0 / u0
        x = exp ( - fq * ( t1 - t0 ) )
	z = lamda * lamda - fq * fq
	fkb = z
	if ( ib .le. mbs ) then
           alfa = gamma3
           beta = gamma4
           fw = 3.1415927 * f0 * w * exp ( - fq * t0 )
           cadd0 = fw * ( ( gamma1 - fq ) * alfa +
     1             beta * gamma2 ) / z
	   cmin0 = fw * ( ( gamma1 + fq ) * beta +
     1             alfa * gamma2 ) / z
        else
           fw = 3.1415927 * f0
           cadd0 = fw * ( g1g2 - fq ) / z
	   cmin0 = fw * ( g1g2 + fq ) / z
	endif
	cadd1 = cadd0 * x
	cmin1 = cmin0 * x
	return
	end


	subroutine qccgts ( ib, asbs, ee )
c **********************************************************************
c In the solar band  asbs is the surface albedo, while in the infrared
c band asbs is  blackbody intensity emitted at the surface temperature
c times surface emissivity.  In this subroutine,  the generalized two-
c stream is applied to nonhomogeneous atmospheres. ee is the IR surface
c emissivity. 
c **********************************************************************
	include 'para.file'
	real lamdan
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0(100), f0(100)
	common /coedfi/ ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n
	common / gtscoe / lamdan(100), gamman(100), caddn(100),
     1                    cminn(100), caddn0(100), cminn0(100),
     1                    aa(100), bb(100), expn(100), g1g2n(100),
     1                    fkbn(100)

	dimension a(ndfs2), b(ndfs2), c(ndfs2), r(ndfs2), u(ndfs2),
     1            gam(ndfs2)
	dimension xn(ndfs), yn(ndfs), zn(ndfs)
	ibn = ib
	do 40 k = 1, ndfs
           wn = w(k)
           w1n = w1(k)
           if ( k .eq. 1 ) then
              t0n = 0.0
           else
              t0n = t(k-1)
           endif
           t1n = t(k)
           u0n = u0(k)
           f0n = f0(k)
           if ( wn .ge. 0.999999 ) then
              wn = 0.999999
           endif
           call cfgts ( lamdan(k), gamman(k), caddn0(k), caddn(k),
     1      cminn0(k),cminn(k),g1g2n(k),fkbn(k) )
           expn(k) = exp ( - lamdan(k) * ( t1n - t0n ) )
	   xn(k) = gamman(k) * expn(k)
	   yn(k) = ( expn(k) - gamman(k) ) / ( xn(k) - 1.0 )
	   zn(k) = ( expn(k) + gamman(k) ) / ( xn(k) + 1.0 )
40	continue
	a(1) = 0.0
        b(1) = xn(1) + 1.0
        c(1) = xn(1) - 1.0
        r(1) = - cminn0(1)
	do 50 k = 1, ndfs - 1
        k1 = k + k
        k2 = k + k + 1
        a(k1) = 1.0 + xn(k) - yn(k+1) * ( gamman(k) + expn(k) )
        b(k1) = 1.0 - xn(k) - yn(k+1) * ( gamman(k) - expn(k) )
	c(k1) = yn(k+1) * ( 1.0 + xn(k+1) ) - expn(k+1) - gamman(k+1)
	r(k1) = caddn0(k+1) - caddn(k) - yn(k+1) *
     1          ( cminn0(k+1) - cminn(k) )
	a(k2) = gamman(k) - expn(k) - zn(k) * ( 1.0 - xn(k) )
	b(k2) = -1.0 - xn(k+1) + zn(k) * ( expn(k+1) + gamman(k+1) )
        c(k2) = zn(k) * ( expn(k+1) - gamman(k+1) ) - xn(k+1) + 1.0
        r(k2) = cminn0(k+1) - cminn(k) - zn(k) *
     1          ( caddn0(k+1) - caddn(k) )
50	continue
	if ( ib .le. mbs ) then
           rsfc = asbs
           ssfc = 3.1415927 * u0(1) * exp(-t(ndfs)/u0(1)) * rsfc *
     1            f0(1)
        else
           rsfc = 1.0 - ee
           ssfc = 3.1415927 * asbs
	endif
	wm1 = 1.0 - rsfc * gamman(ndfs)
	wm2 = xn(ndfs) - rsfc * expn(ndfs)
	a(ndfs2) = wm1 + wm2
	b(ndfs2) = wm1 - wm2
	c(ndfs2) = 0.0
	r(ndfs2) = rsfc * cminn(ndfs) - caddn(ndfs) + ssfc
	call tridag (a,b,c,r,u,gam,ndfs2 )
	do 60 k = 1, ndfs
        k1 = k + k - 1
        k2 = k + k
	aa(k) = u(k1) + u(k2)
	bb(k) = u(k1) - u(k2)
60	continue
	return
	end


	subroutine tridag(a,b,c,r,u,gam,n )
c **********************************************************************
c
c   | b1 c1 0  ...                |   | u1   |   | r1   |                    
c   | a2 b2 c2 ...                |   | u2   |   | r2   |                
c   |          ...                | . | .    | = | .    |
c   |          ... an-1 bn-1 cn-1 |   | un-1 |   | rn-1 |                    
c   |              0    an   bn   |   | un   |   | rn   |                
c
c This  subroutine solves for  a vector U of length N the tridiagonal
c linear set given by above equation. A, B, C and R are input vectors
c and are not modified (Numerical Recipes by Press et al. 1989).
c **********************************************************************
	dimension gam(n), a(n), b(n), c(n), r(n), u(n)
c	if ( b(1) .eq. 0. ) pause
c If this happens then you should rewrite your equations as a set of
c order n-1, with u2 trivially eliminated.
	bet = b(1)
	u(1) = r(1) / bet
c Decomposition and forward substitution
	do 11 j = 2, n
           gam(j) = c(j-1) / bet
	   bet = b(j) - a(j) * gam(j)
c           if ( bet .eq. 0. ) pause
c Algorithm fails; see Numerical Recipes.
	   u(j) = ( r(j) - a(j) * u(j-1) ) / bet
11	continue
c Backsubstitution
	do 12 j = n - 1, 1, -1
           u(j) = u(j) - gam(j+1) * u(j+1)
12	continue
	return
	end


	subroutine qftsts (ib,as,u0,f0 )
c **********************************************************************
c The generalized two stream approximation for nonhomgeneous atmospheres
c in  the  solar  wavelengths.  The  input  parameters are those through
c 'para.file', through argument of 'qftsts' and through common statement
c 'dfsin' and 'gtslog'.
c **********************************************************************
	include 'para.file'
	real lamdan
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dfsin/ ww1(100), ww2(100), ww3(100), ww4(100),
     1	               ww(100), tt(100)
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0a(100), f0a(100)
	common / gtscoe / lamdan(100), gamman(100), caddn(100),
     1                    cminn(100), caddn0(100), cminn0(100),
     1                    aa(100), bb(100), expn(100), g1g2n(100),
     1                    fkbn(100)
	common /dfsout/ ffu(100), ffd(100)

	dimension yy(ndfs)
	n = ndfs
	m = mdfs
	ee = 0.0
	asbs = as
	call adjust2
	do 5 i = 1, n
           u0a(i) = u0
           f0a(i) = f0
5	continue
	call qccgts ( ib, asbs, ee )
	fw3 = u0 * 3.1415927 * f0
	do k = 1, ndfs
           yy(k) = exp(-t(k)/u0)
	enddo
	xx = aa(1) * expn(1)
	ffu(1) = xx + gamman(1) * bb(1) + caddn0(1)
	ffd(1) = gamman(1) * xx + bb(1) + cminn0(1) + fw3
	do 10 i = 2, m
           k = i - 1
	   xx = bb(k) * expn(k)
	   ffu(i) = aa(k) + gamman(k) * xx + caddn(k)
	   ffd(i) = gamman(k) * aa(k) + xx + cminn(k) + fw3 * yy(k)
10	continue
	return
	end


	subroutine qftits ( ib, ee )
c **********************************************************************
c The exponential approximation for the Planck function in optical depth
c is used for the infrared ( Fu, 1991). Since the direct solar radiation
c source has an exponential function form in terms of optical depth, the
c formulation of generalized two stream approximation for infrared  wave
c lengths is the same as that for solar wavelengths. 
c The generalized two stream approximation for nonhomgeneous atmospheres
c in the infrared wavelengths.  The  input  parameters are those through
c 'para.file', through argument of 'qftits' and through common statement
c 'dfsin', 'gtslog', and 'planci'.
c **********************************************************************
	include 'para.file'
	real lamdan
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dfsin/ ww1(100), ww2(100), ww3(100), ww4(100),
     1	               ww(100), tt(100)
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0(100), f0(100)
	common / gtscoe / lamdan(100), gamman(100), caddn(100),
     1                    cminn(100), caddn0(100), cminn0(100),
     1                    aa(100), bb(100), expn(100), g1g2n(100),
     1                    fkbn(100)
	common /dfsout/ ffu(100), ffd(100)
	common /planci/ bf(100), bs
	logical edding, quadra, hemisp, mquadr
	common / gtslog / edding, quadra, hemisp, mquadr
	n = ndfs
	m = mdfs
	asbs = bs * ee
	call adjust2
	t0 = 0.0
	do 3 i = 1, n
	   q1 = alog ( bf(i+1) / bf(i) )
c Paul Stackhouse: July 1, 98
c	   q2 = 1.0 / ( t(i) - t0 )
	   timt0 = t(i) - t0 
	   if (timt0 .lt. 1.0e-10) then
              timt0 = 1.0e-10
	   endif
           q2 = 1.0 / timt0
c Paul Stackhouse: July 1, 98
	   if ( mquadr ) then
	      f0(i) = 1.66 * ( 1.0 - w(i) ) * bf(i)
	   else
	      f0(i) = 2.0 * ( 1.0 - w(i) ) * bf(i)
	   endif
	   if ( abs(q1) .le. 1.0e-10 ) then
	     u0(i) = - 1.0e+10 / q2
           else
	     u0(i) = - 1.0 / ( q1 * q2 )
           endif
	   t0 = t(i)
3	continue
	call qccgts ( ib, asbs, ee ) 
	xx = aa(1) * expn(1)
	ffu(1) = xx + gamman(1) * bb(1) + caddn0(1)
	ffd(1) = gamman(1) * xx + bb(1) + cminn0(1) 
	do 10 i = 2, m
           k = i - 1
	   xx = bb(k) * expn(k)
	   ffu(i) = aa(k) + gamman(k) * xx + caddn(k)
	   ffd(i) = gamman(k) * aa(k) + xx + cminn(k)
10	continue
	return
	end


c 6-24-98 (8)	
c	subroutine qftisf ( ib, ee)
	subroutine qftisf ( ib, ee, ur )
c 6-24-98 (8)	
c **********************************************************************
c In this subroutine, the two- and four- stream combination  scheme  
c (Fu et al. 1997) or the source function technique (Toon et al. 1989)
c is used to calculate the IR radiative fluxes. The exponential approxi-
c mation for the Planck c function in optical depth is used ( Fu, 1991).
c At IR wavelengths, the two-stream results are not exact in the limit 
c of no scattering. It also introduces large error in the case of sca-
c ttering. Since the no-scattering limit is of considerable significance
c at IR wavelengths, we have used  the source function technique  that
c would be exact in the limit of the pure absorption and would also en-
c hance the accuracy of the two-stream approach when scattering occurs
c in the IR wavelengths.
c Here, we use nq Gauss points to obtain the fluxes: when nq=2, we use
c double Gaussian quadrature as in Fu and Liou (1993) for  four-stream
c approximation.
c The two- and four- stream combination  technique  is  only  used for
c "hemisp" and "mquadr" in the infrared. 8/20/96
c The upward radiance at the cosine of zenith angle, ur, in the IR
c atmospheric window 800-980 cm**-1, 980-1100 cm**-1, 1100-1250 cm**-1
c is calculated based on equation (2.25) in Fu et al. (1997). 6-24-98.
c
c **********************************************************************
	include 'para.file'
	parameter ( nq = 2 )
	real lamdan
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dfsin/ ww1(100), ww2(100), ww3(100), ww4(100),
     1	               ww(100), tt(100)
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0(100), f0(100)
	common / gtscoe / lamdan(100), gamman(100), caddn(100),
     1                    cminn(100), caddn0(100), cminn0(100),
     1                    aa(100), bb(100), expn(100), g1g2n(100),
     1                    fkbn(100)
	common /dfsout/ ffu(100), ffd(100)
c 6-24-98 (9)	
  	common /radiance/ fiur(100)  
c 6-24-98 (9)	
	common /planci/ bf(100), bs
	logical edding, quadra, hemisp, mquadr
	common / gtslog / edding, quadra, hemisp, mquadr

	dimension fg(ndfs), fh(ndfs), fj(ndfs), fk(ndfs)
	dimension alfa(ndfs+1), beta(ndfs)
	dimension fiu(mdfs,nq), fid(mdfs,nq)
	dimension fx(ndfs,nq), fy(ndfs), fz1(ndfs,nq), fz2(ndfs,nq),
     1            ub(ndfs,nq)
c 6-24-98 (10)	
	dimension fxr(ndfs), fz1r(ndfs), fz2r(ndfs), ubr(ndfs)
c 6-24-98 (10)	
	dimension fuq1(ndfs), fuq2(ndfs)
	dimension ug(nq), wg(nq), ugwg(nq)
 	data ug / 0.2113248, 0.7886752 /
 	data wg / 0.5, 0.5 /
 	data ugwg / 0.105662, 0.394338 /
	if ( mquadr ) then
	   ugts1 = 0.60241
	else
           ugts1 = 0.5
	endif
	n = ndfs
	m = mdfs
	asbs = bs * ee
	call adjust2
	t0 = 0.0
	do i = 1, n
	   q1 = alog ( bf(i+1) / bf(i) )
c Paul Stackhouse: July 1, 98
c	   q2 = 1.0 / ( t(i) - t0 )
	   timt0 = t(i) - t0 
	   if (timt0 .lt. 1.0e-10) then
              timt0 = 1.0e-10
	   endif
           q2 = 1.0 / timt0
c Paul Stackhouse: July 1, 98
	   if ( mquadr ) then
	      f0(i) = 1.66 * ( 1.0 - w(i) ) * bf(i)
	   else
	      f0(i) = 2.0 * ( 1.0 - w(i) ) * bf(i)
	   endif
	   if ( abs(q1) .le. 1.0e-10 ) then
	     u0(i) = - 1.0e+10 / q2
           else
	     u0(i) = - 1.0 / ( q1 * q2 )
           endif
	   t0 = t(i)
	   beta(i) = - 1.0 / u0(i)
	enddo
	call qccgts ( ib, asbs, ee ) 
	do i = 1, n
	   x = ( 1.0 - w(i) ) * w(i) / fkbn(i) / ugts1
	   y1 = w1(i) / 3.0
	   y = g1g2n(i)
	   z = -y1 * beta(i)
	   fuq1(i) = x * ( y - z ) + 1.0 - w(i)
	   fuq2(i) = x * ( y + z ) + 1.0 - w(i)
	enddo
	do i = 1, n + 1
	   alfa(i) = 6.2832 * bf(i)
	enddo
	do i = 1, n
	   x = lamdan(i) * ugts1
           y = gamman(i) * ( 1.0 + x )
           z = 1.0 - x
           fg(i) = aa(i) * ( z + z )
	   fh(i) = bb(i) * ( y + y )
	   fj(i) = aa(i) * ( y + y )
	   fk(i) = bb(i) * ( z + z )
	enddo
	do j = 1, nq
	   fid(1,j) = 0.0
	enddo
	do j = 1, nq
	t0 = 0.0
	do i = 2, mdfs
	   i1 = i - 1
	   fx(i1,j) = exp ( - ( t(i1) - t0 ) / ug(j) )
	   fy(i1) = expn(i1)
	   xx = lamdan(i1) * ug(j)
           fz1(i1,j) = ( 1.0 - fx(i1,j) * fy(i1) ) / ( xx + 1.0 )
c Stackhouse and Fu: July 16, 2002
           cri = abs ( xx - 1.0 )
           if ( cri .le. 1.0e-10 ) then
              fz2(i1,j) = fx(i1,j) * ( t(i1) - t0 ) / ug(j)
           else
              fz2(i1,j) = ( fx(i1,j) - fy(i1) ) / ( xx - 1.0 )
           endif
c Stackhouse and Fu: July 16, 2002
	   ub(i1,j) = ug(j) * beta(i1)
	   fid(i,j) = fid(i1,j) * fx(i1,j) + fj(i1) * fz1(i1,j) +
     1                fk(i1) * fz2(i1,j) + 
     1                fuq2(i1) / ( ub(i1,j) + 1.0 ) *
     1                ( alfa(i) - alfa(i1) * fx(i1,j) ) 
	   t0 = t(i1)
	enddo
	enddo
	yy = 0.0
	do j = 1, nq
	   yy = yy + ugwg(j) * fid(mdfs,j) 
	enddo
        xx = yy * ( 1.0 - ee ) * 2.0 + 6.2831854 * ee * bs 
	do j = 1, nq
	   fiu(mdfs,j) = xx
	enddo
c 6-24-98 (11)	
	fiur(mdfs) = xx
c 6-24-98 (11)	
	do j = 1, nq
	do i = mdfs - 1, 1, -1
           fiu(i,j) = fiu(i+1,j) * fx(i,j) + fg(i) * fz2(i,j) +
     1                fh(i) * fz1(i,j) + 
     1                fuq1(i) / ( ub(i,j) - 1.0 ) *
     1                ( alfa(i+1) * fx(i,j) - alfa(i) ) 
	enddo
	enddo
	do i = 1, mdfs
	   ffu(i) = 0.0
	   ffd(i) = 0.0
	enddo
	do i = 1, mdfs
	   do j = 1, nq
	        ffu(i) = ffu(i) + ugwg(j) * fiu(i,j) 
	        ffd(i) = ffd(i) + ugwg(j) * fid(i,j)
	   enddo
	enddo
c 6-24-98 (12) 
	if ( ib .eq. 11 .or. ib .eq. 12 .or. ib .eq. 13 ) then
	if ( ur .eq. 0.5 ) then
           ur = 0.500001
	endif
	t0 = 0.0
	do i = 2, mdfs
	   i1 = i - 1
	   fxr(i1) = exp ( - ( t(i1) - t0 ) / ur )
	   xx = lamdan(i1) * ur
           fz1r(i1) = ( 1.0 - fxr(i1) * fy(i1) ) / ( xx + 1.0 )
           fz2r(i1) = ( fxr(i1) - fy(i1) ) / ( xx - 1.0 )
	   ubr(i1) = ur * beta(i1)
	   t0 = t(i1)
	enddo
	do i = mdfs - 1, 1, -1
           fiur(i) = fiur(i+1) * fxr(i) + fg(i) * fz2r(i) +
     1               fh(i) * fz1r(i) + 
     1               fuq1(i) / ( ubr(i) - 1.0 ) *
     1               ( alfa(i+1) * fxr(i) - alfa(i) ) 
	enddo
	else
	do i = 1, mdfs
           fiur(i) = 0.0
	enddo
	endif
c 6-24-98 (12) 
	return
	end
c 8/1/95 (end)


c **********************************************************************
c In this subroutine, we incorporate a delta-function adjustment to
c account for the  forward  diffraction  peak in the context of the 
c two-stream approximation ( Liou, Fu and Ackerman, 1988 ).  w1(n),
c w(n), and t(n) are the adjusted parameters.
c **********************************************************************
	subroutine adjust2 
	include 'para.file'
        common /zdim/ nv1,nv,ndfs,mdfs,ndfs4,ndfs2
	common /dfsin/ ww1(100), ww2(100), ww3(100), ww4(100),
     1	               ww(100), tt(100)
	common /qccfei/ w1(100), w2(100), w3(100), w(100), 
     1	                t(100), u0a(100), f0a(100)

	dimension dtt(ndfs), dt(ndfs)
	n = ndfs
	tt0 = 0.0
	do 10 i = 1, n
	   f = ww2(i) / 5.0
	   fw = 1.0 - f * ww(i) 
	   w1(i) = ( ww1(i) - 3.0 * f ) / ( 1.0 - f )
	   w(i) = ( 1.0 - f ) * ww(i) / fw
	   dtt(i) = tt(i) - tt0
	   tt0 = tt(i)
	   dt(i) = dtt(i) * fw
10	continue
	t(1) = dt(1)
	do 20 i = 2, n
	   t(i) = dt(i) + t(i-1)
20	continue
	return
	end
c 11/4/95 (end) 8/22/96


	block data ckd1
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 to 
c one. fko3 is the corresponding ozone absorption coefficient in units
c of (cm-atm)**-1 (Fu, 1991). The spectral region is from 50000 cm**-1 
c to 14500 cm**-1.
c *********************************************************************
	common /band1/ hk(10), fko3(10)
	data hk / .24, .16, .24, .28, .03, 
     1            .016, .01, .008, .008, .008 /
	data fko3 / .2204e-08,.1207e-01,.4537e-01,.1032e+00,.1740e+00,
     1	            .1210e+01,.7367e+01,.2050e+02,.8100e+02,.2410e+03 /
	end

	block data ckd2
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
c pressures,  and eight cumulative probabilities  ( Fu,  1991 ). The
c spectral region is from 14500 to 7700 cm**-1.
c *********************************************************************
	common /band2/ hk(8), coeh2o(3,11,8)
	data hk / .71, .11, .06, .06, .04, .016, .0034, .0006 /
c   .343849E+03    .532724E+02    .290577E+02    .290577E+02    .193718E+02
c   .774872E+01    .164660E+01    .290577E+00
	data ( ( ( coeh2o(k,j,i), i = 1, 8 ), j = 1, 11 ), k = 1, 3 ) /
     +-.1735E+02,-.1407E+02,-.1268E+02,-.1131E+02,-.9261E+01,-.6666E+01,
     +-.3937E+01,-.5448E+00,-.1690E+02,-.1365E+02,-.1232E+02,-.1101E+02,
     +-.9058E+01,-.6574E+01,-.3914E+01,-.5529E+00,-.1643E+02,-.1323E+02,
     +-.1195E+02,-.1068E+02,-.8840E+01,-.6475E+01,-.3889E+01,-.6143E+00,
     +-.1598E+02,-.1282E+02,-.1157E+02,-.1035E+02,-.8598E+01,-.6339E+01,
     +-.3848E+01,-.6636E+00,-.1551E+02,-.1241E+02,-.1119E+02,-.1001E+02,
     +-.8342E+01,-.6178E+01,-.3788E+01,-.8181E+00,-.1506E+02,-.1201E+02,
     +-.1082E+02,-.9692E+01,-.8073E+01,-.6017E+01,-.3703E+01,-.9003E+00,
     +-.1446E+02,-.1154E+02,-.1042E+02,-.9332E+01,-.7810E+01,-.5846E+01,
     +-.3576E+01,-.1083E+01,-.1394E+02,-.1112E+02,-.1005E+02,-.8992E+01,
     +-.7548E+01,-.5674E+01,-.3477E+01,-.1266E+01,-.1351E+02,-.1076E+02,
     +-.9722E+01,-.8702E+01,-.7334E+01,-.5531E+01,-.3401E+01,-.1524E+01,
     +-.1311E+02,-.1044E+02,-.9422E+01,-.8423E+01,-.7117E+01,-.5383E+01,
     +-.3410E+01,-.1785E+01,-.1274E+02,-.1015E+02,-.9162E+01,-.8190E+01,
     +-.6949E+01,-.5236E+01,-.3477E+01,-.2082E+01, .2407E-02, .2847E-02,
     + .3768E-02, .4626E-02, .5631E-02, .4542E-02, .3475E-02,-.3085E-02,
     + .2428E-02, .2805E-02, .3412E-02, .3893E-02, .4773E-02, .3998E-02,
     + .2742E-02,-.2556E-02, .2428E-02, .2721E-02, .3077E-02, .3161E-02,
     + .4019E-02, .3224E-02, .2512E-02,-.1884E-02, .2449E-02, .2617E-02,
     + .2763E-02, .2658E-02, .3286E-02, .2617E-02, .1989E-02,-.1740E-02,
     + .2512E-02, .2470E-02, .2470E-02, .2282E-02, .2512E-02, .1926E-02,
     + .1465E-02,-.2612E-02, .2554E-02, .2303E-02, .2303E-02, .1842E-02,
     + .2030E-02, .1340E-02, .1068E-02,-.1413E-02, .2449E-02, .2198E-02,
     + .2030E-02, .1465E-02, .1528E-02, .9838E-03, .1005E-02,-.1099E-02,
     + .2868E-02, .2198E-02, .1968E-02, .1382E-02, .1172E-02, .5652E-03,
     + .6070E-03,-.1662E-02, .3077E-02, .2219E-02, .1800E-02, .1277E-02,
     + .1005E-02, .3349E-03, .2512E-03,-.1195E-02, .3182E-02, .2219E-02,
     + .1758E-02, .1172E-02, .7326E-03, .4815E-03, .6280E-04,-.1880E-02,
     + .3265E-02, .2114E-02, .1696E-02, .1298E-02, .4187E-03, .4187E-03,
     +-.3768E-03,-.1467E-02,-.1180E-04,-.1294E-04,-.1142E-04,-.7232E-05,
     +-.8754E-05,-.1484E-04,-.8373E-05, .1028E-04,-.1218E-04,-.1142E-04,
     +-.9515E-05,-.1522E-05,-.9134E-05,-.1484E-04,-.3425E-05, .1142E-06,
     +-.1294E-04,-.9895E-05,-.7231E-05,-.4187E-05,-.7612E-05,-.3806E-05,
     + .1522E-05,-.3882E-05,-.1256E-04,-.8754E-05,-.7612E-05,-.6470E-05,
     +-.4948E-05,-.3425E-05, .4948E-05,-.1054E-04,-.1370E-04,-.6089E-05,
     +-.8373E-05,-.5709E-05,-.3045E-05,-.3806E-05, .5328E-05, .8678E-05,
     +-.1370E-04,-.6851E-05,-.8373E-05,-.1522E-05,-.3425E-05, .0000E+00,
     + .1256E-04,-.1572E-04,-.1484E-04,-.7231E-05,-.7992E-05,-.4567E-05,
     +-.2664E-05,-.3807E-06,-.1522E-05, .2169E-05,-.1713E-04,-.9515E-05,
     +-.6089E-05,-.6851E-05,-.3045E-05,-.1142E-05, .1903E-05, .9363E-05,
     +-.1560E-04,-.9134E-05,-.5328E-05,-.4948E-05, .0000E+00, .7611E-06,
     +-.6851E-05, .1252E-04,-.1522E-04,-.8373E-05,-.6089E-05,-.6089E-05,
     +-.3805E-06,-.1142E-05,-.3807E-06, .2512E-05,-.1599E-04,-.7231E-05,
     +-.5709E-05,-.4567E-05, .1522E-05,-.2284E-05,-.3941E-10, .5290E-05/
	end
      
	block data ckd3
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
c pressures,  and twelve cumulative probabilities ( Fu,  1991 ). The
c spectral region is from 7700 to 5250 cm**-1.
c *********************************************************************
	common /band3/ hk(12), coeh2o(3,11,12)
	data hk / .34, .11, .1, .09, .12, .1,
     1            .06, .04, .026, .01, .0035, .0005 /
c   .509474E+02    .164830E+02    .149845E+02    .134861E+02    .179814E+02
c   .149845E+02    .899071E+01    .599381E+01    .389597E+01    .149845E+01
c   .524458E+00    .749226E-01
	data ( ( ( coeh2o(k,j,i), i = 1, 12 ), j = 1, 11 ), k = 1, 3 ) /
     +-.1900E+02,-.1515E+02,-.1344E+02,-.1224E+02,-.1081E+02,-.9337E+01,
     +-.7965E+01,-.6585E+01,-.4578E+01,-.2247E+01, .1747E+00, .3083E+01,
     +-.1854E+02,-.1471E+02,-.1300E+02,-.1181E+02,-.1039E+02,-.8927E+01,
     +-.7576E+01,-.6238E+01,-.4317E+01,-.2119E+01, .1888E+00, .3033E+01,
     +-.1808E+02,-.1426E+02,-.1257E+02,-.1137E+02,-.9966E+01,-.8513E+01,
     +-.7177E+01,-.5885E+01,-.4053E+01,-.1977E+01, .2245E+00, .3005E+01,
     +-.1763E+02,-.1381E+02,-.1213E+02,-.1094E+02,-.9542E+01,-.8094E+01,
     +-.6779E+01,-.5524E+01,-.3788E+01,-.1796E+01, .2961E+00, .2828E+01,
     +-.1716E+02,-.1337E+02,-.1170E+02,-.1051E+02,-.9116E+01,-.7677E+01,
     +-.6381E+01,-.5153E+01,-.3493E+01,-.1607E+01, .3850E+00, .2660E+01,
     +-.1670E+02,-.1295E+02,-.1127E+02,-.1008E+02,-.8690E+01,-.7265E+01,
     +-.5991E+01,-.4799E+01,-.3212E+01,-.1438E+01, .4582E+00, .2588E+01,
     +-.1596E+02,-.1231E+02,-.1067E+02,-.9501E+01,-.8151E+01,-.6793E+01,
     +-.5588E+01,-.4458E+01,-.2940E+01,-.1257E+01, .4888E+00, .2260E+01,
     +-.1530E+02,-.1184E+02,-.1017E+02,-.8992E+01,-.7661E+01,-.6369E+01,
     +-.5213E+01,-.4145E+01,-.2701E+01,-.1108E+01, .4239E+00, .1974E+01,
     +-.1481E+02,-.1144E+02,-.9756E+01,-.8573E+01,-.7255E+01,-.5994E+01,
     +-.4868E+01,-.3829E+01,-.2485E+01,-.9738E+00, .3343E+00, .1667E+01,
     +-.1439E+02,-.1108E+02,-.9360E+01,-.8183E+01,-.6885E+01,-.5646E+01,
     +-.4559E+01,-.3555E+01,-.2314E+01,-.8904E+00, .2169E+00, .1289E+01,
     +-.1402E+02,-.1073E+02,-.8987E+01,-.7817E+01,-.6551E+01,-.5335E+01,
     +-.4278E+01,-.3316E+01,-.2147E+01,-.8695E+00, .1587E-01, .8658E+00,
     + .1132E-01, .8855E-02, .6698E-02, .5296E-02, .4396E-02, .3370E-02,
     + .3245E-02, .4145E-02, .4731E-02, .4756E-02, .3116E-02,-.2763E-02,
     + .1135E-01, .8917E-02, .6657E-02, .5170E-02, .4207E-02, .3056E-02,
     + .2868E-02, .3433E-02, .3726E-02, .4109E-02, .2836E-02,-.3119E-02,
     + .1135E-01, .8980E-02, .6615E-02, .5045E-02, .4061E-02, .2847E-02,
     + .2491E-02, .2847E-02, .2910E-02, .2671E-02, .2396E-02,-.3245E-02,
     + .1135E-01, .9043E-02, .6594E-02, .4940E-02, .3914E-02, .2638E-02,
     + .2156E-02, .2261E-02, .2051E-02, .1978E-02, .1566E-02,-.3203E-02,
     + .1139E-01, .9085E-02, .6531E-02, .4835E-02, .3768E-02, .2428E-02,
     + .1842E-02, .1612E-02, .1591E-02, .1279E-02, .7201E-03,-.2763E-02,
     + .1143E-01, .9085E-02, .6447E-02, .4752E-02, .3684E-02, .2261E-02,
     + .1570E-02, .1235E-02, .1151E-02, .7243E-03, .6489E-04,-.2240E-02,
     + .1135E-01, .9001E-02, .5694E-02, .4438E-02, .3412E-02, .1968E-02,
     + .1235E-02, .9420E-03, .8792E-03, .5045E-03,-.1821E-03,-.1936E-02,
     + .1174E-01, .9273E-02, .5882E-02, .4689E-02, .3454E-02, .1947E-02,
     + .1151E-02, .6070E-03, .6698E-03, .9420E-04,-.6740E-03,-.2707E-02,
     + .1218E-01, .9336E-02, .6050E-02, .4731E-02, .3475E-02, .1863E-02,
     + .1151E-02, .4605E-03, .3768E-03,-.1214E-03,-.4396E-03,-.1903E-02,
     + .1235E-01, .9294E-02, .6029E-02, .4584E-02, .3370E-02, .1800E-02,
     + .1068E-02, .2303E-03, .1675E-03,-.4501E-03,-.7571E-03,-.1149E-02,
     + .1233E-01, .9315E-02, .6029E-02, .4438E-02, .3203E-02, .1842E-02,
     + .9629E-03, .0000E+00,-.2198E-03,-.5338E-03,-.9721E-03,-.7661E-03,
     +-.3692E-04,-.3844E-04,-.2588E-04,-.1180E-04,-.1066E-04,-.3426E-05,
     +-.2664E-05, .7611E-06, .6089E-05,-.4568E-06,-.2077E-04,-.1142E-04,
     +-.3730E-04,-.3806E-04,-.2360E-04,-.1256E-04,-.1180E-04,-.4567E-05,
     +-.3425E-05,-.2284E-05,-.1522E-05,-.4225E-05,-.9940E-05,-.4187E-05,
     +-.3501E-04,-.3844E-04,-.2131E-04,-.1256E-04,-.9896E-05,-.3806E-05,
     +-.4186E-05, .7612E-06,-.1903E-05, .4110E-05, .1789E-05,-.2169E-04,
     +-.3425E-04,-.3882E-04,-.1941E-04,-.1294E-04,-.9515E-05,-.4567E-05,
     +-.4186E-05, .1522E-05,-.4187E-10, .4605E-05,-.2588E-05, .6470E-05,
     +-.3501E-04,-.3730E-04,-.1751E-04,-.1332E-04,-.1066E-04,-.3806E-05,
     +-.4567E-05,-.1142E-05,-.3045E-05, .1104E-05,-.1058E-04, .2816E-04,
     +-.3578E-04,-.3501E-04,-.1751E-04,-.1332E-04,-.1218E-04,-.3806E-05,
     +-.3425E-05,-.3806E-06,-.4187E-05,-.6090E-06,-.6965E-05,-.3463E-04,
     +-.3578E-04,-.3349E-04,-.1675E-04,-.9895E-05,-.9515E-05,-.6090E-05,
     +-.6470E-05,-.3807E-06,-.5328E-05,-.4186E-06,-.3996E-05, .2074E-04,
     +-.3540E-04,-.3083E-04,-.1789E-04,-.9896E-05,-.1104E-04,-.6470E-05,
     +-.5709E-05, .3425E-05,-.4567E-05, .3463E-05, .5633E-05,-.3159E-05,
     +-.3730E-04,-.2740E-04,-.1484E-04,-.1066E-04,-.1142E-04,-.6470E-05,
     +-.6470E-05, .1522E-05,-.1522E-05,-.3045E-05, .3197E-05,-.1039E-04,
     +-.3425E-04,-.2284E-04,-.1370E-04,-.1028E-04,-.1104E-04,-.8373E-05,
     +-.4948E-05, .1903E-05,-.7612E-06,-.1104E-05, .2455E-05,-.3805E-07,
     +-.3235E-04,-.2093E-04,-.1294E-04,-.1142E-04,-.1180E-04,-.6851E-05,
     +-.3045E-05,-.7611E-06, .1256E-05,-.7231E-06, .9924E-05, .3578E-05/
	end
      
	block data ckd4
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
c pressures,  and seven cumulative probabilities ( Fu,  1991 ). The
c spectral region is from 5250 to 4000 cm**-1.
c *********************************************************************
	common /band4/ hk(7), coeh2o(3,11,7)
	data hk / .52, .21, .11, .1, .04, .015, .005 /
c   .253397E+02    .102333E+02    .536032E+01    .487302E+01    .194921E+01
c   .730953E+00    .243651E+00
	data ( ( ( coeh2o(k,j,i), i = 1, 7 ), j = 1, 11 ), k = 1, 3 ) /
     +-.1722E+02,-.1402E+02,-.1202E+02,-.1001E+02,-.7702E+01,-.5273E+01,
     +-.6530E+00,-.1677E+02,-.1359E+02,-.1164E+02,-.9662E+01,-.7419E+01,
     +-.5001E+01,-.6040E+00,-.1630E+02,-.1316E+02,-.1125E+02,-.9303E+01,
     +-.7092E+01,-.4750E+01,-.5715E+00,-.1584E+02,-.1274E+02,-.1086E+02,
     +-.8939E+01,-.6751E+01,-.4458E+01,-.4928E+00,-.1538E+02,-.1232E+02,
     +-.1048E+02,-.8579E+01,-.6399E+01,-.4191E+01,-.4683E+00,-.1493E+02,
     +-.1192E+02,-.1011E+02,-.8241E+01,-.6065E+01,-.3910E+01,-.4310E+00,
     +-.1440E+02,-.1145E+02,-.9643E+01,-.7873E+01,-.5710E+01,-.3668E+01,
     +-.3304E+00,-.1391E+02,-.1104E+02,-.9238E+01,-.7479E+01,-.5367E+01,
     +-.3387E+01,-.3604E+00,-.1348E+02,-.1069E+02,-.8918E+01,-.7122E+01,
     +-.5086E+01,-.3152E+01,-.3030E+00,-.1310E+02,-.1037E+02,-.8626E+01,
     +-.6790E+01,-.4815E+01,-.2945E+01,-.4789E+00,-.1275E+02,-.1011E+02,
     +-.8347E+01,-.6484E+01,-.4584E+01,-.2788E+01,-.5807E+00, .7934E-02,
     + .9231E-02, .1005E-01, .9043E-02, .8164E-02, .8980E-02, .6403E-02,
     + .7954E-02, .9169E-02, .9797E-02, .8687E-02, .7724E-02, .7954E-02,
     + .6652E-02, .7954E-02, .9043E-02, .9608E-02, .8499E-02, .7347E-02,
     + .7473E-02, .6382E-02, .7996E-02, .8980E-02, .9378E-02, .8289E-02,
     + .7264E-02, .6594E-02, .6674E-02, .8059E-02, .8938E-02, .9294E-02,
     + .8227E-02, .7201E-02, .6678E-02, .7032E-02, .8122E-02, .8896E-02,
     + .9189E-02, .8038E-02, .7033E-02, .5987E-02, .5475E-02, .8268E-02,
     + .9064E-02, .8792E-02, .7975E-02, .6573E-02, .5087E-02, .4657E-02,
     + .8541E-02, .8980E-02, .9085E-02, .7996E-02, .6133E-02, .4501E-02,
     + .3860E-02, .8813E-02, .9043E-02, .9294E-02, .8122E-02, .5861E-02,
     + .4354E-02, .3964E-02, .8875E-02, .8834E-02, .9797E-02, .8164E-02,
     + .5463E-02, .4417E-02, .3270E-02, .8938E-02, .8771E-02, .1005E-01,
     + .8247E-02, .5589E-02, .4835E-02, .3033E-02,-.1484E-04,-.2169E-04,
     +-.2436E-04,-.2588E-04,-.1142E-04,-.1142E-05,-.1519E-04,-.1522E-04,
     +-.2055E-04,-.2131E-04,-.2398E-04,-.4948E-05,-.1675E-04,-.3593E-04,
     +-.1522E-04,-.2055E-04,-.1865E-04,-.2207E-04,-.4948E-05,-.1180E-04,
     +-.1237E-04,-.1598E-04,-.2017E-04,-.1903E-04,-.2284E-04,-.1028E-04,
     +-.1865E-04,-.2381E-04,-.1713E-04,-.2017E-04,-.1827E-04,-.2169E-04,
     +-.1218E-04,-.9515E-05,-.2415E-04,-.1827E-04,-.2093E-04,-.1637E-04,
     +-.1827E-04,-.9134E-05,-.8373E-05,-.1243E-04,-.1560E-04,-.1865E-04,
     +-.1599E-04,-.1256E-04,-.1066E-04,-.1142E-05,-.2181E-04,-.1675E-04,
     +-.1560E-04,-.1522E-04,-.1675E-04,-.1865E-04,-.1865E-04,-.9522E-05,
     +-.1332E-04,-.1370E-04,-.1446E-04,-.2055E-04,-.1142E-04,-.2512E-04,
     +-.3343E-04,-.1294E-04,-.1294E-04,-.1751E-04,-.2512E-04,-.1560E-04,
     +-.2854E-04,-.7003E-05,-.8753E-05,-.1028E-04,-.1751E-04,-.2512E-04,
     +-.1713E-04,-.1713E-04,-.1245E-04 /
	end
      
	block data ckd5
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
c pressures,  and twelve cumulative probabilities ( Fu,  1991 ). The
c spectral region is from 4000 to 2850 cm**-1.
c *********************************************************************
	common /band5/ hk(12), coeh2o(3,11,12)
	data hk / .13, .14, .13, .16, .18, .14, 
     1	          .07, .02, .016, .008, .004, .002 /
c   .411549E+01    .443207E+01    .411549E+01    .506522E+01    .569837E+01
c   .443207E+01    .221603E+01    .633153E+00    .506522E+00    .253261E+00
c   .126631E+00    .633153E-01
	data ( ( ( coeh2o(k,j,i), i = 1, 12 ), j = 1, 11 ), k = 1, 3 ) /
     +-.1499E+02,-.1267E+02,-.1118E+02,-.9696E+01,-.7992E+01,-.6323E+01,
     +-.4414E+01,-.2961E+01,-.1715E+01,-.1406E+00, .1612E+01, .3689E+01,
     +-.1454E+02,-.1223E+02,-.1075E+02,-.9277E+01,-.7576E+01,-.5915E+01,
     +-.4043E+01,-.2630E+01,-.1449E+01, .2314E-01, .1708E+01, .3744E+01,
     +-.1408E+02,-.1178E+02,-.1031E+02,-.8851E+01,-.7154E+01,-.5503E+01,
     +-.3666E+01,-.2288E+01,-.1141E+01, .2772E+00, .1819E+01, .3788E+01,
     +-.1363E+02,-.1134E+02,-.9876E+01,-.8423E+01,-.6733E+01,-.5091E+01,
     +-.3286E+01,-.1938E+01,-.8649E+00, .5349E+00, .1969E+01, .3795E+01,
     +-.1318E+02,-.1091E+02,-.9452E+01,-.8004E+01,-.6309E+01,-.4677E+01,
     +-.2904E+01,-.1595E+01,-.5641E+00, .7592E+00, .2109E+01, .3783E+01,
     +-.1275E+02,-.1048E+02,-.9028E+01,-.7585E+01,-.5892E+01,-.4267E+01,
     +-.2524E+01,-.1274E+01,-.2782E+00, .9376E+00, .2257E+01, .3714E+01,
     +-.1180E+02,-.9887E+01,-.8492E+01,-.7014E+01,-.5390E+01,-.3834E+01,
     +-.2156E+01,-.9775E+00,-.3129E-01, .1151E+01, .2330E+01, .3592E+01,
     +-.1114E+02,-.9367E+01,-.8002E+01,-.6514E+01,-.4928E+01,-.3435E+01,
     +-.1835E+01,-.7064E+00, .2153E+00, .1309E+01, .2422E+01, .3488E+01,
     +-.1074E+02,-.8941E+01,-.7582E+01,-.6116E+01,-.4536E+01,-.3072E+01,
     +-.1521E+01,-.4651E+00, .4053E+00, .1465E+01, .2374E+01, .3260E+01,
     +-.1041E+02,-.8545E+01,-.7180E+01,-.5745E+01,-.4177E+01,-.2735E+01,
     +-.1245E+01,-.2356E+00, .5786E+00, .1516E+01, .2263E+01, .3074E+01,
     +-.1008E+02,-.8149E+01,-.6804E+01,-.5409E+01,-.3855E+01,-.2427E+01,
     +-.9857E+00,-.4939E-01, .7060E+00, .1483E+01, .2159E+01, .2745E+01,
     + .9985E-02, .8373E-02, .7431E-02, .6866E-02, .4584E-02, .2952E-02,
     + .3098E-02, .3768E-02, .4013E-02, .3960E-02, .3228E-02, .3203E-02,
     + .1007E-01, .8436E-02, .7368E-02, .6657E-02, .4375E-02, .2617E-02,
     + .2742E-02, .3286E-02, .3192E-02, .2992E-02, .2612E-02, .1968E-02,
     + .1019E-01, .8457E-02, .7264E-02, .6426E-02, .4187E-02, .2365E-02,
     + .2324E-02, .2614E-02, .2736E-02, .2068E-02, .2085E-02, .1005E-02,
     + .1028E-01, .8478E-02, .7138E-02, .6259E-02, .3998E-02, .2156E-02,
     + .1926E-02, .1953E-02, .2250E-02, .1844E-02, .1869E-02,-.6489E-03,
     + .1030E-01, .8478E-02, .7033E-02, .6112E-02, .3852E-02, .1989E-02,
     + .1716E-02, .1763E-02, .1432E-02, .1193E-02, .1306E-02,-.5861E-03,
     + .1042E-01, .8499E-02, .6887E-02, .5987E-02, .3768E-02, .1800E-02,
     + .1549E-02, .1712E-02, .1287E-02, .7389E-03, .7222E-03,-.1130E-02,
     + .8227E-02, .7201E-02, .6866E-02, .5903E-02, .3412E-02, .1591E-02,
     + .1402E-02, .1346E-02, .1041E-02, .8185E-03, .3349E-03,-.4815E-03,
     + .8268E-02, .6992E-02, .7159E-02, .6384E-02, .3286E-02, .1591E-02,
     + .1271E-02, .1202E-02, .9187E-03, .6531E-03,-.4187E-03,-.7954E-03,
     + .8478E-02, .7159E-02, .7117E-02, .6447E-02, .3349E-02, .1528E-02,
     + .9964E-03, .9210E-03, .6112E-03, .6259E-03,-.3768E-03,-.1298E-02,
     + .8520E-02, .7075E-02, .7096E-02, .6405E-02, .3245E-02, .1528E-02,
     + .1011E-02, .7877E-03, .7536E-03, .9001E-04,-.6719E-03,-.1026E-02,
     + .8561E-02, .6950E-02, .7033E-02, .6280E-02, .2993E-02, .1528E-02,
     + .6698E-03, .5847E-03, .2847E-03,-.6280E-04,-.9420E-03,-.1444E-02,
     +-.1408E-04,-.2664E-04,-.1180E-04,-.1903E-04,-.9515E-05, .3806E-06,
     +-.6851E-05,-.3806E-05,-.4834E-05,-.3239E-05,-.2284E-05,-.1028E-04,
     +-.1484E-04,-.2550E-04,-.1142E-04,-.1827E-04,-.9515E-05, .3805E-06,
     +-.4948E-05, .3806E-06,-.2664E-06, .1058E-04,-.1012E-04,-.1142E-04,
     +-.1560E-04,-.2512E-04,-.1256E-04,-.1865E-04,-.9134E-05, .1142E-05,
     +-.3425E-05, .2474E-05,-.9781E-05,-.1519E-05,-.7916E-05,-.1294E-04,
     +-.1560E-04,-.2474E-04,-.1180E-04,-.2017E-04,-.7992E-05, .3805E-06,
     +-.2283E-05,-.4453E-05,-.1180E-05,-.5138E-05,-.4453E-05,-.3425E-05,
     +-.1522E-04,-.2550E-04,-.9896E-05,-.1903E-04,-.9134E-05,-.1142E-05,
     +-.7611E-06,-.5252E-05,-.4567E-06,-.4643E-05,-.4567E-06,-.4567E-05,
     +-.1294E-04,-.2512E-04,-.1028E-04,-.2055E-04,-.9896E-05,-.4567E-05,
     +-.2284E-05,-.5100E-05,-.4339E-06,-.9515E-06,-.1252E-04,-.7612E-06,
     +-.2246E-04,-.1370E-04,-.1066E-04,-.1598E-04,-.8754E-05,-.5328E-05,
     +-.6622E-05,-.5138E-05,-.8754E-07,-.9515E-06, .6090E-05, .4187E-05,
     +-.3463E-04,-.1599E-04,-.1218E-04,-.2093E-04,-.9515E-05,-.4567E-05,
     +-.1104E-05,-.1903E-05,-.1488E-05,-.3730E-05,-.4567E-05, .3045E-05,
     +-.3463E-04,-.1675E-04,-.1294E-04,-.1979E-04,-.1066E-04,-.4187E-05,
     +-.4034E-05,-.2893E-05,-.2588E-05,-.9401E-05, .2284E-05, .3045E-05,
     +-.2778E-04,-.1522E-04,-.1560E-04,-.1751E-04,-.1256E-04,-.5709E-05,
     +-.2474E-05,-.2577E-05,-.2284E-05,-.4187E-06, .7650E-05,-.3425E-05,
     +-.3083E-04,-.1827E-04,-.1370E-04,-.1751E-04,-.1104E-04,-.9515E-05,
     +-.6318E-05,-.4358E-05,-.7613E-07, .4643E-05, .4415E-05, .1028E-04/
	end
      
	block data ckd6
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, eleven 
c pressures,  and  five  cumulative probabilities ( Fu,  1991 ). The
c spectral region is from 2850 to 2500 cm**-1.
c *********************************************************************
	common /band6/ hk(5), coeh2o(3,11,5)
	data hk / .3, .2, .2, .2, .1 /
c   .173978E+01    .115985E+01    .115985E+01    .115985E+01    .579927E+00
	data ( ( ( coeh2o(k,j,i), i = 1, 5 ), j = 1, 11 ), k = 1, 3 ) /
     +-.1905E+02,-.1602E+02,-.1472E+02,-.1307E+02,-.1024E+02,-.1823E+02,
     +-.1555E+02,-.1427E+02,-.1266E+02,-.9938E+01,-.1749E+02,-.1508E+02,
     +-.1381E+02,-.1225E+02,-.9641E+01,-.1684E+02,-.1462E+02,-.1337E+02,
     +-.1185E+02,-.9367E+01,-.1630E+02,-.1417E+02,-.1294E+02,-.1145E+02,
     +-.9123E+01,-.1578E+02,-.1373E+02,-.1251E+02,-.1108E+02,-.8881E+01,
     +-.1517E+02,-.1327E+02,-.1209E+02,-.1072E+02,-.8653E+01,-.1463E+02,
     +-.1284E+02,-.1169E+02,-.1040E+02,-.8453E+01,-.1421E+02,-.1244E+02,
     +-.1133E+02,-.1014E+02,-.8312E+01,-.1382E+02,-.1207E+02,-.1100E+02,
     +-.9887E+01,-.8220E+01,-.1348E+02,-.1173E+02,-.1071E+02,-.9685E+01,
     +-.8220E+01, .1024E-01, .1842E-02, .6908E-03, .1737E-02, .3517E-02,
     + .8394E-02, .2072E-02, .8164E-03, .1716E-02, .2805E-02, .8143E-02,
     + .2240E-02, .9001E-03, .1570E-02, .1800E-02, .8227E-02, .2386E-02,
     + .9420E-03, .1486E-02, .1068E-02, .8373E-02, .2533E-02, .9210E-03,
     + .1319E-02, .9420E-03, .8394E-02, .2700E-02, .9629E-03, .1026E-02,
     + .5233E-03, .8917E-02, .2575E-02, .8792E-03, .7536E-03, .4187E-03,
     + .9378E-02, .2617E-02, .7955E-03, .6070E-03, .4815E-03, .9797E-02,
     + .2638E-02, .6908E-03, .5233E-03, .6280E-03, .1009E-01, .2638E-02,
     + .4815E-03, .2931E-03, .4815E-03, .1036E-01, .2428E-02, .3140E-03,
     + .3977E-03, .2093E-03,-.5366E-04,-.1522E-04,-.5709E-05,-.2664E-05,
     + .3806E-05,-.4301E-04,-.1484E-04,-.4948E-05,-.7610E-06, .7610E-06,
     +-.3920E-04,-.1484E-04,-.4948E-05, .3804E-06,-.3806E-05,-.3920E-04,
     +-.1522E-04,-.4948E-05, .3425E-05, .1903E-05,-.3806E-04,-.1484E-04,
     +-.3045E-05, .2664E-05, .7993E-05,-.4148E-04,-.1408E-04,-.3806E-05,
     + .4187E-05, .7993E-05,-.5481E-04,-.1180E-04,-.3045E-05, .3045E-05,
     + .2284E-05,-.5709E-04,-.1104E-04,-.2283E-05,-.2664E-05,-.1142E-05,
     +-.6090E-04,-.1218E-04,-.2664E-05, .3804E-06, .3045E-05,-.6698E-04,
     +-.1218E-04,-.2664E-05, .1523E-05,-.1142E-05,-.6508E-04,-.1218E-04,
     +-.3425E-05, .1903E-05, .7612E-06 /
	end

	block data ckd7
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and  two  cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 2200 to 1900 cm**-1.
c *********************************************************************
	common /band7/ hk(2), coeh2o(3,19,2)
	data hk / 0.7, 0.3 /
	data ( ( ( coeh2o(k,j,i), i = 1, 2 ), j = 1, 19 ), k = 1, 3 ) /
     +-.2008E+02,-.1467E+02,-.2004E+02,-.1426E+02,-.2001E+02,-.1386E+02,
     +-.1998E+02,-.1345E+02,-.1995E+02,-.1304E+02,-.1992E+02,-.1263E+02,
     +-.1989E+02,-.1223E+02,-.1986E+02,-.1183E+02,-.1984E+02,-.1143E+02,
     +-.1758E+02,-.1038E+02,-.1602E+02,-.9480E+01,-.1469E+02,-.8752E+01,
     +-.1349E+02,-.8218E+01,-.1255E+02,-.7677E+01,-.1174E+02,-.7184E+01,
     +-.1110E+02,-.6735E+01,-.1056E+02,-.6332E+01,-.1019E+02,-.5975E+01,
     +-.9874E+01,-.5644E+01, .2533E-02, .2269E-01, .2575E-02, .2263E-01,
     + .2554E-02, .2267E-01, .2491E-02, .2250E-01, .2449E-02, .2244E-01,
     + .2344E-02, .2234E-01, .2219E-02, .2208E-01, .5694E-02, .2190E-01,
     + .9650E-02, .2162E-01, .3286E-01, .1848E-01, .2987E-01, .1578E-01,
     + .2527E-01, .1465E-01, .2175E-01, .1386E-01, .2056E-01, .1235E-01,
     + .1963E-01, .1116E-01, .1926E-01, .1040E-01, .2014E-01, .1040E-01,
     + .2024E-01, .1042E-01, .1972E-01, .1080E-01,-.8754E-05,-.6698E-04,
     +-.1104E-04,-.6432E-04,-.1142E-04,-.6051E-04,-.1180E-04,-.6128E-04,
     +-.1180E-04,-.6242E-04,-.1218E-04,-.6280E-04,-.1218E-04,-.6204E-04,
     + .5328E-04,-.5709E-04, .1275E-03,-.5214E-04,-.1370E-03,-.4148E-04,
     +-.1100E-03,-.3045E-04,-.9248E-04,-.3197E-04,-.7346E-04,-.2436E-04,
     +-.5100E-04,-.2131E-04,-.5861E-04,-.2550E-04,-.5328E-04,-.3311E-04,
     +-.6090E-04,-.4225E-04,-.5443E-04,-.4415E-04,-.4034E-04,-.4339E-04/
	end

	block data ckd8
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and  three cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 1900 to 1700 cm**-1.
c *********************************************************************
	common /band8/ hk(3), coeh2o(3,19,3)
 	data hk / 0.2, 0.7, 0.1 /
	data ( ( ( coeh2o(k,j,i), i = 1, 3 ), j = 1, 19 ), k = 1, 3 ) /
     +-.2283E+02,-.1639E+02,-.6155E+01,-.2237E+02,-.1595E+02,-.5775E+01,
     +-.2191E+02,-.1551E+02,-.5381E+01,-.2145E+02,-.1507E+02,-.5004E+01,
     +-.2099E+02,-.1463E+02,-.4617E+01,-.2053E+02,-.1419E+02,-.4218E+01,
     +-.2025E+02,-.1375E+02,-.3806E+01,-.2021E+02,-.1330E+02,-.3403E+01,
     +-.2018E+02,-.1287E+02,-.2993E+01,-.1998E+02,-.1091E+02,-.2586E+01,
     +-.1744E+02,-.9171E+01,-.2162E+01,-.1490E+02,-.7642E+01,-.1763E+01,
     +-.1303E+02,-.6526E+01,-.1373E+01,-.1113E+02,-.5846E+01,-.9699E+00,
     +-.9814E+01,-.5280E+01,-.5955E+00,-.8582E+01,-.4787E+01,-.2510E+00,
     +-.8020E+01,-.4350E+01, .2770E-01,-.7571E+01,-.3942E+01, .2406E+00,
     +-.7140E+01,-.3537E+01, .3567E+00, .3722E-01, .1505E-01, .6615E-02,
     + .3722E-01, .1518E-01, .5840E-02, .3720E-01, .1526E-01, .5170E-02,
     + .3399E-01, .1530E-01, .4773E-02, .3012E-01, .1551E-01, .4333E-02,
     + .2625E-01, .1553E-01, .3956E-02, .2240E-01, .1562E-01, .3454E-02,
     + .1846E-01, .1574E-01, .3161E-02, .1446E-01, .1572E-01, .3098E-02,
     + .5924E-02, .8875E-02, .2658E-02, .2204E-01, .7096E-02, .2504E-02,
     + .1591E-01, .5233E-02, .2292E-02, .8855E-02, .4249E-02, .2190E-02,
     + .5422E-02, .3496E-02, .2041E-02, .4919E-02, .3621E-02, .2200E-02,
     + .6657E-02, .3663E-02, .2248E-02, .8645E-02, .3852E-02, .2118E-02,
     + .8771E-02, .3873E-02, .2176E-02, .9043E-02, .3747E-02, .2079E-02,
     +-.1568E-03,-.4681E-04, .4567E-05,-.1568E-03,-.4605E-04,-.3425E-05,
     +-.1572E-03,-.4605E-04,-.1104E-04,-.2154E-03,-.4453E-04,-.6851E-05,
     +-.2843E-03,-.4225E-04,-.7231E-05,-.3562E-03,-.4110E-04,-.7231E-05,
     +-.3692E-03,-.4110E-04,-.1028E-04,-.3007E-03,-.4263E-04,-.6470E-05,
     +-.2325E-03,-.3996E-04,-.8373E-05,-.5290E-04,-.7612E-05,-.4948E-05,
     +-.7422E-04,-.1256E-04,-.8449E-05,-.3501E-04,-.1446E-04,-.4834E-05,
     + .4529E-04,-.2246E-04,-.2893E-05, .6470E-05,-.1789E-04,-.7498E-05,
     +-.4948E-05,-.1713E-04,-.8183E-05,-.5481E-04,-.1713E-04,-.1447E-04,
     +-.4986E-04,-.1903E-04,-.1353E-04,-.5138E-04,-.1484E-04,-.1147E-04,
     +-.5328E-04,-.1560E-04,-.6588E-05/
	end

	block data ckd9
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and  four cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 1700 to 1400 cm**-1.
c *********************************************************************
	common /band9/ hk(4), coeh2o(3,19,4)
 	data hk / 0.22, 0.51, 0.22, 0.05 /
	data ( ( ( coeh2o(k,j,i), i = 1, 4 ), j = 1, 19 ), k = 1, 3 ) /
     +-.2066E+02,-.1464E+02,-.8301E+01,-.3548E+01,-.2025E+02,-.1419E+02,
     +-.7905E+01,-.3260E+01,-.2019E+02,-.1374E+02,-.7495E+01,-.2927E+01,
     +-.2013E+02,-.1329E+02,-.7078E+01,-.2584E+01,-.2007E+02,-.1284E+02,
     +-.6675E+01,-.2247E+01,-.2001E+02,-.1239E+02,-.6268E+01,-.1890E+01,
     +-.1996E+02,-.1194E+02,-.5853E+01,-.1530E+01,-.1991E+02,-.1150E+02,
     +-.5441E+01,-.1133E+01,-.1987E+02,-.1105E+02,-.5022E+01,-.7447E+00,
     +-.1575E+02,-.9657E+01,-.4191E+01,-.3728E+00,-.1329E+02,-.8133E+01,
     +-.3638E+01, .1616E-01,-.1181E+02,-.6675E+01,-.3178E+01, .4083E+00,
     +-.1036E+02,-.5655E+01,-.2731E+01, .7953E+00,-.8628E+01,-.4990E+01,
     +-.2303E+01, .1153E+01,-.7223E+01,-.4453E+01,-.1877E+01, .1454E+01,
     +-.6567E+01,-.3974E+01,-.1461E+01, .1663E+01,-.6077E+01,-.3551E+01,
     +-.1071E+01, .1800E+01,-.5651E+01,-.3136E+01,-.7005E+00, .1809E+01,
     +-.5241E+01,-.2726E+01,-.3859E+00, .1781E+01, .1315E-01, .4542E-02,
     + .3496E-02, .4877E-02, .9650E-02, .4542E-02, .3098E-02, .3956E-02,
     + .6154E-02, .4626E-02, .2763E-02, .3077E-02, .2658E-02, .4626E-02,
     + .2512E-02, .2261E-02, .2658E-02, .4689E-02, .2219E-02, .1405E-02,
     + .2700E-02, .4752E-02, .1926E-02, .7473E-03, .2658E-02, .4773E-02,
     + .1737E-02, .5066E-03, .4668E-02, .4815E-02, .1507E-02, .1842E-03,
     + .8541E-02, .4794E-02, .1382E-02,-.2156E-03, .1022E-01, .2198E-02,
     + .3977E-03,-.2910E-03, .5484E-02, .6698E-03, .0000E+00,-.2339E-03,
     + .3349E-02, .1068E-02,-.2512E-03,-.4228E-03, .1884E-02, .2093E-03,
     +-.3977E-03,-.6405E-03,-.8373E-04,-.5233E-03,-.4124E-03,-.5945E-03,
     + .7536E-03,-.6698E-03,-.4919E-03,-.4794E-03, .3600E-02,-.4605E-03,
     +-.4375E-03,-.3517E-03, .3873E-02,-.5861E-03,-.3203E-03,-.4689E-03,
     + .3935E-02,-.7326E-03,-.2072E-03,-.4228E-03, .4124E-02,-.8582E-03,
     +-.4187E-04,-.5945E-03,-.8525E-04, .1865E-04,-.1142E-05, .2664E-05,
     +-.1313E-03, .1865E-04, .0000E+00, .1256E-04,-.6470E-04, .1865E-04,
     +-.3045E-05, .8754E-05, .3805E-06, .1789E-04,-.6851E-05, .5328E-05,
     + .1142E-05, .1827E-04,-.6090E-05, .4148E-05, .1142E-05, .1865E-04,
     +-.3806E-05,-.3768E-05,-.1903E-05, .1751E-04,-.4948E-05, .3121E-05,
     + .3159E-04, .1979E-04,-.3045E-05,-.9896E-06, .1005E-03, .1789E-04,
     +-.6089E-05,-.1865E-05,-.2207E-04, .1941E-04, .1903E-05, .2322E-05,
     +-.1675E-04, .6090E-05,-.7611E-06, .4397E-05, .3425E-04, .3806E-06,
     + .1522E-05, .3806E-05, .4796E-04, .1522E-05,-.3806E-06, .3654E-05,
     +-.6851E-05, .2664E-05,-.3920E-05,-.6850E-06,-.1370E-04, .5328E-05,
     +-.6584E-05,-.8716E-05,-.8374E-10, .1522E-05,-.6356E-05, .1294E-05,
     +-.9515E-05, .7612E-06,-.3235E-05,-.1066E-05,-.7612E-05, .1142E-05,
     +-.4529E-05, .3730E-05,-.2664E-05,-.3806E-06,-.3501E-05,-.5328E-06/
	end

	block data ckd10
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and  four cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 1400 to 1250 cm**-1. coech4 and coen2o
c are the coefficients to calculate the CH4 and N2O absorption coe-
c fficients in units of (cm-atm)**-1 at three temperature, nineteen
c pressures, and one cumulative probability (Fu, 1991), respectively.
c *********************************************************************
	common /band10/hk(4), coeh2o(3,19,4), coech4(3,19), coen2o(3,19)
 	data hk / 0.28, 0.42, 0.25, 0.05 /
	data ( ( ( coeh2o(k,j,i), i = 1, 4 ), j = 1, 19 ), k = 1, 3 ) /
     +-.2023E+02,-.1641E+02,-.1171E+02,-.6090E+01,-.2016E+02,-.1595E+02,
     +-.1133E+02,-.5867E+01,-.2011E+02,-.1550E+02,-.1095E+02,-.5660E+01,
     +-.2005E+02,-.1504E+02,-.1055E+02,-.5407E+01,-.2001E+02,-.1459E+02,
     +-.1015E+02,-.5137E+01,-.1997E+02,-.1413E+02,-.9749E+01,-.4852E+01,
     +-.1993E+02,-.1367E+02,-.9337E+01,-.4534E+01,-.1990E+02,-.1321E+02,
     +-.8920E+01,-.4211E+01,-.1987E+02,-.1276E+02,-.8506E+01,-.3889E+01,
     +-.1645E+02,-.1179E+02,-.7711E+01,-.3613E+01,-.1442E+02,-.1081E+02,
     +-.6942E+01,-.3316E+01,-.1308E+02,-.9950E+01,-.6344E+01,-.2950E+01,
     +-.1212E+02,-.9217E+01,-.5904E+01,-.2577E+01,-.1131E+02,-.8559E+01,
     +-.5519E+01,-.2256E+01,-.1064E+02,-.7962E+01,-.5183E+01,-.1929E+01,
     +-.1013E+02,-.7447E+01,-.4833E+01,-.1643E+01,-.9712E+01,-.7071E+01,
     +-.4485E+01,-.1410E+01,-.9305E+01,-.6760E+01,-.4145E+01,-.1249E+01,
     +-.8966E+01,-.6477E+01,-.3820E+01,-.1114E+01, .7913E-02, .8206E-02,
     + .1509E-01, .1869E-01, .4228E-02, .8247E-02, .1467E-01, .1783E-01,
     + .2010E-02, .8227E-02, .1442E-01, .1687E-01, .1947E-02, .8289E-02,
     + .1394E-01, .1568E-01, .1863E-02, .8289E-02, .1346E-01, .1484E-01,
     + .1842E-02, .8415E-02, .1310E-01, .1400E-01, .1800E-02, .8457E-02,
     + .1275E-01, .1377E-01, .1696E-02, .8478E-02, .1220E-01, .1321E-01,
     + .1842E-02, .8478E-02, .1189E-01, .1250E-01, .1409E-01, .8624E-02,
     + .1254E-01, .1214E-01, .9043E-02, .1045E-01, .1225E-01, .1260E-01,
     + .8561E-02, .1202E-01, .1181E-01, .1296E-01, .1114E-01, .1235E-01,
     + .1191E-01, .1330E-01, .1199E-01, .1271E-01, .1195E-01, .1371E-01,
     + .1415E-01, .1315E-01, .1218E-01, .1361E-01, .1478E-01, .1338E-01,
     + .1296E-01, .1306E-01, .1518E-01, .1375E-01, .1365E-01, .1334E-01,
     + .1530E-01, .1411E-01, .1392E-01, .1327E-01, .1547E-01, .1507E-01,
     + .1390E-01, .1264E-01,-.1089E-03,-.2740E-04,-.2017E-04,-.5519E-04,
     +-.4491E-04,-.2740E-04,-.1408E-04,-.5937E-04,-.6090E-05,-.2702E-04,
     +-.6470E-05,-.4719E-04,-.7232E-05,-.2740E-04,-.6089E-05,-.4910E-04,
     +-.7231E-05,-.2969E-04,-.4186E-05,-.5366E-04,-.6090E-05,-.3045E-04,
     +-.2284E-05,-.4986E-04,-.4568E-05,-.3121E-04,-.4948E-05,-.5100E-04,
     +-.3426E-05,-.3007E-04,-.7993E-05,-.4910E-04, .1522E-05,-.2931E-04,
     +-.9896E-05,-.5366E-04,-.5823E-04,-.1599E-04,-.1713E-04,-.4110E-04,
     +-.3121E-04,-.1713E-04,-.3159E-04,-.3578E-04,-.3996E-04,-.1598E-04,
     +-.3958E-04,-.4605E-04,-.3349E-04,-.1751E-04,-.3844E-04,-.5576E-04,
     +-.2626E-04,-.2474E-04,-.3920E-04,-.4464E-04,-.1979E-04,-.3045E-04,
     +-.3958E-04,-.5336E-04,-.2893E-04,-.3616E-04,-.3996E-04,-.4754E-04,
     +-.2398E-04,-.3083E-04,-.4415E-04,-.5119E-04,-.2702E-04,-.2664E-04,
     +-.4605E-04,-.4038E-04,-.2398E-04,-.2360E-04,-.4948E-04,-.5149E-04/
	data ( ( coech4(k,j), j = 1, 19 ), k = 1, 3 ) /
     +-.8909E+01,-.8464E+01,-.8018E+01,-.7573E+01,-.7133E+01,-.6687E+01,
     +-.6240E+01,-.5803E+01,-.5377E+01,-.4534E+01,-.3983E+01,-.3502E+01,
     +-.3062E+01,-.2648E+01,-.2265E+01,-.1896E+01,-.1568E+01,-.1234E+01,
     +-.9298E+00, .9629E-03, .9838E-03, .1088E-02, .1172E-02, .1256E-02,
     + .1402E-02, .1528E-02, .1633E-02, .1716E-02, .4815E-03,-.3977E-03,
     +-.5652E-03,-.5024E-03,-.4605E-03,-.4563E-03,-.4438E-03,-.4521E-03,
     +-.4312E-03,-.3789E-03,-.1294E-04,-.1408E-04,-.1522E-04,-.1675E-04,
     +-.1751E-04,-.1941E-04,-.2246E-04,-.2207E-04,-.1827E-04,-.1256E-04,
     +-.9515E-05,-.6470E-05,-.3045E-05,-.3806E-05,-.2055E-05,-.3730E-05,
     +-.7612E-06,-.3806E-05, .1256E-05/
	data ( ( coen2o(k,j), j = 1, 19 ), k = 1, 3 ) /
     +-.7863E+01,-.7412E+01,-.6963E+01,-.6514E+01,-.6065E+01,-.5611E+01,
     +-.5167E+01,-.4720E+01,-.4283E+01,-.3454E+01,-.2858E+01,-.2404E+01,
     +-.1922E+01,-.1491E+01,-.1097E+01,-.7177E+00,-.3548E+00, .1218E-01,
     + .3088E+00, .4459E-02, .4542E-02, .4668E-02, .4752E-02, .4815E-02,
     + .4919E-02, .5087E-02, .5254E-02, .5296E-02, .2324E-02, .2093E-02,
     + .2294E-02, .2125E-02, .2058E-02, .1920E-02, .1786E-02, .1689E-02,
     + .1788E-02, .2144E-02,-.7231E-05,-.7231E-05,-.7231E-05,-.6470E-05,
     +-.6851E-05,-.7231E-05,-.5709E-05,-.6470E-05,-.4186E-05, .8754E-05,
     +-.7612E-05,-.9134E-06,-.8640E-05,-.8487E-05,-.8259E-05,-.9553E-05,
     +-.8107E-05,-.1654E-04,-.1858E-04/
	end

	block data ckd11
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and three cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 1250 to 1100 cm**-1. coech4 and coen2o
c are the coefficients to calculate the CH4 and N2O absorption coe-
c fficients in units of (cm-atm)**-1 at three temperature, nineteen
c pressures, and one cumulative probability (Fu, 1991), respectively.
c *********************************************************************
	common /band11/hk(3), coeh2o(3,19,3), coech4(3,19), coen2o(3,19)
 	data hk / 0.80, 0.15, 0.05 /
	data ( ( ( coeh2o(k,j,i), i = 1, 3 ), j = 1, 19 ), k = 1, 3 ) /
     +-.2005E+02,-.1548E+02,-.1021E+02,-.2001E+02,-.1504E+02,-.1001E+02,
     +-.1997E+02,-.1459E+02,-.9814E+01,-.1993E+02,-.1416E+02,-.9595E+01,
     +-.1989E+02,-.1373E+02,-.9349E+01,-.1985E+02,-.1328E+02,-.9072E+01,
     +-.1982E+02,-.1286E+02,-.8833E+01,-.1957E+02,-.1243E+02,-.8566E+01,
     +-.1911E+02,-.1200E+02,-.8276E+01,-.1743E+02,-.1134E+02,-.7958E+01,
     +-.1625E+02,-.1078E+02,-.7629E+01,-.1524E+02,-.1036E+02,-.7334E+01,
     +-.1429E+02,-.9970E+01,-.7051E+01,-.1348E+02,-.9620E+01,-.6749E+01,
     +-.1282E+02,-.9270E+01,-.6505E+01,-.1229E+02,-.8932E+01,-.6277E+01,
     +-.1186E+02,-.8628E+01,-.6120E+01,-.1148E+02,-.8345E+01,-.6049E+01,
     +-.1112E+02,-.8066E+01,-.5906E+01, .1842E-02, .2131E-01, .3033E-01,
     + .1905E-02, .2137E-01, .2841E-01, .1926E-02, .2135E-01, .2696E-01,
     + .1926E-02, .2133E-01, .2514E-01, .1884E-02, .2154E-01, .2401E-01,
     + .5589E-02, .2156E-01, .2321E-01, .9483E-02, .2156E-01, .2210E-01,
     + .1333E-01, .2150E-01, .2133E-01, .1725E-01, .2154E-01, .2074E-01,
     + .2254E-01, .1999E-01, .2005E-01, .2118E-01, .1926E-01, .1978E-01,
     + .1936E-01, .1920E-01, .1963E-01, .1905E-01, .1911E-01, .1934E-01,
     + .1909E-01, .1903E-01, .1920E-01, .1922E-01, .1901E-01, .1899E-01,
     + .1934E-01, .1930E-01, .1974E-01, .1966E-01, .1909E-01, .2014E-01,
     + .1976E-01, .1905E-01, .1984E-01, .1963E-01, .1940E-01, .1897E-01,
     +-.1522E-05,-.6013E-04,-.5062E-04,-.2665E-05,-.6204E-04,-.5519E-04,
     +-.3806E-05,-.6394E-04,-.5633E-04,-.4567E-05,-.6280E-04,-.5214E-04,
     +-.6090E-05,-.6128E-04,-.5290E-04, .6051E-04,-.6242E-04,-.5823E-04,
     + .1313E-03,-.6013E-04,-.5176E-04, .1336E-03,-.5747E-04,-.4072E-04,
     + .6318E-04,-.5671E-04,-.3996E-04,-.5595E-04,-.3996E-04,-.4263E-04,
     +-.3958E-04,-.4719E-04,-.4453E-04,-.3387E-04,-.5138E-04,-.5100E-04,
     +-.5252E-04,-.4986E-04,-.4491E-04,-.5100E-04,-.4453E-04,-.4529E-04,
     +-.5176E-04,-.4795E-04,-.4453E-04,-.5557E-04,-.5176E-04,-.5062E-04,
     +-.5747E-04,-.4795E-04,-.5633E-04,-.5709E-04,-.4643E-04,-.3806E-04,
     +-.5481E-04,-.5671E-04,-.4948E-04/
	data ( ( coech4(k,j), j = 1, 19 ), k = 1, 3 ) /
     +-.1207E+02,-.1162E+02,-.1116E+02,-.1070E+02,-.1024E+02,-.9777E+01,
     +-.9319E+01,-.8858E+01,-.8398E+01,-.7384E+01,-.6643E+01,-.6081E+01,
     +-.5602E+01,-.5188E+01,-.4822E+01,-.4479E+01,-.4184E+01,-.3884E+01,
     +-.3627E+01, .1036E-01, .1036E-01, .1040E-01, .1040E-01, .1045E-01,
     + .1047E-01, .1049E-01, .1055E-01, .1059E-01, .1059E-01, .1026E-01,
     + .1011E-01, .1024E-01, .1049E-01, .1072E-01, .1089E-01, .1109E-01,
     + .1153E-01, .1191E-01,-.4910E-04,-.4834E-04,-.4910E-04,-.4910E-04,
     +-.4910E-04,-.4872E-04,-.4834E-04,-.4948E-04,-.5100E-04,-.5633E-04,
     +-.6166E-04,-.5595E-04,-.5366E-04,-.5366E-04,-.5328E-04,-.5328E-04,
     +-.4948E-04,-.5519E-04,-.5595E-04/
	data ( ( coen2o(k,j), j = 1, 19 ), k = 1, 3 ) /
     +-.9461E+01,-.9003E+01,-.8543E+01,-.8084E+01,-.7629E+01,-.7166E+01,
     +-.6707E+01,-.6249E+01,-.5793E+01,-.5312E+01,-.4847E+01,-.4393E+01,
     +-.3974E+01,-.3587E+01,-.3231E+01,-.2885E+01,-.2602E+01,-.2358E+01,
     +-.2108E+01, .4710E-02, .4752E-02, .4773E-02, .4773E-02, .4815E-02,
     + .4877E-02, .4898E-02, .4982E-02, .5066E-02, .5296E-02, .5149E-02,
     + .5129E-02, .5024E-02, .4752E-02, .4501E-02, .4270E-02, .4019E-02,
     + .3646E-02, .2759E-02,-.1484E-04,-.1408E-04,-.1446E-04,-.1446E-04,
     +-.1522E-04,-.1560E-04,-.1522E-04,-.1522E-04,-.1598E-04,-.1484E-04,
     +-.9895E-05,-.1028E-04,-.7612E-05,-.1903E-05, .1903E-05, .0000E+00,
     + .2283E-05, .6166E-05,-.2740E-05/
	end

	block data ckd12
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeo3 is the coefficient to calculate the ozone absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and  five cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 1100 to  980 cm**-1.    coeh2o is the
c coefficient to calculate the H2O absorption coefficient in units
c of (cm-atm)**-1 at three temperature, nineteen pressures, and one
c cumulative probability ( Fu, 1991 ).
c *********************************************************************
	common /band12/ hk(5), coeo3(3,19,5), coeh2o(3,19)
	data hk / 0.45, 0.30, 0.2, 0.04, 0.01 /
	data ( ( ( coeo3(k,j,i), i = 1, 5 ), j = 1, 19 ), k = 1, 3 ) /
     +-.6590E+01,-.3912E+01,-.8513E+00, .2731E+01, .5515E+01,-.6157E+01,
     +-.3583E+01,-.7292E+00, .2740E+01, .5508E+01,-.5731E+01,-.3242E+01,
     +-.5800E+00, .2782E+01, .5485E+01,-.5301E+01,-.2901E+01,-.4131E+00,
     + .2805E+01, .5455E+01,-.4879E+01,-.2551E+01,-.2288E+00, .2878E+01,
     + .5416E+01,-.4449E+01,-.2201E+01,-.2228E-01, .3000E+01, .5374E+01,
     +-.4018E+01,-.1843E+01, .2055E+00, .3143E+01, .5342E+01,-.3615E+01,
     +-.1502E+01, .4561E+00, .3288E+01, .5204E+01,-.3228E+01,-.1172E+01,
     + .7099E+00, .3396E+01, .5077E+01,-.2828E+01,-.8499E+00, .9664E+00,
     + .3463E+01, .4893E+01,-.2480E+01,-.5393E+00, .1229E+01, .3493E+01,
     + .4656E+01,-.2181E+01,-.2653E+00, .1504E+01, .3456E+01, .4398E+01,
     +-.1950E+01,-.1469E-01, .1735E+01, .3387E+01, .4115E+01,-.1788E+01,
     + .2517E+00, .1919E+01, .3251E+01, .3832E+01,-.1677E+01, .5027E+00,
     + .2032E+01, .3088E+01, .3581E+01,-.1637E+01, .7373E+00, .2100E+01,
     + .2910E+01, .3364E+01,-.1650E+01, .9383E+00, .2123E+01, .2793E+01,
     + .3150E+01,-.1658E+01, .1091E+01, .2112E+01, .2683E+01, .3021E+01,
     +-.1654E+01, .1163E+01, .2099E+01, .2602E+01, .2871E+01, .9498E-02,
     + .8894E-02, .1161E-01, .8828E-02,-.1669E-02, .9613E-02, .8347E-02,
     + .1053E-01, .8462E-02,-.1612E-02, .9700E-02, .7829E-02, .9101E-02,
     + .7915E-02,-.1439E-02, .9815E-02, .7167E-02, .7981E-02, .7282E-02,
     +-.1094E-02, .9671E-02, .6764E-02, .6930E-02, .5613E-02,-.8347E-03,
     + .9613E-02, .6312E-02, .6225E-02, .4145E-02,-.1295E-02, .9728E-02,
     + .6099E-02, .5293E-02, .2965E-02,-.1756E-02, .9844E-02, .5915E-02,
     + .4496E-02, .1871E-02,-.2044E-02, .9930E-02, .5817E-02, .3509E-02,
     + .1324E-02,-.2044E-02, .9988E-02, .5535E-02, .2711E-02, .6620E-03,
     +-.1813E-02, .1034E-01, .5247E-02, .1926E-02,-.2303E-03,-.1842E-02,
     + .1058E-01, .4795E-02, .1197E-02,-.9498E-03,-.2216E-02, .1084E-01,
     + .4414E-02, .6188E-03,-.1123E-02,-.2303E-02, .1079E-01, .3926E-02,
     + .1756E-03,-.1497E-02,-.2274E-02, .1039E-01, .3425E-02,-.1900E-03,
     +-.1353E-02,-.2389E-02, .9815E-02, .2769E-02,-.6620E-03,-.1756E-02,
     +-.1785E-02, .9818E-02, .2444E-02,-.1016E-02,-.1410E-02,-.1698E-02,
     + .1074E-01, .3218E-02,-.1235E-02,-.1900E-02,-.2533E-02, .1145E-01,
     + .3684E-02,-.1364E-02,-.1353E-02,-.1957E-02,-.4030E-04,-.2375E-04,
     +-.3814E-05,-.4943E-04,-.3166E-04,-.3742E-04,-.1871E-04,-.1137E-04,
     +-.4317E-04,-.2878E-04,-.3526E-04,-.2015E-04,-.1295E-04,-.4821E-04,
     +-.2303E-04,-.3382E-04,-.2087E-04,-.1519E-04,-.2231E-04,-.1871E-04,
     +-.3454E-04,-.2087E-04,-.8109E-05,-.6476E-05,-.1511E-04,-.3454E-04,
     +-.1820E-04,-.1269E-05,-.1439E-04,-.5037E-05,-.4173E-04,-.2598E-04,
     + .6645E-05,-.1943E-04,-.2087E-04,-.3454E-04,-.2267E-04, .2159E-05,
     +-.2231E-04,-.2159E-05,-.2950E-04,-.2080E-04, .2159E-06,-.4317E-05,
     + .1799E-04,-.3670E-04,-.1590E-04,-.4461E-05,-.9354E-05,-.3598E-05,
     +-.3216E-04,-.1475E-04,-.2231E-05,-.1295E-04,-.2878E-05,-.3576E-04,
     +-.7347E-05,-.1022E-04,-.2159E-05,-.7915E-05,-.3015E-04,-.5230E-05,
     +-.5109E-05,-.6476E-05,-.7196E-05,-.2331E-04,-.1079E-04,-.4102E-05,
     + .1439E-05,-.1223E-04,-.2216E-04,-.1094E-04,-.5325E-05,-.7196E-06,
     +-.1655E-04,-.1036E-04,-.7627E-05,-.2878E-05, .5037E-05,-.1295E-04,
     + .1029E-04,-.1346E-04,-.4821E-05,-.7915E-05, .7915E-05, .2835E-04,
     +-.2893E-04,-.1367E-05,-.7196E-05,-.1871E-04, .3965E-04,-.3310E-04,
     +-.3310E-05,-.7195E-06, .2303E-04/
	data ( ( coeh2o(k,j), j = 1, 19 ), k = 1, 3 ) /
     +-.1984E+02,-.1983E+02,-.1982E+02,-.1981E+02,-.1963E+02,-.1917E+02,
     +-.1871E+02,-.1825E+02,-.1779E+02,-.1639E+02,-.1545E+02,-.1484E+02,
     +-.1433E+02,-.1387E+02,-.1345E+02,-.1305E+02,-.1268E+02,-.1231E+02,
     +-.1196E+02, .6071E-03, .2072E-02, .6196E-02, .1030E-01, .1436E-01,
     + .1846E-01, .2259E-01, .2667E-01, .2993E-01, .2878E-01, .2803E-01,
     + .2851E-01, .2864E-01, .2874E-01, .2862E-01, .2859E-01, .2853E-01,
     + .2868E-01, .2887E-01,-.3808E-06, .2474E-04, .9895E-04, .1728E-03,
     + .1911E-03, .1165E-03, .4225E-04,-.3121E-04,-.8982E-04,-.9553E-04,
     +-.9705E-04,-.9591E-04,-.9287E-04,-.9172E-04,-.9096E-04,-.9134E-04,
     +-.9248E-04,-.1050E-03,-.1031E-03/
	end

	block data ckd13
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and  two  cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 980 to 800 cm**-1.
c *********************************************************************
	common /band13/ hk(2), coeh2o(3,19,2)
 	data hk / 0.95, 0.05 /
	data ( ( ( coeh2o(k,j,i), i = 1, 2 ), j = 1, 19 ), k = 1, 3 ) /
     +-.1992E+02,-.1446E+02,-.1992E+02,-.1405E+02,-.1991E+02,-.1363E+02,
     +-.1990E+02,-.1322E+02,-.1989E+02,-.1282E+02,-.1989E+02,-.1242E+02,
     +-.1988E+02,-.1201E+02,-.1987E+02,-.1159E+02,-.1986E+02,-.1119E+02,
     +-.1982E+02,-.1079E+02,-.1817E+02,-.1039E+02,-.1659E+02,-.1000E+02,
     +-.1537E+02,-.9623E+01,-.1460E+02,-.9266E+01,-.1406E+02,-.8959E+01,
     +-.1354E+02,-.8676E+01,-.1309E+02,-.8411E+01,-.1267E+02,-.8232E+01,
     +-.1229E+02,-.8094E+01, .5024E-03, .3199E-01, .5652E-03, .3199E-01,
     + .6071E-03, .3211E-01, .6489E-03, .3199E-01, .6699E-03, .3178E-01,
     + .6908E-03, .3157E-01, .6908E-03, .3109E-01, .6698E-03, .3075E-01,
     + .6698E-03, .3054E-01, .1474E-01, .3000E-01, .3085E-01, .2960E-01,
     + .3659E-01, .2935E-01, .3016E-01, .2920E-01, .2834E-01, .2895E-01,
     + .2780E-01, .2870E-01, .2753E-01, .2843E-01, .2755E-01, .2820E-01,
     + .2765E-01, .2732E-01, .2769E-01, .2705E-01, .6299E-09,-.7993E-04,
     +-.3802E-06,-.7992E-04,-.3802E-06,-.8525E-04,-.3808E-06,-.8449E-04,
     +-.7610E-06,-.7764E-04,-.1142E-05,-.7231E-04,-.1142E-05,-.7345E-04,
     +-.2284E-05,-.8259E-04,-.2284E-05,-.8031E-04, .2436E-03,-.7878E-04,
     + .7612E-05,-.8525E-04,-.1248E-03,-.9439E-04,-.9477E-04,-.9172E-04,
     +-.8982E-04,-.8640E-04,-.7916E-04,-.6813E-04,-.7574E-04,-.6090E-04,
     +-.7612E-04,-.7117E-04,-.7498E-04,-.7041E-04,-.7269E-04,-.7992E-04/
	end

	block data ckd14
c **********************************************************************
c hk is the interval in the g (cumulative probability) space from 0
c to one. coehca and coehcb are the coefficients to calculate the
c H2O and CO2 overlapping absorption coefficients in units of (cm-
c atm)**-1 at three temperature, nineteen pressures, and ten cumu-
c lative probabilities (Fu, 1991). The spectral region is from 800
c to 670 cm**-1.
c **********************************************************************
	common /band14/ hk(10), coehca(3,19,10), coehcb(3,19,10)
 	data hk / .3,.3,.2,.12,.06,.012,.004,.0025,.0011,.0004 /
	data ( ( ( coehca(k,j,i), i = 1, 10 ), j = 1, 19 ), k = 1, 3 ) /
     +-.1847E+02,-.1399E+02,-.1106E+02,-.8539E+01,-.5852E+01,-.3295E+01,
     +-.1208E+01,-.6272E-01, .2055E+01, .6071E+01,-.1801E+02,-.1357E+02,
     +-.1067E+02,-.8171E+01,-.5562E+01,-.3071E+01,-.1073E+01, .1033E+00,
     + .2055E+01, .6071E+01,-.1755E+02,-.1314E+02,-.1027E+02,-.7798E+01,
     +-.5224E+01,-.2823E+01,-.9280E+00, .2723E+00, .2165E+01, .5969E+01,
     +-.1709E+02,-.1272E+02,-.9868E+01,-.7404E+01,-.4880E+01,-.2569E+01,
     +-.6908E+00, .4453E+00, .2241E+01, .5969E+01,-.1663E+02,-.1230E+02,
     +-.9467E+01,-.7013E+01,-.4535E+01,-.2297E+01,-.4408E+00, .6353E+00,
     + .2359E+01, .5969E+01,-.1617E+02,-.1188E+02,-.9050E+01,-.6619E+01,
     +-.4160E+01,-.1967E+01,-.1687E+00, .8213E+00, .2421E+01, .5969E+01,
     +-.1571E+02,-.1147E+02,-.8629E+01,-.6230E+01,-.3771E+01,-.1648E+01,
     + .1573E+00, .1019E+01, .2511E+01, .5884E+01,-.1525E+02,-.1106E+02,
     +-.8215E+01,-.5841E+01,-.3393E+01,-.1331E+01, .4013E+00, .1198E+01,
     + .2654E+01, .5794E+01,-.1480E+02,-.1066E+02,-.7800E+01,-.5454E+01,
     +-.3032E+01,-.9870E+00, .6323E+00, .1373E+01, .2905E+01, .5647E+01,
     +-.1402E+02,-.9693E+01,-.7206E+01,-.4846E+01,-.2656E+01,-.6540E+00,
     + .8323E+00, .1530E+01, .3211E+01, .5355E+01,-.1343E+02,-.9060E+01,
     +-.6596E+01,-.4399E+01,-.2294E+01,-.3519E+00, .9823E+00, .1673E+01,
     + .3420E+01, .5083E+01,-.1279E+02,-.8611E+01,-.5785E+01,-.4010E+01,
     +-.1936E+01,-.1177E+00, .1134E+01, .1974E+01, .3591E+01, .4770E+01,
     +-.1230E+02,-.8174E+01,-.5298E+01,-.3611E+01,-.1607E+01, .3636E-01,
     + .1433E+01, .2260E+01, .3539E+01, .4439E+01,-.1192E+02,-.7763E+01,
     +-.4946E+01,-.3228E+01,-.1321E+01, .1991E+00, .1720E+01, .2420E+01,
     + .3383E+01, .4041E+01,-.1154E+02,-.7377E+01,-.4576E+01,-.2851E+01,
     +-.1093E+01, .4430E+00, .1896E+01, .2462E+01, .3122E+01, .3620E+01,
     +-.1118E+02,-.7003E+01,-.4210E+01,-.2524E+01,-.8973E+00, .7490E+00,
     + .1966E+01, .2363E+01, .2818E+01, .3182E+01,-.1080E+02,-.6677E+01,
     +-.3872E+01,-.2264E+01,-.6846E+00, .9392E+00, .1867E+01, .2138E+01,
     + .2505E+01, .2738E+01,-.1031E+02,-.6353E+01,-.3596E+01,-.1938E+01,
     +-.4537E+00, .1015E+01, .1659E+01, .1830E+01, .2142E+01, .2287E+01,
     +-.9695E+01,-.5977E+01,-.3427E+01,-.1596E+01,-.1979E+00, .9458E+00,
     + .1363E+01, .1545E+01, .1743E+01, .1832E+01, .3628E-01, .2728E-01,
     + .2213E-01, .1656E-01, .1507E-01, .1564E-01, .1623E-01, .1419E-01,
     + .1455E-01, .1089E-02, .3632E-01, .2740E-01, .2164E-01, .1606E-01,
     + .1369E-01, .1418E-01, .1444E-01, .1275E-01, .1331E-01, .9210E-03,
     + .3636E-01, .2746E-01, .2114E-01, .1557E-01, .1239E-01, .1285E-01,
     + .1237E-01, .1141E-01, .1141E-01, .9210E-03, .3640E-01, .2748E-01,
     + .2064E-01, .1516E-01, .1141E-01, .1125E-01, .1092E-01, .1026E-01,
     + .1011E-01,-.5652E-03, .3646E-01, .2746E-01, .2024E-01, .1478E-01,
     + .1036E-01, .9688E-02, .9610E-02, .9305E-02, .9399E-02,-.6489E-03,
     + .3651E-01, .2734E-01, .1984E-01, .1438E-01, .9436E-02, .8486E-02,
     + .8214E-02, .8995E-02, .7892E-02,-.8582E-03, .3655E-01, .2723E-01,
     + .1951E-01, .1402E-01, .8716E-02, .7433E-02, .7169E-02, .8072E-02,
     + .5443E-02,-.1172E-02, .3659E-01, .2709E-01, .1911E-01, .1379E-01,
     + .8107E-02, .6818E-02, .6818E-02, .7033E-02, .3056E-02,-.1047E-02,
     + .3670E-01, .2698E-01, .1890E-01, .1363E-01, .7502E-02, .6371E-02,
     + .6558E-02, .6489E-02,-.5652E-03,-.1340E-02, .3592E-01, .2238E-01,
     + .1804E-01, .1007E-01, .6730E-02, .5512E-02, .6194E-02, .4375E-02,
     +-.1109E-02,-.3559E-03, .3609E-01, .2242E-01, .1526E-01, .8582E-02,
     + .6284E-02, .5809E-02, .4501E-02, .9420E-03,-.9001E-03,-.1005E-02,
     + .3703E-01, .2196E-01, .1281E-01, .7860E-02, .5861E-02, .5842E-02,
     + .1800E-02,-.1591E-02,-.1235E-02,-.9420E-03, .3728E-01, .2114E-01,
     + .1347E-01, .6678E-02, .5449E-02, .4837E-02,-.1084E-02,-.1361E-02,
     +-.6699E-03,-.1256E-03, .3683E-01, .2061E-01, .1350E-01, .6133E-02,
     + .5449E-02, .2111E-02,-.1386E-02,-.1235E-02,-.5652E-03,-.8373E-04,
     + .3656E-01, .1988E-01, .1348E-01, .5441E-02, .5149E-02,-.8813E-03,
     +-.1116E-02,-.8373E-03,-.3140E-03,-.6280E-04, .3669E-01, .1934E-01,
     + .1363E-01, .5035E-02, .3585E-02,-.1250E-02,-.9357E-03,-.8227E-03,
     +-.3140E-03,-.4187E-04, .3618E-01, .1856E-01, .1390E-01, .3836E-02,
     + .1470E-02,-.1096E-02,-.8080E-03,-.4480E-03,-.2093E-03,-.2093E-04,
     + .3416E-01, .1741E-01, .1431E-01, .1951E-02,-.2923E-04,-.9422E-03,
     +-.4576E-03,-.2395E-03,-.1565E-03,-.2799E-04, .3219E-01, .1674E-01,
     + .1516E-01, .6652E-03,-.5051E-03,-.7052E-03,-.2002E-03,-.2135E-03,
     +-.7633E-04,-.7300E-04,-.1290E-03,-.9934E-04,-.5595E-04,-.3996E-04,
     + .1294E-04,-.9134E-05, .1294E-05,-.3121E-05,-.4757E-04,-.1979E-04,
     +-.1305E-03,-.9629E-04,-.5481E-04,-.4301E-04, .1827E-04,-.9363E-05,
     + .1777E-04,-.2185E-04,-.1903E-04,-.1675E-04,-.1313E-03,-.9439E-04,
     +-.5404E-04,-.4263E-04, .9134E-05,-.1020E-04, .3524E-04,-.2599E-04,
     +-.2093E-04, .1675E-04,-.1313E-03,-.9172E-04,-.5252E-04,-.4567E-04,
     + .4186E-05,-.3920E-05, .2552E-04,-.2059E-04,-.2246E-04,-.1028E-04,
     +-.1324E-03,-.9210E-04,-.5138E-04,-.4491E-04, .6470E-05,-.2131E-05,
     + .1496E-04,-.1572E-04,-.3311E-04,-.8754E-05,-.1324E-03,-.9058E-04,
     +-.5328E-04,-.4225E-04, .1827E-05,-.8411E-06, .4719E-05,-.6813E-05,
     +-.2474E-04,-.1256E-04,-.1340E-03,-.8868E-04,-.5633E-04,-.4187E-04,
     +-.4415E-05, .6055E-05,-.1648E-04,-.1507E-04, .1979E-04,-.2131E-04,
     +-.1340E-03,-.8373E-04,-.5899E-04,-.3920E-04,-.4072E-05, .1491E-04,
     +-.9781E-05,-.5328E-05, .3578E-04,-.1979E-04,-.1321E-03,-.7954E-04,
     +-.5899E-04,-.4072E-04, .1066E-05, .5728E-05,-.5138E-05,-.8373E-05,
     + .2626E-04,-.2436E-04,-.1363E-03,-.6432E-04,-.5176E-04,-.3083E-04,
     + .2169E-05,-.8944E-05, .3159E-05, .6470E-05,-.4187E-05, .4948E-05,
     +-.1302E-03,-.7802E-04,-.3311E-04,-.1903E-04, .5328E-05,-.1884E-04,
     + .1408E-04, .3311E-04, .1142E-05,-.7613E-06,-.1473E-03,-.6737E-04,
     +-.7536E-04,-.1085E-04,-.1903E-05,-.1458E-04, .4034E-04,-.3941E-10,
     +-.7992E-05, .2664E-05,-.1361E-03,-.5709E-04,-.8550E-04,-.5709E-05,
     +-.8640E-05, .6523E-05, .1903E-05,-.8221E-05,-.3045E-05,-.9134E-05,
     +-.1329E-03,-.5529E-04,-.7107E-04, .2664E-05,-.9020E-05, .3320E-04,
     +-.2131E-05,-.4187E-05,-.7231E-05,-.3806E-05,-.1278E-03,-.5247E-04,
     +-.6465E-04, .3806E-05,-.6091E-05, .1245E-04,-.3844E-05,-.6090E-05,
     +-.8754E-05,-.2664E-05,-.1321E-03,-.5632E-04,-.5897E-04, .1012E-04,
     + .1168E-04,-.4196E-06,-.8411E-05,-.8868E-05,-.1484E-04,-.1522E-05,
     +-.1252E-03,-.4907E-04,-.5932E-04, .3245E-04, .1996E-04,-.3325E-05,
     +-.5785E-05,-.6394E-05,-.6851E-05,-.1142E-05,-.1093E-03,-.4731E-04,
     +-.6761E-04, .1808E-04, .1754E-04,-.5079E-05,-.5809E-05,-.5649E-05,
     +-.3988E-05,-.5849E-06,-.1151E-03,-.4965E-04,-.7163E-04, .7839E-05,
     + .5505E-05,-.6084E-05,-.3344E-05,-.3894E-05,-.1391E-05,-.1327E-05/
	data ( ( ( coehcb(k,j,i), i = 1, 10 ), j = 1, 19 ), k = 1, 3 ) /
     +-.9398E+01,-.5678E+01,-.3606E+01,-.2192E+01, .2104E+01, .3044E+01,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.9094E+01,-.5422E+01,
     +-.3448E+01,-.1650E+01, .2046E+01, .2749E+01,-.4587E+02,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.8760E+01,-.5270E+01,-.3329E+01,-.1147E+01,
     + .2112E+01, .2709E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.8537E+01,-.5152E+01,-.3129E+01,-.9544E+00, .2254E+01, .2771E+01,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.8176E+01,-.4936E+01,
     +-.2680E+01,-.9259E+00, .2247E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.7836E+01,-.4676E+01,-.2378E+01,-.3550E+00,
     + .1396E+01, .1976E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.7419E+01,-.4122E+01,-.2407E+01,-.1204E-01, .1744E+01,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.7124E+01,-.3727E+01,
     +-.2160E+01, .6158E+00, .1953E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.6823E+01,-.3324E+01,-.1748E+01,-.9806E-01,
     + .2319E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.5957E+01,-.3017E+01,-.1647E+01, .1398E+01,-.4587E+02,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.5115E+01,-.2290E+01,
     +-.5273E+00, .5662E+00, .1459E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.4162E+01,-.1453E+01, .1116E+00,-.4587E+02,
     + .9569E+00,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.3611E+01,-.9744E+00,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.3075E+01,-.4176E+00,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.3469E+01,-.9395E+00, .5092E+00, .6200E+00,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.3808E+01,-.1505E+01, .3901E+00, .6264E+00,-.1155E+01,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4058E+01,-.1818E+01,
     + .2693E+00, .7087E+00, .3820E+00,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.4587E+02,-.4587E+02,-.4262E+01,-.2097E+01,-.5711E-01, .5681E+00,
     + .1310E+01, .7371E+00,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.3997E+01,-.1784E+01, .4388E-01, .5167E+00, .6930E+00,-.6906E+00,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02, .2944E-01, .2723E-01,
     + .1854E-01, .2023E-01, .2254E-01, .3059E-02, .4788E+00, .3059E-02,
     + .3059E-02, .3059E-02, .3080E-01, .2549E-01, .1547E-01, .2225E-01,
     + .2107E-01, .3059E-02, .4737E+00, .3059E-02, .3059E-02, .3059E-02,
     + .3269E-01, .2656E-01, .2125E-01, .2179E-01, .2162E-01, .4589E+00,
     + .4643E+00, .3059E-02, .3059E-02, .3059E-02, .3322E-01, .2476E-01,
     + .2075E-01, .2139E-01, .1907E-01, .4501E+00, .4441E+00, .3059E-02,
     + .3059E-02, .3059E-02, .3387E-01, .2182E-01, .2665E-01, .1841E-01,
     + .2506E-01, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .3532E-01, .2091E-01, .1995E-01, .2067E-01, .1949E-01, .4491E+00,
     + .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3468E-01, .2075E-01,
     + .2587E-01, .1401E-01, .8646E-02, .3059E-02, .3059E-02, .3059E-02,
     + .3059E-02, .3059E-02, .3666E-01, .2430E-01, .1919E-01, .2007E-01,
     + .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .3613E-01, .2147E-01, .1892E-01, .1361E-01, .3059E-02, .4506E+00,
     + .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3129E-01, .1954E-01,
     + .2442E-01, .1011E-01, .4420E+00, .3059E-02, .3059E-02, .3059E-02,
     + .3059E-02, .3059E-02, .3177E-01, .2101E-01, .1526E-01, .4376E+00,
     + .4379E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .2887E-01, .2044E-01, .1285E-01, .3059E-02,-.4862E-03, .3059E-02,
     + .3059E-02, .3059E-02, .3059E-02, .3059E-02, .2759E-01, .2114E-01,
     + .4303E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .3059E-02, .3059E-02, .2880E-01, .1690E-01,-.4187E+00, .3059E-02,
     + .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .2852E-01, .2255E-01, .2184E-01, .4334E+00, .4217E+00, .3059E-02,
     + .3059E-02, .3059E-02, .3059E-02, .3059E-02, .2840E-01, .2136E-01,
     + .1644E-01, .2812E-01, .4358E+00, .4288E+00, .3059E-02, .3059E-02,
     + .3059E-02, .3059E-02, .2809E-01, .2173E-01, .1708E-01, .3346E-01,
     + .4225E-01, .4419E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .2702E-01, .2260E-01, .1607E-01, .2720E-01, .3982E-01, .4452E+00,
     + .4365E+00, .4345E+00, .4432E+00, .4623E+00, .2684E-01, .2328E-01,
     + .2099E-01, .3040E-01, .3867E-01, .4389E+00, .3132E-01, .3158E-01,
     + .4083E-01, .4580E+00,-.1581E-03,-.9707E-04,-.1250E-03, .2580E-03,
     + .7378E-04,-.1617E-01, .8646E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.1319E-03,-.9528E-04,-.1710E-03, .7118E-04, .2076E-04,-.1608E-01,
     + .8552E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.1721E-03,-.4680E-04,
     +-.5522E-04,-.6242E-04, .4517E-04,-.7777E-02, .8382E-02,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.1482E-03,-.4208E-04,-.5216E-04,-.6514E-04,
     +-.8378E-04,-.7956E-02, .8013E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.1501E-03,-.4002E-04,-.1664E-03, .2272E-04,-.1888E-03,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1201E-03,-.4709E-04,
     +-.5371E-04,-.1574E-03, .1854E-03,-.7712E-02,-.4656E-05,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.1333E-03,-.1062E-03, .5785E-04,-.4150E-04,
     +-.5717E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.1212E-03,-.8524E-04,-.5895E-04,-.2884E-03,-.1581E-01,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.8148E-04,-.9361E-04,
     +-.2873E-03, .1883E-03,-.1594E-01, .8133E-02,-.4656E-05,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.1221E-03,-.1430E-04, .6335E-04,-.2581E-03,
     + .7977E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.9257E-04,-.5008E-04, .6389E-04,-.7455E-02,-.7745E-02,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1186E-03,-.9037E-04,
     +-.7461E-04,-.4656E-05, .1168E-03,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.8513E-04,-.5708E-04, .7763E-02,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.1124E-03,-.1228E-03, .7663E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.1015E-03,-.8369E-04,
     +-.2167E-03,-.7548E-02, .7608E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.4656E-05,-.4656E-05,-.1049E-03,-.6414E-04,-.1384E-03,-.1644E-03,
     +-.6919E-02, .7736E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.1008E-03,-.7047E-04,-.1276E-03,-.2445E-03,-.1860E-03, .7975E-02,
     +-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.9629E-04,-.1007E-03,
     +-.1127E-03,-.1527E-03,-.3238E-03,-.7373E-02, .7877E-02, .7840E-02,
     + .7997E-02, .8345E-02,-.8800E-04,-.1072E-03,-.1046E-03,-.1777E-03,
     +-.2146E-03,-.7016E-02, .1516E-01, .1532E-01, .1509E-01, .8268E-02/
	end

	block data ckd15
c **********************************************************************
c hk is the interval in the g (cumulative probability) space from 0
c to one. coehca and coehcb are the coefficients to calculate the
c H2O and CO2 overlapping absorption coefficients in units of (cm-
c atm)**-1 at three temperatures, nineteen pressures, and 12 cumu-
c lative probabilities (Fu, 1991). The spectral region is from 670
c to 540 cm**-1.
c **********************************************************************
	common /band15/ hk(12), coehca(3,19,12), coehcb(3,19,12)
 	data hk /.24,.36,.18,.1,.05,.02,.016,.012,.01,.006,.0039,.0021/
	data ( ( ( coehca(k,j,i), i = 1, 12 ), j = 1, 19 ), k = 1, 2 ) /
     +-.1921E+02,-.1363E+02,-.1080E+02,-.8392E+01,-.6776E+01,-.5696E+01,
     +-.4572E+01,-.3752E+01,-.2382E+01,-.1110E+01, .6803E+00, .3259E+01,
     +-.1875E+02,-.1321E+02,-.1040E+02,-.8026E+01,-.6449E+01,-.5401E+01,
     +-.4316E+01,-.3498E+01,-.2141E+01,-.9439E+00, .8103E+00, .3314E+01,
     +-.1829E+02,-.1278E+02,-.1000E+02,-.7646E+01,-.6089E+01,-.5085E+01,
     +-.4047E+01,-.3217E+01,-.1872E+01,-.7106E+00, .9573E+00, .3390E+01,
     +-.1783E+02,-.1236E+02,-.9596E+01,-.7264E+01,-.5735E+01,-.4740E+01,
     +-.3743E+01,-.2882E+01,-.1587E+01,-.4714E+00, .1120E+01, .3425E+01,
     +-.1737E+02,-.1195E+02,-.9193E+01,-.6877E+01,-.5371E+01,-.4404E+01,
     +-.3405E+01,-.2574E+01,-.1298E+01,-.1747E+00, .1327E+01, .3547E+01,
     +-.1691E+02,-.1153E+02,-.8776E+01,-.6490E+01,-.4993E+01,-.4049E+01,
     +-.3039E+01,-.2256E+01,-.1012E+01, .1103E+00, .1530E+01, .3651E+01,
     +-.1644E+02,-.1112E+02,-.8360E+01,-.6105E+01,-.4623E+01,-.3688E+01,
     +-.2694E+01,-.1915E+01,-.6855E+00, .3993E+00, .1714E+01, .3950E+01,
     +-.1598E+02,-.1073E+02,-.7943E+01,-.5723E+01,-.4236E+01,-.3314E+01,
     +-.2338E+01,-.1596E+01,-.3583E+00, .6963E+00, .1868E+01, .4127E+01,
     +-.1553E+02,-.1034E+02,-.7542E+01,-.5357E+01,-.3856E+01,-.2942E+01,
     +-.1986E+01,-.1299E+01,-.5472E-01, .9443E+00, .2149E+01, .4261E+01,
     +-.1485E+02,-.9661E+01,-.7008E+01,-.4830E+01,-.3458E+01,-.2566E+01,
     +-.1658E+01,-.9639E+00, .2083E+00, .1182E+01, .2458E+01, .4452E+01,
     +-.1427E+02,-.9166E+01,-.6373E+01,-.4404E+01,-.3073E+01,-.2209E+01,
     +-.1349E+01,-.6648E+00, .4023E+00, .1452E+01, .2739E+01, .4466E+01,
     +-.1380E+02,-.8726E+01,-.5772E+01,-.3982E+01,-.2732E+01,-.1874E+01,
     +-.1052E+01,-.4403E+00, .5763E+00, .1792E+01, .2999E+01, .4335E+01,
     +-.1305E+02,-.8270E+01,-.5304E+01,-.3586E+01,-.2392E+01,-.1568E+01,
     +-.8299E+00,-.2650E+00, .8584E+00, .2062E+01, .3141E+01, .4168E+01,
     +-.1269E+02,-.7900E+01,-.4956E+01,-.3205E+01,-.2065E+01,-.1332E+01,
     +-.6415E+00,-.7921E-01, .1170E+01, .2269E+01, .3198E+01, .4066E+01,
     +-.1227E+02,-.7536E+01,-.4576E+01,-.2859E+01,-.1815E+01,-.1139E+01,
     +-.4520E+00, .2272E+00, .1371E+01, .2351E+01, .3150E+01, .3935E+01,
     +-.1186E+02,-.7159E+01,-.4223E+01,-.2538E+01,-.1619E+01,-.9324E+00,
     +-.1566E+00, .5151E+00, .1520E+01, .2339E+01, .3132E+01, .3880E+01,
     +-.1120E+02,-.6777E+01,-.3919E+01,-.2330E+01,-.1387E+01,-.6737E+00,
     + .1108E+00, .6991E+00, .1531E+01, .2163E+01, .3150E+01, .3767E+01,
     +-.9973E+01,-.6279E+01,-.3638E+01,-.2048E+01,-.1098E+01,-.4407E+00,
     + .3043E+00, .7797E+00, .1424E+01, .2002E+01, .3122E+01, .3611E+01,
     +-.8483E+01,-.5607E+01,-.3357E+01,-.1744E+01,-.8884E+00,-.2264E+00,
     + .3800E+00, .7504E+00, .1245E+01, .2032E+01, .3097E+01, .3546E+01,
     + .3762E-01, .2372E-01, .1643E-01, .1208E-01, .1170E-01, .1164E-01,
     + .1214E-01, .1161E-01, .1028E-01, .9185E-02, .7712E-02, .1001E-01,
     + .3762E-01, .2382E-01, .1593E-01, .1145E-01, .1059E-01, .1049E-01,
     + .1080E-01, .1057E-01, .8894E-02, .7807E-02, .7132E-02, .1032E-01,
     + .3764E-01, .2386E-01, .1555E-01, .1080E-01, .9692E-02, .9231E-02,
     + .9585E-02, .9644E-02, .7711E-02, .6443E-02, .6223E-02, .9922E-02,
     + .3764E-01, .2395E-01, .1516E-01, .1028E-01, .8917E-02, .8415E-02,
     + .8457E-02, .8777E-02, .6436E-02, .5428E-02, .5499E-02, .8017E-02,
     + .3768E-01, .2399E-01, .1482E-01, .9692E-02, .8247E-02, .7640E-02,
     + .7582E-02, .7783E-02, .5432E-02, .4482E-02, .4919E-02, .5903E-02,
     + .3770E-01, .2401E-01, .1449E-01, .9252E-02, .7620E-02, .6678E-02,
     + .6845E-02, .6925E-02, .4939E-02, .3471E-02, .4124E-02, .3873E-02,
     + .3776E-01, .2395E-01, .1419E-01, .8959E-02, .7096E-02, .6184E-02,
     + .6110E-02, .6075E-02, .4419E-02, .2891E-02, .3056E-02, .1214E-02,
     + .3780E-01, .2391E-01, .1392E-01, .8687E-02, .6573E-02, .5733E-02,
     + .5359E-02, .5009E-02, .4034E-02, .2755E-02, .1968E-02,-.4187E-04,
     + .3791E-01, .2382E-01, .1373E-01, .8561E-02, .6060E-02, .5120E-02,
     + .4618E-02, .4713E-02, .3965E-02, .2481E-02, .8164E-03,-.1088E-02,
     + .3843E-01, .2148E-01, .1302E-01, .6384E-02, .5256E-02, .4260E-02,
     + .4077E-02, .4181E-02, .4132E-02, .2135E-02,-.2931E-03,-.1151E-02,
     + .3896E-01, .2081E-01, .1097E-01, .5568E-02, .4475E-02, .3795E-02,
     + .3828E-02, .3996E-02, .3766E-02, .1193E-02,-.1089E-02,-.9420E-03,
     + .3973E-01, .2024E-01, .9943E-02, .4815E-02, .3820E-02, .3663E-02,
     + .3568E-02, .3881E-02, .2859E-02, .6698E-03,-.1549E-02,-.6280E-03,
     + .3635E-01, .1963E-01, .1061E-01, .3812E-02, .3509E-02, .3429E-02,
     + .3693E-02, .3316E-02, .1120E-02, .6552E-03,-.1193E-02,-.1109E-02,
     + .3631E-01, .1893E-01, .1056E-01, .3172E-02, .3378E-02, .3164E-02,
     + .2751E-02, .1722E-02, .1112E-02, .4354E-03,-.7327E-03,-.1319E-02,
     + .3500E-01, .1828E-01, .1050E-01, .2831E-02, .2784E-02, .2564E-02,
     + .1469E-02, .7739E-03, .1209E-02, .7913E-03,-.2512E-03,-.1758E-02,
     + .3352E-01, .1763E-01, .1045E-01, .2401E-02, .1928E-02, .1340E-02,
     + .3753E-03, .5794E-03, .9060E-03, .1042E-02, .1465E-03,-.2533E-02,
     + .2880E-01, .1729E-01, .1077E-01, .1347E-02, .1194E-02,-.1191E-03,
     + .2828E-03, .6606E-03, .9743E-03, .1002E-02, .0000E+00,-.3140E-02,
     + .2040E-01, .1585E-01, .1165E-01, .3871E-05, .1509E-04,-.1046E-02,
     + .2444E-03, .4359E-03, .1041E-02, .2429E-02,-.1721E-03,-.2786E-02,
     + .1737E-01, .1560E-01, .1240E-01,-.2139E-03,-.1025E-02,-.1248E-02,
     +-.6934E-04, .1649E-03, .4062E-03, .1554E-02,-.4179E-03,-.7795E-03/
	data ( ( ( coehca(k,j,i), i = 1, 12 ), j = 1, 19 ), k = 3, 3 ) /
     +-.1488E-03,-.9248E-04,-.2322E-04,-.4187E-05, .1104E-04, .9895E-05,
     +-.2283E-05, .2512E-05,-.9058E-05, .8449E-05, .8297E-05,-.3882E-04,
     +-.1488E-03,-.9058E-04,-.2398E-04,-.5709E-05, .1218E-04, .1180E-04,
     + .1522E-05, .6927E-05,-.1161E-04, .1714E-04,-.4948E-06,-.3540E-04,
     +-.1500E-03,-.8830E-04,-.2474E-04,-.8373E-05, .6470E-05, .7992E-05,
     + .9096E-05, .6737E-05,-.1485E-04, .1873E-04,-.4948E-06,-.4491E-04,
     +-.1500E-03,-.8601E-04,-.2664E-04,-.1028E-04, .6851E-05, .6851E-05,
     + .1294E-04,-.2550E-05,-.1520E-04, .2310E-04, .4948E-06,-.2017E-04,
     +-.1507E-03,-.8373E-04,-.2664E-04,-.1256E-04, .4567E-05, .1028E-04,
     + .9210E-05,-.2131E-05,-.6995E-05, .7498E-05,-.1104E-04,-.2284E-05,
     +-.1519E-03,-.8183E-04,-.2816E-04,-.1142E-04, .7611E-06, .7231E-05,
     + .1751E-05,-.7612E-06, .8312E-05, .2436E-05,-.7231E-05, .2398E-04,
     +-.1530E-03,-.7992E-04,-.2893E-04,-.9896E-05, .3806E-06, .8906E-05,
     + .3159E-05,-.5328E-05, .3692E-05,-.2093E-05,-.6851E-05,-.3045E-05,
     +-.1538E-03,-.7536E-04,-.3007E-04,-.8754E-05,-.3045E-05, .5138E-05,
     + .9134E-06,-.1979E-06, .1560E-05,-.1507E-04, .2284E-04, .9895E-05,
     +-.1541E-03,-.7688E-04,-.2969E-04,-.5709E-05,-.3996E-05, .1142E-05,
     +-.8373E-06, .1235E-04,-.7079E-05,-.6737E-05, .1028E-04, .3578E-04,
     +-.1560E-03,-.6851E-04,-.1903E-04,-.4187E-05,-.4605E-05,-.1142E-06,
     + .3878E-05, .3597E-05,-.9591E-05, .5328E-05, .7612E-05,-.4948E-05,
     +-.1587E-03,-.6546E-04,-.2740E-04,-.7612E-06,-.3578E-05, .1713E-05,
     + .6064E-05,-.9781E-05, .1408E-05, .5709E-05, .8373E-05,-.1256E-04,
     +-.1484E-03,-.5823E-04,-.4301E-04,-.1522E-05, .7498E-05,-.5328E-06,
     +-.7855E-05,-.1599E-05, .1964E-04,-.2284E-05, .7882E-10, .5328E-05,
     +-.1238E-03,-.5700E-04,-.5266E-04, .3286E-05, .4910E-05,-.8602E-05,
     + .6090E-06, .8454E-05, .1256E-05,-.4072E-05,-.1903E-05, .6470E-05,
     +-.1155E-03,-.5231E-04,-.4396E-04, .3626E-05,-.7051E-05,-.1743E-05,
     + .9667E-05, .2064E-04,-.2778E-05,-.6546E-05,-.4948E-05, .1903E-05,
     +-.1024E-03,-.5129E-04,-.4506E-04, .7943E-06, .3074E-06, .3243E-05,
     + .2754E-04,-.1479E-05, .1661E-05,-.2969E-05,-.1066E-04, .7612E-06,
     +-.8473E-04,-.5418E-04,-.4674E-04,-.3418E-05, .9460E-05, .1151E-04,
     + .5714E-05,-.1069E-04,-.2022E-05,-.9061E-05,-.1104E-04,-.3083E-04,
     +-.4283E-04,-.5037E-04,-.4476E-04, .1951E-04, .8922E-05, .1296E-04,
     +-.4053E-05,-.4355E-05,-.2355E-05,-.5004E-05,-.1218E-04,-.1522E-04,
     + .6411E-05,-.5937E-04,-.5331E-04, .1934E-04, .5284E-05, .1129E-04,
     +-.2166E-05,-.1484E-06,-.5407E-05,-.1364E-04,-.3115E-05, .3004E-04,
     +-.5074E-04,-.6256E-04,-.5097E-04, .2218E-04, .1228E-04,-.1160E-05,
     +-.1105E-05, .1618E-06,-.6089E-05,-.4216E-06,-.5314E-05, .7903E-05/
	data ( ( ( coehcb(k,j,i), i = 1, 12 ), j = 1, 19 ), k = 1, 2 ) /
     +-.9593E+01,-.4078E+01,-.2812E+01,-.6506E+00,-.4123E+00, .2055E+01,
     + .4097E+01, .4671E+01, .4639E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.9276E+01,-.3757E+01,-.2467E+01,-.5784E+00, .8833E-01, .2232E+01,
     + .3826E+01, .4723E+01, .4942E+01, .5135E+01,-.4587E+02,-.4587E+02,
     +-.8968E+01,-.3508E+01,-.2116E+01,-.1363E+00, .1662E+00, .2424E+01,
     + .4220E+01, .4513E+01, .1375E+01, .4601E+01,-.4587E+02,-.4587E+02,
     +-.8662E+01,-.3164E+01,-.1722E+01, .5178E-01, .7288E+00, .2411E+01,
     + .3805E+01, .4766E+01, .4342E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.8292E+01,-.2799E+01,-.1359E+01, .3271E+00, .1650E+01, .2395E+01,
     + .4192E+01, .4758E+01, .2470E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.7812E+01,-.2404E+01,-.1085E+01, .7167E+00, .2202E+01, .2922E+01,
     + .4322E+01, .4591E+01, .4186E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.7441E+01,-.2066E+01,-.7142E+00, .1057E+01, .2524E+01, .2946E+01,
     + .4220E+01, .3607E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.7191E+01,-.1745E+01,-.3487E+00, .1453E+01, .2739E+01, .3660E+01,
     + .4114E+01, .3245E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.6895E+01,-.1326E+01,-.3500E+00, .1647E+01, .2899E+01, .4023E+01,
     + .3361E+01, .3360E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.5876E+01,-.9573E+00, .2014E+00, .2130E+01, .3493E+01, .4088E+01,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.4429E+01,-.3417E+00, .1204E+01, .2780E+01, .3843E+01, .3099E+01,
     +-.4587E+02, .3605E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.3122E+01, .2697E+00, .1866E+01, .3526E+01, .3569E+01, .1025E+01,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.2284E+01, .8186E+00, .2754E+01, .3206E+01, .3704E+01,-.4587E+02,
     +-.4587E+02, .4625E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.1711E+01, .1220E+01, .3248E+01,-.4587E+02, .2565E+01, .3297E+01,
     +-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.1758E+01, .7970E+00, .2758E+01, .2926E+01, .2613E+01, .1974E+01,
     +-.4587E+02, .2310E+01,-.4587E+02,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.1737E+01, .3499E+00, .2246E+01, .2673E+01, .3308E+01, .3463E+01,
     + .3103E+01, .2611E+01, .2178E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.1559E+01, .2215E+00, .1875E+01, .2500E+01, .3346E+01, .3585E+01,
     + .3946E+01, .3533E+01, .3205E+01,-.4587E+02,-.4587E+02,-.4587E+02,
     +-.1601E+01, .5060E-01, .1275E+01, .2176E+01, .3081E+01, .3649E+01,
     + .3940E+01, .4106E+01, .4112E+01, .4349E+01, .2292E+01,-.4587E+02,
     +-.1222E+01, .3199E+00, .1642E+01, .2380E+01, .3254E+01, .3534E+01,
     + .3687E+01, .3717E+01, .3402E+01, .3868E+01,-.4587E+02,-.4587E+02,
     + .2967E-01, .1697E-01, .1795E-01, .1387E-01, .2032E-01, .1187E-01,
     + .2560E-01, .1044E-01,-.4560E+00, .3059E-02, .3059E-02, .3059E-02,
     + .2998E-01, .1586E-01, .1786E-01, .1521E-01, .1710E-01, .1061E-01,
     + .2030E-01, .1158E-01, .4452E+00, .3059E-02, .3059E-02, .3059E-02,
     + .2993E-01, .1551E-01, .1481E-01, .9846E-02, .2443E-01, .1150E-01,
     + .1865E-01, .1376E-01, .4617E+00, .3059E-02, .3059E-02, .3059E-02,
     + .3035E-01, .1417E-01, .1438E-01, .1511E-01, .1901E-01, .8582E-02,
     + .1746E-01, .1450E-01, .4523E+00, .3059E-02, .3059E-02, .3059E-02,
     + .2970E-01, .1347E-01, .1322E-01, .1252E-01, .1665E-01, .1037E-01,
     + .1320E-01, .1199E-01, .4436E+00, .3059E-02, .3059E-02, .3059E-02,
     + .2949E-01, .1291E-01, .1671E-01, .1111E-01, .1400E-01, .1318E-01,
     + .1060E-01, .1046E-01, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .3004E-01, .1300E-01, .1413E-01, .9085E-02, .9764E-02, .2260E-01,
     + .9778E-02, .4671E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .3086E-01, .1436E-01, .1205E-01, .1081E-01, .4681E-02, .1479E-01,
     + .1888E-01, .3494E-01, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .3094E-01, .1500E-01, .1457E-01, .1060E-01, .8319E-02, .8983E-02,
     + .3791E-01, .2232E-01, .4631E+00, .3059E-02, .3059E-02, .3059E-02,
     + .3158E-01, .1585E-01, .1292E-01, .6531E-02, .1383E-01, .4605E+00,
     + .4662E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .3182E-01, .1586E-01, .8724E-02, .5798E-02, .2454E-01, .4607E+00,
     + .4560E+00, .4511E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .2369E-01, .1606E-01, .5477E-02, .1228E-01, .4579E+00, .4561E+00,
     + .4497E+00, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .2190E-01, .1779E-01, .6267E-02, .4535E+00, .4533E+00, .3059E-02,
     + .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .2100E-01, .1653E-01, .7449E-02, .4543E+00, .4472E+00, .4439E+00,
     + .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02, .3059E-02,
     + .1864E-01, .1771E-01, .7040E-02, .2877E-01, .3381E-01, .2691E-01,
     + .4466E+00, .3059E-02, .4613E+00, .3059E-02, .3059E-02, .3059E-02,
     + .1637E-01, .1641E-01, .8424E-02, .1318E-01, .2060E-01, .3426E-01,
     + .4122E-01, .4621E+00, .4555E+00, .4525E+00, .3059E-02, .3059E-02,
     + .1607E-01, .1452E-01, .8013E-02, .1213E-01, .1482E-01, .2125E-01,
     + .3379E-01, .3562E-01, .4619E+00, .4569E+00, .3059E-02, .3059E-02,
     + .1698E-01, .1538E-01, .6616E-02, .1147E-01, .1217E-01, .1696E-01,
     + .1871E-01, .2273E-01, .4513E-01, .4702E+00, .4617E+00, .4553E+00,
     + .1700E-01, .1547E-01, .6456E-02, .1324E-01, .1502E-01, .2095E-01,
     + .2547E-01, .2823E-01, .4107E-01, .4676E+00, .4583E+00, .4498E+00/
	data ( ( ( coehcb(k,j,i), i = 1, 12 ), j = 1, 19 ), k = 3, 3 ) /
     +-.6747E-05,-.2483E-04, .6575E-04, .1026E-03, .3888E-03,-.8519E-04,
     +-.1629E-03,-.1808E-04,-.8355E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.2270E-04,-.3427E-04, .5118E-04, .1218E-03, .1245E-03,-.1245E-03,
     + .3841E-05,-.4151E-04,-.8763E-02,-.1687E-01,-.4656E-05,-.4656E-05,
     +-.4557E-04,-.3023E-04, .2286E-04, .5656E-04, .4113E-04,-.1407E-03,
     +-.1301E-03, .8503E-04,-.7284E-02,-.1669E-01,-.4656E-05,-.4656E-05,
     +-.5325E-04,-.5309E-04,-.1246E-04, .2244E-04, .5136E-04,-.1272E-03,
     + .4217E-04,-.1749E-04,-.8435E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.6857E-04,-.7217E-04, .1740E-05, .3653E-04,-.1490E-03,-.4090E-04,
     +-.2376E-04, .2047E-04,-.7974E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.1232E-03,-.9826E-04,-.2849E-04, .1703E-04,-.1895E-03,-.3363E-03,
     + .7102E-04,-.1838E-05,-.1655E-01,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.9896E-04,-.5127E-04,-.2704E-04,-.1218E-04,-.1207E-03,-.5883E-04,
     + .6893E-04,-.7924E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.7837E-04,-.4980E-04, .6902E-05,-.1072E-03,-.4051E-04,-.1991E-05,
     +-.1173E-03,-.5195E-04,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.8136E-04,-.8102E-04, .1254E-03,-.4658E-04, .3173E-04,-.4461E-05,
     +-.1558E-03,-.2036E-03, .8360E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.2232E-04,-.6411E-04, .9486E-04,-.2322E-03,-.8282E-04,-.8202E-02,
     + .8416E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.1398E-03,-.7165E-04,-.4258E-04,-.3970E-04,-.2839E-03,-.7873E-02,
     + .8231E-02,-.8213E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.6754E-04,-.7469E-04,-.6898E-04,-.1702E-03,-.8079E-02,-.7270E-02,
     + .8116E-02,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.2396E-04,-.2361E-04,-.8664E-04,-.8038E-02,-.8207E-02,-.4656E-05,
     +-.4656E-05,-.1670E-01,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.5479E-04,-.7593E-04,-.1005E-03, .8199E-02,-.7942E-02,-.8244E-02,
     +-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.3806E-04,-.5825E-04,-.1003E-03,-.2925E-03,-.1506E-03, .3148E-04,
     + .8060E-02,-.1593E-01, .8327E-02,-.4656E-05,-.4656E-05,-.4656E-05,
     +-.4706E-04,-.3630E-04,-.7811E-04,-.6881E-04,-.1822E-03,-.3091E-03,
     +-.3033E-03,-.7684E-02,-.7663E-02, .8167E-02,-.4656E-05,-.4656E-05,
     +-.7669E-04,-.4610E-04,-.8063E-04,-.7250E-04,-.1094E-03,-.1241E-03,
     +-.2944E-03,-.1736E-03,-.7886E-02, .8248E-02,-.4656E-05,-.4656E-05,
     +-.7138E-04,-.4545E-04,-.3653E-04,-.6075E-04,-.4528E-04,-.1077E-03,
     +-.1119E-03,-.1657E-03,-.4695E-03,-.8112E-02,-.7587E-02, .8217E-02,
     +-.6812E-04,-.4558E-04,-.6739E-04,-.8861E-04,-.9386E-04,-.1334E-03,
     +-.2007E-03,-.2179E-03,-.1650E-03,-.8001E-02, .8273E-02, .8118E-02/
	end

	block data ckd16
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and  seven cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 540 to 400 cm**-1.
c *********************************************************************
	common /band16/ hk(7), coeh2o(3,19,7)
 	data hk / .12, .24, .24, .20, .12, .06, .02 /
	data ( ( ( coeh2o(k,j,i), i = 1, 7 ), j = 1, 19 ), k = 1, 3 ) /
     +-.2344E+02,-.2016E+02,-.1986E+02,-.1655E+02,-.1243E+02,-.8437E+01,
     +-.4858E+01,-.2298E+02,-.2014E+02,-.1984E+02,-.1609E+02,-.1198E+02,
     +-.8020E+01,-.4548E+01,-.2252E+02,-.2012E+02,-.1981E+02,-.1564E+02,
     +-.1153E+02,-.7596E+01,-.4239E+01,-.2206E+02,-.2009E+02,-.1957E+02,
     +-.1517E+02,-.1111E+02,-.7161E+01,-.3871E+01,-.2160E+02,-.2007E+02,
     +-.1911E+02,-.1472E+02,-.1065E+02,-.6721E+01,-.3479E+01,-.2113E+02,
     +-.2005E+02,-.1865E+02,-.1426E+02,-.1021E+02,-.6302E+01,-.3081E+01,
     +-.2067E+02,-.2003E+02,-.1819E+02,-.1379E+02,-.9765E+01,-.5883E+01,
     +-.2678E+01,-.2026E+02,-.2001E+02,-.1773E+02,-.1333E+02,-.9332E+01,
     +-.5443E+01,-.2253E+01,-.2024E+02,-.1999E+02,-.1727E+02,-.1288E+02,
     +-.8897E+01,-.5029E+01,-.1858E+01,-.2026E+02,-.1959E+02,-.1481E+02,
     +-.1147E+02,-.7477E+01,-.4555E+01,-.1464E+01,-.2022E+02,-.1632E+02,
     +-.1305E+02,-.9885E+01,-.6689E+01,-.4108E+01,-.1068E+01,-.1936E+02,
     +-.1438E+02,-.1163E+02,-.8499E+01,-.6146E+01,-.3673E+01,-.6816E+00,
     +-.1675E+02,-.1281E+02,-.1020E+02,-.7716E+01,-.5678E+01,-.3256E+01,
     +-.3125E+00,-.1510E+02,-.1124E+02,-.8821E+01,-.7140E+01,-.5243E+01,
     +-.2851E+01,-.2560E-01,-.1334E+02,-.9708E+01,-.8061E+01,-.6611E+01,
     +-.4842E+01,-.2459E+01, .1711E+00,-.1155E+02,-.8798E+01,-.7440E+01,
     +-.6123E+01,-.4439E+01,-.2089E+01, .2480E+00,-.1020E+02,-.8154E+01,
     +-.6945E+01,-.5681E+01,-.4055E+01,-.1737E+01, .2390E+00,-.9464E+01,
     +-.7677E+01,-.6512E+01,-.5284E+01,-.3707E+01,-.1453E+01, .2015E+00,
     +-.9033E+01,-.7246E+01,-.6093E+01,-.4882E+01,-.3346E+01,-.1264E+01,
     + .1033E+00, .4658E-01, .5840E-02, .4626E-02, .2688E-01, .2395E-01,
     + .1804E-01, .2074E-01, .4660E-01, .1884E-02, .8561E-02, .2690E-01,
     + .2403E-01, .1788E-01, .1934E-01, .4660E-01, .1800E-02, .1252E-01,
     + .2694E-01, .2393E-01, .1786E-01, .1825E-01, .4660E-01, .1779E-02,
     + .1649E-01, .2696E-01, .2397E-01, .1779E-01, .1765E-01, .4348E-01,
     + .1758E-02, .2043E-01, .2696E-01, .2393E-01, .1748E-01, .1675E-01,
     + .3944E-01, .1737E-02, .2445E-01, .2698E-01, .2384E-01, .1752E-01,
     + .1549E-01, .3538E-01, .1654E-02, .2847E-01, .2702E-01, .2384E-01,
     + .1714E-01, .1565E-01, .3127E-01, .1570E-02, .3245E-01, .2705E-01,
     + .2374E-01, .1712E-01, .1514E-01, .2715E-01, .1444E-02, .3540E-01,
     + .2711E-01, .2363E-01, .1702E-01, .1446E-01, .2960E-01, .1760E-01,
     + .2977E-01, .2397E-01, .2087E-01, .1618E-01, .1445E-01, .2466E-01,
     + .3039E-01, .2428E-01, .2217E-01, .1821E-01, .1593E-01, .1463E-01,
     + .2640E-01, .2545E-01, .2231E-01, .2060E-01, .1773E-01, .1555E-01,
     + .1473E-01, .3456E-01, .2135E-01, .2030E-01, .1844E-01, .1740E-01,
     + .1559E-01, .1428E-01, .3203E-01, .2047E-01, .1809E-01, .1760E-01,
     + .1725E-01, .1545E-01, .1541E-01, .2137E-01, .1857E-01, .1616E-01,
     + .1698E-01, .1700E-01, .1537E-01, .1636E-01, .1338E-01, .1518E-01,
     + .1580E-01, .1658E-01, .1710E-01, .1518E-01, .1513E-01, .1570E-01,
     + .1614E-01, .1603E-01, .1673E-01, .1706E-01, .1497E-01, .1439E-01,
     + .1987E-01, .1731E-01, .1601E-01, .1675E-01, .1681E-01, .1535E-01,
     + .1425E-01, .2018E-01, .1723E-01, .1597E-01, .1691E-01, .1666E-01,
     + .1509E-01, .1446E-01,-.2873E-03,-.8031E-04, .4225E-04,-.9287E-04,
     +-.6013E-04,-.4339E-04,-.2474E-04,-.2862E-03,-.8372E-05, .1146E-03,
     +-.9248E-04,-.6166E-04,-.3882E-04,-.1827E-04,-.2870E-03,-.6851E-05,
     + .1865E-03,-.9172E-04,-.6128E-04,-.3616E-04,-.7612E-05,-.2877E-03,
     +-.7231E-05, .1880E-03,-.9287E-04,-.5671E-04,-.4110E-04,-.1104E-04,
     +-.3429E-03,-.7612E-05, .1149E-03,-.9287E-04,-.6356E-04,-.4529E-04,
     +-.2436E-04,-.4187E-03,-.7992E-05, .4339E-04,-.9325E-04,-.6280E-04,
     +-.4225E-04,-.3197E-04,-.4925E-03,-.8754E-05,-.2740E-04,-.9477E-04,
     +-.6432E-04,-.3768E-04,-.3361E-04,-.5511E-03,-.8753E-05,-.9972E-04,
     +-.9515E-04,-.6394E-04,-.3806E-04,-.3787E-04,-.4792E-03,-.1028E-04,
     +-.1534E-03,-.9477E-04,-.6356E-04,-.3616E-04,-.2923E-04,-.5070E-03,
     + .1922E-03,-.1028E-03,-.5823E-04,-.7954E-04,-.2550E-04,-.3893E-04,
     +-.3776E-03,-.1043E-03,-.7993E-04,-.7422E-04,-.4948E-04,-.3007E-04,
     +-.3863E-04, .8335E-04,-.5709E-04,-.6090E-04,-.7840E-04,-.3692E-04,
     +-.3007E-04,-.4251E-04,-.6204E-04,-.4872E-04,-.3806E-04,-.4681E-04,
     +-.3463E-04,-.3007E-04,-.4312E-04,-.1142E-04,-.5176E-04,-.5024E-04,
     +-.3007E-04,-.3730E-04,-.3037E-04,-.3888E-04, .2550E-04,-.6508E-04,
     +-.2512E-04,-.3083E-04,-.3197E-04,-.3041E-04,-.3750E-04, .1484E-04,
     +-.1941E-04,-.2626E-04,-.3349E-04,-.3463E-04,-.2896E-04,-.1716E-04,
     +-.7231E-04,-.3920E-04,-.2893E-04,-.3540E-04,-.3311E-04,-.3734E-04,
     +-.2550E-05,-.7650E-04,-.3159E-04,-.2778E-04,-.3121E-04,-.2169E-04,
     +-.4365E-04,-.1546E-04,-.7916E-04,-.2931E-04,-.2854E-04,-.3654E-04,
     +-.1979E-04,-.4811E-04,-.1435E-04/
	end

	block data ckd17
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and  seven cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 400 to 280 cm**-1.
c *********************************************************************
	common /band17/ hk(7), coeh2o(3,19,7)
 	data hk / .12, .26, .22, .20, .10, .085, .015 /
	data ( ( ( coeh2o(k,j,i), i = 1, 7 ), j = 1, 19 ), k = 1, 3 ) /
     +-.2255E+02,-.2000E+02,-.1703E+02,-.1282E+02,-.9215E+01,-.5938E+01,
     +-.2009E+01,-.2209E+02,-.1997E+02,-.1657E+02,-.1236E+02,-.8764E+01,
     +-.5499E+01,-.1582E+01,-.2163E+02,-.1993E+02,-.1611E+02,-.1191E+02,
     +-.8324E+01,-.5061E+01,-.1170E+01,-.2117E+02,-.1990E+02,-.1565E+02,
     +-.1146E+02,-.7889E+01,-.4631E+01,-.7737E+00,-.2071E+02,-.1987E+02,
     +-.1519E+02,-.1100E+02,-.7440E+01,-.4179E+01,-.3719E+00,-.2026E+02,
     +-.1985E+02,-.1473E+02,-.1054E+02,-.6995E+01,-.3721E+01, .0000E+00,
     +-.2024E+02,-.1982E+02,-.1426E+02,-.1009E+02,-.6549E+01,-.3284E+01,
     + .4053E+00,-.2022E+02,-.1980E+02,-.1381E+02,-.9639E+01,-.6097E+01,
     +-.2821E+01, .8375E+00,-.2021E+02,-.1933E+02,-.1335E+02,-.9187E+01,
     +-.5653E+01,-.2379E+01, .1272E+01,-.2010E+02,-.1503E+02,-.1125E+02,
     +-.7665E+01,-.4492E+01,-.1893E+01, .1642E+01,-.1747E+02,-.1278E+02,
     +-.9547E+01,-.6120E+01,-.3756E+01,-.1443E+01, .1995E+01,-.1529E+02,
     +-.1095E+02,-.8107E+01,-.5036E+01,-.3182E+01,-.1032E+01, .2429E+01,
     +-.1370E+02,-.9303E+01,-.6691E+01,-.4357E+01,-.2683E+01,-.6173E+00,
     + .2805E+01,-.1150E+02,-.7859E+01,-.5618E+01,-.3843E+01,-.2234E+01,
     +-.2171E+00, .2973E+01,-.9590E+01,-.6537E+01,-.4886E+01,-.3355E+01,
     +-.1805E+01, .1615E+00, .3157E+01,-.7530E+01,-.5699E+01,-.4306E+01,
     +-.2892E+01,-.1388E+01, .5448E+00, .3155E+01,-.6758E+01,-.5112E+01,
     +-.3809E+01,-.2464E+01,-.9947E+00, .8713E+00, .3203E+01,-.6245E+01,
     +-.4610E+01,-.3376E+01,-.2058E+01,-.6166E+00, .1073E+01, .3109E+01,
     +-.5777E+01,-.4175E+01,-.2963E+01,-.1671E+01,-.2556E+00, .1241E+01,
     + .3014E+01, .4264E-01, .1968E-02, .1863E-01, .1436E-01, .1101E-01,
     + .1055E-01, .1281E-01, .4264E-01, .1989E-02, .1861E-01, .1438E-01,
     + .1095E-01, .1030E-01, .1211E-01, .3996E-01, .1968E-02, .1861E-01,
     + .1434E-01, .1103E-01, .1019E-01, .1160E-01, .3600E-01, .1947E-02,
     + .1861E-01, .1442E-01, .1086E-01, .1003E-01, .1157E-01, .3203E-01,
     + .5756E-02, .1861E-01, .1444E-01, .1080E-01, .9922E-02, .1151E-01,
     + .2801E-01, .9713E-02, .1859E-01, .1446E-01, .1070E-01, .9880E-02,
     + .1066E-01, .2393E-01, .1369E-01, .1859E-01, .1451E-01, .1057E-01,
     + .9880E-02, .1072E-01, .1987E-01, .1767E-01, .1863E-01, .1451E-01,
     + .1040E-01, .9880E-02, .1057E-01, .1572E-01, .2169E-01, .1863E-01,
     + .1442E-01, .1022E-01, .9742E-02, .1036E-01, .3391E-02, .1884E-01,
     + .1566E-01, .1105E-01, .1011E-01, .1001E-01, .1017E-01, .1982E-01,
     + .1444E-01, .1189E-01, .1030E-01, .9859E-02, .9861E-02, .1038E-01,
     + .1748E-01, .1321E-01, .9922E-02, .1068E-01, .1013E-01, .9937E-02,
     + .9958E-02, .1346E-01, .9943E-02, .9566E-02, .1097E-01, .9815E-02,
     + .9964E-02, .1059E-01, .9817E-02, .7159E-02, .8687E-02, .1114E-01,
     + .1007E-01, .1014E-01, .1058E-01, .3370E-02, .7264E-02, .9378E-02,
     + .1112E-01, .9767E-02, .1016E-01, .1101E-01, .2993E-02, .8017E-02,
     + .9566E-02, .1116E-01, .9738E-02, .1025E-01, .1086E-01, .8331E-02,
     + .8771E-02, .1001E-01, .1117E-01, .9847E-02, .1076E-01, .1084E-01,
     + .7850E-02, .9378E-02, .1001E-01, .1105E-01, .9964E-02, .1113E-01,
     + .1168E-01, .8038E-02, .9336E-02, .9817E-02, .1096E-01, .1024E-01,
     + .1175E-01, .1107E-01,-.2188E-03,-.2283E-05,-.8069E-04,-.4415E-04,
     +-.2284E-04,-.4491E-04,-.4518E-04,-.2196E-03,-.2665E-05,-.8107E-04,
     +-.4301E-04,-.2398E-04,-.4795E-04,-.4693E-04,-.2683E-03,-.3045E-05,
     +-.8107E-04,-.4301E-04,-.2246E-04,-.4757E-04,-.4152E-04,-.3403E-03,
     +-.4187E-05,-.8031E-04,-.3996E-04,-.1865E-04,-.4301E-04,-.4350E-04,
     +-.4118E-03, .6584E-04,-.8107E-04,-.4034E-04,-.1903E-04,-.4643E-04,
     +-.4834E-04,-.4803E-03, .1378E-03,-.8069E-04,-.4072E-04,-.1713E-04,
     +-.5176E-04,-.3460E-04,-.4099E-03, .2101E-03,-.8069E-04,-.3920E-04,
     +-.1713E-04,-.5024E-04,-.3524E-04,-.3391E-03, .2809E-03,-.7992E-04,
     +-.3616E-04,-.2017E-04,-.5633E-04,-.4886E-04,-.2668E-03, .2078E-03,
     +-.8069E-04,-.3768E-04,-.2131E-04,-.5580E-04,-.5454E-04,-.2207E-04,
     +-.8601E-04,-.4643E-04,-.2436E-04,-.4148E-04,-.5458E-04,-.4579E-04,
     +-.5138E-04,-.2893E-04,-.3273E-04,-.3882E-04,-.3920E-04,-.5035E-04,
     +-.3170E-04,-.2169E-04,-.3007E-04,-.2740E-04,-.5328E-04,-.4491E-04,
     +-.4403E-04,-.6383E-04, .4834E-04,-.2702E-04,-.4453E-04,-.4339E-04,
     +-.4457E-04,-.4551E-04,-.8133E-04, .3768E-04,-.7611E-06,-.2626E-04,
     +-.4643E-04,-.4305E-04,-.4840E-04,-.5149E-04, .7193E-04,-.2169E-04,
     +-.4491E-04,-.3996E-04,-.4483E-04,-.4487E-04,-.6698E-04,-.4834E-04,
     +-.3463E-04,-.4986E-04,-.4377E-04,-.4514E-04,-.5377E-04,-.2626E-04,
     +-.4187E-04,-.3692E-04,-.5100E-04,-.4651E-04,-.4392E-04,-.5386E-04,
     +-.4643E-04,-.4301E-04,-.3578E-04,-.5176E-04,-.4594E-04,-.4551E-04,
     +-.3920E-04,-.3425E-04,-.4491E-04,-.3654E-04,-.5138E-04,-.4377E-04,
     +-.5614E-04,-.5758E-04,-.3600E-04/
	end

	block data ckd18
c *********************************************************************
c hk is the interval in the g (cumulative probability) space from 0 
c to one. coeh2o is the coefficient to calculate the H2O absorption
c coefficient in units of (cm-atm)**-1 at there temperatures, nine-
c teen pressures, and eight cumulative probabilities ( Fu,  1991 ).
c The spectral region is from 280 to 0 cm**-1.
c *********************************************************************
	common /band18/ hk(8), coeh2o(3,19,8)
c	data hk / .07, .1, .2, .25, .2, .1, .03, .02 /
 	data hk / .10, .1, .2, .25, .2, .1, .03, .02 /
	data ( ( ( coeh2o(k,j,i), i = 1, 8 ), j = 1, 19 ), k = 1, 3 ) /
     +-.2121E+02,-.2002E+02,-.1676E+02,-.1274E+02,-.8780E+01,-.5167E+01,
     +-.2692E+01,-.6275E+00,-.2075E+02,-.1996E+02,-.1630E+02,-.1228E+02,
     +-.8324E+01,-.4718E+01,-.2260E+01,-.2303E+00,-.2029E+02,-.1990E+02,
     +-.1584E+02,-.1182E+02,-.7868E+01,-.4269E+01,-.1806E+01, .1645E+00,
     +-.2022E+02,-.1985E+02,-.1538E+02,-.1136E+02,-.7417E+01,-.3820E+01,
     +-.1373E+01, .5657E+00,-.2018E+02,-.1981E+02,-.1492E+02,-.1090E+02,
     +-.6965E+01,-.3369E+01,-.9319E+00, .9577E+00,-.2013E+02,-.1937E+02,
     +-.1446E+02,-.1044E+02,-.6512E+01,-.2917E+01,-.4928E+00, .1376E+01,
     +-.2009E+02,-.1891E+02,-.1400E+02,-.9984E+01,-.6063E+01,-.2466E+01,
     +-.6887E-01, .1768E+01,-.2006E+02,-.1845E+02,-.1354E+02,-.9530E+01,
     +-.5618E+01,-.2024E+01, .3615E+00, .2196E+01,-.2003E+02,-.1800E+02,
     +-.1308E+02,-.9075E+01,-.5174E+01,-.1593E+01, .7820E+00, .2600E+01,
     +-.1827E+02,-.1464E+02,-.1097E+02,-.7525E+01,-.3733E+01,-.1077E+01,
     + .1204E+01, .3014E+01,-.1525E+02,-.1210E+02,-.9275E+01,-.5876E+01,
     +-.2768E+01,-.6286E+00, .1622E+01, .3394E+01,-.1298E+02,-.1060E+02,
     +-.7764E+01,-.4462E+01,-.2154E+01,-.2001E+00, .2034E+01, .3756E+01,
     +-.1157E+02,-.8941E+01,-.5984E+01,-.3509E+01,-.1651E+01, .2279E+00,
     + .2422E+01, .4066E+01,-.9986E+01,-.7062E+01,-.4794E+01,-.2818E+01,
     +-.1196E+01, .6394E+00, .2791E+01, .4283E+01,-.8064E+01,-.5512E+01,
     +-.3933E+01,-.2274E+01,-.7559E+00, .1036E+01, .3085E+01, .4444E+01,
     +-.6440E+01,-.4863E+01,-.3219E+01,-.1791E+01,-.3279E+00, .1427E+01,
     + .3304E+01, .4527E+01,-.5902E+01,-.4207E+01,-.2756E+01,-.1350E+01,
     + .7686E-01, .1776E+01, .3475E+01, .4550E+01,-.5439E+01,-.3739E+01,
     +-.2330E+01,-.9233E+00, .4612E+00, .2066E+01, .3564E+01, .4502E+01,
     +-.5006E+01,-.3316E+01,-.1906E+01,-.5066E+00, .8352E+00, .2272E+01,
     + .3587E+01, .4419E+01, .2338E-01, .1968E-02, .9503E-02, .3412E-02,
     + .6280E-03,-.1109E-02,-.1089E-02,-.1026E-02, .1972E-01, .2093E-02,
     + .9503E-02, .3391E-02, .6489E-03,-.1172E-02,-.1164E-02,-.1158E-02,
     + .1603E-01, .3328E-02, .9524E-02, .3391E-02, .6489E-03,-.1277E-02,
     +-.1229E-02,-.1296E-02, .1229E-01, .7138E-02, .9524E-02, .3370E-02,
     + .6070E-03,-.1319E-02,-.1264E-02,-.1610E-02, .8478E-02, .1095E-01,
     + .9566E-02, .3412E-02, .5652E-03,-.1382E-02,-.1266E-02,-.1566E-02,
     + .4563E-02, .1480E-01, .9566E-02, .3412E-02, .5443E-03,-.1423E-02,
     +-.1199E-02,-.1679E-02, .2261E-02, .1865E-01, .9608E-02, .3454E-02,
     + .4815E-03,-.1423E-02,-.1296E-02,-.1555E-02, .2198E-02, .2250E-01,
     + .9671E-02, .3412E-02, .4187E-03,-.1426E-02,-.1472E-02,-.1800E-02,
     + .2072E-02, .2600E-01, .9734E-02, .3433E-02, .3977E-03,-.1428E-02,
     +-.1541E-02,-.1591E-02, .1987E-01, .8645E-02, .6280E-02, .1298E-02,
     +-.1151E-02,-.1509E-02,-.1662E-02,-.1570E-02, .4668E-02, .8373E-02,
     + .3956E-02,-.4187E-04,-.1968E-02,-.1624E-02,-.1700E-02,-.1947E-02,
     + .9231E-02, .5694E-02, .1444E-02,-.2512E-03,-.1827E-02,-.1662E-02,
     +-.1576E-02,-.1633E-02, .8666E-02, .3077E-02,-.1737E-02,-.1277E-02,
     +-.1507E-02,-.1757E-02,-.1612E-02,-.1612E-02, .8164E-03,-.4375E-02,
     +-.1884E-02,-.1277E-02,-.1564E-02,-.1853E-02,-.1591E-02,-.1486E-02,
     +-.1486E-02,-.2596E-02,-.1633E-02,-.1539E-02,-.1662E-02,-.1846E-02,
     +-.1423E-02,-.1277E-02,-.1423E-02,-.2617E-02,-.1005E-02,-.1379E-02,
     +-.1687E-02,-.1905E-02,-.1528E-02,-.1298E-02,-.1675E-03,-.1947E-02,
     +-.5024E-03,-.1325E-02,-.1696E-02,-.1698E-02,-.1486E-02,-.1277E-02,
     + .1047E-03,-.1109E-02,-.5861E-03,-.1363E-02,-.1620E-02,-.1666E-02,
     +-.1507E-02,-.9210E-03, .1047E-03,-.1047E-02,-.8394E-03,-.1342E-02,
     +-.1591E-02,-.1323E-02,-.1340E-02,-.9420E-03,-.1085E-03, .2283E-05,
     +-.4719E-04,-.3807E-06,-.1522E-05,-.3425E-05,-.7612E-06, .1751E-05,
     +-.1766E-03, .1523E-05,-.4719E-04,-.7609E-06,-.3807E-06,-.3045E-05,
     + .1599E-05, .8723E-05,-.2443E-03, .1941E-04,-.4757E-04,-.1522E-05,
     +-.3806E-06,-.1903E-05,-.2778E-05, .1294E-04,-.1838E-03, .8563E-04,
     +-.4757E-04,-.1903E-05, .1142E-05,-.2664E-05,-.6090E-06, .1321E-04,
     +-.1161E-03, .1526E-03,-.4757E-04,-.2664E-05,-.3805E-06,-.3806E-05,
     +-.2093E-05, .2253E-04,-.4795E-04, .9248E-04,-.4757E-04,-.1903E-05,
     + .0000E+00,-.3045E-05,-.7992E-06, .1393E-04,-.9134E-05, .2246E-04,
     +-.4834E-04,-.2664E-05, .3804E-06,-.5328E-05,-.1510E-05, .1465E-04,
     +-.1028E-04,-.4757E-04,-.4948E-04,-.1142E-05, .7614E-06,-.4910E-05,
     +-.5709E-06, .1477E-04,-.1256E-04,-.1066E-03,-.4910E-04,-.1523E-05,
     +-.3805E-06,-.3121E-05,-.2512E-05, .1142E-04,-.7878E-04,-.2664E-05,
     +-.8373E-05,-.7612E-06, .1104E-04,-.3311E-05,-.1979E-05, .5709E-05,
     +-.2626E-04,-.4872E-04,-.3808E-06,-.2283E-05, .2284E-05,-.3349E-05,
     +-.4034E-05, .7231E-05,-.4910E-04, .1599E-04, .1256E-04,-.7612E-05,
     + .1180E-05,-.1815E-05,-.7193E-05, .3045E-05, .1576E-09, .6470E-05,
     +-.1408E-04,-.1903E-05, .1522E-05,-.4746E-05,-.4948E-05, .3806E-06,
     + .9020E-04, .5214E-04, .6090E-05,-.1104E-04, .1180E-05,-.2778E-05,
     +-.6090E-05,-.2664E-05,-.6737E-04,-.1218E-04,-.3806E-05,-.5214E-05,
     +-.1066E-05,-.1294E-05,-.3045E-05,-.2664E-05,-.4643E-04, .1713E-04,
     +-.1218E-04,-.6204E-05,-.2360E-05,-.1979E-05,-.1903E-05,-.3806E-05,
     +-.3045E-04,-.1256E-04,-.9134E-05,-.6508E-05,-.1027E-05,-.7993E-06,
     +-.1142E-05,-.7992E-05,-.3616E-04,-.1028E-04,-.1066E-04,-.6051E-05,
     + .1066E-05,-.1751E-05,-.2284E-05,-.2284E-05,-.3920E-04,-.9895E-05,
     +-.1321E-04,-.3844E-05,-.2055E-05,-.2512E-05,-.3806E-05,-.3425E-05/
	end
