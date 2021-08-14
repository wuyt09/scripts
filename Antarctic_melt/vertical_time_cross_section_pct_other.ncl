load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

 yrStrt = 201512
 yrLast = 201601
 nvar = 7
 njob = 62
;************************************************
; read air temperature in 201512-201601
;************************************************
 f1     = addfile("/home/yangsong3/huxm/antarctic_melt/201512_201601_t_4time_daily.nc","r")
 TIME   = cd_calendar(f1->time,-3)
 time   = TIME(0::4)/100
 level  = f1->level
 lat    = f1->latitude
 lon    = f1->longitude
 delete(f1)

 nt = dimsizes(time)
 nlev = dimsizes(level)
 nlat = dimsizes(lat)
 nlon = dimsizes(lon)

 lev = new((/nlev+1/), typeof(level))
 lev(0:nlev-1) = level
 lev(nlev) = 1001
;************************************************
; read climatological mean air temperature
;************************************************
 case    = "9192"
 dirin   = "/home/yangsong3/wuyt/sysu/scripts/Antarctic_melt/"+case+"/"
 ; dirin = "/home/yangsong3/huxm/antarctic_melt/"
 fname1  = "ta_warm.dat"
 fname2  = "ta_base.dat"
 twarm   = new ( (/nt,nlev,nlat,nlon/),"float")
 tbase   = new ( (/nt,nlev,nlat,nlon/),"float")

 twarm!0 = "time"
 twarm!1 = "lev"
 twarm!2 = "lat"
 twarm!3 = "lon"
 twarm&time = time
 twarm&lev  = level
 twarm&lat  = lat
 twarm&lon  = lon
 twarm&lev@units = "hPa"
 twarm&lat@units = "degrees_north"
 twarm&lon@units = "degrees_east"
 copy_VarMeta(twarm, tbase)

 setfileoption("bin","ReadByteOrder","Native")

 irec = 0
 do tt = 0, nt-1
   do ll = 0, nlev-1
      twarm(tt,ll,:,:) = fbindirread(dirin+fname1, irec, (/nlat, nlon/), "float")
      tbase(tt,ll,:,:) = fbindirread(dirin+fname2, irec, (/nlat, nlon/), "float")
      irec = irec+1
   end do
 end do

 rad  = 4.0*atan(1.0)/180.0
 clat = cos(lat*rad)
 clat!0 = "lat"
 clat&lat = lat

 twais = wgt_areaave_Wrap(twarm(:,:,{-85:-75},{180:240}),clat({-85:-75}),1.0,0)

 twarm = twarm-tbase
 twaisano = wgt_areaave_Wrap(twarm(:,:,{-85:-75},{180:240}),clat({-85:-75}),1.0,0)

 x = ispan(1,62,1)
 x(31:61) = x(0:30)
 y = ispan(1,62,1)
 y(0:30) = 1512
 y(31:61) = 1601
 date = y*100+x
 ; print(date)

 ;---------------------------CFRAM data input-----------------------------
   UNDEF = -999.
   ptc  = new ( (/njob,nvar,nlev+1,nlat,nlon/), float, UNDEF)

   ptc@units = "W/m**2"
   ptc!0 = "time"
   ptc!1 = "var"
   ptc!2 = "lev"
   ptc!3 = "lat"
   ptc!4 = "lon"
   ptc&lev  = lev
   ptc&lat  = lat
   ptc&lon  = lon
   ptc&lat@units = "degrees_north"
   ptc&lon@units = "degrees_east"

   do job = 1,62
      fName1 = dirin+"partial_T_"+job+".grd"

     setfileoption("bin","ReadByteOrder","Native")

     do nv = 0,nvar-1
       do nl = 0, nlev
        irec = nl+nv*38
        ptc(job-1,nv,nl,:,:) = fbindirread(fName1,irec,(/nlat,nlon/), "float")
       end do
     end do

     ; print("job="+job)
   end do
   print("end of data input")
      print(ptc(12,1,10,{-85:-75},{180:240}))
      exit
   if (any(isnan_ieee(ptc))) then
        value = -999.
        replace_ieeenan (ptc, value, 0)
        ptc@_FillValue = value
   end if

   ptc = where(abs(ptc).gt.400, ptc@_FillValue, ptc)

   rad  = 4.0*atan(1.0)/180.0
   clat = cos(lat*rad)
   clat!0 = "lat"
   clat&lat = lat



;original value is {-85:-70},{180:240}
   PTC_temp = wgt_areaave_Wrap(ptc(:,:,:,{-85:-75},{180:240}),clat({-85:-75}),1.0,0)
   printVarSummary(PTC_temp)
   PTC = PTC_temp(var|:, lev|:, time|:)
   PTC = runave_n_Wrap(PTC, 3, 1, 2)
;************************************************
; create plot 2
;************************************************
    wks = gsn_open_wks("pdf","vertical_time_cross_ptc"+case)	; open a pdf file,"X11"disable at once
;    wks = gsn_open_wks("pdf","./vertical_time_cross_ptc")	; open a pdf file,"X11"disable at once
    gsn_define_colormap(wks,"temp_19lev")
    plot = new(3,graphic)
    plotover = new(3, graphic)

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnSpreadColors  = True

    res@cnFillOn         = True
    res@cnLinesOn        = False
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
    res@lbLabelBarOn     = False

    res@tmXTOn     = False
    res@trYLog     = False
    res@tmYRMode   = "Automatic"
    res@tmYLMode   = "Explicit"
    res@tmYLValues = (/"1000","800","600","400","200"/)
    res@tmYLLabels = (/"1000","800","600","400","200"/)
    res@tmYUseLeft = False

    res@tmXTOn            = False
    res@tmXBMode          = "Explicit"
    res@tmXBValues        = ispan(1,60,10)
    res@tmXBLabels        = date(0::10)
    res@tmXBMinorOn       = True
    res@tmXBMinorValues   = ispan(1,62,2)
    res@tmXBLabelFontHeightF = 0.026
    res@tmYLLabelFontHeightF = 0.026
    res@gsnStringFontHeightF = 0.025

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels          = (/-20,-15,-12,-10,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,10,12,15,20/)

    res@gsnRightString    = "K"
    res@gsnLeftString     = "(b) PTC due to WV"
    plot(1) = gsn_csm_pres_hgt(wks, PTC(1,{50:1000},:), res)

    res@gsnLeftString     = "(a) PTC due to Cloud"
    plot(0) = gsn_csm_pres_hgt(wks, PTC(2,{50:1000},:), res)

    res@gsnLeftString     = "(c) PTC due to Dynamics"
    plot(2) = gsn_csm_pres_hgt(wks, PTC(6,{50:1000},:), res)

;    res@cnFillOn         = False
;    res@cnLinesOn        = True
;    res@cnLineLabelsOn   = True
;    total = PTC(6,:,:)
;    total = PTC(6,:,:)+PTC(4,:,:)+PTC(2,:,:)

;    delete([/res@cnLevels,res@gsnLeftString/])
;    res@gsnContourNegLineDashPattern = 1
;    res@cnLevelSelectionMode = "ExplicitLevels"
;    res@cnLevels = (/-40,-30,-20,-10,-5,0,5,10,20,30,40/)
;    plotover(0) = gsn_csm_pres_hgt(wks, total({50:1000},:), res)
;    plotover(0) = ColorNegDashZeroPosContour(plotover(0),"red","magenta","blue")
;    plotover(1) = gsn_csm_pres_hgt(wks, total({50:1000},:), res)
;    plotover(1) = ColorNegDashZeroPosContour(plotover(1),"red","magenta","blue")
;    plotover(2) = gsn_csm_pres_hgt(wks, total({50:1000},:), res)
;    plotover(2) = ColorNegDashZeroPosContour(plotover(2),"red","magenta","blue")
;    overlay(plot(0),plotover(0))
;    overlay(plot(1),plotover(1))
;    overlay(plot(2),plotover(2))
    ;************************************************
    ; create panel
    ;************************************************
      resP                     = True
      resP@gsnFrame            = False
      resP@gsnMaximize         = True
      resP@gsnPanelLabelBar    = True
      resP@lbOrientation       = "Horizontal"
      resP@pmLabelBarOrthogonalPosF = 0.02
      resP@gsngsnPaperOrientation = "landscape"
      gsn_panel(wks,plot,(/1,3/),resP)

 end
