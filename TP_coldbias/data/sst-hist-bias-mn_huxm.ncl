;*************************************************
; corel_3.ncl
;
; Concepts illustrated:
;   - Calculating a two-dimensional correlation in time
;   - Reordering an array
;   - Copying attributes from one tariable to another
;   - Copying coordinate arrays from one tariable to another
;
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

ntime = 600
nlat = 91
nlon = 180
nmodel = 28
;************************************************
; open file and read in tariable
;************************************************
  yrStrtmodel1 = 1950          ;定义
  yrLastmodel1 = 1999
  yrStrtmodel2 = 2050          ;定义
  yrLastmodel2 = 2099
  dirh="ls /home/yangsong3/jyt/ensocro/sst/ssth/"
  dirf="ls /home/yangsong3/jyt/ensocro/sst/sstf/"

  ;read the name of all of models
  dirhname = asciiread("/home/yangsong3/jyt/name/dirssth28.txt",(/28/),"string")
  dirfname = asciiread("/home/yangsong3/jyt/name/dirsstf28.txt",(/28/),"string")

  ssthi= new((/600,91,180/),float)   ;historical SST 2D
  sstfi= new((/600,91,180/),float)   ;future SST 2D

  ssthclm=new((/28,12,91,180/),float) ;historical SST mean of 28 models
  sstfclm=new((/28,12,91,180/),float) ;historical SST mean of 28 models

  time3=ispan(1,600,1) ;months

;===============================================
;models
;===============================================
   do i=0,27
	  ;==============read SST==============================
	  filsh = systemfunc (dirh+dirhname(i))
    filsf = systemfunc (dirf+dirfname(i))

	  ssthi = fbindirread(filsh,0,(/600,91,180/),"float")
	  sstfi = fbindirread(filsf,0,(/600,91,180/),"float")
	  ssthi!0 = "time"
    ssthi!1 = "lat"
	  ssthi!2 = "lon"
	  ssthi&time = time3
    ssthi&lat  = fspan(-90,90,91)
    ssthi&lon  = fspan(0,358,180)
    ssthi&lat@units = "degrees_north"
    ssthi&lon@units = "degrees_east"

	  sstfi!0 = "time"
    sstfi!1 = "lat"
	  sstfi!2 = "lon"
	  sstfi&time = time3
    sstfi&lat  = fspan(-90,90,91)
    sstfi&lon  = fspan(0,358,180)
    sstfi&lat@units = "degrees_north"
    sstfi&lon@units = "degrees_east"

	  if (any(isnan_ieee(ssthi))) then
       value = 1.e20
       replace_ieeenan (ssthi, value, 0)
       ssthi@_FillValue = value
     end if

	  if (any(isnan_ieee(sstfi))) then
       value = 1.e20
       replace_ieeenan (sstfi, value, 0)
       sstfi@_FillValue = value
    end if

   ;------detrend---------
   ssthdt  =  dtrend_msg_n(time3,ssthi,False,False,0);去趋势
	 sstfdt  =  dtrend_msg_n(time3,sstfi,False,False,0)
	 copy_VarAtts(ssthi,ssthdt)
   copy_VarCoords(ssthi,ssthdt)
	 copy_VarAtts(sstfi,sstfdt)
   copy_VarCoords(sstfi,sstfdt)

   ;------anomaly---------
	 ssthclm(i,:,:,:) = clmMonTLL(ssthdt)
	 sstfclm(i,:,:,:) = clmMonTLL(sstfdt)
  end do

;===============================================
;bias of models
;===============================================
;----------------read reanalysis----------------
  dir1 = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
  fName = addfile(dir1+"sst.mnmean.v5.nc","r")

  TIME      = fName->time
  YYYY      = cd_calendar(TIME,-1)/100
  iYYYY     = ind(YYYY.ge.yrStrtmodel1 .and. YYYY.le.yrLastmodel1)
  Time      = fName->time(iYYYY)
  time      = ut_calendar(Time,-1)
  lon       = fName->lon
  lat       = fName->lat(::-1)

  sstobs    = fName->sst(iYYYY,::-1,:)
  sstobs!0  = "time"
  sstobs!1  = "lat"
  sstobs!2  = "lon"
  sstobs&time = time

 ;------interpolate sstobs to same resolution to SST from CMIP5---------
  LON = fspan(0,358,180)
  LAT = fspan(-90,90,91)
  copy_VarAtts(lon, LON)
  copy_VarAtts(lat, LAT)
  SSTOBS = linint2_Wrap(lon,lat,sstobs,True,LON,LAT,0)

  ;------detrend---------
  sstodt  =  dtrend_msg_n(time,SSTOBS,False,False,0);去趋势
  copy_VarAtts(SSTOBS,sstodt)
  copy_VarCoords(SSTOBS,sstodt)

  ;------anomaly---------
  sstoclm = clmMonTLL(sstodt)

  ;-------model bias------
  biasclm = ssthclm-conform(ssthclm, sstoclm, (/1,2,3/))-273.15
  shiftclm = sstfclm-ssthclm
  copy_VarMeta(ssthclm, biasclm)
  copy_VarMeta(ssthclm, shiftclm)

  biasclm!0 = "time"
  shiftclm!0 = "time"
  ;------------ensemble mean--------------------------
  biasclmem = dim_avg_n_Wrap(biasclm, 0)
  shiftclmem = dim_avg_n_Wrap(shiftclm, 0)

  ;-----------data output-----------------------------
  ;---Remove file just in case
    file1 = "cmip5_hist_sst_mn_clm_bias_mme.dat"
    system("rm -f " + file1)
  ;---Write first variable to file
    fbindirwrite(file1,biasclmem)
    system("ls -l " + file1)

  ;---Remove file just in case
    file2 = "cmip5_futr_sst_mn_clm_shift_mme.dat"
    system("rm -f " + file2)
  ;---Write first variable to file
    fbindirwrite(file2,shiftclmem)
    system("ls -l " + file2)

  ;-------------eof pts-------------------------------
  biasclmeof = biasclmem(:,{-30:30},{110:290})
  shiftclmeof = shiftclmem(:,{-30:30},{110:290})

  neof   = 1        ; number of EOFs
  optEOF = True
  optEOF@jopt = 0   ; This is the default; most commonly used; no need to specify.

  do i = 0,11
    in = biasclm(:,i,{-30:30},{110:290})
    input = in(lat|:,lon|:,time|:)
    biasclmeof(i,:,:) = eofunc(input,neof,optEOF)
  end do
  copy_VarMeta(biasclmem(:,{-30:30},{110:290}),biasclmeof)
  printVarSummary(biasclmeof)
;=======================================================
;figure output
;=======================================================
  wks   = gsn_open_wks("pdf" ,"cmip5_sst_clm_bias_gl")
  gsn_define_colormap(wks,"BlueDarkRed18")
  plot  = new(12,graphic)

  res                      = True            ; plot mods desired

  res@gsnDraw              = False
  res@gsnFrame             = False
  res@cnFillOn             = True            ; turn on color

  ;---This resource not needed in NCL V6.1.0
  res@gsnSpreadColors      = True            ; use full color table
  res@cnLinesOn            = False           ; no contour lines
  res@lbLabelBarOn         = True

  resh                     = res

  resh@mpMinLonF         = 0;110
  resh@mpMaxLonF         = 360;290
  resh@mpMinLatF         = -90
  resh@mpMaxLatF         = 90
  resh@mpCenterLonF      = 180;200
  resh@mpGeophysicalLineThicknessF=1.

  resh@cnLevelSelectionMode = "ManualLevels"
  resh@cnMinLevelValF       = -5
  resh@cnMaxLevelValF       = 5
  resh@cnLevelSpacingF      = 1.

  do i = 0, 11
  resh@gsnLeftString     = "Month = "+(i+1)
  plot(i) = gsn_csm_contour_map(wks,biasclmem(i,:,:),resh)
  end do

  resP                     = True
  resP@gsnMaximize         = True
  resP@txString            = "Climatology Bias (MME)"
  gsn_panel(wks,plot(0:5),(/3,2/),resP)
  gsn_panel(wks,plot(6:11),(/3,2/),resP)

  resh@cnLevelSelectionMode = "ManualLevels"
  resh@cnMinLevelValF       = -5
  resh@cnMaxLevelValF       = 5
  resh@cnLevelSpacingF      = 1.

  do i = 0, 11
  resh@gsnLeftString     = "Month = "+(i+1)
  plot(i) = gsn_csm_contour_map(wks,shiftclmem(i,:,:),resh)
  end do

  resP@txString            = "Future Climatology Shift (MME)"
  gsn_panel(wks,plot(0:5),(/3,2/),resP)
  gsn_panel(wks,plot(6:11),(/3,2/),resP)
end
