load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  yrStrt = 1981
  yrLast = 2014

;************************************************
; read data from *_base/*_warm.dat
;************************************************
vars  = (/"albedo","water vapor","cloud","cloud_sw","cloud_lw","co2","o3"\
          ,"solar","dynamic","atmos_dyn","ocean","shflx","lhflx"/)
nvar  = dimsizes(vars)
;----------------------get coordinate----------------------------
  f1 = addfile("../hus_era_mn_197901_201507_4time.nc","r")

  lat  = f1->latitude
  lon  = f1->longitude
  level= f1->level

  nlat = dimsizes(lat)
  mlon = dimsizes(lon)
  nlev = 38

  lev  = new (nlev, integer,"No_FillValue")
  lev(0:nlev-2) = level
  lev(nlev-1) = 1001
  delete(f1)
;---------------------read data----------------------------------
fili1  = "forcing_1_ERA_SRF_res.grd"
fili2  = "forcing_by_temp_fb.grd"

  UNDEF   = -999.                           ; UNDEF
  dir     = new ( (/nvar,nlev,nlat,mlon/), float, UNDEF)
  indir   = new ( (/nvar,nlat,mlon/), float, UNDEF)
  dir@long_name = "Energy Perturbation"
  dir@units     = "W m~S~-2~N~"

  setfileoption("bin","ReadByteOrder","Native")

  nl = 37

  do nv = 0,nvar-1
    do nl = 0,nlev-1
     dir(nv,nl,:,:) = fbindirread(fili1,nv*nlev+nl, (/nlat,mlon/), "float")
    end do
     indir(nv,:,:) = fbindirread(fili2,nv, (/nlat,mlon/), "float")
  end do

  dir!0 = "var"
  dir!1 = "lev"
  dir!2 = "lat"
  dir!3 = "lon"
  dir&lat  = lat
  dir&lon  = lon
  dir&lev  = lev
  dir&lat@units = "degrees_north"
  dir&lon@units = "degrees_east"

  if (any(isnan_ieee(dir))) then
         value = -999.
         replace_ieeenan (dir, value, 0)
         dir@_FillValue = value
  end if
  if (any(isnan_ieee(indir))) then
         value = -999.
         replace_ieeenan (indir, value, 0)
         indir@_FillValue = value
  end if

  copy_VarMeta(dir(:,0,:,:), indir)
  print ("read data end")

  indir(9,:,:)=where(abs(indir(9,:,:)).lt.250,indir(9,:,:),0)
  dir(9,:,:,:)=where(abs(dir(9,:,:,:)).lt.250,dir(9,:,:,:),0)
  dir_zn = dim_avg_n_Wrap(dir, 3)
  total = indir
  total = indir+dir(:,37,:,:)
;********************************************************
; Create  Plot page1
;********************************************************
  i = 5
  vname = "CO2"

  wks   = gsn_open_wks ("ps","qxxb-Fig7-4figs")             ; open workstation
;  gsn_define_colormap(wks,"BlueWhiteOrangeRed")
  plot  = new(4,graphic)

  res                  = True                      ; plot mods desired
  res@gsnFrame         = False
  res@gsnDraw          = False
;  res@gsnSpreadColors  = True                      ; use full colormap
  res@cnFillOn         = True                      ; turn on color
  res@cnLinesOn        = False
  res@cnLineLabelsOn   = False                     ; no contour line labels
  res@cnInfoLabelOn    = False                     ; no info label
  res@lbLabelBarOn     = True


  res@gsnStringFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.016                 ; res2ize tick labels
  res@tmYLLabelFontHeightF = 0.016
  res@lbLabelFontHeightF   = 0.02

  res@tmXTOn = False
  res@tmYROn = False

  sres = res

  res@mpOutlineOn      = True
  res@mpFillOn         = True
  res@mpMinLatF        = -90.	; zoom in on a subregion
  res@mpMaxLatF        = 90.
  res@mpMinLonF        = 0.
  res@mpCenterLonF     = 180.
  res@mpMaxLonF        = 360.
  res@mpGeophysicalLineThicknessF=1.3
  res@pmLabelBarWidthF = 0.7
  res@pmLabelBarOrthogonalPosF = 0.12

  res1 = res
  res1@cnFillPalette        = "BlWhRe"
  res1@cnLevelSelectionMode = "ManualLevels"

  res1@cnMinLevelValF = -1.
  res1@cnMaxLevelValF = 1.
  res1@cnLevelSpacingF = .1

  res1@gsnLeftString = "(a) "+vname+" induced direct change in downward rad. flux"
  res1@gsnRightString = "W m~S~-2~N~"
  plot(0) = gsn_csm_contour_map(wks,dir(i,37,:,:),res1)

; this controls the size and location of the second plot
;  sres@vpXF                 = 0.15
;  sres@vpYF                 = 0.4
  sres@vpWidthF             = 0.6
  sres@vpHeightF            = 0.4

  sres@cnFillPalette        = "BlWhRe"
  sres@gsnLeftString = "(d) "+vname+" induced direct energy perturb."
  sres@gsnRightString = "W m~S~-2~N~"
  sres@cnLevelSelectionMode = "ManualLevels"

  sres@cnMinLevelValF = -.1
  sres@cnMaxLevelValF = .1
  sres@cnLevelSpacingF = .01

  sres@tmYLLabelFontHeightF = 0.017
  sres@tmYLLabelDeltaF = -0.4
  sres@tmYRMode = "Automatic"
  sres@tmYROn = False

  sres@tiYAxisString        = "Pressure (hPa)"
  sres@tiYAxisOffsetXF      = 0.017
  sres@lbOrientation        = "Vertical"
  plot(3)  = gsn_csm_pres_hgt(wks, dir_zn(5,18:36,:), sres)   ; place holder

  res2 = res
  res2@cnFillPalette        = "BlWhRe"
  res2@cnLevelSelectionMode = "ManualLevels"

  res2@cnMinLevelValF = -2.
  res2@cnMaxLevelValF = 2.
  res2@cnLevelSpacingF = .2

  res2@gsnRightString    = "W m~S~-2~N~"
  res2@gsnLeftString     = "(b) "+vname+" induced indirect change in downward rad. flux"
  plot(1) = gsn_csm_contour_map_ce(wks,indir(i,:,:),res2)

  res2@gsnRightString    = "W m~S~-2~N~"
  res2@gsnLeftString     = "(c) Total "+vname+" induced change in downward rad. flux"
  plot(2) = gsn_csm_contour_map_ce(wks,total(i,:,:),res2)

;************************************************
; create panel
;************************************************
  resP                     = True
  resP@gsnFrame            = False
  resP@gsnMaximize         = True
  resP@gsnPanelLabelBar    = False                ; add common colorbar
  resP@pmLabelBarHeightF        = 0.1
  resP@pmLabelBarWidthF         = 0.85
  resP@lbLabelFontHeightF       = 0.013

  gsn_panel(wks,plot,(/2,2/),resP)
  frame(wks)
end
