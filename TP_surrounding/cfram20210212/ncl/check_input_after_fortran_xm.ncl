;************************************************************
;draw all the input data case by case baseline_ran3_ncl.f
;************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;------------------------------------------------------------
  fName = "/home/yangsong3/data-observation/ERA5-monthly/pressure/t_mon_1x1_1979_2019.nc"
  f1    = addfile(fName,"r")
  lat  = f1->latitude
  lon  = f1->longitude
  level= f1->level
  delete(fName)
  print(level)
  print(lat)

  vname = (/"solarin","ssrd","ssru","ts","huss","ps","o3","ta","hus","cc","clwc","ciwc"/)

  nlat  = dimsizes(lat)
  mlon  = dimsizes(lon)
  nlev  = dimsizes(level)
  nvar  = dimsizes(vname)
  nrec  = 6+nlev*6
;------------------------------------------------------------
  UNDEF = -999.
  in_3d = new ( (/nrec,nlat,mlon/), float, UNDEF)

  setfileoption("bin","ReadByteOrder","Native")
  fName1 = "../baseline_input_1.dat"

  do i = 0,5
    in_3d(i,:,:) = fbindirread(fName1,i,(/nlat,mlon/), "float")
  end do

  irec = 6
  do i = 0,5
    do j = 0,nlev-1
      in_3d(irec,:,:) = fbindirread(fName1,irec,(/nlat,mlon/), "float")
      irec = irec+1
    end do
  end do

  in_3d!0 = "lev"
  in_3d!1 = "lat"
  in_3d!2 = "lon"
  in_3d&lat  = lat
  in_3d&lon  = lon
  in_3d&lat@units = "degrees_north"
  in_3d&lon@units = "degrees_east"

  print("end of input data")

  print(in_3d(:,60,160))
;-------------------------------------------------------------
  figpath = "/home/yangsong3/wuyt/sysu/figures/TP_surrounding/cfram20210212/ncl/"
  figname = "check_input_after_fortran_xm"
  wks  = gsn_open_wks("pdf",figpath+figname)
  gsn_define_colormap(wks,"BlueDarkRed18")

  plot = new(6,graphic)

  res                   = True
  res@mpFillOn          = True
  res@mpCenterLonF      = 180
  res@mpGeophysicalLineThicknessF=1.

  res@gsnDraw           = False
  res@gsnFrame          = False
  res@gsnSpreadColors   = True

  res@cnFillOn          = True
  res@cnLinesOn         = False

  pres                  = True
  pres@gsnMaximize      = True

  do k = 0,5
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = min(in_3d(k,:,:))
  res@cnMaxLevelValF       = max(in_3d(k,:,:))
  res@cnLevelSpacingF      = (max(in_3d(k,:,:))-min(in_3d(k,:,:)))/10.

  res@gsnLeftString        = vname(k)
  plot(k) = gsn_csm_contour_map(wks,in_3d(k,:,:),res)
  end do
  gsn_panel(wks,plot,(/3,2/),pres)

  l = (/36,26,11/)
  do i = 0,5
    do j = 0,2
      ll = l(j)+nlev*i+6
      res@cnLevelSelectionMode = "ManualLevels"
      res@cnMinLevelValF       = min(in_3d(ll,:,:))
      res@cnMaxLevelValF       = max(in_3d(ll,:,:))
      res@cnLevelSpacingF      = (max(in_3d(ll,:,:))-min(in_3d(ll,:,:)))/10.
      res@gsnLeftString        = level(l(j))+" hPa"
      plot(j) = gsn_csm_contour_map(wks,in_3d(ll,:,:),res)
    end do
    pres@txString         = vname(i+6)
    gsn_panel(wks,plot(0:2),(/3,1/),pres)
  end do

end
