
begin

yrStrt = 1950
yrLast = 1999
season = "DJF"

;===============================================
; OBS 
;===============================================
;----------------read reanalysis----------------
  diri = "/home/yangsong3/data-observation/SST/NOAA_ERSST/"
  fName = addfile(diri+"sst.mnmean.v5.nc","r")

  TIME      = fName->time
  YYYY      = cd_calendar(TIME,-1)/100
  iYYYY     = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
  Time      = fName->time(iYYYY)
  time      = ut_calendar(Time,-1)
  lon       = fName->lon
  lat       = fName->lat(::-1)
  nlat      = dimsizes(lat)
  nlon      = dimsizes(lon)

  sstobs    = fName->sst(iYYYY,::-1,:)
  sstobs    = sstobs+273.15
  sstobs!0  = "time"
  sstobs!1  = "lat"
  sstobs!2  = "lon"
  sstobs&time = time
  sstobs@units = "K"
  printVarSummary(sstobs)


;===================================================================
; CMIP5 28 GCMs sst data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/"

models = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
bias_sst5 = new((/nm,600,nlat,nlon/), float)
bias_sst5!0 = "model"
bias_sst5&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+models(i)+"/ts_Amon_"+models(i)+"_historical_r1i1p1_18*-20*.nc")
    print(filigcm)
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    TS_GCMs := fgcm->ts(iYYYY,:,:)
    TS_GCMs!1 = "lat"
    TS_GCMs!2 = "lon"

    TS_GCM := linint2_Wrap(TS_GCMs&lon, TS_GCMs&lat, TS_GCMs, True, lon, lat, 0)
    TS_GCM = mask(TS_GCM, ismissing(sstobs), False)
    bias_sst5(i,:,:,:) = TS_GCM - sstobs
    print("end of "+models(i))
end do
copy_VarCoords(sstobs, bias_sst5(0,:,:,:))
bias_sst5A = dim_avg_n_Wrap(bias_sst5, 0)
printVarSummary(bias_sst5A)
;===================================================================
; CMIP6 20 GCMs sst data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP6/"

models := (/"BCC-CSM2-MR","BCC-ESM1","CanESM5","CAS-ESM2-0","CESM2","CESM2-FV2","CESM2-WACCM",\
          "CESM2-WACCM-FV2","CIESM","E3SM-1-0","E3SM-1-1","E3SM-1-1-ECA","FGOALS-f3-L",\
          "FGOALS-g3","KIOST-ESM","MCM-UA-1-0","MIROC6","MRI-ESM2-0","NESM3","SAM0-UNICON"/)
; models = (/"B" /)

nm := dimsizes(models)
bias_sst6 = new((/nm,600,nlat,nlon/), float)
bias_sst6!0 = "model"
bias_sst6&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"ts_Amon_"+models(i)+"_historical_r1i1p1f1_*_18*-20*.nc")
    print(filigcm)
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    TS_GCMs := fgcm->ts(iYYYY,:,:)
    TS_GCMs!1 = "lat"
    TS_GCMs!2 = "lon"

    TS_GCM := linint2_Wrap(TS_GCMs&lon, TS_GCMs&lat, TS_GCMs, True, lon, lat, 0)
    printVarSummary(TS_GCM)
    TS_GCM = mask(TS_GCM, ismissing(sstobs), False)
    bias_sst6(i,:,:,:) = TS_GCM - sstobs
    print("end of "+models(i))
end do
copy_VarCoords(sstobs, bias_sst6(0,:,:,:))
bias_sst6A = dim_avg_n_Wrap(bias_sst6, 0)
printVarSummary(bias_sst6A)

;==================================================================
; DJF map
;==================================================================
bias_sst5_seax = month_to_season(bias_sst5A, season)
bias_sst5_sea = dim_avg_n_Wrap(bias_sst5_seax, 0)
printVarSummary(bias_sst5_sea)

bias_sst6_seax = month_to_season(bias_sst6A, season)
bias_sst6_sea = dim_avg_n_Wrap(bias_sst6_seax, 0)
printVarSummary(bias_sst6_sea)

diff = bias_sst5_sea
diff = bias_sst6_sea - bias_sst5_sea

;==================================================================
; 12 month time series
;==================================================================
bias_sst5_clm = clmMonTLL(bias_sst5A)
printVarSummary(bias_sst5_clm)

bias_sst6_clm = clmMonTLL(bias_sst6A)
printVarSummary(bias_sst6_clm)

fili_landsea = "/home/yangsong3/soft/ncl650/lib/ncarg/data/cdf/landsea.nc"
flandsea = addfile(fili_landsea, "r")
LANDSEA = flandsea->LSMASK
landsea = linint2_Wrap(LANDSEA&lon, LANDSEA&lat, LANDSEA, True, bias_sst5_clm&lon, bias_sst5_clm&lat, 0)
printVarSummary(landsea)

;CMIP5
bias_sst5_io = bias_sst5_clm
bias_sst5_io = bias_sst5_clm@_FillValue
bias_sst5_po = bias_sst5_io
bias_sst5_ao = bias_sst5_io
bias_sst5_pio = bias_sst5_io
bias_sst5_iao = bias_sst5_io
bias_sst5_pao = bias_sst5_io
    ;==============================================================
    ; Indian Ocean
    ;==============================================================
    bias_sst5_io(:,{-30:30},{30:105}) = bias_sst5_clm(:,{-30:30},{30:105})
    bias_sst5_io(:,{-30:-10},{105:130}) = bias_sst5_clm(:,{-30:-10},{105:130})

    bias_sst5_pio(:,{-30:30},{30:105}) = bias_sst5_clm(:,{-30:30},{30:105})
    bias_sst5_pio(:,{-30:-10},{105:130}) = bias_sst5_clm(:,{-30:-10},{105:130})

    bias_sst5_iao(:,{-30:30},{30:105}) = bias_sst5_clm(:,{-30:30},{30:105})
    bias_sst5_iao(:,{-30:-10},{105:130}) = bias_sst5_clm(:,{-30:-10},{105:130})

    ;==============================================================
    ; Pacific Ocean
    ;==============================================================
    bias_sst5_po(:,{-10:30},{105:130}) = bias_sst5_clm(:,{-10:30},{105:130})
    bias_sst5_po(:,{-30:30},{130:240}) = bias_sst5_clm(:,{-30:30},{130:240})

    bias_sst5_pio(:,{-10:30},{105:130}) = bias_sst5_clm(:,{-10:30},{105:130})
    bias_sst5_pio(:,{-30:30},{130:240}) = bias_sst5_clm(:,{-30:30},{130:240})

    bias_sst5_pao(:,{-10:30},{105:130}) = bias_sst5_clm(:,{-10:30},{105:130})
    bias_sst5_pao(:,{-30:30},{130:240}) = bias_sst5_clm(:,{-30:30},{130:240})

    do j = 32, 63 ; -29.36 to 29.36 56, 57;
        count = 0
        do i = 120, 145 ; 240 to 290
            if (landsea(j,i).lt.1 .and. count.eq.0) then
                bias_sst5_po(:,j,i) = bias_sst5_clm(:,j,i)
                bias_sst5_pio(:,j,i) = bias_sst5_clm(:,j,i)
                bias_sst5_pao(:,j,i) = bias_sst5_clm(:,j,i)
                print("lat "+lat(j)+"   lon "+(360-lon(i))) ;lon(i)
            else
                count = count + 1
            end if
        end do
    end do
    print("=====================================================")
    ;==============================================================
    ; Atlantic Ocean
    ;==============================================================
    bias_sst5_ao(:,{-30:30},{0:30}) = bias_sst5_clm(:,{-30:30},{0:30})
    bias_sst5_pao(:,{-30:30},{0:30}) = bias_sst5_clm(:,{-30:30},{0:30})
    bias_sst5_iao(:,{-30:30},{0:30}) = bias_sst5_clm(:,{-30:30},{0:30})
    do j = 32, 63 ; -29.36 to 29.36
        count = 0
        do i = 130, 179 ; 260 to 357.5
            if (landsea(j,i).eq.1) then
                count = count + 1
            end if
            if (landsea(j,i).eq.0 .and. count.ge.1) then ;
                bias_sst5_ao(:,j,i) = bias_sst5_clm(:,j,i)
                bias_sst5_pao(:,j,i) = bias_sst5_clm(:,j,i)
                bias_sst5_iao(:,j,i) = bias_sst5_clm(:,j,i)
                ;print("lat "+lat(j)+"   lon "+(360-lon(i))) ;lon(i)
            end if
        end do
    end do
    print("=====================================================")


;CMIP6
bias_sst6_io = bias_sst6_clm
bias_sst6_io = bias_sst6_clm@_FillValue
bias_sst6_po = bias_sst6_io
bias_sst6_ao = bias_sst6_io
bias_sst6_pio = bias_sst6_io
bias_sst6_iao = bias_sst6_io
bias_sst6_pao = bias_sst6_io
    ;==============================================================
    ; Indian Ocean
    ;==============================================================
    bias_sst6_io(:,{-30:30},{30:105}) = bias_sst6_clm(:,{-30:30},{30:105})
    bias_sst6_io(:,{-30:-10},{105:130}) = bias_sst6_clm(:,{-30:-10},{105:130})

    bias_sst6_pio(:,{-30:30},{30:105}) = bias_sst6_clm(:,{-30:30},{30:105})
    bias_sst6_pio(:,{-30:-10},{105:130}) = bias_sst6_clm(:,{-30:-10},{105:130})

    bias_sst6_iao(:,{-30:30},{30:105}) = bias_sst6_clm(:,{-30:30},{30:105})
    bias_sst6_iao(:,{-30:-10},{105:130}) = bias_sst6_clm(:,{-30:-10},{105:130})

    ;==============================================================
    ; Pacific Ocean
    ;==============================================================
    bias_sst6_po(:,{-10:30},{105:130}) = bias_sst6_clm(:,{-10:30},{105:130})
    bias_sst6_po(:,{-30:30},{130:240}) = bias_sst6_clm(:,{-30:30},{130:240})

    bias_sst6_pio(:,{-10:30},{105:130}) = bias_sst6_clm(:,{-10:30},{105:130})
    bias_sst6_pio(:,{-30:30},{130:240}) = bias_sst6_clm(:,{-30:30},{130:240})

    bias_sst6_pao(:,{-10:30},{105:130}) = bias_sst6_clm(:,{-10:30},{105:130})
    bias_sst6_pao(:,{-30:30},{130:240}) = bias_sst6_clm(:,{-30:30},{130:240})

    do j = 32, 63 ; -29.36 to 29.36 56, 57;
        count = 0
        do i = 120, 145 ; 240 to 290
            if (landsea(j,i).lt.1 .and. count.eq.0) then
                bias_sst6_po(:,j,i) = bias_sst6_clm(:,j,i)
                bias_sst6_pio(:,j,i) = bias_sst6_clm(:,j,i)
                bias_sst6_pao(:,j,i) = bias_sst6_clm(:,j,i)
                print("lat "+lat(j)+"   lon "+(360-lon(i))) ;lon(i)
            else
                count = count + 1
            end if
        end do
    end do
    print("=====================================================")
    ;==============================================================
    ; Atlantic Ocean
    ;==============================================================
    bias_sst6_ao(:,{-30:30},{0:30}) = bias_sst6_clm(:,{-30:30},{0:30})
    bias_sst6_pao(:,{-30:30},{0:30}) = bias_sst6_clm(:,{-30:30},{0:30})
    bias_sst6_iao(:,{-30:30},{0:30}) = bias_sst6_clm(:,{-30:30},{0:30})
    do j = 32, 63 ; -29.36 to 29.36
        count = 0
        do i = 130, 179 ; 260 to 357.5
            if (landsea(j,i).eq.1) then
                count = count + 1
            end if
            if (landsea(j,i).eq.0 .and. count.ge.1) then ;
                bias_sst6_ao(:,j,i) = bias_sst6_clm(:,j,i)
                bias_sst6_pao(:,j,i) = bias_sst6_clm(:,j,i)
                bias_sst6_iao(:,j,i) = bias_sst6_clm(:,j,i)
                ;print("lat "+lat(j)+"   lon "+(360-lon(i))) ;lon(i)
            end if
        end do
    end do
    print("=====================================================")


;===================================================================
leftsti = (/"All Basin","Pacific Only","Atlantic Only","Indian Only","Pacific & Atlantic","Pacific & Indian","Atlantic & Indian"/)
line5 = new((/7,12/), float)
line5(0,:) = wgt_areaave_Wrap(bias_sst5_clm(:,{-30:30},:), 1., 1., 0) ;
line5(1,:) = wgt_areaave(bias_sst5_po(:,{-30:30},:), 1., 1., 0)
line5(2,:) = wgt_areaave(bias_sst5_ao(:,{-30:30},:), 1., 1., 0)
line5(3,:) = wgt_areaave(bias_sst5_io(:,{-30:30},:), 1., 1., 0)
line5(4,:) = wgt_areaave(bias_sst5_pao(:,{-30:30},:), 1., 1., 0)
line5(5,:) = wgt_areaave(bias_sst5_pio(:,{-30:30},:), 1., 1., 0)
line5(6,:) = wgt_areaave(bias_sst5_iao(:,{-30:30},:), 1., 1., 0)

line6 = new((/7,12/), float)
line6(0,:) = wgt_areaave_Wrap(bias_sst6_clm(:,{-30:30},:), 1., 1., 0) ;
line6(1,:) = wgt_areaave(bias_sst6_po(:,{-30:30},:), 1., 1., 0)
line6(2,:) = wgt_areaave(bias_sst6_ao(:,{-30:30},:), 1., 1., 0)
line6(3,:) = wgt_areaave(bias_sst6_io(:,{-30:30},:), 1., 1., 0)
line6(4,:) = wgt_areaave(bias_sst6_pao(:,{-30:30},:), 1., 1., 0)
line6(5,:) = wgt_areaave(bias_sst6_pio(:,{-30:30},:), 1., 1., 0)
line6(6,:) = wgt_areaave(bias_sst6_iao(:,{-30:30},:), 1., 1., 0)
;=========================================================
; plot 
;=========================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/model/a20200910/"
figname = season+"_sst_bias_CMIP5_CMIP6_OBS"

wks = gsn_open_wks("pdf", figpath+figname)
plot = new(5, graphic)


res                       = True     ; plot mods desired
res@gsnDraw               = False
res@gsnFrame              = False

res@gsnMaximize           = True     ; Make plot fit the frame
res@cnFillOn              = True     ; turn on color fill
res@cnLinesOn             = False    ; turn of contour lines
;res@cnLevelSpacingF       = 0.5      ; contour spacing
res@cnFillPalette         = "hotcold_18lev"
res@lbLabelBarOn          = False
res@cnInfoLabelOn         = False
res@cnLineLabelsOn        = False

res@mpMinLatF            = -30      ; range to zoom in on
res@mpMaxLatF            =  30.
;res@mpMinLonF            =  30.
;res@mpMaxLonF            = 120.
res@mpCenterLonF         = 180.

res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02
res@gsnStringFontHeightF = 0.025

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -2
res@cnMaxLevelValF = 2
res@cnLevelSpacingF = 0.2


res@gsnLeftString = "(a) CMIP5 sst bias"
plot(0) = gsn_csm_contour_map(wks, bias_sst5_sea, res);
res@gsnLeftString = "(b) CMIP6 sst bias"
plot(1) = gsn_csm_contour_map(wks, bias_sst6_sea, res) ;
res@gsnLeftString = "(c) b-a"
plot(2) = gsn_csm_contour_map(wks, diff, res) ;

pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot(0:2), (/3,1/), pres)

delete(res)
Month = ispan(1, 12, 1)
res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = 1
res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = Month
res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@tmXTOn = False
res@tiXAxisString = "Month"

res@trYMinF = -1
res@trYMaxF = 1
res@gsnYRefLine = 0
res@gsnYRefLineColor = "gray40"
res@gsnYRefLineDashPattern = 2
res@gsnYRefLineThicknessF = 5
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-1, 1, 0.1)
; res@tmYLLabels = ispan(-1, 1, 0.1)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"
;res@gsnLeftString = "Tas: grid2station"

Colors = (/"black","red","blue","forestgreen","purple","orange2","deeppink2"/);
; Colors = (//)
; Patterns = (/0,0,0/)
res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPattern = 0

plot(3) = gsn_csm_xy(wks, Month, line5, res)


res@xyDashPattern = 16
plot(4) = gsn_csm_xy(wks, Month, line6, res)
overlay(plot(3), plot(4))
;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 8                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = leftsti

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    interval = 0.15
    gy1 = (/0.8,0.8/)
    gy2 = gy1-interval
    gy3 = gy2-interval
    gy4 = gy3-interval
    
    gx2 = (/5,5.7/)
    text_pos2 = 5.8

    res_lines@gsLineColor = Colors(0)
    ; res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(3), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(3), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ; res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(3), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(3), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    ; res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(3), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(3), Texts(2), text_pos, gy3, res_text)
    res_lines@gsLineColor = Colors(3)
    dum41 = gsn_add_polyline(wks, plot(3), gx, gy4, res_lines)
    dum42 = gsn_add_text(wks, plot(3), Texts(3), text_pos, gy4, res_text)

    res_lines@gsLineColor = Colors(4)
    dum51 = gsn_add_polyline(wks, plot(3), gx2, gy1, res_lines)
    dum52 = gsn_add_text(wks, plot(3), Texts(4), text_pos2, gy1, res_text)
    res_lines@gsLineColor = Colors(5)
    dum61 = gsn_add_polyline(wks, plot(3), gx2, gy2, res_lines)
    dum62 = gsn_add_text(wks, plot(3), Texts(5), text_pos2, gy2, res_text)
    res_lines@gsLineColor = Colors(6)
    dum71 = gsn_add_polyline(wks, plot(3), gx2, gy3, res_lines)
    dum72 = gsn_add_text(wks, plot(3), Texts(6), text_pos2, gy3, res_text)

draw(plot(3))

end