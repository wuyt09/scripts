
begin

yrStrt = 1980
yrLast = 2005
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/) ;ispan(mnStrt, mnLast, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
season = (/"DJF","JJA"/)
nseason = dimsizes(season)

lat = ispan(-90, 90, 1)
lat@units = "degree_north"
lon = ispan(0, 359, 1)
lon@units = "degree_east"
; print(lat)
; print(lon)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
LEV = (/1000,925,850,700,600,500,400,300,250,200,150,100,50,10/)
LEV@units = "hPa"
LEV2 = LEV*100
LEV2@units = "Pa"
nlev = dimsizes(LEV)


lonstr = 75 ; 72 ;75
lonend = 105 ;136 ;105
latstr = 25 ;18 ;25
latend = 40 ;54 ; 40

;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo(::-1,:)
printVarSummary(topo)

;===================================================================
; OBS grid data
;===================================================================

t = new((/nyear*nmonth,nlev,nlat,nlon/), float)

diri = "/home/yangsong3/data-observation/MERRA2/monthly/t/"

do i = 0, nyear-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        fili = systemfunc("ls "+diri+"MERRA2_*.instM_3d_asm_Np."+yr+mn+".SUB.nc4")
        ; print(fili)
        f = addfile(fili, "r")
        T = f->T(0,{LEV},:,:)

        t(i*12+j,:,:,:) = linint2_Wrap(T&lon, T&lat, T, True, lon, lat, 0)
    end do
    print("end of year "+year(i))
end do
t!0 = "time"
t&time = ispan(1, nyear*nmonth, 1)
t!1 = "lev"
t&lev = LEV
t!2 = "lat"
t&lat = lat
t!3 = "lon"
t&lon = lon
printVarSummary(t)

tp_obs = t(:,:,{latstr:latend},{lonstr:lonend})
TP_obsx = month_to_seasonN(tp_obs, season) ;clmMonTLL(t)
TP_obs = dim_avg_n_Wrap(TP_obsx, 1)
printVarSummary(TP_obs)

; topox = linint2_Wrap(topo&lon, topo&lat, topo, False, t&lon, t&lat, 0)
; print(topox({35},{90}))

; TS_obs = mask(TS_obs, topox.ge.1500, True)
; TP_obs = wgt_areaave_Wrap(TS_obs(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
; printVarSummary(TP_obs)

;===================================================================
; GCM 20 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/ta/"

models1 = (/"ACCESS1-0","bcc-csm1-1","BNU-ESM","CanCM4","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CESM","CNRM-CM5","GFDL-ESM2M","GISS-E2-H","CMCC-CM","CMCC-CMS","GFDL-CM3","GFDL-ESM2G",\
        "IPSL-CM5A-MR","NorESM1-ME","GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR",\
        "IPSL-CM5B-LR","MIROC5","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M"/)
;"CESM1-WACCM","FIO-ESM","MPI-ESM-P","CSIRO-Mk3-6-0",
;

models2 = (/"ACCESS1-0","ACCESS1-3","bcc-csm1-1","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",\
        "GISS-E2-H","GISS-E2-R","HadGEM2-CC","HadGEM2-ES","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

models = venn2_intersection(models1,models2)

nm = dimsizes(models)
t_gcm = new((/2,latend-latstr+1,lonend-lonstr+1/), float)

TP_final = new((/nm,2/), float)
TP_final!0 = "model"
TP_final&model = models

diff = new((/nlev/), float)
diff = 1
last = diff


do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"/ta_Amon_"+models(i)+"_historical_r1i1p1_1*-20*12.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    T_GCMs := fgcm->ta(iYYYY,{LEV2},:,:)
    ; printVarSummary(T_GCMs)
    T_GCMs!1 = "lev"
    T_GCMs!2 = "lat"
    T_GCMs!3 = "lon"

    T_GCM := linint2_Wrap(T_GCMs&lon, T_GCMs&lat, T_GCMs, True, lon, lat, 0)
    printVarSummary(T_GCM)


    plev := tofloat(T_GCM&lev)
    ptop = min(plev)
    psfc = max(plev)
    dP := dpres_plevel(plev, psfc, ptop, 0)
    dP!0 = "lev"
    dP&lev = T_GCM&lev


    TP_GCM = T_GCM(:,:,{latstr:latend},{lonstr:lonend})
    TP_GCM_seax = month_to_seasonN(TP_GCM, season)
    TP_GCM_sea = dim_avg_n_Wrap(TP_GCM_seax, 1)

    TP_diff = TP_GCM_sea
    TP_diff = TP_GCM_sea - TP_obs

    dp = conform_dims(dimsizes(TP_diff), dP, 1)
    tgcmdp = TP_diff*dp
    copy_VarCoords(TP_diff, tgcmdp)


    do isea = 0, nseason-1
        do ilat = 0, dimsizes(TP_obs&lat)-1
            do ilon = 0, dimsizes(TP_obs&lon)-1
                var = TP_diff(isea,:,ilat,ilon)
                do ilev = 1, nlev-2
                    diff(ilev) = var(ilev) * var(ilev-1)
                    last(ilev) = var(ilev) * var(ilev+1)
                end do
                a := ind((diff.ge.0 .and. var.le.0 .and. last.le.0))
                lev_zero = (/a(0)/)
                if (ismissing(lev_zero)) then
                    lev_zero = nlev-1
                end if
                print(LEV2(lev_zero))
                x_med = dim_sum_n_Wrap(tgcmdp(isea,{100000:LEV2(lev_zero)},ilat,ilon), 0)/sum(dP({100000:LEV2(lev_zero)}))
                t_gcm(isea,ilat,ilon) = x_med
            end do
        end do
        print(lev_zero)
    end do
    copy_VarCoords(TP_GCM_sea(:,0,:,:), t_gcm)

    t_gcm_mask = mask(t_gcm, topo({latstr:latend},{lonstr:lonend}).ge.1500, True)

    TP_final(i,:) = wgt_areaave_Wrap(t_gcm_mask, 1., 1., 0)

    delete([/T_GCMs,T_GCM,plev,dP,dp,tgcmdp,TP_GCM,TP_GCM_seax,TP_GCM_sea,t_gcm_mask/])
end do
printVarSummary(t_gcm)

;========================================================
Solid = new((/2,nm/), float)
Solid(0,:) = TP_final(:,0)
Solid(1,:) = TP_final(:,1)


write_table("write_Tropo_T_bias.txt", "w", [/"DJF","JJA"/], "%4s %4s")
write_table("write_Tropo_T_bias.txt", "a", [/Solid(0,:),Solid(1,:)/], "%4.2f %4.2f")

Month = ispan(1, 12, 1)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/20201228/"
figname = "write_Tropo_T_bias"
pwks = "png"
pwks@wkHeight              = 3048
pwks@wkWidth               = 6048
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

; res@trXMinF = 1
; res@trXMaxF = 12
res@tmXBMode = "Explicit"
res@tmXBValues = ispan(1, nm, 1)
res@tmXBLabels = models
res@tmXTOn = False
res@tmXBLabelAngleF = 45
res@tmXBLabelFontHeightF = 0.007
; res@tiXAxisString = "Month"

; res@trYMinF = -8
; res@trYMaxF = 4
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-5, 4, 1)
res@tmYLPrecision = 2
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYROn = False
res@tiYAxisString = "Diff Tropo T"

res@tmYROn = False


;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","red","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, ispan(1, nm, 1), Solid, res)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = season

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/257,257/)
    interval = 2
    gy2 = gy1 - interval
    gy3 = gy2 -interval
    gy4 = gy3 - interval


    res_lines@gsLineThicknessF = 15
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)


draw(plot(0))
frame(wks)

end