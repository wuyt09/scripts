
begin

yrStrt = 1980
yrLast = 2005
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/) ;ispan(mnStrt, mnLast, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)

nlat = 361
nlon = 576
nlev = 42


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

t = new((/nyear*nmonth,nlat,nlon/), float)

diri = "/home/yangsong3/data-observation/MERRA2/monthly/t/"

do i = 0, nyear-1
    yr = yrStrt+i
    do j = 0, 11
        mn = sprinti("%0.2i", j+1)
        fili = systemfunc("ls "+diri+"MERRA2_*.instM_3d_asm_Np."+yr+mn+".SUB.nc4")
        ; print(fili)
        f = addfile(fili, "r")
        T = f->T(0,:,:,:)

        plev = tofloat(T&lev)
        plev = plev*100
        plev@units = "Pa"
        psfc = max(plev)
        ptop = min(plev)

        dP = dpres_plevel(plev, psfc, ptop, 0)
        dP!0 = "lev"
        dP&lev = T&lev
        dp = conform_dims(dimsizes(T), dP, 0)
        tdp = T*dp
        copy_VarCoords(T, tdp)
        t(i*12+j,:,:) = dim_sum_n_Wrap(tdp({500:200},:,:), 0)/sum(dP({500:200}))
        copy_VarCoords(tdp(0,:,:), t(0,:,:))
    end do
    print("end of year "+year(i))
end do
printVarSummary(t)
t!0 = "time"
t&time = ispan(1, nyear*nmonth, 1)

TS_obs = clmMonTLL(t)

topox = linint2_Wrap(topo&lon, topo&lat, topo, False, t&lon, t&lat, 0)
print(topox({35},{90}))

TS_obs = mask(TS_obs, topox.ge.1500, True)
TP_obs = wgt_areaave_Wrap(TS_obs(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
printVarSummary(TP_obs)

;===================================================================
; GCM 20 data
;===================================================================
dirigcm = "/home/yangsong3/data-model/wuyt/CMIP5/aaa-post/ta/"

models = (/"ACCESS1-0","bcc-csm1-1","BNU-ESM","CanCM4","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5",\
        "CESM1-WACCM","CMCC-CESM","CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","FIO-ESM","GFDL-CM3",\
        "GFDL-ESM2G","GFDL-ESM2M","GISS-E2-H","GISS-E2-R","HadCM3","inmcm4","IPSL-CM5A-LR","IPSL-CM5A-MR",\
        "IPSL-CM5B-LR","MIROC5","MIROC-ESM","MPI-ESM-LR","MPI-ESM-MR","MPI-ESM-P","MRI-CGCM3","NorESM1-M","NorESM1-ME"/)

nm = dimsizes(models)
TP_GCM5 = new((/nm,12/), float)
TP_GCM5!0 = "model"
TP_GCM5&model = models
do i = 0, nm-1
    filigcm = systemfunc("ls "+dirigcm+"/ta_Amon_"+models(i)+"_historical_r1i1p1_1*-20*12.nc")
    fgcm = addfile(filigcm, "r")
    TIME := fgcm->time
    YYYY := cd_calendar(TIME, -1)/100
    iYYYY := ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    T_GCMs := fgcm->ta(iYYYY,:,:,:)
    ; printVarSummary(T_GCMs)
    T_GCMs!1 = "lev"
    T_GCMs!2 = "lat"
    T_GCMs!3 = "lon"

    T_GCM := linint2_Wrap(T_GCMs&lon, T_GCMs&lat, T_GCMs, True, t&lon, t&lat, 0)
    printVarSummary(T_GCM)

    plev := tofloat(T_GCM&lev)
    ptop = min(plev)
    psfc = max(plev)
    dP := dpres_plevel(plev, psfc, ptop, 0)
    dP!0 = "lev"
    dP&lev = T_GCM&lev

    dp := conform_dims(dimsizes(T_GCM), dP, 1)
    tgcmdp := T_GCM*dp
    copy_VarCoords(T_GCM, tgcmdp)
    t_gcm := dim_sum_n_Wrap(tgcmdp(:,{50000:20000},:,:), 1)/sum(dP({50000:20000}))
    copy_VarCoords(T_GCM(:,0,:,:), t_gcm)

    t_gcm_clm := clmMonTLL(t_gcm)
    copy_VarCoords(T_GCM(0,0,:,:), t_gcm_clm(0,:,:))
    t_gcm_clm := mask(t_gcm_clm, topox.ge.1500, True)
    TP_GCM5(i,:) := wgt_areaave_Wrap(t_gcm_clm(:,{latstr:latend},{lonstr:lonend}), 1, 1, 0)
    print("end of "+models(i))
end do
printVarSummary(TP_GCM5)

TP_GCMA5 = dim_avg_n_Wrap(TP_GCM5, 0)
printVarSummary(TP_GCMA5)

print(TP_GCMA5)
;========================================================
Solid = new((/3,12/), float)
Solid(0,:) = (/TP_obs/)
Solid(1,:) = (/TP_GCMA5/)


Month = ispan(1, 12, 1)
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_coldbias/20201228/"
figname = "Tropo_T_bias_MME"
pwks = "png"
pwks@wkHeight              = 3048
pwks@wkWidth               = 6048
wks = gsn_open_wks("pdf", figpath+figname)

plot = new(2, graphic)

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
res@tmXBLabels = month
res@tmXTOn = False
; res@tiXAxisString = "Month"

; res@trYMinF = -20
; res@trYMaxF = 20
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-20, 20, 5)
; res@tmYLLabels = ispan(-20, 20, 5)
res@tmYROn = False
res@tiYAxisString = "Temperature (~S~o~N~C)"
res@gsnLeftString = "Tropo T"


;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","red","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, Month, Solid, res)

delete(res@xyLineColors)
delete(res@xyDashPatterns)
res@xyLineThicknessF = 7

gsn_define_colormap(wks, "WhiteYellowOrangeRed")
res@xyLineColors = ispan(25, 25+28*8, 8)
res@xyDashPattern = 16
plot(1) = gsn_csm_xy(wks, Month, TP_GCM5, res)
overlay(plot(1), plot(0))


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = (/"MERRA2","CMIP5-MME"/)

;============================================================
    gx = (/1.5,2.2/)
    text_pos = 2.3
    gy1 = (/254,254/)
    interval = 2
    gy2 = gy1 - interval
    gy3 = gy2 -interval
    gy4 = gy3 - interval


    res_lines@gsLineThicknessF = 15
    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(1), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(1), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(1), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(1), Texts(1), text_pos, gy2, res_text)
    ; res_lines@gsLineColor = Colors(2)
    ; res_lines@gsLineDashPattern = Patterns(2)
    ; dum31 = gsn_add_polyline(wks, plot(1), gx, gy3, res_lines)
    ; dum32 = gsn_add_text(wks, plot(1), Texts(2), text_pos, gy3, res_text)

    res_lines@gsLineColor = "red3"
    res_lines@gsLineThicknessF = 7
    res_lines@gsLineDashPattern = 16
    dum41 = gsn_add_polyline(wks, plot(1), gx, gy3, res_lines)
    dum42 = gsn_add_text(wks, plot(1), "CMIP5 GCMs", text_pos, gy4, res_text)


draw(plot(1))
;draw(plot(1))
frame(wks)

end