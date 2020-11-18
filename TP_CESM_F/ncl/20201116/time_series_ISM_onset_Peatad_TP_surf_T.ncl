
begin

yrStrt = 1979
yrLast = 2017
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)

nlat = 181
nlon = 360
nlev = 19

svarlist = (/"t2m"/)
mvarlist = (/"temp","uwind","q"/)

sdiri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/"
mdiri = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/"

t2m = new((/73,nlat,nlon/), float)
t = new((/73,nlev,nlat,nlon/), float)
u = t
q = t
ttgi = new((/nyear,73/), float)
wyi = ttgi
howi = ttgi
wfi = ttgi
ts = ttgi

do i = 0, nyear-1
    filist = systemfunc("ls "+sdiri+svarlist(0)+"/*."+year(i)+"*.daily.nc")
    print(filist)
    fst = addfiles(filist, "r")
    T2MX = short2flt(fst[:]->t2m(:,::-1,:))

    if (dimsizes(T2MX&time).eq.366) then
        T2M = T2MX(0:364,:,:)
        T2M(0:58,:,:) = T2MX(0:58,:,:)
        T2M(59:364,:,:) = T2MX(60:365,:,:)
    else
        T2M = T2MX
    end if
    printVarSummary(T2M)
    delete([/T2MX/])

    ; filips = systemfunc("ls "+sdiri+svarlist(1)+"/*."+year(i)+"*.daily.nc")
    ; print(filips)
    ; fps = addfiles(filips, "r")
    ; PSX = short2flt(fst[:]->$svarlist(1)$)

    ; if (dimsizes(PSX&time).eq.366) then
    ;     PS = PSX(0:364,:,:)
    ;     PS(0:58,:,:) = PSX(0:58,:,:)
    ;     PS(59:364,:,:) = PSX(60:365,:,:)
    ; else
    ;     PS = PSX
    ; end if
    ; printVarSummary(PS)
    ; delete([/PSX/])

    filit = systemfunc("ls "+mdiri+mvarlist(0)+"/*."+year(i)+"*.daily.nc")
    ft = addfiles(filit, "r")    
    TX = short2flt(ft[:]->t(:,:,::-1,:))

    if (dimsizes(TX&time).eq.366) then
        T = TX(0:364,:,:,:)
        T(0:58,:,:,:) = TX(0:58,:,:,:)
        T(59:364,:,:,:) = TX(60:365,:,:,:)
    else
        T = TX
    end if
    printVarSummary(T)
    delete([/TX/])

    filiu = systemfunc("ls "+mdiri+mvarlist(1)+"/*."+year(i)+"*.daily.nc")
    fu = addfiles(filiu, "r")    
    UX = short2flt(fu[:]->u(:,:,::-1,:))

    if (dimsizes(UX&time).eq.366) then
        U = UX(0:364,:,:,:)
        U(0:58,:,:,:) = UX(0:58,:,:,:)
        U(59:364,:,:,:) = UX(60:365,:,:,:)
    else
        U = UX
    end if
    printVarSummary(U)
    delete([/UX/])

    filiq = systemfunc("ls "+mdiri+mvarlist(2)+"/*."+year(i)+"*.daily.nc")
    fq = addfiles(filiq, "r")    
    QX = short2flt(fq[:]->q(:,:,::-1,:))

    if (dimsizes(QX&time).eq.366) then
        Q = QX(0:364,:,:,:)
        Q(0:58,:,:,:) = QX(0:58,:,:,:)
        Q(59:364,:,:,:) = QX(60:365,:,:,:)
    else
        Q = QX
    end if
    printVarSummary(Q)
    delete([/QX/])

    do p = 0, 72
        ; t2m(i,p,:,:) = dim_sum_n_Wrap(T2M(5*p:5*p+4,:,:), 0)/5.0
        ; t(i,p,:,:,:) = dim_sum_n_Wrap(T(5*p:5*p+4,:,:,:), 0)/5.0
        ; u(i,p,:,:,:) = dim_sum_n_Wrap(U(5*p:5*p+4,:,:,:), 0)/5.0
        ; q(i,p,:,:,:) = dim_sum_n_Wrap(Q(5*p:5*p+4,:,:,:), 0)/5.0

        t2m(p,:,:) = dim_sum_n_Wrap(T2M(5*p:5*p+4,:,:), 0)/5.0
        t(p,:,:,:) = dim_sum_n_Wrap(T(5*p:5*p+4,:,:,:), 0)/5.0
        u(p,:,:,:) = dim_sum_n_Wrap(U(5*p:5*p+4,:,:,:), 0)/5.0
        q(p,:,:,:) = dim_sum_n_Wrap(Q(5*p:5*p+4,:,:,:), 0)/5.0
    end do
    copy_VarCoords(T2M(0,:,:), t2m(0,:,:))
    copy_VarCoords(T(0,:,:,:), t(0,:,:,:))
    copy_VarCoords(U(0,:,:,:), u(0,:,:,:))
    copy_VarCoords(Q(0,:,:,:), q(0,:,:,:))


    ;===============================================================
    ; T2m
    ;===============================================================
    ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
    topo = ftopo->topo(::-1,:)
    t2m = mask(t2m, topo.ge.1500, True)
    ts(i,:) = wgt_areaave(t2m(:,{25:40},{75:105}), 1, 1, 0)
    copy_VarCoords(t2m(:,0,0), ts(i,:))




    lev = tofloat(t&level)
    plev = lev*100
    plev@units = "Pa"
    psfc = max(plev)
    ptop = min(plev)

    ;===============================================================
    ; TTGI
    ;===============================================================
    dP = dpres_plevel(plev, psfc, ptop, 0)
    dP!0 = "lev"
    dP&lev = lev
    dp = conform_dims(dimsizes(t), dP, 1)
    tdp = t*dp
    copy_VarCoords(t, tdp)
    t_tropo = dim_sum_n_Wrap(tdp(:,{600:200},:,:), 1)/sum(dP({600:200}))
    copy_VarCoords(tdp(:,0,:,:), t_tropo)
    tN = wgt_areaave_Wrap(t_tropo(:,{5:35},{40:100}), 1, 1, 0)
    tS = wgt_areaave_Wrap(t_tropo(:,{-15:5},{40:100}), 1, 1, 0)
    ttgi(i,:) = tN - tS
    copy_VarCoords(tN, ttgi(i,:))

    ;===============================================================
    ; WYI
    ;===============================================================
    u_wgt = wgt_areaave_Wrap(u(:,:,{0:20},{40:110}), 1, 1, 0)
    wyix = u_wgt(:,{850}) - u_wgt(:,{200})
    u_May30 = wgt_areaave_Wrap(U(149,:,{0:20},{40:110}), 1, 1, 0)
    u_wgt_May30 = u_May30({850}) - u_May30({200})
    wyi(i,:) = wyix - u_wgt_May30
    copy_VarCoords(u_wgt(:,0), wyi(i,:))

    ; ;===============================================================
    ; ; HOWI
    ; ;===============================================================
    ; qudp = q*u*dp
    ; copy_VarCoords(q, qudp)
    ; qdp = dim_sum_n_Wrap(qudp(:,{1000:300},:,:), 1)
    ; qdp_wgt = wgt_areaave_Wrap(qdp(:,{5:20},{45:80}), 1, 1, 0)

    ; qdpmax = max(qdp_wgt)
    ; qdpmin = min(qdp_wgt)

    ; howi(i,:) = 2*((qdp_wgt - qdpmin)/(qdpmax - qdpmin))-1
    ; copy_VarCoords(qdp_wgt, howi(i,:))

    ;===============================================================
    ; WFI
    ;===============================================================
    u_850 = u(:,{850},:,:)
    uS = wgt_areaave_Wrap(u_850(:,{5:15},{40:80}), 1, 1, 0)
    uN = wgt_areaave_Wrap(u_850(:,{20:30},{70:90}), 1, 1, 0)
    wfi(i,:) = uS - uN
    copy_VarCoords(uS, wfi(i,:))
print("end of "+year(i))
end do
printVarSummary(wyi)

inx_name = (/"TP_t2m","TTGI","WYI","WFI"/) ;,"HOWI"
pp = new((/4,nyear,73/), float)
pp(0,:,:) = ts
pp(1,:,:) = ttgi
pp(2,:,:) = wyi
; pp(3,:,:) = howi
pp(3,:,:) = wfi
pp!0 = "index"
pp&index = inx_name
pp!1 = "year"
pp&year = year
pp!2 = "pentad"
pp&pentad = ispan(1, 73, 1)


ptime = new((/4,nyear/), integer)

diff = pp
do i = 1, 72
    diff(:,:,i) = pp(:,:,i)*pp(:,:,i-1)    
end do

do ll = 1,3
    do mm = 0, nyear-1
        a = ind((diff(ll,mm,:).le.0) .and. (pp(ll,mm,:).gt.0))
        print(ll+", "+mm+", "+a)
        ptime(ll,mm) = (/a(0)/)
        delete(a)
    end do
end do
copy_VarCoords(pp(:,:,0), ptime)

ts_spr = dim_avg_n_Wrap(ts(:,12:29), 1) ; spring
copy_VarCoords(ts(:,0), ts_spr)

;===============================================================================
; wtite data to txt file
;===============================================================================
filipath = "./"
filiname = "ISM_onset_pentad.txt"

write_table(filipath+filiname, "w", [/"TP_t2m","TTGI","WYI","WFI"/], "%10s %10s %10s %10s")
write_table(filipath+filiname, "a", [/ts_spr, ptime(1,:),ptime(2,:),ptime(3,:)/], "%10.2f %0.4i %0.4i %0.4i")

corr = new(3, float)

do kk = 1,3
    corr(kk-1) = escorc(ts_spr, ptime(kk,:))
end do
print(corr)
;=============================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20201116/"
figname = "time_series_ISM_onset_Pentad_TP_surf_T"



wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

res@trXMinF = yrStrt
res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = Month
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
; res@tmXTOn = False
; res@tiXAxisString = "Month"

res@trYMinF = 22
res@trYMaxF = 36
res@tmYLMode = "Explicit"
res@tmYLValues = ispan(22, 36, 1)
res@tmYLLabels = "P"+res@tmYLValues
; res@tmYROn = False
; res@tiYAxisString = "Temperature (~S~o~N~C)"

;Colors = (/"black","deeppink2", "orange2", "green4", "black" /);
Colors = (/"red","blue","forestgreen","purple"/)
Patterns = (/0,0,0,0/)
res@xyLineThicknessF = 10
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

; plot(0) = gsn_csm_xy(wks, year, ptime, res)

res2 = True
res2@gsnDraw  = False
res2@gsnFrame = False
res2@gsnMaximize = True

res2@vpWidthF = 0.7
res2@vpHeightF = 0.4

res2@trXMinF = yrStrt
res2@trXMaxF = yrLast

res2@xyLineThicknessF = 10
res2@xyLineColor = "black"
res2@xyDashPattern = 0

plot(0) = gsn_csm_xy2(wks, year, ptime(1:3,:), ts_spr, res, res2)

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.018                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

    Texts  = inx_name(1:dimsizes(inx_name)-1)

;============================================================
    gx = (/1983,1985/)
    text_pos = 1985.5

    interval = 1
    gy1 = (/35,35/)
    gy2 = gy1 - interval
    gy3 = gy2 - interval
    ; gy4 = gy3 - interval

    gx2 = (/1991,1993/)
    text_pos2 = 1993.5
    gy5 = gy1
    ; gy6 = (/8,8/)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)
    ; res_lines@gsLineColor = Colors(3)
    ; res_lines@gsLineDashPattern = Patterns(3)
    ; dum41 = gsn_add_polyline(wks, plot(0), gx, gy4, res_lines)
    ; dum42 = gsn_add_text(wks, plot(0), Texts(3), text_pos, gy4, res_text)
    res_lines@gsLineColor = "black"
    res_lines@gsLineDashPattern = 0
    dum51 = gsn_add_polyline(wks, plot(0), gx2, gy5, res_lines)
    dum52 = gsn_add_text(wks, plot(0), "TP_t2m", text_pos2, gy5, res_text)

    draw(plot(0))
end