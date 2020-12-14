begin

yrStrt = 1979
yrLast = 2015
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)

nlat = 121
nlon = 240
nlev = 37


f3  = addfile("/home/yangsong3/data-observation/wuyt/ERAI/monthly/1.5x1.5/sp_era_mn_197901_201507.nc","r")
PS0 = short2flt(f3->sp(0:431,:,:))

ps0 = clmMonTLL(PS0)


mdiri = "/home/yangsong3/data-observation/ERA-interim/pressure/air/"

t = new((/npenta,nlev,nlat,nlon/), float)
mtg = new((/nyear,npenta,nlat,nlon/), float)

do i = 0, 1 ;nyear-1
    filit = "air.interim."+year(i)+".nc"
    ft = addfile(mdiri+filit, "r")    
    TX = short2flt(ft->air(:,::-1,::-1,:))

    if (dimsizes(TX&time).eq.366) then
        T = TX(0:364,:,:,:)
        T(0:58,:,:,:) = TX(0:58,:,:,:)
        T(59:364,:,:,:) = TX(60:365,:,:,:)
    else
        T = TX
    end if
    ; printVarSummary(T)
    delete([/TX/])

    do p = 0, 72
        t(p,:,:,:) = dim_sum_n_Wrap(T(5*p:5*p+4,:,:,:), 0)/5.0
    end do
    copy_VarCoords(T(0,:,:,:), t(0,:,:,:))
    t!0 = "pentad"
    t&pentad = penta
    copy_VarCoords(t(:,0,:,:), mtg(0,:,:,:))
    ;-------------------------------------------------------------------
    ; partial T / partial y
    rearth = 6371229.0 ; meters

    lon = ft->lon
    dlon = (lon(2) - lon(1))*0.0174533
    dY = tofloat(rearth*dlon) 
    ;dX = cos(lat*0.0174533)*dlon
    printVarSummary(dY)

    MTG = t
    MTG = center_finite_diff_n(t, dY, True, 0, 2)


    lev = tofloat(t&level)
    plev = lev*100
    plev@units = "Pa"
    psfc = max(plev)
    ptop = min(plev)

    ; dP = dpres_plevel(plev, psfc, ptop, 0)
    ; dP!0 = "lev"
    ; dP&lev = lev
    ; dp = conform_dims(dimsizes(MTG), dP, 1)
    ; mtgdp = MTG*dp
    ; copy_VarCoords(MTG, mtgdp)
    ; mtg(i,:,:,:) = dim_sum_n_Wrap(mtgdp(:,{500:200},:,:), 1)/sum(dP({500:200}))
    ; copy_VarCoords(mtgdp(:,0,:,:), mtg(i,:,:,:))
    linlog =1
    pbot   =50000
    ptop   =20000

    ppppp = ps0(6,:,:)
    do n=0, npenta-1
        kkkkk=MTG(pentad|n,lat|:,lon|:,level|:)
        ;ppppp=p_ClmPnt(26,:,:)
;       printVarSummary(kkkkk)
;       printVarSummary(ppppp)
        mtg(i,n,:,:) = vibeta(plev,kkkkk,linlog,ppppp,pbot,ptop)
    end do

end do
mtg!0 = "year"
mtg&year = year
mtg@units = "K"
mtg@long_name = "troposphere (500hPa-200hPa) temperature gradient derived from ERA5 daily data"
printVarSummary(mtg)



outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/ERAI/"
outname = "1_MTG_pentad_ERAI_1979_2015.nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->MTG = mtg

;================================================================
; plot
;================================================================

bob_lat1 = 5
bob_lat2 = 15
bob_lon1 = 85
bob_lon2 = 100

scs_lat1 = 5
scs_lat2 = 20
scs_lon1 = 110
scs_lon2 = 120

ism_lat1 = 10
ism_lat2 = 20
ism_lon1 = 60
ism_lon2 = 75

mtg_bob = wgt_areaave_Wrap(mtg(:,:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
mtg_scs = wgt_areaave_Wrap(mtg(:,:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
mtg_ism = wgt_areaave_Wrap(mtg(:,:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)

x = new((/6,npenta/), float)
x(0,:) = mtg_bob(0,:)
x(1,:) = mtg_scs(0,:)
x(2,:) = mtg_ism(0,:)
x(3,:) = dim_avg_n_Wrap(mtg_bob, 0)
x(4,:) = dim_avg_n_Wrap(mtg_scs, 0)
x(5,:) = dim_avg_n_Wrap(mtg_ism, 0)

figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/ERAI/"
figname = "1_calc_MTG_pentad"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)
Colors = (/"black","red","blue","black","red","blue"/)
Patterns = (/0,0,0,2,2,2/)

res = True
res@gsnFrame = False
res@gsnDraw  = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.45

res@xyLineThicknessF = 5
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@trXMinF = min(penta)
res@trXMaxF = max(penta)

res@gsnYRefLine = 0
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

res@gsnRightString = ""
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, penta, x, res)
draw(plot)

end

