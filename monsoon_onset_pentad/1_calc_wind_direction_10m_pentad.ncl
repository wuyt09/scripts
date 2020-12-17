
begin

yrStrt = 1979
yrLast = 2018
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)

nlat = 181
nlon = 360
nlev = 19


diriu = "/home/yangsong3/data-observation/ERA5-daily/Single-level/u10m/"
diriv = "/home/yangsong3/data-observation/ERA5-daily/Single-level/v10m/"

U = new((/365,nlat,nlon/), float)
u = new((/npenta,nlat,nlon/), float)
V = new((/365,nlat,nlon/), float)
v = new((/npenta,nlat,nlon/), float)

theta = new((/nyear,npenta,nlat,nlon/), float)
r2d = 45.0/atan(1.0)

do i = 0, nyear-1
    irecS = 0
    do j = 0, nmonth-1
        mth = sprinti("%0.2i", month(j))
        filiu = "u10m."+year(i)+"-"+mth+".daily.nc"
        filiv = "v10m."+year(i)+"-"+mth+".daily.nc"

        fu = addfile(diriu+filiu, "r")
        UX := short2flt(fu->u10(:,::-1,:))
        fv = addfile(diriv+filiv, "r")
        VX := short2flt(fv->v10(:,::-1,:))

        ndays := dimsizes(UX&time)
        if (ndays.eq.29) then
            irecL = irecS+ndays-2
            U(irecS:irecL,:,:) = UX(0:27,:,:)
            V(irecS:irecL,:,:) = VX(0:27,:,:)
            irecS = irecS+ndays-1
        else
            irecL = irecS+ndays-1
            U(irecS:irecL,:,:) = UX
            V(irecS:irecL,:,:) = VX
            irecS = irecS+ndays
        end if
        ; print(irecL)
    end do

    do p = 0, 72
        u(p,:,:) = dim_sum_n_Wrap(U(5*p:5*p+4,:,:), 0)/5.0
        v(p,:,:) = dim_sum_n_Wrap(V(5*p:5*p+4,:,:), 0)/5.0
    end do
    copy_VarCoords(U(0,:,:), u(0,:,:))
    copy_VarCoords(V(0,:,:), v(0,:,:))

    theta(i,:,:,:) = atan2(u, v) * r2d + 180
    copy_VarCoords(U(0,:,:), theta(0,0,:,:))
    print("end of year "+year(i))
end do
theta!0 = "year"
theta&year = year
theta!1 = "pentad"
theta&pentad = penta
printVarSummary(theta)

theta@units = "degree"
theta@long_name = "10m wind direction derived from ERA5 daily data"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "1_wind_direction_10m_pentad_ERA5_"+yrStrt+"_"+yrLast+".nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->theta = theta

print("end of writing data")
;================================================================
; plot
;================================================================
u := theta
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

u850_bob = wgt_areaave_Wrap(u(:,:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
u850_scs = wgt_areaave_Wrap(u(:,:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
u850_ism = wgt_areaave_Wrap(u(:,:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)

x = new((/6,npenta/), float)
x(0,:) = u850_bob(0,:)
x(1,:) = u850_scs(0,:)
x(2,:) = u850_ism(0,:)
x(3,:) = dim_avg_n_Wrap(u850_bob, 0)
x(4,:) = dim_avg_n_Wrap(u850_scs, 0)
x(5,:) = dim_avg_n_Wrap(u850_ism, 0)

;====================================================================
; plot
;====================================================================
figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "1_calc_wind_direction_10m_pentad"
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

res@gsnYRefLine = (/90,180,270/)
res@gsnYRefLineColor = "grey40"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2

res@gsnRightString = ""
res@gsnStringFontHeightF = 0.025

plot(0) = gsn_csm_xy(wks, penta, x, res)
draw(plot)

end
