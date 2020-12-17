begin
; this pentad data are quite different with GPCP pentad data, careful


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


mdiri = "/home/yangsong3/data-observation/ERA5-daily/Single-level/prect/"
PREC = new((/365,nlat,nlon/), float)
prec = new((/nyear,npenta,nlat,nlon/), float)

do i = 0, nyear-1
    irecS = 0
    do j = 0, nmonth-1
        mth = sprinti("%0.2i", month(j))
        filip = "prect."+year(i)+"-"+mth+".daily.nc"

        fp = addfile(mdiri+filip, "r")    
        UX := short2flt(fp->tp(:,::-1,:))
        UX = UX * 24 * 1000 ; convert from m/hour to mm/day 

        ndays := dimsizes(UX&time)
        if (ndays.eq.29) then
            irecL = irecS+ndays-2
            PREC(irecS:irecL,:,:) = UX(0:27,:,:)
            irecS = irecS+ndays-1
        else
            irecL = irecS+ndays-1
            PREC(irecS:irecL,:,:) = UX
            irecS = irecS+ndays
        end if
        ; print(irecL)
    end do

    do p = 0, 72
        prec(i,p,:,:) = dim_sum_n_Wrap(PREC(5*p:5*p+4,:,:), 0)/5.0
    end do
    copy_VarCoords(PREC(0,:,:), prec(0,0,:,:))
    print("end of year "+year(i))
end do
prec!0 = "year"
prec&year = year
prec!1 = "pentad"
prec&pentad = penta
printVarSummary(prec)

prec@units = "mm/day"
prec@long_name = "precipitation derived from ERA5 daily data"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "1_prec_pentad_ERA5_"+yrStrt+"_"+yrLast+".nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->prec = prec

print("end of writing data")
;================================================================
; plot
;================================================================

u = prec
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
figname = "1_calc_prec_pentad_ERA5"
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
