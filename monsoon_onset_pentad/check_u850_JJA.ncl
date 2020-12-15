
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

mdiri = "/home/yangsong3/data-observation/ERA5-daily/Multi-levels/uwind/"
LEV = 850

U = new((/365,nlat,nlon/), float)
u = new((/nyear,npenta,nlat,nlon/), float)



do i = 0, 2 ;nyear-1
    irecS = 0
    do j = 0, nmonth-1
        mth = sprinti("%0.2i", month(j))
        filiu = "uwind."+year(i)+"-"+mth+".daily.nc"

        fu = addfile(mdiri+filiu, "r")    
        UX := short2flt(fu->u(:,{LEV},::-1,:))

        ndays := dimsizes(UX&time)
        if (ndays.eq.29) then
            irecL = irecS+ndays-2
            U(irecS:irecL,:,:) = UX(0:27,:,:)
            irecS = irecS+ndays-1
        else
            irecL = irecS+ndays-1
            U(irecS:irecL,:,:) = UX
            irecS = irecS+ndays
        end if
        ; print(irecL)
    end do

    do p = 0, 72
        u(i,p,:,:) = dim_sum_n_Wrap(U(5*p:5*p+4,:,:), 0)/5.0
    end do
    copy_VarCoords(U(0,:,:), u(0,0,:,:))
    print("end of year "+year(i))
end do
u!0 = "year"
u&year = year
u!1 = "pentad"
u&pentad = penta
printVarSummary(u)

u_jja = dim_avg_n_Wrap(u(:,30:47,:,:), (/0,1/))


;================================================================
; plot
;================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "check_u850_JJA"
wks = gsn_open_wks("pdf", figpath+figname)
plot = new(1, graphic)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red" 
; res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 1
res@cnMinLevelValF = -15
res@cnMaxLevelValF = 15
plot(0) = gsn_csm_contour_map(wks, u_jja, res)
draw(plot)

end