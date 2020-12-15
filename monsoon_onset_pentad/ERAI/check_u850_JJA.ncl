
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

nlat = 121
nlon = 240
nlev = 37

mdiri = "/home/yangsong3/data-observation/ERA-interim/pressure/uwnd/"
LEV = 850

u = new((/nyear,npenta,nlat,nlon/), float)

do i = 0, 2 ;nyear-1
    filiu = "uwnd.interim."+year(i)+".nc"
    fu = addfile(mdiri+filiu, "r")    
    UX = short2flt(fu->uwnd(:,{LEV},::-1,:))

    if (dimsizes(UX&time).eq.366) then
        U = UX(0:364,:,:)
        U(0:58,:,:) = UX(0:58,:,:)
        U(59:364,:,:) = UX(60:365,:,:)
    else
        U = UX
    end if
    ; printVarSummary(U)
    delete([/UX/])

    do p = 0, 72
        u(i,p,:,:) = dim_sum_n_Wrap(U(5*p:5*p+4,:,:), 0)/5.0
    end do
    copy_VarCoords(U(0,:,:), u(0,0,:,:))
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

figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/ERAI/"
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