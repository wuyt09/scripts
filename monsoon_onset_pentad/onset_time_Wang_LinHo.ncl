
begin

yrStrt = 1979
yrLast = 2017
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)

bob_lat1 = 5
bob_lat2 = 15
bob_lon1 = 90
bob_lon2 = 100

scs_lat1 = 10
scs_lat2 = 20
scs_lon1 = 110
scs_lon2 = 120

ism_lat1 = 10
ism_lat2 = 20
ism_lon1 = 60
ism_lon2 = 85


lat = fspan(-88.75, 88.75, 72)
lon = fspan(1.25, 358.75, 144)
lat@units = "degree_north"
lon@units = "degree_east"
nlat = dimsizes(lat)
nlon = dimsizes(lon)

diri = "/home/yangsong3/data-observation/GPCP/pentad/"
prec = new((/nyear,npenta,nlat,nlon/), float)
prec@_FillValue = -999.0
prec!0 = "year"
prec&year = year
prec!1 = "penta"
prec&penta = penta
prec!2 = "lat"
prec&lat = lat
prec!3 = "lon"
prec&lon = lon
do i = 0, nyear-1
    fili = "gpcp_pen_v2.2.lnx."+year(i)
    do p = 0, npenta-1
        prec(i,p,:,:) = fbindirread(diri+fili, p, (/nlat,nlon/), "float")
    end do
end do
printVarSummary(prec)

prec_Jan = dim_avg_n_Wrap(prec(:,0:5,:,:), 1)
printVarSummary(prec_Jan)


pens = 18
precx = prec(:,pens:npenta-1,:,:)
prec_new = precx
prec_new = precx - conform_dims(dimsizes(precx), prec_Jan, (/0,2,3/))

prec_bob = wgt_areaave_Wrap(prec_new(:,:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
prec_scs = wgt_areaave_Wrap(prec_new(:,:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
prec_ism = wgt_areaave_Wrap(prec_new(:,:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)

diff = prec_new
last = prec_new
diff = 0
last = 0
diff1 = prec_bob
diff1 = 0
diff2 = diff1
diff3 = diff1
last1 = prec_bob
last1 = 0
last2 = last1
last3 = last1

do n = 1, npenta-pens-2
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            diff(:,n,ilat,ilon) = (prec_new(:,n-1,ilat,ilon)-5)*(prec_new(:,n,ilat,ilon)-5)
            last(:,n,ilat,ilon) = (prec_new(:,n+1,ilat,ilon)-5)*(prec_new(:,n,ilat,ilon)-5)

            diff1(:,n) = (prec_bob(:,n-1)-5)*(prec_bob(:,n)-5)
            diff2(:,n) = (prec_scs(:,n-1)-5)*(prec_scs(:,n)-5)
            diff3(:,n) = (prec_ism(:,n-1)-5)*(prec_ism(:,n)-5)
            
            last1(:,n) = (prec_bob(:,n+1)-5)*(prec_bob(:,n)-5)
            last2(:,n) = (prec_scs(:,n+1)-5)*(prec_scs(:,n)-5)
            last3(:,n) = (prec_ism(:,n+1)-5)*(prec_ism(:,n)-5)
        end do
    end do
end do
printVarSummary(diff)
peak = prec_new(:,0,:,:)
BOB = new(nyear, integer)
SCS = BOB
ISM = BOB

do yy = 0, nyear-1
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            a = ind(prec_new(yy,:,ilat,ilon).gt.5. .and. diff(yy,:,ilat,ilon).lt.0 .and. last(yy,:,ilat,ilon).gt.0)
            peak(yy,ilat,ilon) = (/a(0)/) + pens
            delete(a)
        end do
    end do

    a = ind(prec_bob(yy,:).ge.5. .and. (diff1(yy,:).lt.0) .and. (last1(yy,:).gt.0))
    BOB(yy) =(/a(0)/) + pens
    delete(a)

    a = ind(prec_scs(yy,:).ge.5. .and. (diff2(yy,:).lt.0) .and. (last2(yy,:).gt.0))
    SCS(yy) =(/a(0)/) + pens
    delete(a)

    a = ind(prec_ism(yy,:).ge.5. .and. (diff3(yy,:).lt.0) .and. (last3(yy,:).gt.0))
    ISM(yy) =(/a(0)/) + pens
    delete(a)
    print("end of year "+year(yy))
end do

print(year+"    "+BOB+" "+SCS+" "+ISM)

write_table("onset_time_Wang_LinHo_BOB_SCS_ISM.txt", "w", [/"YEAR","BOB","SCS","ISM"/], "%10s %10s %10s %10s")
write_table("onset_time_Wang_LinHo_BOB_SCS_ISM.txt", "a", [/year,BOB,SCS,ISM/], "%0.4i %0.2i %0.2i %0.2i")


; bob = wgt_areaave_Wrap(peak(:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
; scs = wgt_areaave_Wrap(peak(:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
; ism = wgt_areaave_Wrap(peak(:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)
; print(bob+" "+scs+" "+ism)

peak@units = "penta"
peak@long_name = "the onset pentad define by prec"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "onset_time_Wang_LinHo_"+yrStrt+"_"+yrLast+".nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->prec_peak = peak

;================================================================
; plot
;================================================================
lonstr = 50
lonend = 150
latstr = 0
latend = 30

peak_clm = dim_avg_n_Wrap(peak, 0)

figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "onset_time_Wang_LinHo"
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
; res@mpCenterLonF = 180
res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend


res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.02
res@tmYLLabelFontHeightF = 0.02

res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = 20
res@cnMaxLevelValF = 40
res@cnLevelSpacingF = 2

plot(0) = gsn_csm_contour_map(wks, peak_clm, res)
draw(plot)
end