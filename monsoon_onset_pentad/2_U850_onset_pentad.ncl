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

nlat = 181
nlon = 360
nlev = 19

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


pens = 18
f = addfile("1_U850_pentad_ERA5_"+yrStrt+"_"+yrLast+".nc", "r")
ux = f->u850
u = ux(:,pens:npenta-1,:,:)
printVarSummary(u)

u_bob = wgt_areaave_Wrap(u(:,:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
u_scs = wgt_areaave_Wrap(u(:,:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
u_ism = wgt_areaave_Wrap(u(:,:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)

diff = u
last = u
diff1 = u_bob
diff2 = diff1
diff3 = diff1
last1 = u_bob
last2 = last1
last3 = last1

do n = 1, npenta-pens-2
    diff(:,n,:,:) = u(:,n-1,:,:)*u(:,n,:,:)
    diff1(:,n) = u_bob(:,n-1)*u_bob(:,n)
    diff2(:,n) = u_scs(:,n-1)*u_scs(:,n)
    diff3(:,n) = u_ism(:,n-1)*u_ism(:,n)

    last(:,n,:,:) = u(:,n+1,:,:)*u(:,n,:,:)
    last1(:,n) = u_bob(:,n+1)*u_bob(:,n)
    last2(:,n) = u_scs(:,n+1)*u_scs(:,n)
    last3(:,n) = u_ism(:,n+1)*u_ism(:,n)
end do
printVarSummary(diff)

peak = u(:,0,:,:)
BOB = new(nyear, integer)
SCS = BOB
ISM = BOB

do yy = 0, nyear-1
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            a = ind((u(yy,:,ilat,ilon).gt.0) .and. (diff(yy,:,ilat,ilon).lt.0) .and. (last(yy,:,ilat,ilon).gt.0))
            peak(yy,ilat,ilon) = (/a(0)/) + pens +1
            delete(a)
        end do
    end do

    a= ind((diff1(yy,:).lt.0) .and. (u_bob(yy,:).gt.0) .and. (last1(yy,:).gt.0))
    BOB(yy) =(/a(0)/) + pens +1
    delete(a)

    a= ind((diff2(yy,:).lt.0) .and. (u_scs(yy,:).gt.0) .and. (last2(yy,:).gt.0))
    SCS(yy) =(/a(0)/) + pens +1
    delete(a)

    a= ind((diff3(yy,:).lt.0) .and. (u_ism(yy,:).gt.0) .and. (last3(yy,:).gt.0))
    ISM(yy) =(/a(0)/) + pens +1
    delete(a)
end do
printVarSummary(peak)

print(BOB+" "+SCS+" "+ISM)

write_table("U850_onset_BOB_SCS_ISM.txt", "w", [/"YEAR","BOB","SCS","ISM"/], "%10s %10s %10s %10s")
write_table("U850_onset_BOB_SCS_ISM.txt", "a", [/year,BOB,SCS,ISM/], "%0.4i %0.2i %0.2i %0.2i")


peak@units = "penta"
peak@long_name = "the onset pentad define by U850 reverse"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "2_U850_onset_pentad_"+yrStrt+"_"+yrLast+".nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->u850_peak = peak

;================================================================
; plot
;================================================================

lonstr = 50
lonend = 150
latstr = 0
latend = 30

peak_clm = dim_avg_n_Wrap(peak, 0)

figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "2_U850_onset"
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

;res@mpCenterLonF = 180
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