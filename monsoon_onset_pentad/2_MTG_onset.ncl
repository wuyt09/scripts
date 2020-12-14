
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

f = addfile("./1_MTG_pentad_ERA5_1979_2017.nc", "r")
MTG = f->MTG
printVarSummary(MTG)

diff = MTG
last = MTG
last2 = MTG
do n = 1, npenta-3
    diff(:,n,:,:) = MTG(:,n-1,:,:)*MTG(:,n,:,:)
    last(:,n,:,:) = MTG(:,n+1,:,:)*MTG(:,n,:,:)
    last2(:,n,:,:) = MTG(:,n+1,:,:)*MTG(:,n+2,:,:)
end do
printVarSummary(diff)

peak = MTG(:,0,:,:)
do yy = 0, nyear-1
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            a = ind(MTG(yy,:,ilat,ilon).gt.0 .and. diff(yy,:,ilat,ilon).lt.0 .and. last(yy,:,ilat,ilon).gt.0 .and. last2(yy,:,ilat,ilon).gt.0)
            peak(yy,ilat,ilon) = (/a(0)/)
            delete(a)
        end do
    end do
end do
printVarSummary(peak)

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

peak_clm = dim_avg_n_Wrap(peak, 0)
print(peak_clm({bob_lat1:bob_lat2},{bob_lon1:bob_lon2}))
print(peak_clm({scs_lat1:scs_lat2},{scs_lon1:scs_lon2}))
print(peak_clm({ism_lat1:ism_lat2},{ism_lon1:ism_lon2}))
printVarSummary(peak_clm)


peak@units = "penta"
peak@long_name = "the onset pentad define by MTG reverse"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "2_MTG_onset_pentad.nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->MTG_peak = peak

;================================================================
; plot
;================================================================

lonstr = 50
lonend = 150
latstr = 0
latend = 30


figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "2_MTG_onset"
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