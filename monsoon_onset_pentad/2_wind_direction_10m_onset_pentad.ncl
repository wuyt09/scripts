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


pens = 18
f = addfile("1_wind_direction_10m_pentad_ERA5_"+yrStrt+"_"+yrLast+".nc", "r")
ux = f->theta
u = ux(:,pens:npenta-1,:,:)
printVarSummary(u)
u_Jan = dim_avg_n_Wrap(ux(:,0:5,:,:), 1)
printVarSummary(u_Jan)
delete(ux)

d_theta = u
barr = 80

do yy = 0, nyear-1
    do pp = 0, npenta-pens-1
        do ilat = 0, nlat-1
            do ilon = 0, nlon-1
                delta = u(yy,pp,ilat,ilon) - u_Jan(yy,ilat,ilon)
                if (delta.gt.180) then
                    d_theta(yy,pp,ilat,ilon) = delta - 360
                else if (delta.lt.-180) then
                    d_theta(yy,pp,ilat,ilon) = abs(delta + 360)
                else
                    d_theta(yy,pp,ilat,ilon) = delta
                end if
                end if
            end do
        end do
    end do
    print("end of year "+year(yy))
end do

u_bob = wgt_areaave_Wrap(d_theta(:,:,{bob_lat1:bob_lat2},{bob_lon1:bob_lon2}), 1., 1., 0)
u_scs = wgt_areaave_Wrap(d_theta(:,:,{scs_lat1:scs_lat2},{scs_lon1:scs_lon2}), 1., 1., 0)
u_ism = wgt_areaave_Wrap(d_theta(:,:,{ism_lat1:ism_lat2},{ism_lon1:ism_lon2}), 1., 1., 0)

diff = u
last = u

diff1 = u_bob
diff1 = 0
diff2 = diff1
diff3 = diff1
last1 = u_bob
last1 = 0
last2 = last1
last3 = last1


do n = 1, npenta-pens-2
    diff1(:,n) = (u_bob(:,n-1)-barr)*(u_bob(:,n)-barr)
    diff2(:,n) = (u_scs(:,n-1)-barr)*(u_scs(:,n)-barr)
    diff3(:,n) = (u_ism(:,n-1)-barr)*(u_ism(:,n)-barr)

    last1(:,n) = (u_bob(:,n+1)-barr)*(u_bob(:,n)-barr)
    last2(:,n) = (u_scs(:,n+1)-barr)*(u_scs(:,n)-barr)
    last3(:,n) = (u_ism(:,n+1)-barr)*(u_ism(:,n)-barr)
end do
printVarSummary(diff)

peak = u(:,0,:,:)
BOB = new(nyear, integer)
SCS = BOB
ISM = BOB

do yy = 0, nyear-1
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            a = ind((d_theta(yy,:,ilat,ilon).gt.barr))
            peak(yy,ilat,ilon) = (/a(0)/) + pens +1
            delete(a)
        end do
    end do

    a= ind((diff1(yy,:).lt.0))
    ; a= ind((u_bob(yy,:).gt.barr))
    BOB(yy) =(/a(0)/) + pens +1
    delete(a)

    a= ind((diff2(yy,:).lt.0))
    ; a= ind((u_scs(yy,:).gt.barr))
    SCS(yy) =(/a(0)/) + pens +1
    delete(a)

    a= ind((diff3(yy,:).lt.0))
    ; a= ind((u_ism(yy,:).gt.barr))
    ISM(yy) =(/a(0)/) + pens +1
    delete(a)
end do
printVarSummary(peak)
delete([/u,u_Jan,delta,d_theta,diff,last,diff1,diff2,diff3,last1,last2,last3,u_bob,u_scs,u_ism/])

print(BOB+" "+SCS+" "+ISM)

write_table("wind_direction_10m_onset_BOB_SCS_ISM.txt", "w", [/"YEAR","BOB","SCS","ISM"/], "%10s %10s %10s %10s")
write_table("wind_direction_10m_onset_BOB_SCS_ISM.txt", "a", [/year,BOB,SCS,ISM/], "%0.4i %0.2i %0.2i %0.2i")


peak@units = "penta"
peak@long_name = "the onset pentad define by 10m wind direction"

outpath = "/home/yangsong3/wuyt/sysu/scripts/monsoon_onset_pentad/"
outname = "2_wind_direction_10m_onset_pentad_"+yrStrt+"_"+yrLast+".nc"

system("rm -rf "+outpath+outname)

fout = addfile(outpath+outname, "c")
setfileoption(fout, "DefineMode", True)

fout->wind_peak = peak

;================================================================
; plot
;================================================================

lonstr = 50
lonend = 150
latstr = 0
latend = 30

peak_clm = dim_avg_n_Wrap(peak, 0)

figpath = "/home/yangsong3/wuyt/sysu/figures/monsoon_onset_pentad/"
figname = "2_wind_direction_10m_onset"
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