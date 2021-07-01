
begin

case = "BCC-CSM2-MR"

yrStrt = 1979
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
season = (/"JJA"/)

outpath = "./hist_79-14/"
outname = "P_3D.dat"

diri = "/home/yangsong3/CMIP6/"+case+"/historical/"

filip = systemfunc("ls "+diri+"ps_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
print(filip)
nfili = dimsizes(filip)
if (nfili.ge.2) then
    f1 = addfiles(filip, "r")
    TIME = f1[:]->time
    YYYY = cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    ps = f1[:]->ps(iYYYY,:,:)
else
    f1 = addfile(filip, "r")
    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
    ps = f1->ps(iYYYY,:,:)
end if
printVarSummary(ps)
delete([/f1,TIME,YYYY,iYYYY/])

filic = systemfunc("ls "+diri+"cl_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
print(filic)
nfilic = dimsizes(filic)
if (nfilic.ge.2) then
    f1 = addfiles(filic, "r")
    hyam = f1[0]->a
    hybm = f1[0]->b
    level = f1[0]->lev
else
    f1 = addfile(filic, "r")
    hyam = f1->a
    hybm = f1->b
    level = f1->lev
end if
printVarSummary(hyam)

ntime = dimsizes(ps&time)
nlev = dimsizes(level)
nlat = dimsizes(ps&lat)
nlon = dimsizes(ps&lon)

P0 = 1.e+5
P0@units = "Pa"

P = new((/nlev,nlat,nlon/), float) ;, "1e+36"
hyao = tofloat(hyam)
hybo = tofloat(hybm)

ps_season = month_to_season(ps, season)
ps_clm = dim_avg_n_Wrap(ps_season, 0)

do ilev = 0, nlev-1
    do ilat = 0, nlat-1
        do ilon = 0, nlon-1
            print(ilev+" "+ilat+" "+ilon)
            P(ilev,ilat,ilon) = hyao(ilev)*P0+hybo(ilev)*ps_clm(ilat,ilon)
            print(hyao(ilev)+" "+hybo(ilev)+" "+ps_clm(ilat,ilon)+" "+P(ilev,ilat,ilon))            
        end do
    end do
end do
copy_VarCoords(ps_clm, P(0,:,:))
copy_VarCoords(hyam, P(:,0,0))
printVarSummary(P)

; P_season = month_to_season(P, season)
; P_clm = dim_avg_n_Wrap(P_season, 0)
; printVarSummary(P_clm)

system("rm -f "+outpath+outname)
fbindirwrite(outpath+outname,P)

;============================================================
lonstr = 75
lonend = 105
latstr = 25
latend = 40

wks = gsn_open_wks("pdf","write_P_3D")
gsn_define_colormap(wks,"BlueWhiteOrangeRed")
plot = new(1,graphic)

 res = True
 res@gsnDraw          = False
 res@gsnFrame         = False
 res@gsnSpreadColors  = True
 
 res@cnFillOn         = True
 res@cnLineLabelsOn   = False
 res@cnInfoLabelOn    = False
 
 res@lbLabelAutoStride = True


 res@mpFillOn  = True
 res@mpOutlineOn = True
 ; res@mpMinLatF            = latstr      ; range to zoom in on
 ; res@mpMaxLatF            = latend
 ; res@mpMinLonF            = lonstr
 ; res@mpMaxLonF            = lonend

 res@gsnLeftString   = season
 res@tiMainString = "Ps "
 plot(0) = gsn_csm_contour_map(wks, P(0,:,:), res)
 
 draw(plot)
; ;--------------------------------
; ; creat panel
; ;--------------------------------

;  resP = True
;  resP@gsnFrame = False
;  resP@gsnMaximize = True
; ; resP@gsnPanelLabelBar = True

;  gsn_panel(wks,plot,(/3,1/),resP)

end