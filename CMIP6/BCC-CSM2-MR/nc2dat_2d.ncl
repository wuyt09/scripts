
begin

case = "BCC-CSM2-MR"

yrStrt = 1979
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
season = (/"JJA"/)

varname = (/"tas","hfls","hfss","ps","rsds","rsus","rlds","rlus","rsdt","huss"/)
nvar = dimsizes(varname)

outname = (/"t2m","slhf","sshf","sp","ssrd","ssru","strd","stru","solarin","huss"/) ;
outpath = "./hist_79-14/"

diri = "/home/yangsong3/CMIP6/"+case+"/historical/"

base_s = ispan(1979, 2004, 1)
base_e = ispan(1988, 2013, 1)
warm_s = ispan(1980, 2005, 1)
warm_e = ispan(1989, 2014, 1)
n = dimsizes(base_s)        
print("n = "+n)

do i = 0, nvar-1
    Title  = (/outname(i)+"_base",outname(i)+"_warm","diff_"+outname(i)/)
    basename = outname(i)+"_base.dat"
    warmname = outname(i)+"_warm.dat"
    
    fili = systemfunc("ls "+diri+varname(i)+"_Amon_"+case+"_historical_r1i1p1f1_g*.nc")
    nfili = dimsizes(fili)

    if (nfili.ge.2) then
        f1 = addfiles(fili, "r")
        TIME = f1[:]->time
        YYYY = cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        X = short2flt(f1[:]->$varname(i)$(iYYYY,:,:))
    else
        f1 = addfile(fili, "r")
        TIME = f1->time
        YYYY = cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        X = short2flt(f1->$varname(i)$(iYYYY,:,:))
    end if
    printVarSummary(X)

    x = month_to_season(X, season)
    printVarSummary(x)
    x!0 = "year"
    x&year = year

    nlat = dimsizes(x&lat)
    nlon = dimsizes(x&lon)

    X_b = new((/n,10,nlat,nlon/), float)
    X_w = X_b

    do ll = 0, n-1 
        X_b(ll,:,:,:) = x({year|base_s(ll):base_e(ll)},lat|:,lon|:)
        X_w(ll,:,:,:) = x({year|warm_s(ll):warm_e(ll)},lat|:,lon|:)
    end do

    xb = dim_avg_n_Wrap(X_b, 1)
    xw = dim_avg_n_Wrap(X_w, 1)
    copy_VarCoords(x(0:n-1,:,:), xb)
    copy_VarCoords(x(0:n-1,:,:), xw)
    printVarSummary(xb)

    system("rm -f "+outpath+basename)
    fbindirwrite(outpath+basename,xb)
    system("rm -f "+outpath+warmname)
    fbindirwrite(outpath+warmname,xw)


;===================================================================
x_diff = xw(20,:,:)
x_diff = xw(20,:,:)-xb(20,:,:)
lonstr = 75
lonend = 105
latstr = 25
latend = 40

wks = gsn_open_wks("pdf",outname(i))
gsn_define_colormap(wks,"BlueWhiteOrangeRed")
plot = new(3,graphic)
plots = plot 

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
 res@mpMinLatF            = latstr      ; range to zoom in on
 res@mpMaxLatF            = latend
 res@mpMinLonF            = lonstr
 res@mpMaxLonF            = lonend

 res@gsnLeftString   = season
 res@tiMainString = Title(0)
 plot(0) = gsn_csm_contour_map(wks, xb(20,:,:), res)
 
 res@tiMainString = Title(1)
 plot(1) = gsn_csm_contour_map(wks, xw(20,:,:), res)
 
 res@tiMainString = Title(2)
 plot(2) = gsn_csm_contour_map(wks, x_diff, res)

;--------------------------------
; creat panel
;--------------------------------

 resP = True
 resP@gsnFrame = False
 resP@gsnMaximize = True
; resP@gsnPanelLabelBar = True

 gsn_panel(wks,plot,(/3,1/),resP)

end do


end