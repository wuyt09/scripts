
begin

filiname = (/"cloud_cover", "cloud_liquid_water_content","cloud_ice_water_content"/)
varname = (/"cc","clwc","ciwc"/)
nvar = dimsizes(varname)

yrStrt = 1979
yrLast = 1981
year = ispan(yrStrt, yrLast, 1)
season = (/"JJA"/)

outname = (/"cc_clear","clwc_clear","ciwc_clear"/)
outpath = "./"

diri = "/home/yangsong3/data-observation/ERA5-monthly/pressure/"

base_s = ispan(1979, 2009, 1)
base_e = ispan(1988, 2018, 1)
warm_s = ispan(1980, 2010, 1)
warm_e = ispan(1989, 2019, 1)
n = dimsizes(base_s)
print("n = "+n)

do i = 0, nvar-1
    Title  = (/outname(i)+"_base",outname(i)+"_warm","diff_"+outname(i)/)
    basename = outname(i)+"_base.dat"
    warmname = outname(i)+"_warm.dat"
    
    fili = filiname(i)+"_mon_1x1_1979_2019.nc"

    f1 = addfile(diri+fili, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)/100
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    X = short2flt(f1->$varname(i)$(iYYYY,:,:,:))
    x = month_to_season(X, season)
    x!0 = "year"
    x&year = year

    nlat = dimsizes(x&latitude)
    nlon = dimsizes(x&longitude)
    nlev = dimsizes(x&level)

    X_b = new((/n,2,nlev,nlat,nlon/), float)
    X_b = 0
    X_w = X_b

    ; do ll = 0, n-1 
    ;     X_b(ll,:,:,:,:) = x({year|1979:1980},level|:,latitude|:,longitude|:)
    ;     X_w(ll,:,:,:,:) = x({year|1979:1980},level|:,latitude|:,longitude|:)
    ; end do

    xb = dim_avg_n_Wrap(X_b, 1)
    xw = dim_avg_n_Wrap(X_w, 1)
    copy_VarCoords(X(0:n-1,:,:,:), xb)
    copy_VarCoords(X(0:n-1,:,:,:), xw)
    printVarSummary(xb)

    system("rm -f "+outpath+basename)
    fbindirwrite(outpath+basename,xb)
    system("rm -f "+outpath+warmname)
    fbindirwrite(outpath+warmname,xw)
;===================================================================
x_diff = xw
x_diff = xw-xb
lonstr = 75
lonend = 105
latstr = 25
latend = 40
LEV = 500

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
 plot(0) = gsn_csm_contour_map(wks, xb(20,{LEV},:,:), res)
 
 res@tiMainString = Title(1)
 plot(1) = gsn_csm_contour_map(wks, xw(20,{LEV},:,:), res)
 
 res@tiMainString = Title(2)
 plot(2) = gsn_csm_contour_map(wks, x_diff(20,{LEV},:,:), res)

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

