
begin

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
varname = (/"FSDS","FSNS"/)
nvar = dimsizes(varname)

yrStrt = 198102
yrLast = 200601
Season = (/"DJF"/)

outname = (/"ssru"/)
outpath = "../data/"

i = 0

    Title  = (/outname(i)+"_base",outname(i)+"_warm","diff_"+outname(i)/)
    basename = outname(i)+"_base.dat"
    warmname = outname(i)+"_warm.dat"

    diri1 = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(0)+"/a20191206/"
    fili1 = casename(0)+".cam.h0.FSDS.1979-2005.nc"
    fili2 = casename(0)+".cam.h0.FSNS.1979-2005.nc"

    diri2 = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(1)+"/a20191206/"
    fili3 = casename(1)+".cam.h0.FSDS.1979-2005.nc"
    fili4 = casename(1)+".cam.h0.FSNS.1979-2005.nc"

    f1 = addfile(diri1+fili1, "r")
    f2 = addfile(diri1+fili2, "r")

    f3 = addfile(diri2+fili3, "r")
    f4 = addfile(diri2+fili4, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    X1 = f1->FSDS(iYYYY,:,:)
    X2 = f2->FSNS(iYYYY,:,:)
    X3 = f3->FSDS(iYYYY,:,:)
    X4 = f4->FSNS(iYYYY,:,:)

    X_base = X1
    X_base = X1 - X2  ; ssrd - ssr
    X_warm = X3
    X_warm = X3 - X4
    
    x_base = month_to_season(X_base, Season)
    x_warm = month_to_season(X_warm, Season)

    xb = dim_avg_n_Wrap(x_base, 0)
    xw = dim_avg_n_Wrap(x_warm, 0)

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

 res@gsnLeftString   = Season
 res@tiMainString = Title(0)
 plot(0) = gsn_csm_contour_map(wks, xb, res)
 
 res@tiMainString = Title(1)
 plot(1) = gsn_csm_contour_map(wks, xw, res)
 
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


end

