
begin

casename = (/"Ctrl_FAMIP","Hist_FAMIP"/)
varname = (/"FSDTOA","QREFHT","TREFHT","FSDS","PS"/) ;"LHFLX","SHFLX",
nvar = dimsizes(varname)

yrStrt = 198102
yrLast = 200601
Season = (/"DJF"/)

outname = (/"solarin","huss","t2m","ssrd","sp"/) ;"slhf","sshf",
outpath = "../data/"


do i = 0, nvar-1
    Title  = (/outname(i)+"_base",outname(i)+"_warm","diff_"+outname(i)/)
    basename = outname(i)+"_base.dat"
    warmname = outname(i)+"_warm.dat"

    diri1 = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(0)+"/a20191206/"
    fili1 = casename(0)+".cam.h0."+varname(i)+".1979-2005.nc"

    diri2 = "/home/yangsong3/data-model/wuyt/TPbias/TPbias_"+casename(1)+"/a20191206/"
    fili2 = casename(1)+".cam.h0."+varname(i)+".1979-2005.nc"

    f1 = addfile(diri1+fili1, "r")
    f2 = addfile(diri2+fili2, "r")

    TIME = f1->time
    YYYY = cd_calendar(TIME, -1)
    iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

    X_base = f1->$varname(i)$(iYYYY,:,:)
    X_warm = f2->$varname(i)$(iYYYY,:,:)
    
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

end do

end

