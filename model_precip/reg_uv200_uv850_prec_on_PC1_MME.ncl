begin

lonstr = 40
lonend = 150
latstr = -15
latend = 55

yrStrt = 1951
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

LEV = (/20000,85000/)
;=======================================================
diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
lat = fobs->lat({latstr:latend})
lon = fobs->lon({lonstr:lonend})
nlat = dimsizes(lat)
nlon = dimsizes(lon)

;===================================================================
; dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"
; models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
;            "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
;            "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
models = (/"BCC-CSM2-MR","BCC-ESM1","CAMS-CSM1-0","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2",\
          "CESM2-WACCM","CIESM","E3SM-1-0","E3SM-1-1-ECA","E3SM-1-1","EC-Earth3","EC-Earth3-Veg",\ ;
          "FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G","GISS-E2-1-H","IPSL-CM6A-LR",\
          "KACE-1-0-G","MCM-UA-1-0","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","SAM0-UNICON"/)
nm = dimsizes(models)

multi_u = new((/nm,nyear,2,nlat,nlon/), float)
multi_u!0 = "model"
multi_u&model = models
multi_v = multi_u
multi_p = new((/nm,nyear,nlat,nlon/), float)

sig = 0.05
do i = 0, nm-1
    ftxt = "./EOF_PC1_"+models(i)+"_"+yrStrt+"-"+yrLast+"_std.txt"
    x = asciiread(ftxt, (/nyear/), "float")

    ;-----------------------------------------------------------------------------------
    dirigcm = "/home/yangsong3/CMIP6/"+models(i)+"/historical/"
    filiu := systemfunc("ls "+dirigcm+"ua_Amon_*"+models(i)+"*_r1i1p1f1_g*.nc")
    filiv := systemfunc("ls "+dirigcm+"va_Amon_*"+models(i)+"*_r1i1p1f1_g*.nc")
    filip := systemfunc("ls "+dirigcm+"pr_Amon_*"+models(i)+"*_r1i1p1f1_g*.nc")
    numu = dimsizes(filiu)
    print(numu)

    if (numu.lt.2) then
        fu := addfile(filiu, "r")
        fv := addfile(filiv, "r")
        fp := addfile(filip, "r")

        TIME := fu->time
        YYYY := cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        U := fu->ua(iYYYY,{LEV},{latstr:latend},{lonstr:lonend})
        V := fv->va(iYYYY,{LEV},{latstr:latend},{lonstr:lonend})
        P := fp->pr(iYYYY,{latstr:latend},{lonstr:lonend})
        P = P*86400. ; units mm/day

        u := (U(6::12,:,:,:)+U(7::12,:,:,:))/2
        v := (V(6::12,:,:,:)+V(7::12,:,:,:))/2
        p := (P(6::12,:,:)+P(7::12,:,:))/2
        copy_VarCoords(U(0,:,:,:), u(0,:,:,:))
        copy_VarCoords(V(0,:,:,:), v(0,:,:,:))
        copy_VarCoords(P(0,:,:), p(0,:,:))
        u!2 = "lat"
        u!3 = "lon"
        v!2 = "lat"
        v!3 = "lon"
        p!1 = "lat"
        p!2 = "lon"

    else
        fu := addfiles(filiu, "r")
        fv := addfiles(filiv, "r")
        fp := addfiles(filip, "r")

        TIME := fu[:]->time
        YYYY := cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        U := fu[:]->ua(iYYYY,{LEV},{latstr:latend},{lonstr:lonend})
        V := fv[:]->va(iYYYY,{LEV},{latstr:latend},{lonstr:lonend})
        P := fp[:]->pr(iYYYY,{latstr:latend},{lonstr:lonend})
        P = P*86400. ; units mm/day

        u := (U(6::12,:,:,:)+U(7::12,:,:,:))/2
        v := (V(6::12,:,:,:)+V(7::12,:,:,:))/2
        p := (P(6::12,:,:)+P(7::12,:,:))/2
        copy_VarCoords(U(0,:,:,:), u(0,:,:,:))
        copy_VarCoords(V(0,:,:,:), v(0,:,:,:))
        copy_VarCoords(P(0,:,:), p(0,:,:))
        u!2 = "lat"
        u!3 = "lon"
        v!2 = "lat"
        v!3 = "lon"
        p!1 = "lat"
        p!2 = "lon"
    end if
    ; printVarSummary(u)

    u_med = linint2_Wrap(u&lon, u&lat, u, False, lon, lat, 0)
    v_med = linint2_Wrap(v&lon, v&lat, v, False, lon, lat, 0)
    p_med = linint2_Wrap(p&lon, p&lat, p, False, lon, lat, 0)
    multi_u(i,:,:,:,:) = u_med
    multi_v(i,:,:,:,:) = v_med
    multi_p(i,:,:,:) = p_med


    ; rc_u = regCoef_n(x, u, 0, 0)
    ; rc_v = regCoef_n(x, v, 0, 0)
    ; rc_p = regCoef_n(x, p, 0, 0)
    ; copy_VarCoords(u(0,:,:,:), rc_u)
    ; copy_VarCoords(v(0,:,:,:), rc_v)
    ; copy_VarCoords(p(0,:,:,:), rc_p)
    ; ; printVarSummary(rc_v)

    ; tval_u = onedtond(rc_u@tval, dimsizes(rc_u))
    ; df_u = onedtond(rc_u@nptxy, dimsizes(rc_u)) - 2
    ; b_u = tval_u
    ; b_u = 0.5
    ; prob_u = 1 - betainc(df_u/(df_u+tval_u^2), df_u/2.0, b_u)

    ; tval_v = onedtond(rc_v@tval, dimsizes(rc_v))
    ; df_v = onedtond(rc_v@nptxy, dimsizes(rc_v)) - 2
    ; b_v = tval_v
    ; b_v = 0.5
    ; prob_v = 1 - betainc(df_v/(df_v+tval_v^2), df_v/2.0, b_v)

    ; tval_p = onedtond(rc_p@tval, dimsizes(rc_p))
    ; df_p = onedtond(rc_p@nptxy, dimsizes(rc_p)) - 2
    ; b_p = tval_p
    ; b_p = 0.5
    ; prob_p = 1 - betainc(df_p/(df_p+tval_p^2), df_p/2.0, b_p)

    ; rc_u = mask(rc_u, prob_u.le.sig .and. prob_v.le.sig, True)
    ; rc_v = mask(rc_v, prob_u.le.sig .and. prob_v.le.sig, True)
    ; rc_p = mask(rc_p, prob_p.le.sig, True)
    ; rc_u = rc_u*10
    ; rc_v = rc_v*10
    ; rc_p = rc_p*10

    ;=================================================================

    ; res@gsnLeftString = models(i)
    ; plot(i) = gsn_csm_contour_map(wks, rc_p, res)
    ; plotv(i) = gsn_csm_vector(wks, rc_u, rc_v, vecres)
    ; overlay(plot(i), plotv(i))

    ; plot_topo(i) = gsn_csm_contour(wks, topo, res_topo)
    ; overlay(plot(i), plot_topo(i))
    ; print("end of "+models(i))
    ; delete([/ftxt,x/])
    ; delete([/fu,fv,TIME,YYYY,iYYYY,U,V,u,v,rc_u,rc_v,tval_u,tval_v,df_u,df_v,b_u,b_v,prob_u,prob_v/])
    ; delete([/fp,P,p,rc_p,tval_p,df_p,b_p,prob_p/])
end do

ftxt = "./EOF_PC1_MME_"+yrStrt+"-"+yrLast+"_std.txt"
x_mme = asciiread(ftxt, (/nyear/), "float")

mme_u = dim_avg_n_Wrap(multi_u, 0)
mme_v = dim_avg_n_Wrap(multi_v, 0)
mme_p = dim_avg_n_Wrap(multi_p, 0)
rc_u_mme = regCoef_n(x_mme, mme_u, 0, 0)
rc_v_mme = regCoef_n(x_mme, mme_v, 0, 0)
rc_p_mme = regCoef_n(x_mme, mme_p, 0, 0)
copy_VarCoords(mme_u(0,:,:,:), rc_u_mme)
copy_VarCoords(mme_v(0,:,:,:), rc_v_mme)
copy_VarCoords(mme_p(0,:,:), rc_p_mme)
printVarSummary(rc_u_mme)


;===================================================================
; topo data
;===================================================================
ftopo = addfile("~/wuyt/data/topo.sfc.1x1.nc", "r")
topo = ftopo->topo({latstr:latend},{lonstr:lonend})

TOPO = linint2_Wrap(topo&lon, topo&lat, topo, False, lon, lat, 0)
rc_u_mme(1,:,:) = mask(rc_u_mme(1,:,:), TOPO.le.1500, True)
rc_v_mme(1,:,:) = mask(rc_v_mme(1,:,:), TOPO.le.1500, True)


 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = (/"grey40"/)
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/2000/)
 res_topo@gsnLeftString = ""
 res_topo@gsnRightString = ""
;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "reg_uv200_uv850_prec_on_PC1_MME"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(2, graphic)
plotv = plot
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True
res@gsnAddCyclic = False

res@mpMinLatF = latstr
res@mpMaxLatF = latend
res@mpMinLonF = lonstr
res@mpMaxLonF = lonend

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

res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

res@lbLabelBarOn = False
res@cnFillPalette = "MPL_BrBG"

;===========================================================
 vecres  = True
 vecres@gsnDraw  = False
 vecres@gsnFrame = False

 vecres@vcGlyphStyle     = "LineArrow"   ; curly vectors"Fillarrow" 
 vecres@vcRefLengthF     = 0.05            ; define length of vec ref
 vecres@vcMinDistanceF   = 0.02
 vecres@vcLineArrowThicknessF = 3.0
 vecres@tiXAxisString    = " "             ; turn off axis label
 vecres@vcRefAnnoOn      = True            ; draw the reference vector annotation.
 vecres@vcRefAnnoString2On =  False          ;display a string below or to the right
 vecres@vcRefAnnoOrientation = "horizontal"
 vecres@vcRefAnnoString2  = " "
 vecres@vcRefAnnoOrthogonalPosF = -1.0  ; move ref vector into plot
 vecres@vcFillArrowsOn       = True
 vecres@vcLineArrowHeadMinSizeF = 0.008
 vecres@vcLineArrowHeadMaxSizeF = 0.3
 vecres@gsnLeftString  = ""
 ;vecres@gsnRightString  = "g*s/Kg"
 vecres@gsnRightString  = ""
 vecres@gsnRightStringFontHeightF = 0.02
 
 ; vecres@vcLineArrowColor = "gray52"
 vecres@vcRefAnnoFontHeightF = 0.02
 ; vecres@vcRefAnnoString1  = "3 Kg/m/s"
;=======================================================================
res@gsnRightString = "MME"
res@gsnLeftString = "200 hPa"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.5
res@cnMinLevelValF = -0.5
vecres@vcRefMagnitudeF  = 0.8

plot(0) = gsn_csm_contour_map(wks, rc_p_mme, res)
plotv(0) = gsn_csm_vector(wks, rc_u_mme(0,:,:), rc_v_mme(0,:,:), vecres)
overlay(plot(0), plotv(0))

plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(0), plot_topo(0))
;-------------------------------------------
res@gsnLeftString = "850 hPa"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.1
res@cnMaxLevelValF = 0.5
res@cnMinLevelValF = -0.5
vecres@vcRefMagnitudeF  = 0.2

plot(1) = gsn_csm_contour_map(wks, rc_p_mme, res)
plotv(1) = gsn_csm_vector(wks, rc_u_mme(1,:,:), rc_v_mme(1,:,:), vecres)
overlay(plot(1), plotv(1))

plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
overlay(plot(1), plot_topo(1))


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@lbLabelFontHeightF = 0.012
gsn_panel(wks, plot, (/1,2/), pres)

end


