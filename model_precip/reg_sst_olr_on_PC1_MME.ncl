begin

; lonstr = 40
; lonend = 150
; latstr = -15
; latend = 55

yrStrt = 1951
yrLast = 2014
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

;=======================================================
diriobs = "/home/yangsong3/data-observation/GPCP/"
filiobs = "precip.mon.mean.nc"
fobs = addfile(diriobs+filiobs, "r")
lat = fobs->lat ;({latstr:latend})
lon = fobs->lon ;({lonstr:lonend})
nlat = dimsizes(lat)
nlon = dimsizes(lon)

;===================================================================
; dirigcm = "/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/"
; models := (/"CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1",\
;            "E3SM-1-0","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G-CC","GISS-E2-1-G",\
;            "HadGEM3-GC31-LL","IPSL-CM6A-LR","MPI-ESM1-2-HR","NESM3","SAM0-UNICON","UKESM1-0-LL"/)
models = (/"BCC-CSM2-MR","BCC-ESM1","CAMS-CSM1-0","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2",\
          "CESM2-WACCM","CIESM","E3SM-1-0","E3SM-1-1-ECA","E3SM-1-1","EC-Earth3","EC-Earth3-Veg",\
          "FGOALS-f3-L","FGOALS-g3","GISS-E2-1-G","GISS-E2-1-H","IPSL-CM6A-LR",\ ;"GISS-E2-1-G-CC",
          "KACE-1-0-G","MCM-UA-1-0","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","SAM0-UNICON"/)
nm = dimsizes(models)

multi_u = new((/nm,nyear,nlat,nlon/), float)
multi_u!0 = "model"
multi_u&model = models
multi_v = multi_u

sig = 0.05
do i = 0, nm-1
    ; ftxt = "./EOF_PC1_"+models(i)+"_"+yrStrt+"-"+yrLast+"_std.txt"
    ; x = asciiread(ftxt, (/nyear/), "float")

    ;-----------------------------------------------------------------------------------
    dirigcm = "/home/yangsong3/CMIP6/"+models(i)+"/historical/"
    filiu := systemfunc("ls "+dirigcm+"ts_Amon_*"+models(i)+"*_r1i1p1f1_g*.nc")
    filiv := systemfunc("ls "+dirigcm+"rlut_Amon_*"+models(i)+"*_r1i1p1f1_g*.nc")
    numu = dimsizes(filiu)
    numv = dimsizes(filiv)
    print(numu+"    "+numv)

    if (numu.lt.2) then
        fu := addfile(filiu, "r")

        TIME := fu->time
        YYYY := cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        U := fu->ts(iYYYY,:,:)

        u := (U(6::12,:,:)+U(7::12,:,:))/2
        copy_VarCoords(U(0,:,:), u(0,:,:))
        u!1 = "lat"
        u!2 = "lon"

    else
        fu := addfiles(filiu, "r")

        TIME := fu[:]->time
        YYYY := cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        U := fu[:]->ts(iYYYY,:,:)

        u := (U(6::12,:,:)+U(7::12,:,:))/2
        copy_VarCoords(U(0,:,:), u(0,:,:))
        u!1 = "lat"
        u!2 = "lon"
    end if

    if (numv.lt.2) then
        fv := addfile(filiv, "r")

        TIME := fv->time
        YYYY := cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        V := fv->rlut(iYYYY,:,:)

        v := (V(6::12,:,:)+V(7::12,:,:))/2
        copy_VarCoords(V(0,:,:), v(0,:,:))
        v!1 = "lat"
        v!2 = "lon"

    else
        fv := addfiles(filiv, "r")

        TIME := fv[:]->time
        YYYY := cd_calendar(TIME, -1)/100
        iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
        V := fv[:]->rlut(iYYYY,:,:)

        v := (V(6::12,:,:)+V(7::12,:,:))/2
        copy_VarCoords(V(0,:,:), v(0,:,:))
        v!1 = "lat"
        v!2 = "lon"
    end if
    ; printVarSummary(u)

    u_med = linint2_Wrap(u&lon, u&lat, u, False, lon, lat, 0)
    v_med = linint2_Wrap(v&lon, v&lat, v, False, lon, lat, 0)
    multi_u(i,:,:,:) = u_med
    multi_v(i,:,:,:) = v_med

    print("end of "+models(i)+" "+i)
    ; delete([/ftxt,x/])
    delete([/fu,fv,TIME,YYYY,iYYYY,U,V,u,v/])
end do

ftxt = "./EOF_PC1_MME_"+yrStrt+"-"+yrLast+"_std.txt"
x_mme = asciiread(ftxt, (/nyear/), "float")

mme_u = dim_avg_n_Wrap(multi_u, 0)
mme_v = dim_avg_n_Wrap(multi_v, 0)
rc_u_mme = regCoef_n(x_mme, mme_u, 0, 0)
rc_v_mme = regCoef_n(x_mme, mme_v, 0, 0)
copy_VarCoords(mme_u(0,:,:), rc_u_mme)
copy_VarCoords(mme_v(0,:,:), rc_v_mme)
printVarSummary(rc_u_mme)


;===================================================================
; Plot
;===================================================================
filipath = "/home/yangsong3/wuyt/sysu/figures/model_precip/"
filiname = "reg_sst_olr_on_PC1_MME"
wks = gsn_open_wks("pdf", filipath+filiname)
plot = new(2, graphic)
plotc = plot
plot_topo = plot

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True
res@gsnAddCyclic = False

res@mpMinLatF = -30
res@mpMaxLatF = 30
res@mpMinLonF = 40
res@mpMaxLonF = 300

res@mpOutlineOn = True
res@mpFillOn = True
res@mpGeophysicalLineThicknessF = 2
res@mpGeophysicalLineColor = "grey25" ;"red"
res@mpFillDrawOrder = "PostDraw"
res@mpCenterLonF = 180

res@cnFillOn = True
res@cnLinesOn = False
res@cnLineLabelsOn = False
res@cnInfoLabelOn = False

res@gsnStringFontHeightF = 0.025
res@tmXBLabelFontHeightF = 0.018
res@tmYLLabelFontHeightF = 0.018

; res@lbLabelBarOn = False
; res@cnLevelSelectionMode = "ManualLevels"
; res@cnLevelSpacingF = 1.
; res@cnMaxLevelValF = 10
; res@cnMinLevelValF = -10
res@cnFillPalette = "temp_diff_18lev"
res@pmLabelBarOrthogonalPosF = 0.2

res@gsnRightString = "MME"

res@gsnLeftString = "SST reg on PC1"
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 0.2
res@cnMinLevelValF = -0.2
plot(0) = gsn_csm_contour_map(wks, rc_u_mme, res)


res@gsnLeftString = "OLR reg on PC1"
res@cnLevelSpacingF = 0.5
res@cnMaxLevelValF = 2.
res@cnMinLevelValF = -2.
plot(1) = gsn_csm_contour_map(wks, rc_v_mme, res)

pres = True
pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; pres@lbLabelFontHeightF = 0.012

gsn_panel(wks, plot, (/2,1/), pres)


end