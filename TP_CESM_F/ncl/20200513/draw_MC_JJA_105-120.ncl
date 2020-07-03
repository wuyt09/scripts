;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    latS   =  0
    latN   =  60
    
    stdlonL = 105
    stdlonR = 120

;time
    season  = "JJA"
;************************************************
;get data
;************************************************
;V
    V_CON_file  = addfile("/home/yangsong3/data-model/Lu_MM/CTRL/B2000_f19g16_CP_CTRL.V.plvl.mon.0250-0349.nc","r")
    V_CON       = V_CON_file->V(120:599,:,:,:)
    printVarSummary(V_CON)
    lev         = V_CON_file->lev_p
    lev@units   = "hPa"
    V_CON&lev_p = lev

    V_CON_lon       = V_CON(:,:,:,{stdlonL:stdlonR})
    V_CON_lon_avg   = dim_avg_n_Wrap(V_CON_lon,3)
    V_CON_JJA       = month_to_season(V_CON_lon_avg,season)
    printVarSummary(V_CON_JJA)
    V_CON_avg = dim_avg_n_Wrap(V_CON_JJA,0)
    printVarSummary(V_CON_avg)
    V_CON_var = dim_variance_n(V_CON_JJA,0)
    printVarSummary(V_CON_var)

    V_SEN_file  = addfile("/home/yangsong3/data-model/Lu_MM/NorthTP_Alb05/NorthTP_Alb05.V.plvl.mon.0250-0299.nc","r")
    V_SEN       = V_SEN_file->V(120:599,:,:,:)
    printVarSummary(V_SEN)
    V_SEN&lev_p     = lev

    V_SEN_lon       = V_SEN(:,:,:,{stdlonL:stdlonR})
    V_SEN_lon_avg   = dim_avg_n_Wrap(V_SEN_lon,3)
    V_SEN_JJA       = month_to_season(V_SEN_lon_avg,season)
    printVarSummary(V_SEN_JJA)
    V_SEN_avg = dim_avg_n_Wrap(V_SEN_JJA,0)
    printVarSummary(V_SEN_avg)
    V_SEN_var = dim_variance_n(V_SEN_JJA,0)
    printVarSummary(V_SEN_var)

    V_lon_diff = V_SEN_avg
    V_lon_diff = V_SEN_avg - V_CON_avg
    printVarSummary(V_lon_diff)

    W_CON_file  = addfile("/home/yangsong3/data-model/Lu_MM/CTRL/B2000_f19g16_CP_CTRL.OMEGA.plvl.mon.0250-0349.nc","r")
    W_CON       = W_CON_file->OMEGA(120:599,:,:,:)
    printVarSummary(W_CON)
    W_CON       = W_CON*(-1)
    W_CON&lev_p = lev

    W_CON_lon       = W_CON(:,:,:,{stdlonL:stdlonR})
    W_CON_lon_avg   = dim_avg_n_Wrap(W_CON_lon,3)
    W_CONsu         = month_to_season(W_CON_lon_avg,season)
    printVarSummary(W_CONsu)
    W_CON_avg2 = dim_avg_n_Wrap(W_CONsu,0)
    printVarSummary(W_CON_avg2)
    W_CON_var2 = dim_variance_n(W_CONsu,0)
    printVarSummary(W_CON_var2)

    W_CON_avg22 = W_CON_avg2
    W_CON_avg22 = W_CON_avg2*(100)
    
    W_SEN_file      = addfile("/home/yangsong3/data-model/Lu_MM/NorthTP_Alb05/NorthTP_Alb05.OMEGA.mon.0250-0299.nc","r")
    W_SEN           = W_SEN_file->OMEGA(120:599,:,:,:)
    printVarSummary(W_SEN)
    W_SEN           = W_SEN*(-1)
    W_SEN&lev_p     = lev

    W_SEN_lon       = W_SEN(:,:,:,{stdlonL:stdlonR})
    W_SEN_lon_avg   = dim_avg_n_Wrap(W_SEN_lon,3)
    W_SENsu         = month_to_season(W_SEN_lon_avg,season)
    printVarSummary(W_SENsu)
    W_SEN_avg2 = dim_avg_n_Wrap(W_SENsu,0)
    printVarSummary(W_SEN_avg2)
    W_SEN_var2 = dim_variance_n(W_SENsu,0)
    printVarSummary(W_SEN_var2)

    W_lon_diff = W_SEN_avg2
    W_lon_diff = W_SEN_avg2 - W_CON_avg2
    printVarSummary(W_lon_diff)
    
    W_lon_diff1 = W_lon_diff
    W_lon_diff1 = W_lon_diff*(100)
 
;tp
    tp_name =  "TOPO.nc"
    tp_file = "/home/yangsong3/mm/ncl/scripts/study/TP/"+tp_name
    tp_in       = addfile (tp_file, "r")   ; note the "s" of addfile
    tp      = tp_in->TP
    printVarSummary(tp)
    lat         = tp_in->LATIXY(:,0)
    printVarSummary(lat)
    lon         = tp_in->LONGXY(0,:)
    printVarSummary (lon)
    tp!0        = "lsmlat"
    tp&lsmlat   = lat
    tp!1        = "lsmlon"
    tp&lsmlon   = lon   
    tp = lonFlip(tp)
    
    tp_lon      = tp(:,{stdlonL:stdlonR})
    tp_lon_mask = dim_avg_n_Wrap(tp_lon,1)
    printVarSummary(tp_lon_mask)
    elev_lon = 1013.25*(1-tp_lon_mask*0.0065/288.15)^5.25145        ;m转化为hPa
    copy_VarCoords(tp_lon_mask,elev_lon)
    printVarSummary(elev_lon)

;************************************************
;t-test
;************************************************   
    sigr = 0.05
    iflag = False
    
    ntime_CON_V1 = dimsizes(V_CON_JJA&time)
    ntime_SEN_V1 = dimsizes(V_SEN_JJA&time)
    t1_V1 = ntime_CON_V1
    t2_V1 = ntime_SEN_V1
    prob_V1 = ttest(V_CON_avg,V_CON_var,t1_V1,V_SEN_avg,V_SEN_var,t2_V1,iflag,False)
    copy_VarMeta(V_CON_avg,prob_V1)

    ntime_CON_W2 = dimsizes(W_CONsu&time)
    ntime_SEN_W2 = dimsizes(W_SENsu&time)
    t1_W2 = ntime_CON_W2
    t2_W2 = ntime_SEN_W2
    prob_W2 = ttest(W_CON_avg2,W_CON_var2,t1_W2,W_SEN_avg2,W_SEN_var2,t2_W2,iflag,False)
    copy_VarMeta(W_CON_avg2,prob_W2)

    V_lon_diff_nosig = V_lon_diff
    V_lon_diff_nosig = mask(V_lon_diff,((prob_V1 .ge. sigr) .and. (prob_W2 .ge. sigr)),True)
    W_lon_diff_nosig = W_lon_diff1
    
    V_lon_diff_sig = V_lon_diff
    V_lon_diff_sig = mask(V_lon_diff,((prob_V1 .lt. sigr) .or. (prob_W2 .lt. sigr)),True)
    W_lon_diff_sig = W_lon_diff1

    W_lon_diff_sig1 = W_lon_diff
    W_lon_diff_sig1 = mask(W_lon_diff,(prob_W2 .lt. sigr),True)
;************************************************
; plot
;************************************************
    wks = gsn_open_wks("pdf","/home/yangsong3/mm/ncl/figures/study/TP/NorthTP_Alb05/fig_MC_JJA_105-120E")

    plot = new(2,graphic)                ; create graphic array
    plot_tp = new(2,graphic)

    fig_abbr = (/"(a)","(b)"/)
    
    res                = True
    res@gsnDraw        = False
    res@gsnFrame       = False      ;don't advance frame
    res@gsnAddCyclic   = False      ;EOF data is not cyclic
  
    res@vpHeightF = 0.4        ; Changes the aspect ratio
    res@vpWidthF  = 0.8

;  res@cnMissingValFillColor = "gray30"   ; fill the missing value areas with dark gray
    res@cnFillOn       = True         ; turn on color fill
    res@cnLinesOn      = False
    res@cnLineLabelsOn = False      ;turn off cn line labels
    res@cnInfoLabelOn  = False      ;turn off contour information label
    res@cnFillDrawOrder    = "PreDraw"

    res@lbLabelBarOn               = True        ; turn off individual lb's
    res@pmLabelBarHeightF          = 0.04
    res@pmLabelBarWidthF           = 0.65
    res@lbLabelFontHeightF         = 0.015
    res@pmLabelBarOrthogonalPosF   = 0.01   ;move Label Bar down

    res@vcRefLengthF            = 0.07            ;the length of the reference magnitude
    res@vcRefAnnoFontHeightF    = 0.0095
    res@vcMinDistanceF          = 0.018            ;the minimum distance of neighboring vectors
;   res@vcGlyphStyle            = "CurlyVector"   ;the style of vector
    res@vcLineArrowHeadMaxSizeF = 0.011           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.007
    res@vcLineArrowThicknessF   = 1.0
    res@vcMapDirection          = False           ;vectors to be pointing in the right direction
    res@vcRefAnnoOn             = True           ;turn on vector reference box
    res@vcRefAnnoPerimOn        = False
    res@vcRefAnnoString2        = " "
    res@vpClipOn               = False      ;turn off vector clipping
    res@vcRefAnnoOrthogonalPosF = 0          ;move ref vector down
    res@vcRefAnnoParallelPosF   = 1.0
    
    res@tmYLMode = "Explicit" 
    res@tmYLValues = (/1000,700,500,300,200,100/)
    res@tmYLLabels = (/1000,700,500,300,200,100/)

    res@tmXBLabelFontHeightF = 0.018    ;font height of tick labels
    res@tmYLLabelFontHeightF = 0.018
    res@tmXBTickSpacingF     = 10.     ;label spacing

    res@gsnCenterString = ""
    res@gsnLeftStringFontHeightF = 0.018
    res@gsnRightStringFontHeightF = 0.018
    res@gsnCenterStringFontHeightF = 0.018

    res@cnFillPalette        = "BlueWhiteOrangeRed"
    res@vcRefMagnitudeF        = 10.0               ;the magnitude used as the reference magnitude
    res@vcRefAnnoString1       = "10"
    res@vcLineArrowColor      = "black"           ; white vectors
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res@cnLevels=(/-1.0,-0.7,-0.4,-0.1,1.0,2.0,3.0,4.0/)
    res@cnFillColors=(/236,220,190,160,120,90,58,40,2/)
    res@gsnLeftString  = fig_abbr(0)+" Omega & Vertical Cell (CTRL)"
    res@gsnRightString = stdlonL+"~S~o~N~"+"-"+stdlonR+"~S~o~N~E"
    plot(0) = gsn_csm_pres_hgt_vector(wks,W_CON_avg22({1000:100},{latS:latN}),V_CON_avg({1000:100},{latS:latN}),W_CON_avg22({1000:100},{latS:latN}),res)

    res@vcRefAnnoOn           = False 
    res@vcRefMagnitudeF       = 4.0               ;the magnitude used as the reference magnitude
    res@vcRefAnnoString1      = "4"
    res@vcLineArrowColor      = "grey66"           ; white vectors
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res@cnLevels=(/-1.0,-0.7,-0.4,-0.1,0.1,1.0,1.5,2.0/)
    res@cnFillColors=(/236,220,190,160,120,90,58,40,2/)
    res@gsnLeftString  = fig_abbr(1)+" Diffs. in Omega & Vertical Cell (NorthTP_Alb05 - CTRL)"
    plot(1) = gsn_csm_pres_hgt_vector(wks,W_lon_diff1({1000:100},{latS:latN}),V_lon_diff_nosig({1000:100},{latS:latN}),W_lon_diff_nosig({1000:100},{latS:latN}),res)
    
    resv          = True
    resv@tmYLMode = "Explicit" 
    resv@tmYLValues = (/1000,700,500,300,200,100/)
    resv@tmYLLabels = (/1000,700,500,300,200,100/)
    resv@gsnDraw                 = False
    resv@gsnFrame                = False      ;don't advance frame
    resv@vcRefLengthF            = 0.07            ;the length of the reference magnitude
    resv@vcRefAnnoFontHeightF    = 0.0095
    resv@vcMinDistanceF          = 0.015            ;the minimum distance of neighboring vectors
;   res@vcGlyphStyle            = "CurlyVector"   ;the style of vector
    resv@vcLineArrowHeadMaxSizeF = 0.011           ;change the vector line arrows
    resv@vcLineArrowHeadMinSizeF = 0.007
    resv@vcLineArrowThicknessF   = 1.0
    resv@vcMapDirection          = False           ;vectors to be pointing in the right direction
    resv@vcRefAnnoOn             = True           ;turn on vector reference box
    resv@vcRefAnnoPerimOn        = False
    resv@vcRefAnnoString2        = " "
    resv@vpClipOn                = False      ;turn off vector clipping
;    resv@vcVectorDrawOrder       = "PostDraw"      ;draw vectors last
    resv@vcRefAnnoOrthogonalPosF = -0.1          ;move ref vector down
    resv@vcRefAnnoParallelPosF   = 1.0

    resv@vcRefMagnitudeF       = 4.0               ;the magnitude used as the reference magnitude
    resv@vcRefAnnoString1      = "4"
    res@vcLineArrowColor       = "black"           ; white vectors
    resv@gsnLeftString         = " "
    resv@gsnRightString        = " "

    plot_vect = gsn_csm_vector(wks,V_lon_diff_sig({1000:100},{latS:latN}),W_lon_diff_sig({1000:100},{latS:latN}),resv)
    overlay(plot(1),plot_vect)
    
    res_sig                = True       ; use plot options 
    res_sig@cnFillOn       = True       ; Fill contours 
    res_sig@cnLinesOn      = False;True      ; turn off contour lines 
    res_sig@gsnFrame       = False 
    res_sig@gsnDraw        = False 
    res_sig@gsnAddCyclic   = False  
    res_sig@cnLineLabelsOn   = False      ; turn off contour line labels 
    res_sig@cnLineThicknessF = 2 
    res_sig@cnLineColor      = "black" 
    res_sig@cnInfoLabelOn    = False 
    res_sig@cnMonoFillPattern    = True;False 
    res_sig@cnFillPattern = 17 
    res_sig@cnMonoFillScale = True 
    ;res_sig@cnFillScaleF = 2. 
    res_sig@cnMonoFillColor = True 
    ;res_sig@cnFillColor  =  "red" 
    ;res_sig@cnFillDotSizeF = 0.003 
    res_sig@lbLabelBarOn = False;True   ; turn off the label bar  
    res_sig@gsnLeftString      = " " 
    res_sig@gsnRightString      = " "
    
    plot_sig = gsn_csm_pres_hgt(wks,W_lon_diff_sig1({1000:100},{latS:latN}),res_sig)
    overlay(plot(1),plot_sig)
    
    getvalues plot_sig@contour                                                     
    "pmAnnoManagers" : am_ids                                                
    end getvalues                                                              
  
    index = ind(NhlName(am_ids).eq."right_axis")  
    if(.not.ismissing(index)) then                                             
        NhlRemoveAnnotation(plot_sig@contour,am_ids(index))                          
    end if

    res_tp=True
    res_tp@gsnDraw         = False             ; don't draw the plots yet
    res_tp@gsnFrame        = False        ; don't advance the frame  
    res_tp@gsnMaximize     = True
    res_tp@trYMaxF         = 1000.0
    res_tp@gsnYRefLine     = 1000.0          ; create a reference line 
    res_tp@gsnBelowYRefLineColor = "black"              ; above ref line fill grey
;    res_tp@xyCurveDrawOrder       = "PostDraw"      ;draw vectors last

    plot_tp(0) =  gsn_csm_xy(wks,elev_lon&lsmlat,elev_lon,res_tp)
    overlay(plot(0),plot_tp(0))
    plot_tp(1) =  gsn_csm_xy(wks,elev_lon&lsmlat,elev_lon,res_tp)
    overlay(plot(1),plot_tp(1))
;************************************************
; panel plot only resources
;************************************************
    resP                           = True         ; modify the panel plot
    resP@gsnPanelLabelBar          = False        ; add common colorbar
    gsn_panel(wks,plot,(/2,1/),resP)     ; now draw as one plot
    frame(wks)

end