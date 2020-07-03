
begin

nlat = 96
nlon = 144

casename = (/"CTRL","NOQ1","NOTP"/)
ncase = dimsizes(casename)
clim_pentad_u = new((/ncase,73,96,144/), float)
LEV = 850

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
 
    tname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.U.*.daily.nc")
    print(tname)

    ft = addfiles(tname, "r")

    T = ft[:]->U(:,{LEV},:,:)
    t = reshape(T, (/25,365,nlat,nlon/))
    copy_VarCoords(T(0:364,:,:), t(0,:,:,:))

    pentad_t = t(:,0:72,:,:)

    do p = 0,72  ; a year has 73 penta
      pentad_t(:,p,:,:) = dim_sum_n_Wrap(t(:,5*p:5*p+4,:,:), 1)/5.0
    end do

    clim_pentad_T = dim_avg_n_Wrap(pentad_t, 0)
    clim_pentad_u(i,:,:,:) = clim_pentad_T ;(18:35,:,:,:)
    delete([/T,t,pentad_t,clim_pentad_T/])

end do
printVarSummary(clim_pentad_u)
; ;-------------------------------------------------------------------
; ; ISM onset: U850(5-15N, 40-80E) - U850(20-30N, 70-90E) (Wang et al. 2001)

; u1 = dim_avg_n_Wrap(clim_pentad_u(:,:,{5:15},{40:80}), (/2,3/))
; u2 = dim_avg_n_Wrap(clim_pentad_u(:,:,{20:30},{70:90}), (/2,3/))
; ism_index = u1
; ism_index = u1-u2
; printVarSummary(ism_index)
; ;-------------------------------------------------------------------
;-------------------------------------------------------------------
; ISM onset: U850(5-15N, 40-80E) (Wang et al. 2009)

ism_index = dim_avg_n_Wrap(clim_pentad_u(:,:,{5:15},{40:80}), (/2,3/))
printVarSummary(ism_index)
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "ISM_circulation_index_evolution"

penta = ispan(1, 73, 1)
Colors = (/"red","blue","green4"/) ;,"black","purple"
Patterns = (/0,0,0/)
Texts = (/"CTRL","NOHEAT","NOTP"/) ;"CTRL-NOHEAT","NTNH","NT"
pens = 20
pene = pens+45 ; total 46 pentad
print(penta(pens-1:pene-1))

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 1
plot = new(nplot, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpHeightF= 0.4                    ; change aspect ratio of plot
res@vpWidthF = 0.7     
res@gsnStringFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.025
res@tmYLLabelFontHeightF = 0.025

; res@trXMinF = pens
; res@trXMaxF = pene
; res@trYMinF = -1
; res@trYMaxF = 2
res@xyDashPatterns = Patterns
res@xyLineThicknessF = 7
res@xyLineColors = Colors

res@tmXBMode = "Explicit"
;res@tmXBValues = penta(pens-1:pene-1)
res@tmXBLabels = penta(pens:pene)
res@tmXBMinorOn = True
res@tmXTOn = False

res@tmXBPrecision = 2
; res@tmYLPrecision = 1
res@tiYAxisString = "m/s"
res@tiXAxisString = "Penta"

res@gsnYRefLine = 0
res@gsnYRefLineColor = "darkgoldenrod1"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2
res@gsnYRefLine = 0

res@gsnLeftString = "ISM monsoon index" ;: U850 averaged over 5-15N, 110-120E
plot(0) = gsn_csm_xy(wks, penta(pens-1:pene-1), ism_index(:,pens-1:pene-1), res)


;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.025                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

; ;============================================================
    interval = 2
    gx = (/39,42/)
    text_pos = 42.5
    gy1 = (/-3,-3/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    ; ; gy4 = gy3-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    res_lines@gsLineColor = Colors(2)
    res_lines@gsLineDashPattern = Patterns(2)
    dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

draw(plot(0))
end