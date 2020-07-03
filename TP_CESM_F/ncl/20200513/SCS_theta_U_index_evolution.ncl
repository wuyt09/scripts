
begin

casename = (/"CTRL","NOQ1","NOTP"/) ;
ncase = dimsizes(casename)

nlat = 96
nlon = 144
clim_pentad_u = new((/ncase,73,nlat,nlon/), float)
clim_pentad_t = clim_pentad_u
LEV = 850

do i = 0, ncase-1
    diri = "/home/yangsong3/data-model/wuyt/TP_"+casename(i)+"_F2000_F19G16_CAM4/pro/daily/"
    uname = systemfunc("ls "+diri+"F2000_"+casename(i)+".plevel.U.*.daily.nc")
    print(uname)

    fu = addfiles(uname, "r")

    U = fu[:]->U(:,{LEV},:,:)
    u = reshape(U, (/25,365,nlat,nlon/))
    copy_VarCoords(U(0:364,:,:), u(0,:,:,:))

    pentad_u = u(:,0:72,:,:)

    do p = 0,72  ; a year has 73 penta
      pentad_u(:,p,:,:) = dim_sum_n_Wrap(u(:,5*p:5*p+4,:,:), 1)/5.0
    end do

    clim_pentad_U = dim_avg_n_Wrap(pentad_u, 0)
    clim_pentad_u(i,:,:,:) = clim_pentad_U ;(18:35,:,:,:)
    delete([/U,u,pentad_u,clim_pentad_U/])
    ;-------------------------------------------------------------
    tname = "./F2000_"+casename(i)+".plevel.ThetaSe.daily.clm.nc";v2.
    ft = addfile(tname, "r")
    t = ft->thetase(:,{LEV},:,:)

    pentad_t = t(0:72,:,:)
    do p = 0,72  ; a year has 73 penta
      pentad_t(p,:,:) = dim_sum_n_Wrap(t(5*p:5*p+4,:,:), 0)/5.0
    end do
    clim_pentad_t(i,:,:,:) = pentad_t ;(18:35,:,:,:)
    delete([/t,pentad_t/])
end do
printVarSummary(clim_pentad_u)

;-------------------------------------------------------------------
; SCS onset: U850(5-20N, 110-120E) NCC CMA 
lonstr = 110
lonend = 120
latstr = 5
latend = 20

scs_index1 = dim_avg_n_Wrap(clim_pentad_u(:,:,{latstr:latend},{lonstr:lonend}), (/2,3/))
scs_index2 = dim_avg_n_Wrap(clim_pentad_t(:,:,{latstr:latend},{lonstr:lonend}), (/2,3/))

scs_index = new((/ncase,2,73/), float)
scs_index(:,0,:) = scs_index1
scs_index(:,1,:) = scs_index2
printVarSummary(scs_index)
;-------------------------------------------------------------------

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_CESM_F/ncl/20200513/"
figname = "SCS_circulation_index_evolution_ThetaSe"

penta = ispan(1, 73, 1)
; Colors = (/"red","blue","green4"/) ;,"black","purple"
; Patterns = (/0,1,0,1,1,1/)
; Texts = (/"CTRL","NOHEAT","NOTP"/) ;,"CTRL-NOHEAT","NTNH","NT"
pens = 25
pene = pens+30 ; total 31 pentad
print(penta(pens-1:pene-1))

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 4
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
res@xyDashPatterns = 0
res@xyLineThicknessF = 7
res@xyLineColors = "red"

res@tmXBMode = "Explicit"
;res@tmXBValues = penta(pens-1:pene-1)
res@tmXBLabels = penta(pens:pene)
res@tmXBMinorOn = True
res@tmXTOn = False

res@tmXBPrecision = 2
res@tmYLPrecision = 1
res@tiYAxisString = "m/s"
res@tiXAxisString = "Penta"

res@gsnYRefLine = 0
res@gsnYRefLineColor = "darkgoldenrod1"
res@gsnYRefLineThicknessF = 5
res@gsnYRefLineDashPattern = 2
res@gsnYRefLine = 0

res@gsnLeftString = ""

res2 = res
res2@xyLineColors = "blue"
res2@tiYAxisString = "K"
res2@gsnYRefLine = 335
res2@gsnYRefLineColor = "dodgerblue1"
do i = 0, ncase-1
    res2@gsnLeftString = casename(i)
    plot(i) = gsn_csm_xy2(wks, penta(pens-1:pene-1), scs_index1(i,pens-1:pene-1), scs_index2(i,pens-1:pene-1), res, res2)
end do

;==================================================================================
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineThicknessF  = 5                     ; line thicker
    

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.025                 ; change text size
    res_text@txJust             = "CenterLeft"          ; text justification

; ;============================================================
    Texts = (/"U850","ThetaSe"/)
    Colors = (/"red","blue"/)
    interval = 1
    gx = (/39,42/)
    text_pos = 42.5
    gy1 = (/-2.5,-2.5/)
    gy2 = gy1-interval ;(/15,15/)
    gy3 = gy2-interval
    ; ; gy4 = gy3-interval
    print(gy1)

    res_lines@gsLineColor = Colors(0)
    ; res_lines@gsLineDashPattern = Patterns(0)
    dum11 = gsn_add_polyline(wks, plot(0), gx, gy1, res_lines)
    dum12 = gsn_add_text(wks, plot(0), Texts(0), text_pos, gy1, res_text)
    res_lines@gsLineColor = Colors(1)
    ; res_lines@gsLineDashPattern = Patterns(1)
    dum21 = gsn_add_polyline(wks, plot(0), gx, gy2, res_lines)
    dum22 = gsn_add_text(wks, plot(0), Texts(1), text_pos, gy2, res_text)
    ; res_lines@gsLineColor = Colors(2)
    ; res_lines@gsLineDashPattern = Patterns(2)
    ; dum31 = gsn_add_polyline(wks, plot(0), gx, gy3, res_lines)
    ; dum32 = gsn_add_text(wks, plot(0), Texts(2), text_pos, gy3, res_text)

; draw(plot(0))
pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/2,2/), pres)

end