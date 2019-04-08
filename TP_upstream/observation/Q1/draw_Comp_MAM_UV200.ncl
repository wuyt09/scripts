
begin

yrStrt = 1979
yrLast = 2013
year = ispan(yrStrt, yrLast, 1)
ny = dimsizes(year)
Season = (/"MAM","JJA","SON","DJF"/)
ns = dimsizes(Season)

;---------------------------------------------
path = "./"
name = "Q1_TP80E_dtrend_standardize.txt"

x = asciiread(path+name, (/ny,ns/), "float")
x@_FillValue = -32767
x!0 = "year"
x&year = year
x!1 = "season"
x&season = Season
printVarSummary(x)

index = x(:,0) ;MAM
indlow = ind(index .le. -0.8)
indhig = ind(index .ge. 0.8)
numhig = dimsizes(indhig)
numlow = dimsizes(indlow)
yearhig = year(indhig)
yearlow = year(indlow)
;---------------------------------------------

diri = "/home/yangsong3/data-observation/NCEP_DOE/monthly/"
filu = "uwnd.mon.mean.nc"
filv = "vwnd.mon.mean.nc"

fu = addfile(diri+filu, "r")
fv = addfile(diri+filv, "r")

TIME = fu->time
YYYY = cd_calendar(TIME, -1)/100
iYYYY = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)

U = short2flt(fu->uwnd(iYYYY,{200},:,:))
V = short2flt(fv->vwnd(iYYYY,{200},:,:))
printVarSummary(U)

uclm = clmMonTLL(U)
printVarSummary(uclm)
u = calcMonAnomTLL(U, uclm)
vclm = clmMonTLL(V)
v = calcMonAnomTLL(V, vclm)
printVarSummary(u)

u_spr = month_to_season(u, "MAM")
v_spr = month_to_season(v, "MAM")
u_spr!0 = "year"
u_spr&year = year
copy_VarCoords(u_spr, v_spr)
printVarSummary(u_spr)

uhig = u_spr(indhig,:,:)
ulow = u_spr(indlow,:,:)
vhig = v_spr(indhig,:,:)
vlow = v_spr(indlow,:,:)
printVarSummary(uhig)

uhigvar = dim_variance_n_Wrap(uhig, 0)
ulowvar = dim_variance_n_Wrap(ulow, 0)
uhigave = dim_avg_n_Wrap(uhig, 0)
ulowave = dim_avg_n_Wrap(ulow, 0)
probu = ttest(uhigave, uhigvar, numhig, ulowave, ulowvar, numlow, False, False)

vhigvar = dim_variance_n_Wrap(vhig, 0)
vlowvar = dim_variance_n_Wrap(vlow, 0)
vhigave = dim_avg_n_Wrap(vhig, 0)
vlowave = dim_avg_n_Wrap(vlow, 0)
probv = ttest(vhigave, vhigvar, numhig, vlowave, vlowvar, numlow, False, False)

udiff = uhigave - ulowave
vdiff = vhigave - vlowave
copy_VarCoords(u(0,:,:), udiff)
copy_VarCoords(v(0,:,:), vdiff)

;udiff = mask(udiff, (probu.lt.0.05).or.(probv.lt.0.05), True)
;vdiff = mask(vdiff, (probu.lt.0.05).or.(probv.lt.0.05), True)

;==============divergence=======================================================
divhig = uv2dv_cfd(uhigave, vhigave, u&lat, u&lon, 1)
divlow = uv2dv_cfd(ulowave, vlowave, u&lat, u&lon, 1)
divhig = divhig*1.0e+6
divlow = divlow*1.0e+6
divdiff = divhig - divlow
copy_VarCoords(u(0,:,:), divhig(:,:))
copy_VarCoords(u(0,:,:), divlow(:,:))
copy_VarCoords(u(0,:,:), divdiff(:,:))
printVarSummary(divhig)
printVarSummary(divlow)
;==================================================================

figpath = "/home/yangsong3/wuyt/sysu/figures/TP_upstream/climatology/"
figname = "draw_Comp_MAM_UV200"
latstr = 0
latend = 70
lonstr = -30
lonend = 150

wks = gsn_open_wks("pdf", figpath+figname)
nplot = 3
plot = new(nplot, graphic)
vcid = plot
plot_topo = plot

;------------------------------------------------------------
; TOPO
;------------------------------------------------------------
 
 f2 = "~/wuyt/data/TOPO_1500.nc"
 topo_in = addfile(f2,"r")
 topo = topo_in->TP

 xlat = topo_in->LATIXY(:,0)
 xlon = topo_in->LONGXY(0,:)
 topo!0 = "lsmlat"
 topo!1 = "lsmlon"
 topo&lsmlat = xlat
 topo&lsmlon = xlon
 
 res_topo            = True
 res_topo@cnFillOn   = False
 res_topo@cnLinesOn  = True
 res_topo@gsnFrame   = False
 res_topo@gsnDraw    = False
 res_topo@cnLineLabelsOn         = False
 res_topo@cnLineThicknessF       = 3
 res_topo@cnLineColor            = "gray30"
 res_topo@cnInfoLabelOn          = False
 res_topo@cnLevelSelectionMode   ="ExplicitLevels"
 res_topo@cnLevels               = (/1500/)
;------------------------------------------------------------
 res = True
 res@gsnDraw  = False
 res@gsnFrame = False
 res@gsnAddCyclic = True

 res@mpFillOn = True
 res@mpMinLatF = latstr
 res@mpMaxLatF = latend
 res@mpMinLonF = lonstr
 res@mpMaxLonF = lonend
 res@mpGeophysicalLineThicknessF =2
 res@mpGeophysicalLineColor = "grey40"

 res@cnFillOn = True
 res@cnLinesOn = False
 res@cnLineLabelsOn = False
 res@lbLabelBarOn = False
 res@cnFillPalette = "BlueWhiteOrangeRed"

 res@tmYLLabelFontHeightF = 0.02
 res@tmXBLabelFontHeightF = 0.02
 res@gsnStringFontHeightF = 0.03
 res@gsnRightString = "1.0e-6/s"

 res@cnLevelSelectionMode = "ManualLevels"  ;ExplicitLevels
 res@cnMaxLevelValF = 1.
 res@cnMinLevelValF = -1.
 res@cnLevelSpacingF = 0.2

;************************************************
; setting : vector
;************************************************

  vecres                         = True            ; vector only resources
  vecres@gsnDraw                 = False           ; don't draw
  vecres@gsnFrame                = False           ; don't advance frame
  vecres@vcGlyphStyle            = "LineArrow" ;"CurlyVector"   ; curly vectors
  vecres@vcRefMagnitudeF         = 2 ;1.3               ; define vector ref mag
  vecres@vcRefLengthF            = 0.045           ; define length of vec ref
;  vecres@vcRefAnnoOn             = False
  vecres@gsnRightString          = " "             ; turn off right string
  vecres@gsnLeftString           = " "             ; turn off left string
  vecres@tiXAxisString           = " "             ; turn off axis label
  vecres@vcRefAnnoString2On      = False
  vecres@vcRefAnnoOrthogonalPosF = -1.    ; move ref vector into plot
  vecres@vcRefAnnoFontHeightF    = .025    ; move ref vector into plot
  vecres@vcLineArrowHeadMaxSizeF      = 0.01
  vecres@vcLineArrowHeadMinSizeF      = 0.008
  vecres@vcLineArrowThicknessF      = 3

  vecres@vcMonoLineArrowColor    = True

  res@gsnLeftString = "HIGH"
  plot(0) = gsn_csm_contour_map(wks, divhig, res)
  vcid(0) = gsn_csm_vector(wks, uhigave(::2,::2), vhigave(::2,::2), vecres)
  plot_topo(0) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(0), vcid(0))
  overlay(plot(0), plot_topo(0))

  res@gsnLeftString = "LOW"
  plot(1) = gsn_csm_contour_map(wks, divlow, res)
  vcid(1) = gsn_csm_vector(wks, ulowave(::2,::2), vlowave(::2,::2), vecres)
  plot_topo(1) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(1), vcid(1))
  overlay(plot(1), plot_topo(1))

  res@gsnLeftString = "DIFF"
  plot(2) = gsn_csm_contour_map(wks, divdiff, res)
  vcid(2) = gsn_csm_vector(wks, udiff(::2,::2), vdiff(::2,::2), vecres)
  plot_topo(2) = gsn_csm_contour(wks, topo, res_topo)
  overlay(plot(2), vcid(2))
  overlay(plot(2), plot_topo(2))


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@gsnPanelMainString = "Composite MAM UV200 anomaly(A) & Div UV200A"
gsn_panel(wks, plot, (/2,2/), pres)
end