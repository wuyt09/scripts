
;=====================================================================
; Mission: write a heating file based on CTRL output, in terms of 
; different topography, some process are needed.
;
;  1. get the condensation and diffusion heating data.
;
;  2. transform the heating data from hybird to p-level. 
;
;  3. get the TP's average pressure at 26 layers in NoTOPO run.
;    
;  4. Linear interpolation the heating data (p-level) to a new 
;     heating data (p-level but for NoTOPO)
;
;  5. transform the new heating file to hybrid (NoTOPO)
;======================================================================


begin

;=================================================
; 1. read data, get the heating data
;=================================================

diri1 = "/home/yangsong3/data-model/wuyt/TP_CTRL_F2000_F19G16_CAM4/pro/"
filia = "TP_CN.hyam.mon.0601-3012.nc"
filib = "TP_CN.hybm.mon.0601-3012.nc"
filip = "TP_CN.PS.mon.0601-3012.nc"
filih1 = "TP_CN.DTCOND.mon.0601-3012.nc"
filih2 = "TP_CN.DTV.mon.0601-3012.nc"

diri2 = "/home/yangsong3/data-model/wuyt/TP_NoTOPO_F2000_F19G16_CAM4/pro/"
filiA = "TP_CR.hyam.mon.0601-3012.nc" 
filiB = "TP_CR.hybm.mon.0601-3012.nc"
filiP = "TP_CR.PS.mon.0601-3012.nc"
filiphis = "TP_CR.PHIS.mon.0601-3012.nc"

f1a = addfile(diri1+filia, "r")
f1b = addfile(diri1+filib, "r")
f1p = addfile(diri1+filip, "r")
f1h1 = addfile(diri1+filih1, "r")
f1h2 = addfile(diri1+filih2, "r")

f2a = addfile(diri2+filiA, "r")
f2b = addfile(diri2+filiB, "r")
f2p = addfile(diri2+filiP, "r")
f2phis = addfile(diri2+filiphis, "r")

chyam = f1a->hyam
chybm = f1b->hybm
cps = f1p->PS
h1 = f1h1->DTCOND
h2 = f1h2->DTV
printVarSummary(h1)


rhyam = f2a->hyam
rhybm = f2b->hybm
rps = f2p->PS
phis = f2phis->PHIS


;=================================================
; 2. transform from hybrid to p-level
;=================================================

 h = h1 ;copy meta data
 h = h1 + h2
 h@long_name = "total heating"
 printVarSummary(h)

 plevo = (/1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700,\
         650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175,\
         150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1/)
 plevo@units = "hPa"
 p0 = 1000.
 p0@units = "hPa"

 h_plevel = vinth2p(h, chyam, chybm, plevo, cps, 1, p0, 1, False)
 printVarSummary(h_plevel)
 h_plevel!1 = "lev"
 h_plevel&lev = plevo
 h_plevel@_FillValue = -999
 printVarSummary(h_plevel)


;==============================================================
; 3. get the TP's average pressure at 26 layers in NoTOPO run.
;==============================================================
 ; 1. 检查海平面和高原在650hPa左右的Phis位势高度差别，此为误差来源1
 ; 2. 线性插值，按照混合坐标给的参考气压层，此为误差来源2
 ; 3. 误差各自取最大相加，若误差在10%内，则认为该法可信
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;






end