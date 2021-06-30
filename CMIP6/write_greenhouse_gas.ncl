
begin

yrStrt = 1765
yrLast = 2500
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)


fili = "./RCP85_MIDYR_CONC.DAT"
X = asciiread(fili, (/nyear,36/), "float")
printVarSummary(X)
X!0 = "year"
X&year = year
print(X(nyear-1,:))

X(:,4) = X(:,4)/1000. ; units ppm
X(:,5) = X(:,5)/1000. ; units ppm
x = X(:,3:5)
printVarSummary(x)

yrs = 1979
yre = 2014
base_s = 1979 ;ispan(1979, 2004, 1)
base_e = 1988 ;ispan(1988, 2013, 1)
warm_s = ispan(1980, 2005, 1)
warm_e = ispan(1989, 2014, 1)
n = dimsizes(warm_s)

; X_b = x({year|base_s(0):base_e(0)},:)
X_b = x({base_s(0):base_e(0)},:)
X_w = new((/n,10,3/), float)

do ll = 0, n-1
    ; print(warm_s(ll))
    X_w(ll,:,:) = x({warm_s(ll):warm_e(ll)},:)
end do

x_b = dim_avg_n_Wrap(X_b, 0)
x_w = dim_avg_n_Wrap(X_w, 1)
printVarSummary(x_b)
printVarSummary(x_w)

print(sprintf("%4.2f", x_b))

outpath = "./"
filib1 = "co2_warm.dat"
system("rm -f "+outpath+filib1)
fbindirwrite(outpath+filib1,x_w(:,0))

filib2 = "ch4_warm.dat"
system("rm -f "+outpath+filib2)
fbindirwrite(outpath+filib2,x_w(:,1))

filib3 = "n2o_warm.dat"
system("rm -f "+outpath+filib3)
fbindirwrite(outpath+filib3,x_w(:,2))

x_w!0 = "year"
x_w!1 = "var"
x_w&year = ispan(1984, 2009, 1)
lines = x_w(var|:,year|:)
;=====================================================
figpath = "/home/yangsong3/wuyt/sysu/scripts/CMIP6/"
figname = "write_greenhouse_gas"

wks = gsn_open_wks("pdf", figpath+figname)

plot = new(3, graphic)

res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

; res@trXMinF = 1
; res@trXMaxF = 12

Colors = (/"black","purple","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

res@gsnLeftString = "CO2"
plot(0) = gsn_csm_xy(wks, ispan(1984, 2009, 1), lines(0,:), res)
res@gsnLeftString = "CH4"
plot(1) = gsn_csm_xy(wks, ispan(1984, 2009, 1), lines(1,:), res)
res@gsnLeftString = "N2O"
plot(2) = gsn_csm_xy(wks, ispan(1984, 2009, 1), lines(2,:), res)

pres = True
pres@gsnMaximize = True
gsn_panel(wks, plot, (/3,1/), pres)
end