
begin

yrStrt = 1959
yrLast = 2020
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)
season = "JJA"

fili = "/home/yangsong3/wuyt/sysu/scripts/TP_surrounding/cfram20210212/data/co2_mm_mlo.txt/"
x = asciiread(fili, (/nyear*12,8/), "float")

printVarSummary(x)
print(x(nyear*12-1,:))

y = x(:,3)
y!0 = "time"
y&time = ispan(1, nyear*12, 1)
y_sea = month_to_season(y, season)
printVarSummary(y_sea)
y_sea!0 = "year"
y_sea&year = year
printVarSummary(y_sea)

var = y_sea({year|1979:2019})
printVarSummary(var)

base_s = ispan(1979, 2009, 1)
base_e = ispan(1988, 2018, 1)
warm_s = ispan(1980, 2010, 1)
warm_e = ispan(1989, 2019, 1)
n = dimsizes(base_s)
print("n = "+n)

X_b = new((/n,10/), float)
X_w = X_b

do ll = 0, n-1 
    X_b(ll,:) = var({year|base_s(ll):base_e(ll)})
    X_w(ll,:) = var({year|warm_s(ll):warm_e(ll)})
end do

xb = dim_avg_n_Wrap(X_b, 1)
xw = dim_avg_n_Wrap(X_w, 1)
copy_VarCoords(var(0:n-1), xb)
copy_VarCoords(var(0:n-1), xw)
printVarSummary(xb)

i = 0
outname = (/"co2"/)
outpath = "./"
basename = outname(i)+"_base.dat"
warmname = outname(i)+"_warm.dat"

system("rm -f "+outpath+basename)
fbindirwrite(outpath+basename,xb)
system("rm -f "+outpath+warmname)
fbindirwrite(outpath+warmname,xw)

;==================================================================
Solid = new((/3,n/), float)
Solid(0,:) = xb
Solid(1,:) = xw
Solid(2,:) = xw - xb

wks = gsn_open_wks("pdf",outname(i))
gsn_define_colormap(wks,"BlueWhiteOrangeRed")
plot = new(1,graphic)


res = True
res@gsnDraw  = False
res@gsnFrame = False
res@gsnMaximize = True

res@vpWidthF = 0.7
res@vpHeightF = 0.4

; res@trXMinF = yrStrt
; res@trXMaxF = yrLast
; res@tmXBMode = "Explicit"
; res@tmXBValues = year
; res@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
res@tmXTOn = False
; res@tiXAxisString = "Month"

; res@trYMinF = -20
; res@trYMaxF = 20
; res@tmYLMode = "Explicit"
; res@tmYLValues = ispan(-20, 20, 5)
; res@tmYLLabels = ispan(-20, 20, 5)
; res@tmYROn = False
; res@tiYAxisString = "Temperature (~S~o~N~C)"
;res@gsnLeftString = "Tas: grid2station"

;Colors = (/"deeppink2", "blue", "purple","green4","orange2", "black" /);
Colors = (/"black","red","forestgreen"/)
Patterns = (/0,0,0/)
res@xyLineThicknessF = 15
res@xyLineColors = Colors
res@xyDashPatterns = Patterns

plot(0) = gsn_csm_xy(wks, ispan(1, n, 1), Solid, res)
draw(plot(0))

end