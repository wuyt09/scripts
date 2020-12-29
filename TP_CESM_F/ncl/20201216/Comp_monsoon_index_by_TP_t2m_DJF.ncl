
begin

yrStrt = 1979
yrLast = 2018
mnStrt = 1
mnLast = 12

year = ispan(yrStrt, yrLast, 1)
month = ispan(mnStrt, mnLast, 1)
penta = ispan(1, 73, 1)
nyear = dimsizes(year)
nmonth = dimsizes(month)
npenta = dimsizes(penta)
season = (/"JJA"/)

fx = "time_series_TP_t2m_DJF_MAM_OBS.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(0:nyear-1,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))

; hig_ind = hig_ind+1
; low_ind = low_ind+1


index_name = (/"WY","SA-Goswami","SA-Wang","SEA-Wang","EA-LKY","EA-Wang"/)

fy = "calc_monsoon_index_JJA.txt"
y = asciiread(fy, (/nyear,7/), "float")
iny = y(:,1:6)
; print(iny(:,3))

hig_var = dim_avg_n_Wrap(iny(hig_ind,:), 0)
low_var = dim_avg_n_Wrap(iny(low_ind,:), 0)
print(sprintf("%4.2f", hig_var))
print(("=================================="))
print(sprintf("%4.2f", low_var))
end