
begin
yrStrt = 1979
yrLast = 2018
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

fx = "time_series_TP_t2m_DJF_MAM_ERA5.txt"
x = asciiread(fx, (/nyear,3/), "float")
inx = x(:,1)

hig_ind = ind(inx.ge.1.)
low_ind = ind(inx.le.-1.)
print(year(hig_ind))
print(year(low_ind))


inx_name = (/"WY","SEASM","Goswami","EA"/) ;"SA",
fili = "calc_monsoon_index_May_ERA5.txt"
monsoon_inx = asciiread(fili, (/nyear,5/), "float")
print(monsoon_inx(nyear-1,:))

monsoon_inx2 = monsoon_inx(:,1:4)

hig_monsoon_med = monsoon_inx2(hig_ind,:)
low_monsoon_med = monsoon_inx2(low_ind,:)
hig_monsoon = dim_avg_n_Wrap(hig_monsoon_med, 0)
low_monsoon = dim_avg_n_Wrap(low_monsoon_med, 0)

print(inx_name+"    "+sprintf("%4.2f", hig_monsoon))
print(inx_name+"    "+sprintf("%4.2f", low_monsoon))
end