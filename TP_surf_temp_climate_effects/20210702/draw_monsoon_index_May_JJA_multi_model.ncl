
begin

casename = (/"CTRL","NOQ1", "NOTPQ1_APR", "NOTPQ1_MAY"/)
ncase = dimsizes(casename)

x = new((/ncase,2,4/), float)
do i = 0, ncase-1
    fili1 = "May_monsoon_index_25year_"+casename(i)+".txt"
    fili2 = "JJA_monsoon_index_25year_"+casename(i)+".txt"

    x1 = asciiread(fili1, (/25,4/), "float")
    x2 = asciiread(fili2, (/25,4/), "float")
    ; print(x1(24,:))
    x(i,0,:) = dim_avg_n_Wrap(x1, 0)
    x(i,1,:) = dim_avg_n_Wrap(x2, 0)
end do

print(sprintf("%4.2f", x(0,:,:)))
print("theraml "+sprintf("%4.2f", (x(0,:,:)-x(1,:,:))/x(0,:,:)))
print("Apr theraml "+sprintf("%4.2f", (x(0,:,:)-x(2,:,:))/x(0,:,:)))
print("May theraml "+sprintf("%4.2f", (x(0,:,:)-x(3,:,:))/x(0,:,:)))


end