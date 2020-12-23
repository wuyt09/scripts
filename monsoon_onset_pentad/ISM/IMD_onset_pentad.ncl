
begin

yrStrt = 1970
yrLast = 2020
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

x = new((/nyear,3/), integer)
x(:,0) = year


x(0,1) = day_of_year(year(11), 5, 26)
x(1,1) = day_of_year(year(11), 5, 27)
x(2,1) = day_of_year(year(11), 6, 18)
x(3,1) = day_of_year(year(11), 6, 4)
x(4,1) = day_of_year(year(11), 5, 26)
x(5,1) = day_of_year(year(11), 5, 31)
x(6,1) = day_of_year(year(11), 5, 31)
x(7,1) = day_of_year(year(11), 5, 30)
x(8,1) = day_of_year(year(11), 5, 28)
x(9,1) = day_of_year(year(11), 6, 13)

x(10,1) = day_of_year(year(11), 6, 1)
x(11,1) = day_of_year(year(11), 5, 30)
x(12,1) = day_of_year(year(11), 5, 28)
x(13,1) = day_of_year(year(11), 6, 12)
x(14,1) = day_of_year(year(11), 5, 31)
x(15,1) = day_of_year(year(11), 5, 28)
x(16,1) = day_of_year(year(11), 6, 4)
x(17,1) = day_of_year(year(11), 6, 2)
x(18,1) = day_of_year(year(11), 5, 25)
x(19,1) = day_of_year(year(11), 6, 3)

x(20,1) = day_of_year(year(11), 5, 19)
x(21,1) = day_of_year(year(11), 6, 2)
x(22,1) = day_of_year(year(11), 6, 5)
x(23,1) = day_of_year(year(11), 5, 28)
x(24,1) = day_of_year(year(11), 5, 28)
x(25,1) = day_of_year(year(11), 6, 8)
x(26,1) = day_of_year(year(11), 6, 3)
x(27,1) = day_of_year(year(11), 6, 9)
x(28,1) = day_of_year(year(11), 6, 2)
x(29,1) = day_of_year(year(11), 5, 25)

x(30,1) = day_of_year(year(11), 5, 31)
x(31,1) = day_of_year(year(11), 5, 26)
x(32,1) = day_of_year(year(11), 6, 9)
x(33,1) = day_of_year(year(11), 6, 13)
x(34,1) = day_of_year(year(11), 6, 4)
x(35,1) = day_of_year(year(11), 6, 7)
x(36,1) = day_of_year(year(11), 5, 26)
x(37,1) = day_of_year(year(11), 5, 28)
x(38,1) = day_of_year(year(11), 5, 31)
x(39,1) = day_of_year(year(11), 5, 23)

x(40,1) = day_of_year(year(11), 5, 31)
x(41,1) = day_of_year(year(11), 5, 29)
x(42,1) = day_of_year(year(11), 6, 2)
x(43,1) = day_of_year(year(11), 6, 1)
x(44,1) = day_of_year(year(11), 6, 6)
x(45,1) = day_of_year(year(11), 6, 4)
x(46,1) = day_of_year(year(11), 6, 8)
x(47,1) = day_of_year(year(11), 5, 30)
x(48,1) = day_of_year(year(11), 5, 29)
x(49,1) = day_of_year(year(11), 6, 8)
x(50,1) = day_of_year(year(11), 6, 1)


x(:,2) = tointeger(x(:,1)/5) + 1

print(year+"    "+x(:,1)+"    "+x(:,2))

write_table("IMD_onset_pentad_"+yrStrt+"_"+yrLast+".txt", "w", [/"YEAR","day","pentad"/], "%4s %4s %4s")
write_table("IMD_onset_pentad_"+yrStrt+"_"+yrLast+".txt", "a", [/x(:,0),x(:,1),x(:,2)/], "%0.4i %0.2i %0.2i")

end