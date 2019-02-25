begin
 
  
  f1 = addfile("./USGS-gtopo30_1.9x2.5_remap_c050602C_junbin.nc", "r") 
  f2 = addfile("./TPmask.nc", "c")  
  
  M1 = f1->PHIS(:,:) 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;再对原来青藏高原的区域，位势高度统一赋值为10000，地形变率统一赋值为0

;;;;;;;;;;;;;;;;;;;;;;;;;;统一赋值位势高度10000
 do i = 0, 2
    M1(60,40+i)=1
 end do

 do i = 0, 3
    M1(61,39+i)=1
 end do

 do i = 0, 9
    M1(62,33+i)=1
 end do

 do i = 0, 10
    M1(63,32+i)=1
 end do
 do i = 0, 1
    M1(63,25+i)=1
 end do

 do i = 0, 11
    M1(64,31+i)=1
 end do
 do i = 0, 2
    M1(64,25+i)=1
 end do

 do i = 0, 18
    M1(65,25+i)=1
 end do

 do i = 0, 18
    M1(66,25+i)=1
 end do

 do i = 0, 16
    M1(67,27+i)=1
 end do

 do i = 0, 16
    M1(68,27+i)=1
 end do

 do i = 0, 15
    M1(69,28+i)=1
 end do

 do i = 0, 14
    M1(70,29+i)=1
 end do

 do i = 0, 12
    M1(71,31+i)=1
 end do

M1@long_name = "TP mask, TP area=1, others=0"
printVarSummary(M1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;统一赋回原nc文件
  f2->tp=M1

end
