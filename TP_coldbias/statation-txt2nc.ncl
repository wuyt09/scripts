load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

    ;先设置两个常数，一个是缺省值，一个是站点数
    FillValue = -9999.0                                 ;缺省值设置为-9999
    NS=839                                              ;站点数共839个(number of station)

    ;开始设置站点序列，后面用作data文件的一个维度
    sn=ispan(1,839,1)                                   ;给每个站做标记1到839
    sn!0="sn"                                      ;把station最左边维的名字就设为station吧，毕竟本来就是代表站点的
    sn@long_name="Station Sequence Number"              ;给它设置个名字：站点序列个数
    sn@size=839                                         ;设置sn的大小就是839

    file_path="/mnt/d/hong/documents/project/temp/SURF_CLI_CHN_MUL_DAY-TEM-12001-201612.txt"
    nrow=numAsciiRow(file_path)
    ncol=numAsciiCol(file_path)
    data=asciiread(file_path,(/nrow,ncol/),"integer")
    IND=ind(data(:,6) .eq. 1)                           ;找到2016年每个1月所在的位置就是每个新站点的位置
      
    ;找到每个新站点的位置以后，就要开始写这个站点序列，我们从所有文件的最后一个中提取站点名字（应该是839个，越早的文件站点数越少）
    id=int2flt(data(IND,0))                             ;根据刚才找到的IND，设置id为站号
    id!0="sn"
    id&sn=sn                                            ;站号以1-839为坐标排序
    id@units = "None"
    id@long_name = "Current ID of Corresponding Station"
    id@_FillValue = FillValue

    ;站号写完了以后就是要写序列的其他信息了：纬度、经度、高度等
    lat=new(NS,float,FillValue)
    lat=int2flt(data(IND,1)/100)+int2flt(data(IND,1)%100)/60.
    lat!0="sn"
    lat&sn=sn                                           ;纬度也以1-839排序
    lat@units = "degree"
    lat@long_name = "Current Latitude of Corresponding Station"
    lat@_FillValue = FillValue

    lon=new(NS,float,FillValue)
    lon=int2flt(data(IND,2)/100)+int2flt(data(IND,2)%100)/60.
    lon!0="sn"
    lon&sn=sn
    lon@units = "degree"
    lon@long_name = "Current Longitude of Corresponding Station"
    lon@_FillValue = FillValue

    altitude=new(NS,float,FillValue)
    altitude=int2flt(data(IND,3))
    altitude!0="sn"
    altitude&sn=sn
    altitude@units = "0.1 meter"
    altitude@long_name = "Current altitude of Corresponding Station"
    altitude@_FillValue = FillValue

    ;至此已经完成经度、纬度、高度的属性设置了。最后用于nc文件的filevardef
    delete(IND)
    delete(file_path)
    delete(data)

    ;数一下有多少年要读
    year=ispan(1961,2016,1)
    count=0
    do i=0,dimsizes(year)-1
        if (year(i)%4 .eq. 0) then
            count=count+366
        else
            count=count+365
        end if
    end do
    ntime=count

    ;数好有多少年以后，开始制造time1用于和time2匹配写入data
    time   = new(ntime,double)
    Year   = new(ntime,integer)
    Month  = new(ntime,integer)
    Day    = new(ntime,integer)
    Hour   = new(ntime,integer)
    Minute = new(ntime,integer)
    Second = new(ntime,integer)
    Hour = 0
    Minute = 0
    Second = 0
    ;
    count=0;重新记数啦！注意闰年哦
    do i=0,dimsizes(year)-1
        if (year(i)%4 .eq. 0) then
           Year(count:count+365)=year(i)
           Month(count:count+30)=1
           Day(count:count+30)=ispan(1,31,1)
           Month(count+31:count+59)=2
           Day(count+31:count+59)=ispan(1,29,1)
           Month(count+60:count+90)=3
           Day(count+60:count+90)=ispan(1,31,1)
           Month(count+91:count+120)=4
           Day(count+91:count+120)=ispan(1,30,1)
           Month(count+121:count+151)=5
           Day(count+121:count+151)=ispan(1,31,1)
           Month(count+152:count+181)=6
           Day(count+152:count+181)=ispan(1,30,1)
           Month(count+182:count+212)=7
           Day(count+182:count+212)=ispan(1,31,1)
           Month(count+213:count+243)=8
           Day(count+213:count+243)=ispan(1,31,1)
           Month(count+244:count+273)=9
           Day(count+244:count+273)=ispan(1,30,1)
           Month(count+274:count+304)=10
           Day(count+274:count+304)=ispan(1,31,1)
           Month(count+305:count+334)=11
           Day(count+305:count+334)=ispan(1,30,1)
           Month(count+335:count+365)=12
           Day(count+335:count+365)=ispan(1,31,1)
           count=count+366
        else
           Year(count:count+364)=year(i)
           Month(count:count+30)=1
           Day(count:count+30)=ispan(1,31,1)
           Month(count+31:count+58)=2
           Day(count+31:count+58)=ispan(1,28,1)
           Month(count+59:count+89)=3
           Day(count+59:count+89)=ispan(1,31,1)
           Month(count+90:count+119)=4
           Day(count+90:count+119)=ispan(1,30,1)
           Month(count+120:count+150)=5
           Day(count+120:count+150)=ispan(1,31,1)
           Month(count+151:count+180)=6
           Day(count+151:count+180)=ispan(1,30,1)
           Month(count+181:count+211)=7
           Day(count+181:count+211)=ispan(1,31,1)
           Month(count+212:count+242)=8
           Day(count+212:count+242)=ispan(1,31,1)
           Month(count+243:count+272)=9
           Day(count+243:count+272)=ispan(1,30,1)
           Month(count+273:count+303)=10
           Day(count+273:count+303)=ispan(1,31,1)
           Month(count+304:count+333)=11
           Day(count+304:count+333)=ispan(1,30,1)
           Month(count+334:count+364)=12
           Day(count+334:count+364)=ispan(1,31,1)
           count=count+365
        end if
    end do

    time1=Year*10000+Month*100+Day

    ;设置时间维度
    time_units = "days since 1900-1-1 00:00:0.0"
    time = ut_inv_calendar(Year,Month,Day,Hour,Minute,Second,time_units,0);为了写入nc，对time进行了改编
    time!0 = "time"
    time@long_name = "time"
    time@units = time_units
    time@actual_range = (/"19510101-20161231"/)

    ;最后要填数据了！我要的是一天的最高温度，所以设置了Tmax
    Tmax=new((/ntime,NS/),float,FillValue)

    do i=0,dimsizes(year)-1
        do j=1,12
          nrow = numAsciiRow("/mnt/d/hong/documents/project/temp/SURF_CLI_CHN_MUL_DAY-TEM-12001-"+year(i)+sprinti("%0.2i",j)+".txt")
          data=asciiread("/mnt/d/hong/documents/project/temp/SURF_CLI_CHN_MUL_DAY-TEM-12001-"+year(i)+sprinti("%0.2i",j)+".txt",(/nrow,13/),"integer")
          time2=data(:,4)*10000+data(:,5)*100+data(:,6)
          print(i*12+j)
            do k=0,NS-1
                IND=ind(int2flt(data(:,0)) .eq. id(k))
                if (.not. any (ismissing(IND))) then        ;如果IND不是缺省值的话
                  index=get1Dindex(time1,time2(IND))        ;得到time1中和time2(IND)相等的数字所在的位置
                  T(index,k)=int2flt(data(IND,8))           ;最高温度在第九列
                  delete(index)
                end if
                delete(IND)
            end do
        delete(data)
        delete(time2)
        end do
     end do

     Tmax=T/10                                              ;最高温度以0.1°为单位
     ;至此已经完成数据填写了。

    altitude=where(altitude .gt. 100000, altitude-100000,altitude);海拔高度如果为估计值则+100000，我这里还原了，不管它是否估计值一律采用，因为这和我没啥关系啊，哈哈！
    Tmax=where(Tmax .eq. 32766,FillValue,Tmax)              ;数据集里32766是缺测值，如果是32766，那么我们就把他设置为-9999
    Tmax!0="time"
    Tmax!1="sn"
    Tmax&time=time
    Tmax&sn=id
    Tmax@long_name = "Daily Hightest Temperature Observation of 839 Station in China"
    Tmax@units = "0.1 degree"
    Tmax@_FillValue = FillValue
    Tmax@statistic = "Daily Mean"
    Tmax@level_desc = "Surface"
    Tmax@dataset = "China Meteorology Station Observation Data"
    Tmax@var_desc = "Daily Tem. Obs"
    ;至此数据的属性完成了

    file_path="/mnt/c/usrs/hong/desktop/"           
    file_name="Tmax.nc"
    if(isfilepresent(file_path+file_name))then
        system("rm -rf"+file_name)                          ;如果Tmax.nc已经存在 就删除
    end if                                             
    ;创建新的nc文件
    fout=addfile(file_name,"c")

    ;定义文件的模式
    setfileoption(fout,"DefineModel",True)                  
    ;创建文件全球属性(书P74)
    fAtt=True
    fAtt@title="Observed Temperature at meteorology station"
    fAtt@source= "Meteorology Observation Network, CMA"
    fAtt@history = "created July 24th 2018 by Haixu Hong with NCL at CMA "
    fAtt@platform = "Observation"
    fileattdef(fout,fAtt)                                   ;把全球属性放进fout里

    ;重新定义坐标的变量和维
    dimNames = (/"time","sn"/)                              ;最高温度的维度分别时间和站号
    dimSizes = (/-1,NS/)                                    ;时间序列不限制长度，站号序列正好839
    dimUnlim = (/True,False/)                               ;维度时间序列没有限制，维度站号有限制
    filedimdef(fout,dimNames,dimSizes,dimUnlim)           
    ;时间、站点、维度、经度、高度、温度属性进行设置  
    filevardef(fout,"time",typeof(time),getvardims(time))   
    filevardef(fout,"sn",typeof(id),getvardims(id))
    filevardef(fout,"lat",typeof(lat),getvardims(lat))
    filevardef(fout,"lon",typeof(lon),getvardims(lon))
    filevardef(fout,"altitude",typeof(altitude),getvardims(altitude))
    filevardef(fout,"Tmax",typeof(Tmax),getvardims(Tmax))
    filevarattdef(fout,"time",time)
    filevarattdef(fout,"sn",id)
    filevarattdef(fout,"lat",lat)
    filevarattdef(fout,"lon",lon)
    filevarattdef(fout,"altitude",lon)
    filevarattdef(fout,"Tmax",Tmax)

    setfileoption(fout,"DefineModel",False)
    ;退出文件定义模式（可不用）
    fout->time = (/time/)
    fout->sn = (/id/)
    fout->lat = (/lat/)
    fout->lon = (/lon/)
    fout->altitude = (/altitude/)
    fout->Tmax = (/Tmax/)

end
