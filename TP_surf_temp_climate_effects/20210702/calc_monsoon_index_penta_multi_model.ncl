begin

lonstr = 40
lonend = 150
latstr = -15
latend = 50

yrStrt = 6
yrLast = 30
year = ispan(yrStrt, yrLast, 1)
nyear = dimsizes(year)

nlat = 96
nlon = 144
nlev = 37

case=(/"CTRL","NOQ1","NOTPQ1_APR","NOTPQ1_MAY"/)
ncase = dimsizes(case)


do i = 0, ncase-1
    U = new((/nyear,365,nlev,nlat,nlon/), float)
    U!0 = "year"
    U&year = year
    V = U
    u = new((/73,nlev,nlat,nlon/), float)
    v = u

    diri = "/home/yangsong3/data-model/wuyt/TP_"+case(i)+"_F2000_F19G16_CAM4/pro/daily/"
    do j = 0, nyear-1
        yr = sprinti("%0.2i", year(j))
        filiu = "F2000_"+case(i)+".plevel.U."+yr+".daily.nc"
        filiv = "F2000_"+case(i)+".plevel.V."+yr+".daily.nc"
        print(filiu)
        fu = addfile(diri+filiu, "r")
        fv = addfile(diri+filiv, "r")

        U(j,:,:,:,:) = fu->U
        V(j,:,:,:,:) = fv->V
    end do
    printVarSummary(U)

    u_clm = dim_avg_n_Wrap(U, 0)
    v_clm = dim_avg_n_Wrap(V, 0)

    do p = 0, 72
        u(p,:,:,:) = dim_sum_n_Wrap(u_clm(5*p:5*p+4,:,:,:), 0)/5
        v(p,:,:,:) = dim_sum_n_Wrap(v_clm(5*p:5*p+4,:,:,:), 0)/5
    end do
    copy_VarCoords(u_clm(0:72,:,:,:), u)
    copy_VarCoords(v_clm(0:72,:,:,:), v)
    delete([/U,V,u_clm,v_clm/])

    ;===============================================================================
    ; Webster-Yang index (Asian Monsoon) 
    ;(U850-U200 averaged over 0-20N, 40E-110E)
    ;===============================================================================
    u_wy = wgt_areaave_Wrap(u(:,:,{0:20},{40:110}), 1, 1, 0)
    printVarSummary(u_wy)
    index_wy = u_wy(:,{850}) - u_wy(:,{200})
    copy_VarCoords(u_wy(:,0), index_wy)
    printVarSummary(index_wy)

    ;===============================================================================
    ; East Asian - Western North Pacific monsoon index (southeast monsoon)Wang etal (2001)
    ;(U850 (5N-15N, 90E-130E) - U850 (22.5N-32.5N, 110E-140E))
    ;===============================================================================
    u_wnp = u(:,{850},:,:)
    u1_wnp = wgt_areaave_Wrap(u_wnp(:,{5:15},{90:130}), 1, 1, 0)
    u2_wnp = wgt_areaave_Wrap(u_wnp(:,{22.5:32.5},{110:140}), 1, 1, 0)
    index_wnp = u1_wnp - u2_wnp
    copy_VarCoords(u1_wnp, index_wnp)
    printVarSummary(index_wnp)

    ;===============================================================================
    ; Dynamic Indian monsoon index Wang etal (2001)
    ; (U850 (5N-15N, 40E-80E) - (U850 20N-30N, 70E-90E))
    ;===============================================================================
    u_sa = u(:,{850},:,:)
    u1_sa = wgt_areaave_Wrap(u_sa(:,{5:15},{40:80}), 1, 1, 0)
    u2_sa = wgt_areaave_Wrap(u_sa(:,{20:30},{70:90}), 1, 1, 0)
    index_sa = u1_sa - u2_sa
    copy_VarCoords(u1_sa, index_sa)
    printVarSummary(index_sa)

    ;===============================================================================
    ; South Asian monsoon index RM1 Lau Kim Yang (2000)
    ; (V850-V200 averaged over 10N-30N, 70E-110E)
    ;===============================================================================
    v_rm1 = wgt_areaave_Wrap(v(:,:,{10:30},{70:110}), 1, 1, 0)
    index_rm1 = v_rm1(:,{850}) - v_rm1(:,{200})
    copy_VarCoords(v_rm1(:,0), index_rm1)
    printVarSummary(index_rm1)

    ;===============================================================================
    ; LKY monsoon index (For east asian monsoon) RM2
    ; (U200 (40N-50N, 110-150E) - (U200 25N-35N, 110E-150E))
    ;===============================================================================
    u_rm2 = u(:,{200},:,:)
    u1_rm2 = wgt_areaave_Wrap(u_rm2(:,{40:50},{110:150}), 1, 1, 0)
    u2_rm2 = wgt_areaave_Wrap(u_rm2(:,{25:35},{110:150}), 1, 1, 0)
    index_rm2 = u1_rm2 - u2_rm2
    copy_VarCoords(u1_rm2, index_rm2)
    printVarSummary(index_rm2)

    ;===============================================================================
    ; wtite data to txt file
    ;===============================================================================
    filipath = "./"
    filiname = "monsoon_index_pentad_"+case(i)+".txt"

    write_table(filipath+filiname, "w", [/"Pentad","WY","WNP","Gosiwami","EA"/], "%4s %10s %10s %10s %10s")
    write_table(filipath+filiname, "a", [/ispan(1, 73, 1),index_wy,index_wnp,index_rm1,index_rm2/], "%4.2i %10.2f %10.2f %10.2f %10.2f")

end do






end