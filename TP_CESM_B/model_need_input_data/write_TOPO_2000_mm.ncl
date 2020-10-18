;*************************************************************
;Aim: Calculate teleconnectivity of seasonal stream function
;Mission:
;
;   1. Read data from TOPO datasets
;   2. Output the TOPO
;                       Mengmeng
;                       2016-10-22
;*************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;================================================
;Variables List:
;	var--------------TOPO
;================================================
    in_path = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_B/model_need_input_data/"
    out_path = "/home/yangsong3/wuyt/sysu/scripts/TP_CESM_B/model_need_input_data/"
    out_name = "TOPO_1500.nc"
; ==============================================================
; Mission 1. read the TOPO
; ==============================================================
    TOPO_in = addfile (in_path+"surfdata_1.9x2.5_simyr1850_c130421.nc", "r")
    TOPO    = TOPO_in->TOPO
    LATIXY  = TOPO_in->LATIXY
    LONGXY  = TOPO_in->LONGXY
    printVarSummary (TOPO)
    printVarSummary (LATIXY)
    printVarSummary (LONGXY)
	TP1 = TOPO
	dim = dimsizes(TP1)
    printVarSummary (TP1)
	
    TP_0 = where (TP1.eq.0.,0.,0.)
	TP_00 = flt2dble(TP_0)
    copy_VarCoords (TP1,TP_00)
    printVarSummary (TP_00)
    TP_15 = where (TP1.eq.0.,0.,0.)
	TP_1500 = flt2dble(TP_15)
    copy_VarCoords (TP1,TP_1500)
    TP_1500!0 = "lsmlat"
 	TP_1500&lsmlat = LATIXY(:,0)
    TP_1500!1 = "lsmlon"
    TP_1500&lsmlon = LONGXY(0,:)
    printVarSummary (TP_1500)                 
    TP_1500 = where (TP1.ge.2000.,TP1,TP_00)
	TP_1500({43:},{90:}) = 0.0
	TP = TP_1500
	TP = 0.0
	TP({23:45},{62:105}) = TP_1500({23:45},{62:105})
    print(TP({23:45},{62:105}))
; ==============================================================
; Mission 2. output to a new netCDF file
; ==============================================================
system ("rm -f "+out_path+out_name)      ;remove if exists
fout = addfile (out_path+out_name,"c")   ;open output file

setfileoption (fout,"DefineMode",True)   ;declare file definition mode

;predefine the coordinate variables and their dimensionality
dimNames = (/"lsmlat","lsmlon"/)  
dimSizes = (/dim(0),dim(1)/) 
dimUnlim = (/False,False/)   
filedimdef (fout,dimNames,dimSizes,dimUnlim)

;predefine the the dimensionality of the variables to be written out
filevardef (fout,"LATIXY",typeof(LATIXY),getvardims(LATIXY))                          
filevardef (fout,"LONGXY",typeof(LONGXY),getvardims(LONGXY))                          
filevardef (fout,"TP",typeof(TP),getvardims(TP))

;Copy attributes associated with each variable to the file
filevarattdef (fout,"LATIXY",LATIXY)
filevarattdef (fout,"LONGXY",LONGXY)
filevarattdef (fout,"TP",TP)

;output only the data values
fout->LATIXY  = (/LATIXY/)
fout->LONGXY  = (/LONGXY/) 
fout->TP = (/TP/)

end
