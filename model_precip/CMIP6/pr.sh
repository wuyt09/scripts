#!/bin/bash
#cmip6_dir=/home/yangsong3/data-observation/cmip6/
#for Model_Name in `ls $cmip6_dir`
# Model_Name_total=("ACCESS-CM2")   
# Model_Name_total=("CESM2" "CESM2-WACCM", "CNRM-ESM2-1"\
#      "EC-Earth3" "EC-Earth3-Veg" "GFDL-ESM4" "HadGEM3-GC31-LL" "IPSL-CM6A-LR" "MIROC6"  \
#         "MIROC-ES2L" "MPI-ESM1-2-HR" "MRI-ESM2-0" "NESM3" "NorESM2-LM" "UKESM1-0-LL")
Model_Name_total=("CAMS-CSM1-0" "CanESM5" "CESM2" "CESM2-WACCM" "CNRM-ESM2-1" \
  "CNRM-CM6-1" "E3SM-1-0" "EC-Earth3" "EC-Earth3-Veg" "FGOALS-f3-L" "FGOALS-g3" \
  "GISS-E2-1-G-CC" "GISS-E2-1-G" "HadGEM3-GC31-LL" "IPSL-CM6A-LR" "MPI-ESM1-2-HR" \
  "SAM0-UNICON" "NESM3" "UKESM1-0-LL")


for Model_Name in ${Model_Name_total[*]}
do 
  echo "********start  $Model_Name  *************************"
  OUT_DIR=/home/yangsong3/wuyt/sysu/scripts/model_precip/CMIP6/
  # if [ ! -d $OUT_DIR ]; then
  #   mkdir $OUT_DIR

    IN_DIR=/home/yangsong3/data-observation/cmip6/pr/
    cd $IN_DIR
    aa=`ls pr_Amon_${Model_Name}_historical_* | wc -l`

    echo $aa
    cd $OUT_DIR
    if [[ $aa -eq 1 ]]; then
      ln -s $IN_DIR/pr_Amon_${Model_Name}_historical_* ./
    else
      cdo cat $IN_DIR/pr_Amon_${Model_Name}_historical_* ./pr_Amon_${Model_Name}_historical_r1i1p1f1_g_185001-201412.nc
    fi
  #   pwd
  #   for file1 in `ls ua*`
  #   do
  #   cdo -s remapbil,global_2.5 $IN_DIR/$file1 $OUT_DIR/$file1
  #   done

  #   for file2 in `ls va*`
  #   do
  #   cdo -s remapbil,global_2.5 $IN_DIR/$file2 $OUT_DIR/$file2
  #   done
  # # fi
  #  echo "********finish  $Model_Name  *************************"
done

  echo "finish"

