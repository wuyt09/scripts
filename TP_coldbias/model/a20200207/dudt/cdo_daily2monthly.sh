#!bin/bash
casename=("Ctrl_FAMIP" "Hist_FAMIP")
for ((i=1981; i<2006; i++))
do
    cdo monmean ${casename[0]}.dudt.${i}.nc ${casename[0]}.dudt.monmean.${i}.nc
    cdo monmean ${casename[1]}.dudt.${i}.nc ${casename[1]}.dudt.monmean.${i}.nc
    echo "year ${i} is done!"
done

