#!bin/bash

diri="/home/yangsong3/data-observation/APHRO_MA_025deg_V1003R1"

#APHRO_MA_025deg_V1003R1.1989.nc

for (( i = 1951; i < 2007; i++ )); do
	cdo monmean $diri/APHRO_MA_025deg_V1003R1.${i}.nc ./APHRO_MA_025deg_V1003R1.monmean.${i}.nc
done

