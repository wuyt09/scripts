#!bin/bash

# ln -s /home/yangsong3/huxm/antarctic_melt/drdt_ranc_1.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/cc_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/ciwc_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/clwc_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/hus_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/huss_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/o3_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/ps_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/ssrd_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/ssru_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/ta_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/tisr_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/ts_base.dat
# ln -s /home/yangsong3/huxm/antarctic_melt/base_no_cloud_out_1.dat

for (( i = 1; i < 63; i++ )); do
	ln -s /home/yangsong3/huxm/antarctic_melt/baseline_radranc_$i.grd
	ln -s /home/yangsong3/huxm/antarctic_melt/baseline_radranc_$i.grd

done