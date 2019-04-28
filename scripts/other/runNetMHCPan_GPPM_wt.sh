for i in {1..6}
do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_freq.txt | tail -1` -p temp/GPPM_wt.pep -l 9 -xls -xlsfile temp/GPPM_wt_out_$i.xls > /dev/null & 
done;
wait # To avoid running rest already
