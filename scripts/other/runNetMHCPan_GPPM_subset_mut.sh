# for i in {1..5}
# do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_all.txt | tail -1` -p temp/GPPM_subset_mut.pep -l 9 -xls -xlsfile temp/GPPM_subset_mut_out_$i.xls > /dev/null &
# done;
# wait
for i in {6..10}
do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_all.txt | tail -1` -p temp/GPPM_subset_mut.pep -l 9 -xls -xlsfile temp/GPPM_subset_mut_out_$i.xls > /dev/null &
done;
wait
