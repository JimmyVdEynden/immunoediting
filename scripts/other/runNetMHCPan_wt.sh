for i in {1..10} 
do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_all.txt | tail -1` -p temp/wt.pep -l 9 -xls -xlsfile temp/wt_out_$i.xls > /dev/null &
done;
wait # To avoid running rest already
