for i in {1..10} 
do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_all.txt | tail -1` -p temp/wt_sim.pep -l 9 -xls -xlsfile temp/wt_out_sim_$i.xls > /dev/null &
done;
wait # To avoid running rest already
