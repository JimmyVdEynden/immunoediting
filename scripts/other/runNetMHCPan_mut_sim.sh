for i in {1..5} 
do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_all.txt | tail -1` -p temp/mut_sim.pep -l 9 -xls -xlsfile temp/mut_out_sim_$i.xls > /dev/null &
done;
wait
# for i in {6..10}
# do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_all.txt | tail -1` -p temp/mut_sim.pep -l 9 -xls -xlsfile temp/mut_out_sim_$i.xls > /dev/null &
# done;
# wait
