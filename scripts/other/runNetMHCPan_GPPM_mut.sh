# for i in {1..6}
# do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_freq.txt | tail -1` -p temp/GPPM_mut.pep -l 9 -xls -xlsfile temp/GPPM_mut_out_$i.xls > /dev/null & 
# done;
# wait # To avoid running rest already

# # Test: Start at 15h50: +/- 25'--> 2500' ~ 2days for 100 million? ~ 10 days for all?
# for i in {1..6}
# do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_freq.txt | tail -1` -p temp/GPPM_mut_temp1.pep -l 9 -xls -xlsfile temp/GPPM_mut_temp1_out_$i.xls > /dev/null &
# done;
# wait # To avoid running rest already

# # Try
# # Crashes: pep.files too large? 10x larger than previous times!!!
# for i in {1..6}
# do ~/tools/netMHCpan-3.0/netMHCpan -a `head -$i temp/HLA_types_freq.txt | tail -1` -p temp/GPPM_mut_1.pep -l 9 -xls -xlsfile temp/GPPM_mut_1_out_$i.xls > /dev/null &
# done;
# wait # To avoid running rest already

# # Test: with all HLA types at once on 1pep file? 46m --> All by one = 16d!
# time ~/tools/netMHCpan-3.0/netMHCpan -a `head -1 temp/HLA_types_freq_oneLine.txt | tail -1` -p temp/GPPM_mut_temp1.pep -l 9 -xls -xlsfile temp/GPPM_mut_temp1_out_$i.xls > /dev/null &

# Try 1x 50mil pep file? 2d? CHeck RAM!
time ~/tools/netMHCpan-3.0/netMHCpan -a `head -1 temp/HLA_types_freq_oneLine.txt | tail -1` -p temp/GPPM_mut_temp1.pep -l 9 -xls -xlsfile temp/GPPM_mut_temp1_out_$i.xls > /dev/null &
