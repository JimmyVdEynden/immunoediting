#!/bin/bash

# Get all bamfiles
bamfilenames=` ls /downloads/firecloud/bam/*.bam`

# Randomly shuffle so script could be ran in parallel multiple times with minimal risk of doing the same file
bamfilenames=`shuf -e ${bamfilenames}`

# For loop fo all existing bam files
for bamfile in $bamfilenames
do
	# Check whether HLA typing already done
	barcode_temp=${bamfile##.*/} # ## longest match at front of string left out
	barcode_temp=${barcode_temp%_*m} # % shortest match at back of string left out
	TCGA_HLA_type_filename=temp/polysolver/TCGA_HLA_types/${barcode_temp}_HLA.txt
	output_temp=temp/polysolver/output/${barcode_temp}
    if [ -d $output_temp ] # Only continue of folder (-d; -e for files)  not processed yet!
    then
	    echo "Skipping $bamfile"
    	continue
    fi 

	# Run polysolver
	echo "running $bamfile"
	mkdir $output_temp
	source ~/tools/polysolver/scripts/config.sh
	nohup ~/tools/polysolver/scripts/shell_call_hla_type $bamfile Unknown 0 hg19 STDFQ 0 $output_temp > ${output_temp}/polysolver_log.txt

	# Save data
	mv ${output_temp}/winners.hla.nofreq.txt $TCGA_HLA_type_filename	

done


