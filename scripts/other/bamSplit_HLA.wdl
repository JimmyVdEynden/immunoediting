workflow bamSplit_HLA {
	call splitBam_HLA
}

task splitBam_HLA {
	String sampleName
	File inputBAM
	File inputBAI
	command {
		# Create sam for selected regions
		samtools view -H ${inputBAM} > temp_chr6Region.sam
		samtools view ${inputBAM} "6:29909037-29913661" >> temp_chr6Region.sam
		samtools view ${inputBAM} "6:31321649-31324964" >> temp_chr6Region.sam
		samtools view ${inputBAM} "6:31236526-31239869" >> temp_chr6Region.sam
		# Convert to bam
		samtools view -bS temp_chr6Region.sam > temp_chr6Region.bam
		# Sort bam
		samtools sort temp_chr6Region.bam > ${sampleName}_chr6Region.bam
		# Remove temp
		rm temp_*
		# Index
		samtools index ${sampleName}_chr6Region.bam

	}
	output {
		File splittedBAM = "${sampleName}_chr6Region.bam"
		File splittedBAM_IDX = "${sampleName}_chr6Region.bam.bai"		
	}
	runtime {
		docker: "biocontainers/samtools"
		disks: "local-disk 75 SSD"
	}
}