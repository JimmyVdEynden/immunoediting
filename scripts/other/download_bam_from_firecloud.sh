##!/bin/bash

export CLOUDSDK_PYTHON=/usr/local/bin/python

# Get bucket ID from firecloud
bucket_ID="fc-6527f237-73f1-4218-b934-ffe264b9d326"

# Check which analysis are new 
fc_all=`~/tools/google-cloud-sdk/bin/gsutil ls gs://$bucket_ID`
rudy_exists=`ls $bucket_ID`
echo $fc_all | tr " " "\n" | cut -f4 -d "/" > temp.txt
echo $rudy_exists | tr " " "\n" > temp2.txt
# printf '%s\n' "${rudy_exists[@]}" | LC_ALL=C sort > temp2.txt
analyses_new=`comm -23 temp.txt temp2.txt`

# Download new analysis
analyses_length=`echo $analyses_new | wc -w`
for i in `seq $analyses_length`
do 
	echo $i
	analysis_new=`echo $analyses_new | cut -f$i -d " "`
	~/tools/google-cloud-sdk/bin/gsutil cp -r gs://$bucket_ID/$analysis_new ~/downloads/firecloud/$bucket_ID
done

# Copy bam files to project directory (note that hard linking gives "argument too long" issues)
cp -n ~/downloads/firecloud/$bucket_ID/*/bamSplit_HLA/*/call-splitBam_HLA/*_chr6Region.bam* ~/Projects/pub/immunoediting/downloads/firecloud/bam/
# -n prevents overwriting

# remove temp files
rm temp*.txt

# To empty bucket
# ./google-cloud-sdk/bin/gsutil rm -r gs://$bucket_ID/*
