for i in {1..20} 
do Rscript scripts/fuse_GPPM.R $i &
done;
wait # To avoid running rest already

Rscript scripts/fuse_GPPM2.R