##########################################################################
# File Name:        demo.sh
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Thu 14 Apr 2022 04:30:55 PM CST
##########################################################################

## this is a demostate about how to use CTA select pipline##
/mnt/data2/wuzengding/03.biotools/miniconda3/envs/snakemake/bin/snakemake -s /mnt/data2/wuzengding/docker_data/14.Xena_UCSC/TAA_Selection/00.CTAselect.pipline.py --cores 16 --latency-wait 30 --config Cancer=BLCA OutDir=/mnt/data2/wuzengding/docker_data/14.Xena_UCSC/temp -p 
