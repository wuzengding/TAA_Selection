##########################################################################
# File Name:        Fastq2TAA.sh
# Author:           wuzengding
# mail:             wuzengding@163.com                                                         
# Created Time:     Tue. 26 Oct 2022 09:28:15 AM CST
##########################################################################

raw_fq1=$1
outpath=$2

###################################################################
#########         config                              #############
###################################################################

#!/usr/bin/bash
fastp='/mnt/data2/wuzengding/03.biotools/software/fastp'
STAR='/mnt/data2/wuzengding/03.biotools/miniconda3/bin/STAR'
salmon='/mnt/data2/wuzengding/03.biotools/software/salmon-1.6.0_linux_x86_64/bin/salmon'
genome='/mnt/data2/wuzengding/00.database/15.index_mouse_genome/index_GRCm39_ensemble107.star'
transcriptfa='/mnt/data2/wuzengding/00.database/16.annot_mouse/Mus_musculus.GRCm39.transcript.generated.fa'
gtf="/mnt/data2/wuzengding/00.database/16.annot_mouse/Mus_musculus.GRCm39.107.gtf"
#transcriptfa="/mnt/data2/cuifenglei/project/bulkRNAseq/plasimid/egfp/egfp.transcript.fasta"
#transcriptfa="/mnt/data2/cuifenglei/project/bulkRNAseq/plasimid/genome_egfp_luc2_joinVersion/gffread_creat_transcripts.fa"



###################################################################
#########         RUNNING                              #############
###################################################################
echo ${1}
samplename=$(basename $1|sed 's:\.:_:g'|cut -f1 -d'_')
samplepath=$(dirname $1)
pos2=$(basename $1|sed 's:\.:_:g'|cut -f2 -d'_')
pos3=$(basename $1|sed 's:\.:_:g'|cut -f3 -d'_')

#echo $samplename
#echo ${pos2}
#echo ${pos3}
#echo ${samplepath}


if [ ! -d $outpath ];then
	mkdir -p $outpath
fi


if  false;then
    echo $"nothing"
fi
#########         Step1:Trim                            #############

if [ ! -d $outpath/01.Trim ];then
	mkdir -p $outpath/01.Trim
fi

clean_1_fq="$outpath/01.Trim/${samplename}_clean_1.fq"
clean_2_fq="$outpath/01.Trim/${samplename}_clean_2.fq"

if [ "${pos2}" = "1" ] && [ "${pos3}" = "fastq" ] ;then
	raw_fq2=${samplepath}/${samplename}_2.fastq.gz
	echo $raw_fq2
	$fastp -i ${raw_fq1} -o $clean_1_fq -I ${raw_fq2} -O $clean_2_fq --json $outpath/01.Trim/${samplename}.json  --html $outpath/01.Trim/${samplename}.html >$outpath/01.Trim/${samplename}.fastp.log 2>&1
elif [ $# -eq 2 ];then	
	$fastp -i ${raw_fq1} -o $clean_1_fq --json $outpath/01.Trim/${samplename}.json --html $outpath/01.Trim/${samplename}.html >$outpath/01.Trim/${samplename}.fastp.log 2>&1
	clean_2_fq=""
fi


#########         Step2:Align                            #############

if [ ! -d $outpath/02.Align ];then
	mkdir -p $outpath/02.Align
fi
$STAR --runThreadN 20 --readFilesIn $clean_1_fq $clean_2_fq --genomeDir $genome --quantTranscriptomeBan Singleend --quantMode TranscriptomeSAM --outFileNamePrefix $outpath/02.Align/${samplename}. -outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --sjdbGTFfile $gtf --outReadsUnmapped Fastx --outSAMattributes Standard >$outpath/02.Align/${samplename}.STAR.log 2>&1



#########         Step3:Quant                            #############

if [ ! -d $outpath/03.Quant ];then
	mkdir -p $outpath/03.Quant
fi

$salmon quant -l A -a $outpath/02.Align/${samplename}.Aligned.toTranscriptome.out.bam -t $transcriptfa -p 10 -g $gtf  -o $outpath/03.Quant/${samplename}_salmon_quant >$outpath/03.Quant/${samplename}_salmon_quant.log 2>&1

#用于salmon
#--quantTranscriptomeBan Singleend
#--quantMode TranscriptomeSAM
