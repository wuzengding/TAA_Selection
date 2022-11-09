#!/usr/bin/bash

samplefile=$1
if [ ! "$samplefile" ];then
	echo "$0 <samplelistfile>[sametime_oncomput_sample_count]"
fi
#samplelistfile
#SRR20852112	/path/to/SRR20852112_pass_1.fastq.gz	/path/to/SRR20852112_pass_2.fastq.gz

outdir=$2

#同时投递的任务数，默认值为2
same_online_sample_count=$3
if [ ! "$same_online_sample_count" ];then
	same_online_sample_count=2
fi

data_count_supply=0
cat $samplefile|while read line
do
	col_num=`echo "$line"|awk -F "\t" '{print NF}'`
	samplename=`echo  "$line"|awk -F "\t" '{print $1}'`
	raw_fq1=`echo  "$line"|awk -F "\t" '{print $2}'`
	echo ${samplename}
	echo ${raw_fq1}
	/usr/bin/bash /mnt/data2/wuzengding/02.ResDev/05.MC38_TAA_Epitope/Fastq2TAA.sh $raw_fq1 ${outdir}/$samplename > ${outdir}/${samplename}.nohupout &
#	bams="$bam result/$samplename/$${samplename}.Aligned.sortedByCoord.out.bam"
	echo $data_count_supply
	let data_count_supply=$data_count_supply+1
	echo $data_count_supply
	if [[ $data_count_supply -ge $same_online_sample_count ]];then
		data_count_supply=0
		wait;
	fi
done < $samplefile
