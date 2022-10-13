##########################################################################
# File Name:        01.matrix_make.sh
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Thu 07 Apr 2022 02:00:46 PM CST
##########################################################################

#!/bin/bash
# The script that can transfer original data into data matrix which can be used for Different Expression Gene analysis according user's specified requirement such as Cancer_Type,Cohort_Type and so on.

Indir=$1
Outdir=$2
Cancer=$3


declare -A CancerType
declare -A CancerNormal
declare -A CutRegion

CancerType=(
	["LAML"]="Acute Myeloid Leukemia"
	["ACC"]="Adrenocortical Cancer"
	["BLCA"]="Bladder Urothelial Carcinoma"
	["LGG"]="Brain Lower Grade Glioma"
	["BRCA"]="Breast Invasive Carcinoma"
	["CESC"]="Cervical & Endocervical Cancer"
	["CHOL"]="Cholangiocarcinoma"
	["LCML"]="NA"
	["COAD"]="Colon Adenocarcinoma"
	["ESCA"]="Esophageal Carcinoma"
	["GBM"]="Glioblastoma Multiforme"
	["HNSC"]="Head & Neck Squamous Cell Carcinoma"
	["KICH"]="Kidney Chromophobe"
	["KIRC"]="Kidney Clear Cell Carcinoma"
	["KIRP"]="Kidney Papillary Cell Carcinoma"
	["LIHC"]="Liver Hepatocellular Carcinoma"
	["LUAD"]="Lung Adenocarcinoma"
	["LUSC"]="Lung Squamous Cell Carcinoma"
	["DLBC"]="Diffuse Large B-Cell Lymphoma"
	["MESO"]="Mesothelioma"
	["OV"]="Ovarian Serous Cystadenocarcinoma"
	["PAAD"]="Pancreatic Adenocarcinoma"
	["PCPG"]="Pheochromocytoma & Paraganglioma"
	["PRAD"]="Prostate Adenocarcinoma"
	["READ"]="Rectum Adenocarcinoma"
	["SARC"]="Sarcoma"
	["SKCM"]="Skin Cutaneous Melanoma"
	["STAD"]="Stomach Adenocarcinoma"
	["TGCT"]="Testicular Germ Cell Tumor"
	["THYM"]="Thymoma"
	["THCA"]="Thyroid Carcinoma"
	["UCS"]="Uterine Carcinosarcoma"
	["UCEC"]="Uterine Corpus Endometrioid Carcinoma"
	["UVM"]="Uveal Melanoma")
CancerNormal=(
	["LAML"]="Whole Blood"
	["ACC"]="Adrenal Gland"
	["BLCA"]="Bladder"
	["LGG"]="Brain"
	["BRCA"]="Breast"
	["CESC"]="Cervix Uteri"
	["CHOL"]="NA"
	["LCML"]="NA"
	["COAD"]="Colon"
	["ESCA"]="Esophagus"
	["GBM"]="Brain"
	["HNSC"]="NA"
	["KICH"]="Kidney"
	["KIRC"]="Kidney"
	["KIRP"]="Kidney"
	["LIHC"]="Liver"
	["LUAD"]="Lung"
	["LUSC"]="Lung"
	["DLBC"]="Whole Blood"
	["MESO"]="NA"
	["OV"]="Ovary"
	["PAAD"]="Pancreas"
	["PCPG"]="NA"
	["PRAD"]="Prostate"
	["READ"]="Colon"
	["SARC"]="NA"
	["SKCM"]="Skin"
	["STAD"]="Stomach"
	["TGCT"]="Testis"
	["THYM"]="NA"
	["THCA"]="Thyroid"
	["UCS"]="Uterus"
	["UCEC"]="Uterus"
	["UVM"]="NA")
CutRegion=(
	["LAML"]="3"
	["ACC"]="3"
	["BLCA"]="3"
	["LGG"]="4"
	["BRCA"]="4"
	["CESC"]="4"
	["CHOL"]="NA"
	["LCML"]="NA"
	["COAD"]="4"
	["ESCA"]="4"
	["GBM"]="4"
	["HNSC"]="NA"
	["KICH"]="4"
	["KIRC"]="4"
	["KIRP"]="4"
	["LIHC"]="4"
	["LUAD"]="4"
	["LUSC"]="4"
	["DLBC"]="3"
	["MESO"]="NA"
	["OV"]="4"
	["PAAD"]="4"
	["PCPG"]="NA"
	["PRAD"]="4"
	["READ"]="4"
	["SARC"]="NA"
	["SKCM"]="4"
	["STAD"]="4"
	["TGCT"]="4"
	["THYM"]="NA"
	["THCA"]="4"
	["UCS"]="4"
	["UCEC"]="4"
	["UVM"]="NA")

#echo ${CancerType[$Cancer]}
#echo ${CancerNormal[$Cancer]}
#echo ${CutRegion[$Cancer]}

if [[ "${Cancer}" =~ ^(LAML|ACC|BLCA|LGG|BRCA|CESC|CHOL|LCML|COAD|ESCA|GBM|HNSC|KICH|KIRC|KIRP|LIHC|LUAD|LUSC|DLBC|MESO|OV|PAAD|PCPG|PRAD|READ|SARC|SKCM|STAD|TGCT|THYM|THCA|UCS|UCEC|UVM)$ ]];then
	echo "Cancer type in TCGA list"
else
	echo "Input Cancer type not in TCGA list, please input again,"
	echo "TCGA Cancer type: LAML|ACC|BLCA|LGG|BRCA|CESC|CHOL|LCML|COAD|ESCA|GBM|HNSC|KICH|KIRC|KIRP|LIHC|LUAD|LUSC|DLBC|MESO|OV|PAAD|PCPG|PRAD|READ|SARC|SKCM|STAD|TGCT|THYM|THCA|UCS|UCEC|UVM"
	exit 1
fi

CancerName=${CancerType[${Cancer}]}
NormalName=${CancerNormal[${Cancer}]}
Cutfield=${CutRegion[${Cancer}]}


#if [[ $CancerName == "Acute Myeloid Leukemia" ]];then
#	echo "yes"
#fi

if [[ "$CancerName" == "LCML" ]];then 
	echo "LCML Cancer do not contain RNAseq sample data"
	exit 1
else
	if [[ "$NormalName" == "NA" ]];then
		echo "you select ${Cancer} cancer,cannot find paired normal tissue samples"
		exit 1
	else
		echo "Starting select samples from TCGA and GTEX! and then be made to data matrix "
		
		phenotype="${Indir}/TcgaTargetGTEX_phenotype.txt"
		tpm="${Indir}/TcgaTargetGtex_rsem_gene_tpm"
		count="${Indir}/TcgaTargetGtex_RSEM_Hugo_norm_count"
		
		cancersample=$(cat $phenotype|grep "TCGA"|grep -v "Solid Tissue Normal"|awk -F'\t' -v CancerName="$CancerName" '{if ($3 == CancerName) {print $0}}'|cut -f1|tr '\n' ','|sed 's/\(.*\),/\1/')
		pairedsample=$(cat $phenotype|grep "TCGA"|grep "Solid Tissue Normal"|awk -F'\t' -v CancerName="$CancerName" '{if ($3 == CancerName) {print $0}}'|cut -f1|tr '\n' ','|sed 's/\(.*\),/\1/')
		normalsample=$(cat $phenotype|grep "GTEX"|awk -F'\t' -v NormalName="$NormalName" -v Cutfield="$Cutfield" '{if ($Cutfield == NormalName) {print $0}}'|cut -f1|tr '\n' ','|sed 's/\(.*\),/\1/')
		tissuesample=$(cat $phenotype|grep "GTEX"|grep -v "Cell Line"|cut -f1|tr '\n' ','|sed 's/\(.*\),/\1/' )
		cancersamplenum=$(echo ${cancersample}|tr ',' '\n'|wc -l)
		pairedsamplenum=$(echo ${pairedsample}|tr ',' '\n'|wc -l)
		normalsamplenum=$(echo ${normalsample}|tr ',' '\n'|wc -l)
		tissuesamplenum=$(echo ${tissuesample}|tr ',' '\n'|wc -l)
		echo "you select ${Cancer} cancer,whose full name in TCGA is '${CancerName}' with ${cancersamplenum} samples, 
		       and paired Cancerside ${pairedsamplenum}sample,
		       and paired Normal type is '${NormalName}' with ${normalsamplenum} samples,
			   and tissue Normal type with ${tissuesamplenum} samples"
		samplelist="sample,"${cancersample},${normalsample}
		cancertissuelist="sample,"${cancersample},${tissuesample}
		tpmsamplelist=$(head -1 ${tpm})
		position=0
		fieldlist=""
		for sampleid in ${tpmsamplelist};
		do
			let "position++"
			##let position=$((position+1))
			if [[ $(echo "$cancertissuelist" | grep "$sampleid") ]];then
				fieldlist="${fieldlist},$position"
				echo "$position $sampleid in list"
			fi
		done
		fieldlist=$(echo $fieldlist|sed 's:,::')
		echo $fieldlist
		#echo ${cancertissuelist} > "${Outdir}/test0412.txt"
		#echo ${cancertissuelist}
		
		cat ${count}|tr '\t' ','|csvcut -c ${samplelist}|tr ',' '\t' > "${Outdir}/$(echo ${CancerName}|tr ' ' '_')_vs_Paired_RSEM_Hugo_norm_count"
		cat ${tpm}|tr '\t' ','|csvcut -c ${samplelist}|tr ',' '\t' > "${Outdir}/$(echo ${CancerName}|tr ' ' '_')_vs_Paired_rsem_gene_tpm"
		#cat ${tpm}|tr '\t' ','|csvcut -c ${cancertissuelist}|tr ',' '\t' > "${Outdir}/$(echo ${CancerName}|tr ' ' '_')_vs_Tissue_rsem_gene_tpm"
		cat ${tpm}|cut -f ${fieldlist} > "${Outdir}/$(echo ${CancerName}|tr ' ' '_')_vs_Tissue_rsem_gene_tpm"
		echo "Counts and TPM data matrix making were finished!"
	fi
fi
