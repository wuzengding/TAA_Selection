##########################################################################
# File Name:        00.CTAselect.pipline.py
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Thu 14 Apr 2022 09:24:19 AM CST
##########################################################################


configfile: os.path.join("/mnt/data2/wuzengding/05.pipeline_dev/TAA_Selection_human/TAA_Selection","Config.pipline.json")

ConfigDir=config["ConfigDir"]
Cancer=config["Cancer"]
BioType=config["BioType"]
MatrixMakeSh=config["MatrixMakeSh"]
DiffExpGeneSelR=config["DiffExpGeneSelR"]
NegSelR=config["NegSelR"]

try:
    OutDir=os.path.join(config["OutDir"],"DEGselect_"+Cancer)
except:
    config["OutDir"]=os.getcwd()
    OutDir=config["OutDir"]
    
print("DiffExpGeneSelR",DiffExpGeneSelR)
print("NegSelR",NegSelR)
CancerTypeDict={"LAML":"Acute_Myeloid_Leukemia",
	"ACC":"Adrenocortical_Cancer",
	"BLCA":"Bladder_Urothelial_Carcinoma",
	"LGG":"Brain_Lower_Grade_Glioma",
	"BRCA":"Breast_Invasive_Carcinoma",
	"CESC":"Cervical_&_Endocervical_Cancer",
	"CHOL":"Cholangiocarcinoma",
	"LCML":"NA",
	"COAD":"Colon_Adenocarcinoma",
	"ESCA":"Esophageal_Carcinoma",
	"GBM":"Glioblastoma_Multiforme",
	"HNSC":"Head_&_Neck_Squamous_Cell_Carcinoma",
	"KICH":"Kidney_Chromophobe",
	"KIRC":"Kidney_Clear_Cell_Carcinoma",
	"KIRP":"Kidney_Papillary_Cell_Carcinoma",
	"LIHC":"Liver_Hepatocellular_Carcinoma",
	"LUAD":"Lung_Adenocarcinoma",
	"LUSC":"Lung_Squamous_Cell_Carcinoma",
	"DLBC":"Diffuse_Large_B-Cell_Lymphoma",
	"MESO":"Mesothelioma",
	"OV":"Ovarian_Serous_Cystadenocarcinoma",
	"PAAD":"Pancreatic_Adenocarcinoma",
	"PCPG":"Pheochromocytoma_&_Paraganglioma",
	"PRAD":"Prostate_Adenocarcinoma",
	"READ":"Rectum_Adenocarcinoma",
	"SARC":"Sarcoma",
	"SKCM":"Skin_Cutaneous_Melanoma",
	"STAD":"Stomach_Adenocarcinoma",
	"TGCT":"Testicular_Germ_Cell_Tumor",
	"THYM":"Thymoma",
	"THCA":"Thyroid_Carcinoma",
	"UCS":"Uterine_Carcinosarcoma",
	"UCEC":"Uterine_Corpus_Endometrioid_Carcinoma",
	"UVM":"Uveal_Melanoma"}
CancerType = CancerTypeDict[Cancer]
print("ConfigDir",ConfigDir)
print("OutDir",OutDir)

rule all:
    input:
        OutDir+"/"+CancerTypeDict[Cancer]+"_select.gene.csv"

rule matrix_make:
    output:
        PairedTpm = OutDir + "/"+CancerTypeDict[Cancer]+"_vs_Paired_rsem_gene_tpm",
        PairedCount = OutDir + "/"+CancerTypeDict[Cancer]+"_vs_Paired_RSEM_Hugo_norm_count",
        TissueTpm = OutDir+"/"+CancerTypeDict[Cancer]+"_vs_Tissue_rsem_gene_tpm"
    shell:
        "bash {MatrixMakeSh} {ConfigDir} {OutDir} {Cancer} >{OutDir}/matrix_make.log 2>&1"

rule DEGselect:
    input:
        PairedTpm = OutDir+"/"+CancerTypeDict[Cancer]+"_vs_Paired_rsem_gene_tpm",
        PairedCount = OutDir+"/"+CancerTypeDict[Cancer]+"_vs_Paired_RSEM_Hugo_norm_count",
        EnsGeneBio = ConfigDir+"/ensembleID2genename2biotype.txt"
    output:
        Top100Pdf = OutDir+"/DESeq2.DEG."+BioType+".TPM.Heatmap.gene.of.top100.pdf",
        DegF = OutDir+"/DESeq2.DEG."+BioType+".csv"
    shell:
        "/usr/bin/Rscript {DiffExpGeneSelR} {BioType} {input.PairedCount} {input.PairedTpm} {input.EnsGeneBio} {OutDir} > {OutDir}/DEGselect.log 2>&1"

rule negative_select:
    input:
        DegF = os.path.join(OutDir, "DESeq2.DEG."+BioType+".csv"),
        TissueTpm = os.path.join(OutDir, CancerTypeDict[Cancer]+"_vs_Tissue_rsem_gene_tpm"),
        EnsGeneBio = os.path.join(ConfigDir, "ensembleID2genename2biotype.txt"),
        Phenotype = os.path.join(ConfigDir, "TcgaTargetGTEX_phenotype.txt"),
        
    output:
        OutDir+"/"+CancerTypeDict[Cancer]+"_select.gene.csv"
    shell:
        "/usr/bin/Rscript {NegSelR} {input.DegF} {input.TissueTpm} {input.EnsGeneBio} {input.Phenotype} {OutDir} {CancerType} > {OutDir}/negative_select.log 2>&1"