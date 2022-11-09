# TAA_Selection overview
   This is a tool kits that pick up the Cancer-sepecific expession gene through analysis of RNA expression profile data from public websites such as TCGA and GTEx. Those selected gene, usually express the proteins which called Tumor-associated angtigens (TAA). 
   The fisrt clone of  human tumor antigen ,melanoma antigen-1 (MAGE-1) was reported in 1990s. Further studies found that MAGE-1 was expressed in various cancers of different histological origin but not in normal tissues excluding Testis and placenta.So those TAA like MAGE-1 also called Cancer-Testis antigen (CTA)
   Those TAA/CTA were regarded as ideal targets in cancer vaccine therapy, for example  BioNtech's BNT111 of FixVac products include the CTA of MAGE-3.
    
# Install
`git clone git@github.com:wuzengding/TAA_Selection.git`
# Utilise
## Config setting
### Making config files
1.download file of TcgaTargetGtex_RSEM_Hugo_norm_count
`
https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_RSEM_Hugo_norm_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
`

2.download file of TcgaTargetGtex_rsem_gene_tpm
`
https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
`

3.download file of TcgaTargetGTEX_phenotype.txt
`
https://xenabrowser.net/datapages/?dataset=TcgaTargetGTEX_phenotype.txt&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
`

4.make file of ensembleID2genename2biotype.txt
 `grep -v "^#" gencode.v32.chr_patch_hapl_scaff.annotation.gtf |cut -f9 |grep -v "transcript_id"|grep 'gene_name'|sed 's:;:\t:g'|awk  '{print $2"\t"$6"\t"$4}'|sed 's:"::g'|awk -F'\t' '{a[$2]=$0; b[$2]++;} END{for(i in a){if(b[i]==1){print a[i]}}}' >ensembleID2genename2biotype.txt `

5.put all config files to  ConfigDir
 ```
 mkdir -p /path/to/dataset
 cp TcgaTargetGtex_RSEM_Hugo_norm_count   /path/to/dataset
 cp TcgaTargetGtex_rsem_gene_tpm  /path/to/dataset
 cp TcgaTargetGTEX_phenotype.txt  /path/to/dataset
 cp ensembleID2genename2biotype.txt  /path/to/dataset
 ``` 
### Modify config file path
```
cd /path/to/TAA_Selection/
#vim Config.pipline.json
{                                                                                                                                  
  "MatrixMakeSh": "/path/to/TAA_Selection/01.matrix_make.sh",
  "DiffExpGeneSelR":  "/path/to/TAA_Selection/02.different_expression_select.r",   
  "NegSelR": "/path/to/TAA_Selection/03.negative_select_according_tissue.r",
  "ConfigDir": "/path/to/dataset",
  "Cancer": "LAML",
  "BioType": "protein_coding",
  "OutDir": "/output/path"
  }
```
 ## Running
 
```
/path/to/miniconda3/envs/snakemake/bin/snakemake -s /path/to/TAA_Selection/00.CTAselect.pipline.py --cores 16 --latency-wait 30 --config Cancer=BLCA OutDir=/out/path -p
```

## Results
```
Colon_Adenocarcinoma_select.gene.csv
Colon_Adenocarcinoma_vs_Paired_rsem_gene_tpm
Colon_Adenocarcinoma_vs_Paired_RSEM_Hugo_norm_count
Colon_Adenocarcinoma_vs_Tissue_rsem_gene_tpm
DEGselect.log
DESeq2.DEG.protein_coding.csv
DESeq2.DEG.protein_coding.TPM.Heatmap.gene.of.top100.pdf
DESeq2.DEG.protein_coding.TPM.Heatmap.gene.of.top200.pdf
DESeq2.DEG.protein_coding.Volcano.pdf
matrix_make.log
negative_select.log
Expression_level_of_ACAN_in_different_tissues.png
Expression_level_of_ACTL8_in_different_tissues.png
Expression_level_of_ADAM12_in_different_tissues.png
Expression_level_of_ALDH3B2_in_different_tissues.png
Expression_level_of_AQP5_in_different_tissues.png
Expression_level_of_ASCL2_in_different_tissues.png
Expression_level_of_BLACAT1_in_different_tissues.png
Expression_level_of_C17orf78_in_different_tissues.png
Expression_level_of_C1orf105_in_different_tissues.png
Expression_level_of_C2CD4A_in_different_tissues.png
Expression_level_of_C6orf15_in_different_tissues.png
Expression_level_of_CA9_in_different_tissues.png
Expression_level_of_CCL20_in_different_tissues.png
Expression_level_of_CCL24_in_different_tissues.png
Expression_level_of_CCL25_in_different_tissues.png
Expression_level_of_CCR8_in_different_tissues.png
.....
```