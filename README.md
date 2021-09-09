# 1. [cosmic数据库更新](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/)

    prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.normal.vcf > hg19_cosmic94_coding.txt

    prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.normal.vcf > hg19_cosmic94_noncoding.txt

# 2. GnomAD

Download data:
    axel -n20 http://www.openbioinformatics.org/annovar/download/hg19_gnomad211_exome.txt.gz
    axel -n20 http://www.openbioinformatics.org/annovar/download/hg19_gnomad211_genome.txt.gz

Same title information from gnomad211_exome and gnomad211_genome as following:
    
    AF AF_popmax AF_male AF_female AF_raw AF_afr AF_sas AF_amr AF_eas AF_nfe AF_fin AF_asj AF_oth non_topmed_AF_popmax non_neuro_AF_popmax non_cancer_AF_popmax controls_AF_popmax
    
    change the gnomad211_exome table_tile (e.g: from AF_male to exome_AF_male)
    change the gnomad211_genome table_tile (e.g: from AF_male to genome_AF_male)

# 3. 关于预测非同义突变打分：dbNSFP Information

    #SIFT_score                      ≤ 0.05 was regarded as deleterious (D), and a score value >0.05 was regarded as tolerated (T)
    #Polyphen2_HDIV_score            ≥ 0.957, probably damaging (D); 0.453 < Polyphen2 HDIV score < 0.956, possibly damaging (P); Polyphen2 HDIV score ≤ 0.452, benign (B)
    #Polyphen2_HVAR_score            ≥ 0.909, probably damaging (D); 0.447 < Polyphen2 HVAR score < 0.909, possibly damaging (P); Polyphen2 HVAR score ≤ 0.446, benign (B)
    #LRT_score                       D, deleterious; N, neutral
    #MutationTaster_score            A, disease causing automatic; D, disease causing
    #FATHMM_score                    D: Deleterious; T: Tolerated
    #CADD_raw                        higher scores are more deleterious:20 means 1% percentile highest scores, and 30% means 0.1% percentile highest scores.

# 4. annovar.py

    usage: [-h] -v VCF -n SAMPLE_NAME -o OUTDIR -t TRANSCRIPT -r REF

    optional arguments:
      -h, --help            show this help message and exit
      -v VCF, --vcf VCF     format vcf file
      -n SAMPLE_NAME, --sample_name SAMPLE_NAME
                            sample name in vcf
      -o OUTDIR, --outdir OUTDIR
                            output directory
      -t TRANSCRIPT, --transcript TRANSCRIPT
                            Canonical transcript file
      -r REF, --ref REF     annovar directory


流程说明：

+ 输入文件VCF标准化:\<sample name\>.format.vcf
+ annovar注释:\<sample name\>.raw.annovar.tsv
+ annovar结果标准化：\<sample name\>.final.annovar.tsv

# 5.附录

+ 增删数据库对应会增加或减少条目，需要对应修改脚本中变量:**out_name** 、 **par** 
+ annovar软件与数据库文件夹目录结构


    ./
    ├── annotate_variation.pl
    ├── canonical_transcript.txt
    ├── coding_change.pl
    ├── convert2annovar.pl
    ├── example
    │     ├── ex1.avinput
    │     ├── ex2.vcf
    │     ├── example.simple_region
    │     ├── example.tab_region
    │     ├── gene_fullxref.txt
    │     ├── gene_xref.txt
    │     ├── grantham.matrix
    │     ├── README
    │     └── snplist.txt
    ├── humandb
    │     ├── genometrax-sample-files-gff
    │     ├── GRCh37_MT_ensGeneMrna.fa
    │     ├── GRCh37_MT_ensGene.txt
    │     ├── hg19_AFR.sites.2015_08.txt
    │     ├── hg19_AFR.sites.2015_08.txt.idx
    │     ├── hg19_ALL.sites.2015_08.txt
    │     ├── hg19_ALL.sites.2015_08.txt.idx
    │     ├── hg19_AMR.sites.2015_08.txt
    │     ├── hg19_AMR.sites.2015_08.txt.idx
    │     ├── hg19_avsnp150.txt
    │     ├── hg19_avsnp150.txt.idx
    │     ├── hg19_clinvar_20210501.txt
    │     ├── hg19_clinvar_20210501.txt.idx
    │     ├── hg19_cosmic94_coding.txt
    │     ├── hg19_cosmic94_noncoding.txt
    │     ├── hg19_cytoBand.txt
    │     ├── hg19_dbnsfp42a.txt
    │     ├── hg19_dbnsfp42a.txt.idx
    │     ├── hg19_EAS.sites.2015_08.txt
    │     ├── hg19_EAS.sites.2015_08.txt.idx
    │     ├── hg19_esp6500siv2_all.txt
    │     ├── hg19_esp6500siv2_all.txt.idx
    │     ├── hg19_EUR.sites.2015_08.txt
    │     ├── hg19_EUR.sites.2015_08.txt.idx
    │     ├── hg19_exac03.txt
    │     ├── hg19_exac03.txt.idx
    │     ├── hg19_example_db_generic.txt
    │     ├── hg19_example_db_gff3.txt
    │     ├── hg19_gnomad211_exome.txt
    │     ├── hg19_gnomad211_exome.txt.idx
    │     ├── hg19_gnomad211_genome.txt
    │     ├── hg19_gnomad211_genome.txt.idx
    │     ├── hg19_icgc28.txt
    │     ├── hg19_icgc28.txt.idx
    │     ├── hg19_intervar_20180118.txt
    │     ├── hg19_intervar_20180118.txt.idx
    │     ├── hg19_MT_ensGeneMrna.fa
    │     ├── hg19_MT_ensGene.txt
    │     ├── hg19_refGeneMrna.fa
    │     ├── hg19_refGene.txt
    │     ├── hg19_refGeneVersion.txt
    │     ├── hg19_refGeneWithVerMrna.fa
    │     ├── hg19_refGeneWithVer.txt
    │     ├── hg19_SAS.sites.2015_08.txt
    │     ├── hg19_SAS.sites.2015_08.txt.idx
    │     ├── hg19_snp138.txt
    │     └── hg19_snp138.txt.idx
    ├── index_annovar.pl
    ├── prepare_annovar_user.pl
    ├── retrieve_seq_from_fasta.pl
    ├── table_annovar.pl
    └── variants_reduction.pl

#6. 备注
