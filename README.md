# 1. [cosmic数据库更新](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/)

    prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.normal.vcf > hg19_cosmic94_coding.txt

    prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.normal.vcf > hg19_cosmic94_noncoding.txt

# 2. 修改gnomad211_genome与gnomad211_exome的标签

    gnomad211_genome(genome_AF)与gnomad211_exome(exome_AF)两个数据库的标签注意进行修改默认标签都是AF

# 3. 关于预测非同义突变打分：dbNSFP Information

    #SIFT_score                      ≤ 0.05 was regarded as deleterious (D), and a score value >0.05 was regarded as tolerated (T)
    #Polyphen2_HDIV_score            ≥ 0.957, probably damaging (D); 0.453 < Polyphen2 HDIV score < 0.956, possibly damaging (P); Polyphen2 HDIV score ≤ 0.452, benign (B)
    #Polyphen2_HVAR_score            ≥ 0.909, probably damaging (D); 0.447 < Polyphen2 HVAR score < 0.909, possibly damaging (P); Polyphen2 HVAR score ≤ 0.446, benign (B)
    #LRT_score                       D, deleterious; N, neutral
    #MutationTaster_score            A, disease causing automatic; D, disease causing
    #FATHMM_score                    D: Deleterious; T: Tolerated
    #CADD_raw                        higher scores are more deleterious:20 means 1% percentile highest scores, and 30% means 0.1% percentile highest scores.

# 4. annovar.py

    usage: [-h] -v VCF -n SAMPLE_NAME -o OUTDIR
    
    optional arguments:
      -h, --help            show this help message and exit
      -v VCF, --vcf VCF     format vcf file                 #输入要注释的VCF文件
      -n SAMPLE_NAME, --sample_name SAMPLE_NAME             #往往一个VCF文件中对应多个样本，需要输入VCF文件中对应的样本名
                            sample name in vcf              
      -o OUTDIR, --outdir OUTDIR                            #输出文件夹
                            output directory
      -t TRANSCRIPT, --transcript TRANSCRIPT                #经典转录本
                        Canonical transcript file
      -r REF, --ref REF     annovar directory               #annovar软件与数据库文件夹

流程说明：

1:  输入文件标准化，输出<sample name>.format.vcf

2:  使用annovar对标准化后对VCF文件进行注释

3:  增删数据库对应会增加或减少条目，需要对应修改脚本中变量: <out_name>    与   <par>

4:  annovar软件与数据库文件夹目录结构

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

# 5. 输出示例文件

[test.annovar.tsv](./test.annovar.tsv)