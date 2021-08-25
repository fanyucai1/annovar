# 1. [cosmic数据库更新](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/)

    prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.normal.vcf > hg38_cosmic92_coding.txt

    prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.normal.vcf > hg38_cosmic92_noncoding.txt

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
    
流程说明：

    step1   输入文件标准化，输出<sample name>.format.vcf
    step2:  使用annovar对标准化后对VCF文件进行注释

# 5. 备注

增删数据库对应会增加或减少条目，需要对应修改脚本中变量 **out_name**    与   **par**



