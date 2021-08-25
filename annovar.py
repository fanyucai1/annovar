#Email:fanyucai1@126.com
#2021.08.25 version:1.0

import argparse
import os
import subprocess
import re
#################################指定输出标签
out_name=['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene',
          'ExonicFunc.refGene',	'AAChange.refGene','canonical_transcript','cytoBand','avsnp150','snp138','1000g2015aug_all','1000g2015aug_eas','genome_AF','genome_AF_eas','exome_AF','exome_AF_eas','ExAC_ALL','ExAC_EAS','esp6500siv2_all',
          'CLNALLELEID','CLNDN','CLNDISDB',	'CLNREVSTAT','CLNSIG','cosmic94_coding','cosmic94_noncoding','SIFT_score','LRT_score','FATHMM_score','CADD_raw',
          'MutationTaster_score','Polyphen2_HDIV_score','Polyphen2_HVAR_score','InterVar_automated','GT',
          'Ref_Reads','Alt_Reads','Var']
# dbNSFP Information
#SIFT_score                      ≤ 0.05 was regarded as deleterious (D), and a score value >0.05 was regarded as tolerated (T)
#Polyphen2_HDIV_score            ≥ 0.957, probably damaging (D); 0.453 < Polyphen2 HDIV score < 0.956, possibly damaging (P); Polyphen2 HDIV score ≤ 0.452, benign (B)
#Polyphen2_HVAR_score            ≥ 0.909, probably damaging (D); 0.447 < Polyphen2 HVAR score < 0.909, possibly damaging (P); Polyphen2 HVAR score ≤ 0.446, benign (B)
#LRT_score                       D, deleterious; N, neutral
#MutationTaster_score            A, disease causing automatic; D, disease causing
#FATHMM_score                    D: Deleterious; T: Tolerated
#CADD_raw                        higher scores are more deleterious:20 means 1% percentile highest scores, and 30% means 0.1% percentile highest scores.

############################################################vcf annoattion use annovar
def anno(vcf,outdir,prefix,Canonical_transcript_file,annovar):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    out=outdir+"/"+prefix
    #########################get Canonical transcript info
    transcript={}
    infile=open(Canonical_transcript_file,"r")
    for line in infile:
        line=line.strip()
        transcript[line.split(",")[0]]=line.split(",")[1]
    infile.close()
    ##########################run annovar
    par = " -protocol refGene,cytoBand,snp138,avsnp150,exac03,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,gnomad211_exome,gnomad211_genome,cosmic94_coding,cosmic94_noncoding,clinvar_20210501,dbnsfp42a,intervar_20180118"
    #par += ",1000g2015aug_sas,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur "
    #par += " -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f "
    par += " -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f "
    par += " -nastring . -polish "
    subprocess.check_call("perl %s/table_annovar.pl %s %s/humandb -buildver hg19 -out %s -remove %s -vcfinput " % (annovar, vcf, annovar, out, par), shell=True)
    subprocess.check_call("rm -rf %s.hg19_multianno.vcf %s.avinput" % (out, out), shell=True)
    infile = open("%s.hg19_multianno.txt" % (out), "r")
    outfile = open("%s.annovar.tsv" % (out), "w")
    for i in range(len(out_name)):
        if i == 0:
            outfile.write("%s" % (out_name[i]))
        else:
            outfile.write("\t%s" % (out_name[i]))
    outfile.write("\n")
    dict = {}
    for line in infile:
        line = line.strip()
        array = line.split("\t")
        name = []
        if line.startswith("Chr"):
            for i in range(len(array)):
                name.append(array[i])
                dict[array[i]] = i
        else:
            p1 = re.compile(r'Alt_Reads=([0-9]+)')
            p2 = re.compile(r'Ref_Reads=([0-9]+)')
            p3 = re.compile(r'Var=(\d+.\d+)')
            p4 = re.compile(r'GT=(\d+/\d+)')
            Alt_Reads = p1.findall(line)
            Ref_Reads = p2.findall(line)
            Var = p3.findall(line)
            GT = p4.findall(line)
            ##########################format output knownCanonical transcript
            final_nm = "."
            if array[dict['AAChange.refGene']]!=".":
                tmp = array[dict['AAChange.refGene']].split(",")
                final_nm = tmp[0]
                for i in range(len(tmp)):
                    if re.search(r'%s'%(transcript[array[6]]),tmp[i]):
                        final_nm=tmp[i]
            for l in range(len(out_name)):
                if l == 0:
                    outfile.write("%s" % (array[dict[out_name[l]]]))
                elif out_name[l] == "Var":
                    outfile.write("\t%s" % (Var[0]) + "%")
                elif out_name[l] == "Alt_Reads":
                    outfile.write("\t%s" % (Alt_Reads[0]))
                elif out_name[l] == "Ref_Reads":
                    outfile.write("\t%s" % (Ref_Reads[0]))
                elif out_name[l] == "canonical_transcript":
                    outfile.write("\t%s" % (final_nm))
                elif out_name[l] == "GT":
                    outfile.write("\t%s" % (GT[0]))
                else:
                    outfile.write("\t%s" % (array[dict[out_name[l]]]))
            outfile.write("\n")
    infile.close()
    outfile.close()
    if os.path.exists("%s.hg19_multianno.txt" % (out)):
        subprocess.check_call("rm -rf %s.hg19_multianno.txt" % (out), shell=True)

###########################################################format vcf and input annovar
def format_vcf(vcf,sample_name,outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    infile = open(vcf, "r")
    outfile = open("%s/%s.format.vcf" % (outdir, sample_name), "w")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    name,num = "",0
    for line in infile:
        line = line.strip()
        if line.startswith("#CHROM"):
            array = line.split("\t")
            for i in range(len(array)):
                if array[i] == sample_name:
                    name = i
                if array[i] == "FORMAT":
                    num = i
                    continue
        if not line.startswith("#"):
            array = line.split("\t")
            tmp = array[num].split(":")
            info = array[int(name)].split(":")
            GT, a, b, c = "", array[4].split(","), [], []
            for k in range(len(tmp)):
                if tmp[k] == "GT":
                    GT = info[k]  # GT
                elif tmp[k] == "AD":
                    b = info[k].split(",")  # AD
                elif tmp[k] == "AF":
                    c = info[k].split(",")  # AF
                else:
                    pass
            Ref_Reads = b[0]
            if len(a) == 1:
                Var = float(c[0]) * 100
                outfile.write("%s\t%s\t%s\t%s\t%s\t.\t.\tGT=%s;Ref_Reads=%s;Alt_Reads=%s;Var=%.2f"
                              % (array[0], array[1], array[2], array[3], array[4], GT, Ref_Reads, b[1], Var))
                outfile.write("%\n")
            else:
                for i in range(len(a)):
                    ALT = a[i]
                    Alt_Reads = b[i + 1]
                    Var = float(c[i]) * 100
                    outfile.write("%s\t%s\t%s\t%s\t%s\t.\t.\tGT=%s;Ref_Reads=%s;Alt_Reads=%s;Var=%.2f"
                                  % (array[0], array[1], array[2], array[3], ALT, GT, Ref_Reads, Alt_Reads, Var))
                    outfile.write("%\n")
    infile.close()
    outfile.close()


if __name__=="__main__":
    parser=argparse.ArgumentParser("")
    parser.add_argument("-v","--vcf",help="format vcf file",required=True)
    parser.add_argument("-n","--sample_name",help="sample name in vcf",required=True)
    parser.add_argument("-o","--outdir",help="output directory",required=True)
    parser.add_argument("-t","--transcript",help="Canonical transcript file",required=True)
    parser.add_argument("-r","--ref",help="annovar directory",required=True)
    args=parser.parse_args()
    format_vcf(args.vcf,args.sample_name,args.outdir)
    anno("%s/%s.format.vcf"%(args.outdir,args.sample_name),args.outdir,args.sample_name,args.transcript,args.ref)
