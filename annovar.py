#Email:fanyucai1@126.com
#2021.08.26 version:1.0
#2021.09.08 version:2.0
    #   change parameter (annotate_variation.pl) in table_annovar.pl :--hgvs (gene-based annotation))
    #   change 'splicing' annoatation from "NM_030649:exon4:c.226-2A>G" to "NM_030649:intron3:c.226-2A>G"
    #   change delete the filter process using population MAF

import argparse
import os
import subprocess
import re
##########################################################指定输出标签
out_name=['Chr','Start','End','Ref','Alt','Func.refGene','Gene.refGene','GeneDetail.refGene',
          'ExonicFunc.refGene',	'AAChange.refGene','canonical_transcript','cytoBand','avsnp150','snp138','1000g2015aug_all','1000g2015aug_eas','genome_AF','genome_AF_eas','exome_AF','exome_AF_eas','ExAC_ALL','ExAC_EAS','esp6500siv2_all',
          'CLNALLELEID','CLNDN','CLNDISDB',	'CLNREVSTAT','CLNSIG','cosmic94_coding','cosmic94_noncoding','SIFT_pred','LRT_pred','FATHMM_pred','PROVEAN_pred','ClinPred_pred','MetaRNN_pred',
          'MutationTaster_pred','Polyphen2_HDIV_pred','Polyphen2_HVAR_pred','InterVar_automated','DP','Ref_Reads','Alt_Reads','F1R2','F2R1','AF']
AF=['1000g2015aug_all','1000g2015aug_eas','genome_AF','genome_AF_eas','exome_AF','exome_AF_eas','ExAC_ALL','ExAC_EAS','esp6500siv2_all']
########################################################### step1:format vcf and input annovar
def format_vcf(vcf,sample_name,outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    infile = open(vcf, "r")
    outfile = open("%s/%s.format.vcf" % (outdir, sample_name), "w")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    name,num = 0,0
    for line in infile:
        line = line.strip()
        array = line.split("\t")
        if line.startswith("#CHROM"):#############获得对应的样本
            for i in range(len(array)):
                if array[i] == sample_name:
                    name = i
                if array[i] == "FORMAT":
                    num = i
                    continue
        if not line.startswith("#"):
            if array[6]=="PASS" or array[6]==".":#考虑到VCF中如果没有设置过滤条件，则该行为应该都为"."
                tmp = array[num].split(":")#format
                info = array[name].split(":")#sample
                AD, AF, DP, F1R2, F2R1=[],[],"",[],[]
                for k in range(len(tmp)):
                    if tmp[k] == "AD":
                        AD = info[k].split(",")
                    if tmp[k] == "AF":
                        AF = info[k].split(",")
                    if tmp[k] == "F1R2":
                        F1R2 = info[k].split(",")
                    if tmp[k]=="F2R1":
                        F2R1 = info[k].split(",")
                    if tmp[k] == "DP":
                        DP= info[k]
                for j in range(1,len(AD)):
                    outfile.write("%s\t%s\t%s\t%s\t%s\t.\t.\tRef_Reads=%s;Alt_Reads=%s;AF=%.2f;F1R2=%s,%s;F2R1=%s,%s;DP=%s"
                                  % (array[0], array[1], array[2], array[3], array[4].split(",")[j-1],AD[0], AD[j],float(AF[j-1]) * 100,
                                     F1R2[0],F1R2[j],F2R1[0],F2R1[j],DP))
                    outfile.write("\n")
    infile.close()
    outfile.close()
########################################################### step2:vcf annoattion use annovar
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
    par += " -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f "
    par += " -nastring . -polish --thread 20 "
    subprocess.check_call("perl %s/table_annovar.pl %s %s/humandb -buildver hg19 -out %s -remove %s -vcfinput " % (annovar, vcf, annovar, out, par), shell=True)
    subprocess.check_call("rm -rf %s.hg19_multianno.vcf %s.avinput" % (out, out), shell=True)
    infile = open("%s.hg19_multianno.txt" % (out), "r")
    outfile = open("%s.raw.annovar.tsv" % (out), "w")
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
            p3 = re.compile(r'AF=(\d+.\d+)')
            p4 = re.compile(r'DP=([0-9]+)')
            p5 = re.compile(r'F1R2=(\d+,\d+)')
            p6 = re.compile(r'F2R1=(\d+,\d+)')
            Alt_Reads = p1.findall(line)
            Ref_Reads = p2.findall(line)
            AF = p3.findall(line)
            DP = p4.findall(line)
            F1R2=p5.findall(line)
            F2R1=p6.findall(line)
            ##########################format output knownCanonical transcript
            final_nm = "."#定义为点
            if array[dict['AAChange.refGene']]!=".":#如果是定义在编码区
                tmp = array[dict['AAChange.refGene']].split(",")
                final_nm = tmp[0]#默认变为第一个
                for i in range(len(tmp)):
                    if array[6] in transcript:#如果存在经典转录本
                        if re.search(r'%s'%(transcript[array[6]]),tmp[i]):#并且可以在注释中找到
                            final_nm=tmp[i]#重新定义输出
            for l in range(len(out_name)):
                if l == 0:
                    outfile.write("%s" % (array[dict[out_name[l]]]))
                elif out_name[l] == "AF":
                    outfile.write("\t%s" % (AF[0]) + "%")
                elif out_name[l] == "Alt_Reads":
                    outfile.write("\t%s" % (Alt_Reads[0]))
                elif out_name[l] == "Ref_Reads":
                    outfile.write("\t%s" % (Ref_Reads[0]))
                elif out_name[l] == "canonical_transcript":
                    outfile.write("\t%s" % (final_nm))
                elif out_name[l] == "DP":
                    outfile.write("\t%s" % (DP[0]))
                elif out_name[l] == "F1R2":
                    outfile.write("\t%s" % (F1R2[0]))
                elif out_name[l] == "F2R1":
                    outfile.write("\t%s" % (F2R1[0]))
                else:
                    outfile.write("\t%s" % (array[dict[out_name[l]]]))
            outfile.write("\n")
    infile.close()
    outfile.close()
    if os.path.exists("%s.hg19_multianno.txt" % (out)):
        subprocess.check_call("rm -rf %s.hg19_multianno.txt" % (out), shell=True)
########################################################### step3:format annovar
def change_anno(annovar,format_anno_file):
    infile=open(annovar,"r")
    outfile=open(format_anno_file,"w")
    for line in infile:
        line=line.strip()
        array=line.split("\t")
        if array[5]=="splicing" and array[7]!=".":#######判断该列是否是可变剪切列
            a = array[7].split(";")######按照逗号分开该列
            new=""
            for j in range(0,len(a)):
                exon_num=int(re.search(r'exon(\d+)',a[j]).group(1))###获得对应的外显子编号
                if re.search(r'-1', a[j]) or re.search(r'-2', a[j]):
                    new += re.sub(r'exon(\d+)', "intron%s" % (exon_num - 1), a[j])+";"
                elif re.search(r'\+1', a[j]) or re.search(r'\+2', a[j]):
                    new += re.sub(r'exon', 'intron', a[j])+";"
                else:
                    print(line)
                    exit()
            new=new.strip(";")
            #########################################################打印所有更改后的结果
            for k in range(0,len(array)):
                if k==9:
                    outfile.write("%s\t" % (new))
                elif k==7:
                    outfile.write(".\t")
                elif k!=len(array)-1:
                    outfile.write("%s\t"%(array[k]))
                else:
                    outfile.write("%s\n"%(array[k]))
        else:
            outfile.write("%s\n" % (line))
    infile.close()
    outfile.close()
    print("This process end.")

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
    change_anno("%s/%s.raw.annovar.tsv"%(args.outdir,args.sample_name),"%s/%s.final.anno.tsv"%(args.outdir,args.sample_name))