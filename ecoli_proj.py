#!/usr/bin/env python3
import os
from Bio import SeqIO



#Set cwd to the current directory that you are in
cwd = os.getcwd()
print(cwd)
#Create a new directory named 'FirstName_LastName'
if(not os.path.exists('Paul_Risteca')):
    os.mkdir('./Paul_Risteca')
    os.chdir('Paul_Risteca')
else:
    os.chdir('Paul_Risteca')

# Create some storage for each strain
# HM27 = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz','ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz')
# HM46 = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz','ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz')
# HM65 = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz','ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz')
# HM69 = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz','ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz')
#
#
# os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz -O HM27_fasta.fna.gz")
# os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz -O HM27_feature_count.txt.gz")
# os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz -O HM46_fasta.fna.gz")
# os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz -O HM46_feature_count.txt.gz")
# os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz -O HM65_fasta.fna.gz")
# os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz -O HM65_feature_count.txt.gz")
# os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz -O HM69_fasta.fna.gz")
# os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz -O HM69_feature_count.txt.gz")

# os.system("prefetch SRR1278956")
# os.system("prefetch SRR1278960")
# os.system("prefetch SRR1283106")
# os.system("prefetch SRR1278963")
#
# os.system("gunzip -d *.gz")
# SRA File addresses
# HM27 = ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278956/SRR1278956.sra
# HM46 = ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278960/SRR1278960.sra
# HM65 = ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR128/SRR1283106/SRR1283106.sra
# HM69 = ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR127/SRR1278963/SRR1278963.sra
#
HM27 = ("HM27_fasta.fna","HM27_feature_count.txt","HM27", 'hm27_prokka', "./hm27_prokka/hm27_prokka.txt", "SRR1278956","./hm27_prokka/hm27_prokka.gff","./HM27_tophat_out/accepted_hits.bam")
HM46 = ("HM46_fasta.fna","HM46_feature_count.txt","HM46", 'hm46_prokka'," ","SRR1278960","./hm46_prokka/hm46_prokka.gff","./HM46_tophat_out/accepted_hits.bam")
HM65 = ("HM65_fasta.fna","HM65_feature_count.txt","HM65", 'hm65_prokka'," ","SRR1283106","./hm65_prokka/hm65_prokka.gff","./HM65_tophat_out/accepted_hits.bam")
HM69 = ("HM69_fasta.fna","HM69_feature_count.txt","HM69", 'hm69_prokka'," ","SRR1278963","./hm65_prokka/hm65_prokka.gff","./HM69_tophat_out/accepted_hits.bam")

# 2&3
# def recordInfo(strain,log_file):
#     record = list(SeqIO.parse(strain[0],"fasta"))
#     number_of_contigs = len(record)
#     log_file.write("There are "+str(number_of_contigs)+" contigs in the assembly "+strain[2]+". \n")
#     count = 0
#     for read in record:
#         if len(read.seq) > 1000:
#             count += len(read.seq)
#     log_file.write("There are "+str(count)+" bp in the assembly "+strain[2]+". \n")
#     return
#
# # prokka
# def annotaion_prokka(prokka_prefix, strain_fasta,log_file):
#
#     prokka_command = "prokka --usegenus --genus Escherichia --cpus 4 --prefix {} {}".format(prokka_prefix, strain_fasta)
#     log_file.write(prokka_command+"\n")
#     os.system(prokka_command)
#
# def tophat_cufflinks(strain_fasta, read1, read2, log_file):
#     fasta_copy = "cp {} {}.fa".format(strain_fasta[0], strain_fasta[2])
#     bowt2_command = "bowtie2-build --threads 6 -f {} {}".format(strain_fasta[0], strain_fasta[2])
#     tophat_command = "tophat -p 6 -o {} {} {} {}".format(strain_fasta[2]+"_tophat_out",strain_fasta[2],read1,read2)
#     os.system(bowt2_command)
#     os.system(fasta_copy)
#     log_file.write(tophat_command+"\n")
#     os.system(tophat_command)

def cufflinks(gff_file,cuff_out,bam):
    cufflink_command = "cufflinks -p 6 -G {} -o {} {}".format(gff_file, cuff_out, bam)
    os.system(cufflink_command)

def cuffmerge_norm(assembly_file, merged_gtf, strain1, strain2, strain3, strain4):
    cuffmerge_command = "cuffmerge -p 6 -o {} {}".format('merged_ecoli', assembly_file)
    os.system(cuffmerge_command)
    cuffnorm_command = "cuffnorm -o diff_results -p 6 {} {} {} {}".format(merged_gtf, strain1[7], strain2[7], strain3[7], strain4[7])
    os.system(cuffnorm_command)
# def Record(strain):
#     recordSeqIO.parse(strain[0],"fasta"):
def main():
    # output_file = open ("UPEC.log","a")
    # recordInfo(HM27,output_file)
    # recordInfo(HM46,output_file)
    # recordInfo(HM65,output_file)
    # recordInfo(HM69,output_file)
    #
    # # testing with hm27
    # annotaion_prokka(HM27[3],HM27[0],output_file)
    # annotaion_prokka(HM46[3],HM46[0],output_file)
    # annotaion_prokka(HM65[3],HM65[0],output_file)
    # annotaion_prokka(HM69[3],HM69[0],output_file)
    #
    # # try to see if this will add to UPEC
    # os.system("cat ./hm27_prokka/hm27_prokka.txt >> UPEC.log")
    # os.system("cat ./hm46_prokka/hm46_prokka.txt >> UPEC.log")
    # os.system("cat ./hm65_prokka/hm65_prokka.txt >> UPEC.log")
    # os.system("cat ./hm69_prokka/hm69_prokka.txt >> UPEC.log")
    #
    # os.system("prefetch SRR1278956")
    # os.system("prefetch SRR1278960")
    # os.system("prefetch SRR1283106")
    # os.system("prefetch SRR1278963")
    # os.system("fastq-dump -I --split-files ~/ncbi/public/sra/SRR1278956.sra")
    # os.system("fastq-dump -I --split-files ~/ncbi/public/sra/SRR1278960.sra")
    # os.system("fastq-dump -I --split-files ~/ncbi/public/sra/SRR1283106.sra")
    # os.system("fastq-dump -I --split-files ~/ncbi/public/sra/SRR1278963.sra")
    #
    # tophat_cufflinks(HM27,HM27[5]+"_1.fastq", HM27[5]+"_2.fastq",output_file)
    # tophat_cufflinks(HM46,HM46[5]+"_1.fastq", HM46[5]+"_2.fastq",output_file)
    # tophat_cufflinks(HM65,HM65[5]+"_1.fastq", HM65[5]+"_2.fastq",output_file)
    # tophat_cufflinks(HM69,HM69[5]+"_1.fastq", HM69[5]+"_2.fastq",output_file)
    #
    # cufflinks(HM27[6],HM27[2],HM27[7])
    # cufflinks(HM46[6],HM46[2],HM46[7])
    # cufflinks(HM65[6],HM65[2],HM65[7])
    # cufflinks(HM69[6],HM69[2],HM69[7])

    with open('ecoli_assemblies.txt', 'w') as assemble_file:
        assemble_file.write("./HM27/transcripts.gtf\n")
        assemble_file.write("./HM46/transcripts.gtf\n")
        assemble_file.write("./HM65/transcripts.gtf\n")
        assemble_file.write("./HM69/transcripts.gtf\n")

    cuffmerge_norm('ecoli_assemblies.txt', './merged_ecoli/merged.gtf', HM27, HM46, HM65, HM69)
    output_file.close()



if __name__ == '__main__':
    main()
