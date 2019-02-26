#!/usr/bin/env python
import os
import logging
from Bio import SeqIO



#Set cwd to the current directory that you are in
cwd = os.getcwd()
print(cwd)
#Create a new directory named 'FirstName_LastName'
os.mkdir('./Paul_Risteca')
os.chdir('Paul_Risteca')

# Create some storage for each strain
# HM27 = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz','ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz')
# HM46 = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz','ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz')
# HM65 = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz','ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz')
# HM69 = ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz','ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz')


os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz -O HM27_fasta.fna.gz")
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz -O HM27_feature_count.txt.gz")
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz -O HM46_fasta.fna.gz")
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz -O HM46_feature_count.txt.gz")
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz -O HM65_fasta.fna.gz")
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz -O HM65_feature_count.txt.gz")
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz -O HM69_fasta.fna.gz")
os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz -O HM69_feature_count.txt.gz")

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

HM27 = ("HM27_fasta.fna","HM27_feature_count.txt","HM27", 'hm27_prokka', "./hm27_prokka/hm27_prokka.txt")
HM46 = ("HM46_fasta.fna","HM46_feature_count.txt","HM46", 'hm46_prokka')
HM65 = ("HM65_fasta.fna","HM65_feature_count.txt","HM65", 'hm65_prokka')
HM69 = ("HM69_fasta.fna","HM69_feature_count.txt","HM69", 'hm69_prokka')

# 2&3
def recordInfo(strain,log_file):
    record = list(SeqIO.parse(strain[0],"fasta"))
    number_of_contigs = len(record)
    log_file.write("There are "+str(number_of_contigs)+" contigs in the assembly "+strain[2]+". \n")
    count = 0
    for read in record:
        if len(read.seq) > 1000:
            count += len(read.seq)
    log_file.write("There are "+str(count)+" bp in the assembly "+strain[2]+". \n")
    return

# prokka
def annotaion_prokka(prokka_prefix, strain_fasta,log_file):

    prokka_command = "prokka --usegenus --genus Escherichia --cpus [4] --myprefix {} {}".format(prokka_prefix, strain_fasta)
    log_file.write(prokka_command+"\n")
    os.system(prokka_command)


# def Record(strain):
#     recordSeqIO.parse(strain[0],"fasta"):
def main():
    output_file = open ("UPEC.log","a")
    recordInfo(HM27,output_file)
    recordInfo(HM46,output_file)
    recordInfo(HM65,output_file)
    recordInfo(HM69,output_file)

    # testing with hm27
    annotaion_prokka(HM27[3], HM27[0], output_file)
    output_file.close()

    # try to see if this will add to UPEC
    os.system("cat ./hm27_prokka/hm27_prokka.txt >> UPEC.log")

if __name__ == '__main__':
    main()
