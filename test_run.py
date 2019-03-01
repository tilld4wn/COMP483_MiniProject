import os
import sys
from Bio import SeqIO
# Set cwd to the current directory that you are in
cwd = os.getcwd()
print(cwd)
#Create a new directory named 'FirstName_LastName'
if(not os.path.exists('Paul_Risteca')):
    os.mkdir('./Paul_Risteca')
    os.chdir('Paul_Risteca')
else:
    os.chdir('Paul_Risteca')

# ftp address lists
hm27_ftp = ['hm27']
hm46_ftp = ['hm46']
hm65_ftp = ['hm65']
hm69_ftp = ['hm69']

# stores strain name, fasta file, feature count file, and prokka output file prefix
hm27_data = ['hm27', 'hm27_fasta.fna', 'hm27_feature_count.txt', 'hm27_prokka']
hm46_data = ['hm46', 'hm46_fasta.fna', 'hm46_feature_count.txt', 'hm46_prokka']
hm65_data = ['hm65', 'hm65_fasta.fna', 'hm65_feature_count.txt', 'hm65_prokka']
hm69_data = ['hm69', 'hm69_fasta.fna', 'hm69_feature_count.txt', 'hm69_prokka']

# Retrieve ftp addresses from input file "sample_data.txt"
def get_ftp(strain, data_file):
    # with open('sample_data.txt','r') as data_file:
    for line in data_file:
        temp = line.strip('\n').split('\t')
        if temp[0] == strain[0]:
            strain.append(temp[1])
            strain.append(temp[2])
    return

# Pulls data from NCBI using wget and ftp addresses
def download_data(strain):
    wget_command_fasta = "wget {} -O {}".format(strain[1], strain[0]+'_fasta.fna.gz')
    wget_command_feature_content = "wget {} -O {}".format(strain[2], strain[0]+'_feature_count.txt.gz')
    os.system(wget_command_fasta)
    os.system(wget_command_feature_content)
    return

# Reads the fasta file in and determines number of contigs present and length of assembly when all reads longer than 1000 bp
def seqCount(data_lst, log_file):
    for strain in data_lst:
        record = list(SeqIO.parse(strain[1], 'fasta'))
        contigs = str(len(record))
        log_file.write("There are {} contigs in the assembly {}.\n".format(contigs, strain[0]))
        count = 0
        for read in record:
            if len(read.seq) > 1000:
                count += len(read.seq)
        log_file.write("There are {} bp in the assembly {}.\n".format(str(count), strain[0]))
    return


def main():
    log_file = open("UPEC.log", "a")

    data_file = []
    for data in sys.stdin:
        data_file.append(data)

# Getting data into working directory to begin analysis
    ftp_lst = [hm27_ftp, hm46_ftp, hm65_ftp, hm69_ftp]

    for strain in ftp_lst:
        get_ftp(strain, data_file)

    for strain_info in ftp_lst:
        download_data(strain_info)

    os.system("gunzip -d *.gz")

# Finding number of contigs and assembly length
    data_lst = [hm27_data, hm46_data, hm65_data, hm69_data]
    seqCount(data_lst, log_file)

if __name__ == '__main__':
    main()
