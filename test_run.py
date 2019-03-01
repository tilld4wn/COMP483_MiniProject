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

# stores location of prokka .txt output
hm27_prokka = ['./hm27_prokka/hm27_prokka.txt']
hm46_prokka = ['./hm46_prokka/hm46_prokka.txt']
hm65_prokka = ['./hm65_prokka/hm65_prokka.txt']
hm69_prokka = ['./hm69_prokka/hm69_prokka.txt']

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

def prokka_run(data_lst, prokka_lst, log_file):
    log_file.write("Now running prokka.")
    for strain in data_lst:
        prokka_command = "prokka --usegenus --genus Escherichia --cpus 10 --prefix {} {}".format(strain[3], strain[1])
        log_file.write(prokka_command)
        os.system(prokka_command)

    for strain, prokInfo in zip(data_lst, prokka_lst):
        feature_count = {"CDS":0, "tRNA":0}
        prokka_count = {"CDS":0, "tRNA":0}
        with open(strain[2], "r") as in_file:
            line_arr = []
            for line in in_file:
                line_arr.append(line)

            for line in line_arr:
                temp = line.strip('\n').split('\t')

                if temp[0] == "CDS" and temp[1] == "with_protein":
                    feature_count["CDS"] = int(temp[6])

                if temp[0] == "tRNA":
                    feature_count["tRNA"] = int(temp[6])

        with open(prokInfo[0], "r") as in_file:
            line_arr = []
            for line in in_file:
                line_arr.append(line)

            for line in line_arr:
                temp = line.strip('\n').split(':')

                if temp[0] == "CDS":
                    prokka_count["CDS"] = int(temp[1])

                if temp[0] == "tRNA":
                    prokka_count["tRNA"] = int(temp[1])
        diffCDS = feature_count["CDS"]-prokka_count["CDS"]
        diffT = feature_count["tRNA"]-prokka_count["tRNA"]
        if diffCDS > 0 and diffT > 0:
            log_file.write("Prokka found {} less CDS and {} less tRNA than the RefSeq in assembly {}.".format(abs(diffCDS), abs(diffT), strain[0]))
        elif diffCDS > 0 and diffT <= 0:
            log_file.write("Prokka found {} less CDS and {} more tRNA than the RefSeq in assembly {}.".format(abs(diffCDS), abs(diffT), strain[0]))
        elif diffCDS <= 0 and diffT > 0:
            log_file.write("Prokka found {} more CDS and {} less tRNA than the RefSeq in assembly {}.".format(abs(diffCDS), abs(diffT), strain[0]))
        elif diffCDS <= 0 and diffT <= 0:
            log_file.write("Prokka found {} more CDS and {} more tRNA than the RefSeq in assembly {}.".format(abs(diffCDS), abs(diffT), strain[0]))
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

# Running Prokka on data
    prokka_lst = [hm27_prokka, hm46_prokka, hm65_prokka, hm69_prokka]
    prokka_run(data_lst, prokka_lst, log_file)

if __name__ == '__main__':
    main()
