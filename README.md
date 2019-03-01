# COMP483_MiniProject
This repo contains code for the COMP 483 Mini Project for the grad level. It uses RefSeq data from four strains of E. coli that were gathered and analyzed by Subashchandrabose et al.

## Author
* **Paul J. Risteca**

## Version
1.0

## Dependencies
* **Python 3**
* **Biopython**
* **Prokka**
* **Bowtie2**
* **Tophat2**
  * **Samtools**
  * **Also uses Bowtie2**
* **Cufflinks Suite**

## Installing

You must clone the repository in command line to install. Copy the command below.
```
git clone https://github.com/tilld4wn/COMP483_MiniProject.git
```

## Proper Input
The input for this wrapper is a tab delimited text file. Use the supplied test data file titled sample_data.txt
The data contained in the input file should contain strain name, fasta .fna.gz FTP link, and feature_count.txt.gz FTP link. 
**EXAMPLE INPUT FROM TEST DATA FILE**
```
Strain	.fna.gz file ftp address	.txt.gz file ftp address
hm27	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_genomic.fna.gz	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/825/GCF_000387825.2_ASM38782v2/GCF_000387825.2_ASM38782v2_feature_count.txt.gz
hm46	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_genomic.fna.gz	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/845/GCF_000387845.2_ASM38784v2/GCF_000387845.2_ASM38784v2_feature_count.txt.gz
hm65	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_genomic.fna.gz	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/785/GCF_000387785.2_ASM38778v2/GCF_000387785.2_ASM38778v2_feature_count.txt.gz
hm69	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_genomic.fna.gz	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/387/865/GCF_000387865.2_ASM38786v2/GCF_000387865.2_ASM38786v2_feature_count.txt.gz
```

## Running Main 
To run, follow all preceding instructions and enter the following command in your terminal:
```
python3 ecoli_proj.py < sample_data.txt
```
As the pipeline runs, it will produce multiple files. Among these files will be the UPEC.log file. Do not attempt to open this file while running.
The pipeline should run in about 20 hours with the preset cpu input.

## Future Development/Improvements
This wrapper is quite crude and could use some modifications to run more efficiently. 
Additionally, I would like to improve user input functionality.
***IMPORTANT***
CUFFNORM IS GOING TO FAIL. THE FILES USED AS INPUT ARE NOT COMPATIBLE WITH THE SOFTWARE.
This error in the wrapper is another future improvement.

