#Start with imports python will need to run the code
import logging
import wget
import sys
import os
from Bio import Entrez
from Bio import SeqIO

#Part 1 Sample Accession
#Import test data by setting a stock URL variable
url = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/'
#List of ids in order of patient and day 2 or 6 post infection
sra_ids = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

#Create function to see if output folder already exists
def out_folder(folder, ids):
    # If it already exists we do not create the folder again
    if os.path.isdir(folder)==False:
        os.mkdir(folder)
        #Change working directory to the created folder
        os.chdir('PipelineProject_Samantha_Rutherford')
        #Set up an empty list to assign urls and a for loop to create individual urls
        i=0
        sra_urls = []
        for id in ids:
            temp =  url + ids[i] + '/' + ids[i]
            sra_urls.append(temp)
            #Use wget to download sra files and os to use fastq to split the paired-end files
            wget.download(sra_urls[i])
            cmd = 'fastq-dump -I --split-files ' + ids[i]
            os.system(cmd)
            i+=1
    else:
        #Be sure output directory is the correct one
        os.chdir('PipelineProject_Samantha_Rutherford')

#Call output folder function with designated name and our ids list
out_folder('PipelineProject_Samantha_Rutherford', sra_ids)

#Create log file
logging.basicConfig(filename='PipelineProject.log', filemode='w', level=logging.INFO)

#Assign each file name to a variable just in case
d1_2_1 = sra_ids[0] + '_1.fastq'
d1_2_2 = sra_ids[0] + '_2.fastq'
d1_6_1 = sra_ids[1] + '_1.fastq'
d1_6_2 = sra_ids[1] + '_2.fastq'
d3_2_1 = sra_ids[2] + '_1.fastq'
d3_2_2 = sra_ids[2] + '_2.fastq'
d3_6_1 = sra_ids[3] + '_1.fastq'
d3_6_2 = sra_ids[3] + '_2.fastq'


#Step 2 Indexing

#Use Entrez to access the HCMV genome and SeqIO for later parsing
Entrez.email = 'srutherford@luc.edu'
handle = Entrez.efetch(db ='nucleotide', id = 'NC_006273.2', rettype = 'gbwithparts', retmode = 'text')
genome = next(SeqIO.parse(handle, 'genbank'))
#Create an empty list for coding sequences
cds = []
#Use for loop to access all CDS
for feature in genome.features:
    if feature.type == 'CDS':
        cds.append(feature)
#Second loop is created to remove genbank information so the sequence can be put into FASTA format
seqs = []
for feature in cds:
    seq = genome[int(feature.location.start):int(feature.location.end)]
    if feature.location.strand == -1:
        seq.seq = seq.seq.reverse_complement()
    seq.description = feature.qualifiers['product'][0]
    seqs.append(seq)
#Write coding sequences to a fasta file
SeqIO.write(seqs, 'NC_006273.2.fasta', 'fasta')
#Create variable
logging.info('The HCMV genome (NC_006273.2) has %s CDS.', len(seqs))
#Build the HCMV index for Kallisto with coding sequence fasta
os.system('time kallisto index -i NC_006273.2.idx NC_006273.2.fasta --make-unique -k 31')


#Step 3 Quantification


