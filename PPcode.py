#Start with imports python will need to run the code
import wget
import sys
import os
#Create function to see if output folder already exists
def out_folder(folder):
    # If it already exists we do not create the folder again
    if os.path.isdir(folder)==False:
        os.mkdir(folder)
#Call output folder function with designated name and change working directory to that folder
out_folder('PipelineProject_Samantha_Rutherford')
os.chdir('PipelineProject_Samantha_Rutherford')

#Import test data by setting a stock URL variable
url = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/'
#List of ids in order of patient and day 2 or 6 post infection
sra_ids = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
i = 0
sra_urls = []
#Set up an empty list to assign urls and a for loop to create individual urls
for ids in sra_ids:
    temp =  url + sra_ids[i] + '/' + sra_ids[i]
    sra_urls.append(temp)
    #Use wget to download sra files and os to use fastq to split the paired-end files
    #wget.download(sra_urls[i])
    #cmd = 'fastq-dump -I --split-files ' + sra_ids[i]
    #os.system(cmd)
    i+=1
#Assign each file name to a variable for now //break//
d1_2_1 = sra_ids[0] + '_1.fastq'
d1_2_2 = sra_ids[0] + '_2.fastq'
d1_6_1 = sra_ids[1] + '_1.fastq'
d1_6_2 = sra_ids[1] + '_2.fastq'
d3_2_1 = sra_ids[2] + '_1.fastq'
d3_2_2 = sra_ids[2] + '_2.fastq'
d3_6_1 = sra_ids[3] + '_1.fastq'
d3_6_2 = sra_ids[3] + '_2.fastq'
