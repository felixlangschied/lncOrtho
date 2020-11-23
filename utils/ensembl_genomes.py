"""
Download a set of genomes from EnsemblDB.
Should support NCBI-TaxIDs and/or Ensembl species names
"""

from ftplib import FTP

input = 'ensembl_input.txt'

with open(input, 'r') as file:
    ids = file.read().split('\n')
    ids = [eid.lower() for eid in ids]
    print(ids)

ftp = FTP('ftp.ensembl.org')
ftp.login()
ftp.cwd("/pub/current_fasta")

root_dirs = ftp.nlst()
print(root_dirs)
print(ftp.pwd())

for id in ids:
    if id in root_dirs:
        print('/' + id)
        ftp.cwd('/' + id)

        
        ftp.nlst()