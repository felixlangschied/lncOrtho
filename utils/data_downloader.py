"""
Download a set of genomes from EnsemblDB or NCBI
Should support NCBI-TaxIDs and/or Ensembl species names
Makes use of wget which can be installed with:
pip install wget

Created on Fri Nov 20 14:24:19 2020

@author: felixl

List of available flavors when downloading from NCBI:

genomic.fna.gz
genomic.gbff.gz
genomic.gff.gz
protein.faa.gz
protein.gpff.gz
feature_table.txt.gz
wgsmaster.gbff.gz
cds_from_genomic.fna.gz
rna_from_genomic.fna.gz
feature_count.txt.gz
translated_cds.faa.gz
genomic.gtf.gz
genomic_gaps.txt.gz

"""

from ftplib import FTP
import re
import wget
import os


# download genomes as contigs from ensembl FTP based on ensembl ids
def extract_ensembl(ids, output):
    ids = [eid.lower() for eid in ids]
    print(ids)
    ftp = FTP('ftp.ensembl.org')
    ftp.login()
    ftp.cwd("/pub/current_fasta")
    root_dirs = ftp.nlst()

    for id in ids:
        if id in root_dirs:
            ftp.cwd('/pub/current_fasta/' + id + '/dna')
            files = ftp.nlst()
            r = re.compile(".*dna.toplevel.fa.gz")
            target = list(filter(r.match, files))[0]
            try:
                t_location = (
                    'ftp://ftp.ensembl.org/pub/current_fasta/{}/dna/{}'
                        .format(id, target)
                )
                print(t_location)
                os.chdir(output)
                wget.download(t_location)
                print('# Successfully downloaded file for {}'.format(id))
            except:
                print('# {} not found'.format(target))


# download genomes from NCBI based on GCF id
def extract_GCF(ids, flavor, output):
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd("/genomes/all/")

    for id in ids:
        ftp.cwd('/genomes/all/')
        id = id.split('.')[0]
        typ = id.split('_')[0]
        id = id.split('_')[1]
        part1 = id[0:3]
        part2 = id[3:6]
        part3 = id[6:9]

        url = '/'.join([typ, part1, part2, part3])
        ftp.cwd(url)
        assemblies = ftp.nlst()
        # will download the top file
        target = assemblies[0]
        if len(assemblies) > 1:
            print('# Found multiple assemblies for {}.\n# Downloading most recent: {}'
                  .format(id, assemblies[0])
                  )
        ftp.cwd(target)
        files = ftp.nlst()

        # find file types that can be downloaded (flavors)
        # ftypes = [file.split('.')[-3:] for file in files if '.gz' in file]
        # ftypes = ['.'.join(parts) for parts in ftypes]
        # ftypes = [ftype.split('_')[2:] for ftype in ftypes]
        # ftypes = ['_'.join(parts) for parts in ftypes]
        # print('\n'.join(ftypes))

        r = re.compile('.*' + flavor)
        target_file = list(filter(r.match, files))[0]
        download_url = (
            'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}'
                .format(url, target, target_file)
        )
        os.chdir(output)
        try:
            wget.download(download_url)
            print('# Successfully downloaded {}'.format(target))
        except:
            print('# {} not found'.format(target_file))


# testing
flav = 'genomic.gtf.gz'
input = '/home/felixl/genomeDownloader/ncbi_gcf.txt'
# input = 'ensembl_input.txt'
output = output = '/home/felixl/genomeDownloader'
with open(input, 'r') as file:
    id_list = file.read().strip().split('\n')

    # TODO: Test what kind of identifier is in the file
    # ncbi starts with: GCA_ or GCF_

extract_GCF(id_list, flav, output)
# extract_ensembl(id_list, output)