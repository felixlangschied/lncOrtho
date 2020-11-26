"""
Created on Fri Nov 20 14:24:19 2020
@author: felixl

Download a set of genomes from EnsemblDB or NCBI
based on NCBI accession number or Ensembl genus/species name

use genome_downloader.py -h for instructions
"""

from ftplib import FTP
import re
import wget
import os
import subprocess as sp
import sys
import argparse
import textwrap


# download genomes as contigs from ensembl FTP based on ensembl ids
def extract_ensembl(ids, output, raw_flavor):
    ids = [eid.lower() for eid in ids]
    ftp = FTP('ftp.ensembl.org')
    ftp.login()
    ftp.cwd("/pub/current_fasta")
    root_dirs = ftp.nlst()
    flavor = raw_flavor.replace('ensembl_', '')

    for id in ids:
        if id in root_dirs:
            ftp.cwd('/pub/current_fasta/{}/{}'.format(id, flavor))
            files = ftp.nlst()

            if flavor == 'dna':
                r = re.compile(".*dna.toplevel.fa.gz")
                target = list(filter(r.match, files))[0]
            elif flavor == 'pep':
                r = re.compile('.*pep.all.fa.gz')
                target = list(filter(r.match, files))[0]
            try:
                t_location = (
                    'ftp://ftp.ensembl.org/pub/current_fasta/{}/{}/{}'
                        .format(id, flavor, target)
                )
                os.chdir(output)
                wget.download(t_location)
                print('# Successfully downloaded file for {}'.format(id))
            except:
                print('# {} not found'.format(target))
        else:
            print('{} not found on the Ensembl FTP-server'.format(id))


# download genomes from NCBI based on GCF/GCA id
def extract_GCF(ids, output, flavor):
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
        # print(ftypes)
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
            print('# Successfully downloaded {} for: {}'.format(flavor, target))
        except:
            print('# {} not found'.format(target_file))


# allow boolean arguments for parsing
def str2bool(v):
    if isinstance(v, bool):
        return v
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# RUN PROGRAM
def main():
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        prog='python data_downloader.py',
        description='Download a set of genomes from EnsemblDB or NCBI '
                    'based on NCBI accession number or Ensembl genus/species name ',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i', '--input', metavar='<path>', type=str,
        help='Path to list of NCBI-Accessions or ensembl genus/species names',
    )
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='Output path'
    )
    parser.add_argument(
        '-f', '--flavor', metavar='<str>', type=str,
        help=textwrap.dedent('''\
        Type of data you want to download.

        # List of available flavors when downloading from NCBI:
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

        # If you want to download from Ensembl use a flavor from this list:
        ensembl_dna
        ensembl_pep

        # Currently not supported ensembl formats:
        ensembl_cdna
        ensembl_cds
        ensembl_dna_index
        ensembl_ncrna

        ''')
    )
    parser.add_argument(
        '-u', '--unpack', type=str2bool, metavar='True/False', nargs='?', const=True, default=True,
        help="Set to False if you don't want to unpack the downloaded files straight away (Default=True)"
    )
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
    # parse input
    output = args.output
    input = args.input
    flavor = args.flavor
    unpack = args.unpack

    # MAIN BODY
    # Check if output folder exists or create it otherwise
    if not os.path.isdir(output):
        print('# Creating output folder')
        cmd = 'mkdir {}'.format(output)
        try:
            sp.run(cmd, shell=True, check=True, stderr=sp.PIPE)
        except sp.CalledProcessError:
            print('# Could not create output folder at:\n'
                  '{}\n'
                  'Exiting..'.format(output))

    # Verify flavor
    ncbi_flavors = ['genomic.fna.gz', 'genomic.gbff.gz', 'genomic.gff.gz',
                    'protein.faa.gz', 'protein.gpff.gz',
                    'feature_table.txt.gz', 'wgsmaster.gbff.gz',
                    'cds_from_genomic.fna.gz', 'rna_from_genomic.fna.gz',
                    'feature_count.txt.gz', 'translated_cds.faa.gz',
                    'genomic.gtf.gz', 'genomic_gaps.txt.gz']

    ensembl_supported = ['ensembl_dna', 'ensembl_pep']
    ensembl_unsupported = ['ensembl_cdna', 'ensembl_cds',
                           'ensembl_dna_index', 'ensembl_ncrna']

    if not flavor in ncbi_flavors + ensembl_supported:
        print('# Unknown flavor: {}'.format(flavor))
        sys.exit()
    elif flavor in ensembl_unsupported:
        print('# Flavor not yet supported: {}'.format(flavor))
        sys.exit()

    # Read input
    with open(input, 'r') as file:
        id_list = file.read().strip().split('\n')
        # convert spaces to underscores
        id_list = [id.replace(' ', '_') for id in id_list]
    # Download
    if 'ensembl_' in flavor:
        extract_ensembl(id_list, output, flavor)
    else:
        extract_GCF(id_list, output, flavor)
    # Unpack the downloaded files
    if unpack:
        print('# Starting to unpack downloaded files..')
        os.chdir(output)
        cmd = 'gunzip *.gz'
        sp.run(cmd, shell=True)
    print('# Finished')


if __name__ == "__main__":
    main()