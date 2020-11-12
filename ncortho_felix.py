'''
# Re-implementation of ncOrtho in Python
# TODO: include license, author information etc
'''

# Modules import

# Python
import argparse
import multiprocessing as mp
import sys
import os
import glob
import subprocess as sp

# internal modules
from lib.ncortho_main import ncortho


###############################################################################

# Main function
def main():
    # Print header
    print('\n' + '#' * 57)
    print('###' + ' ' * 51 + '###')
    print('###   ncOrtho - ortholog search for non-coding RNAs   ###')
    print('###' + ' ' * 51 + '###')
    print('#' * 57 + '\n')

    # Parse command-line arguments
    # Define global variables
    parser = argparse.ArgumentParser(
        prog='python ncortho.py', description='ncRNA orthology prediction tool'
    )
    # cpu, use maximum number of available cpus unless specified otherwise
    parser.add_argument(
        '-c', '--cpu', metavar='int', type=int,
        help='number of cpu cores ncOrtho should use', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
    )
    # covariance models folder
    parser.add_argument(
        '-m', '--models', metavar='<path>', type=str,
        help='path to your covariance models'
    )
    # mirna data
    parser.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str,
        help='path to your reference micrornas'
    )
    # output folder
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='path for the output folder'
    )
    # query genome
    parser.add_argument(
        '-q', '--query', metavar='<.fa>', type=str,
        help='path to query genome'
    )
    # reference genome blastDB
    parser.add_argument(
        '-r', '--reference', metavar='<basename>', type=str,
        help='Path to reference genome as fasta file or to an existing BlastDB of the reference genome'
    )
    # bit score cutoff for cmsearch hits
    parser.add_argument(
        '-t', '--cutoff', metavar='float', type=float,
        help='cmsearch bit score cutoff', nargs='?', const=0.6, default=0.6
    )
    # length filter to prevent short hits
    parser.add_argument(
        '-l', '--msl', metavar='float', type=float,
        help='hit length filter', nargs='?', const=0.9, default=0.9
    )
    # cleanup
    parser.add_argument(
        '-x', '--cleanup', metavar='True/False', type=bool,
        help='should ncOrtho delete all intermediate files', nargs='?', default=True
    )
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()


    #########################
    # Parse input arguments #
    #########################

    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if args.cpu > available_cpu:
        print(
            '# Error: The provided number of CPU cores is higher than the '
            'number available on this system. Exiting...'
        )
        sys.exit(1)
    else:
        cpu = args.cpu

    mirnas = args.ncrna
    models = args.models
    output = args.output
    query = args.query
    reference = args.reference
    cm_cutoff = args.cutoff
    cleanup = args.cleanup

    # Not in use yet
    msl = args.msl
    # blast_cutoff = args.blastc

    ##################################
    # Checking Validity of Arguments #
    ##################################

    # TODO: remove '/' at the end of path inputs

    # covariance models folder
    if os.path.isdir(models):
        model_filesuff = [suff.split('.')[-1] for suff in glob.glob(models + '/*')]
        if 'cm' in model_filesuff:
            print('Found covariance models.\n')
        else:
            print('No valid covariance models (<.cm>) found in {}.\n'
                  'Please check your input.\n'
                  'Exiting..'
                  .format(models))
            sys.exit()
    else:
        print('Covariance models must be supplied as a directory that contains valid models (<.cm>)\n'
              'Exiting..')
        sys.exit()

    # mirna data
    if os.path.isfile(mirnas):
        print('File with miRNAs found')
        # consistency of file will be checked by ncortho_main.py and mirna_maker.py
    else:
        print('Could not locate file with miRNAs.\n'
              'Exiting..')
        sys.exit()

    # check if output folder exists, create it otherwise
    if not os.path.isdir(output):
        sp.run('mkdir -p {}'.format(output), shell=True)

    # check how the reference is given
    # test if reference DB exists or has to be created
    file_extensions = ['.nhr', '.nin', '.nsq']
    out_data = output + '/data'
    fname = reference.split('/')[-1]
    db_name = fname.replace('.fa', '')
    db_files = []
    if (
            os.path.isfile(reference)
            and fname.split('.')[-1] == 'fa'
    ):
        print(
            'Reference given as FASTA file, testing if BlastDB exists in\n'
            '{}'
            .format(out_data)
        )
        if not os.path.isdir(out_data):
            mkdir_cmd = 'mkdir {}'.format(out_data)
            sp.call(mkdir_cmd, shell=True)
        for fe in file_extensions:
            db_files.append(glob.glob(out_data + '/' + db_name + '.*' + fe))
        if not [] in db_files:
            print("BLAST database for the reference species found.\n")
            ref_blast_db = out_data + '/' + db_name
        else:
            print(
                'BLAST database for the reference species does not exist.\n'
                'Constructing BLAST database.'
            )
            # At least one of the BLAST db files is not existent and has to be
            # created.
            db_command = 'makeblastdb -in {} -dbtype nucl -out {}/{}'.format(reference, out_data, db_name)
            try:
                sp.run(db_command, shell=True, check=True, stderr=sp.PIPE)
                ref_blast_db = out_data + '/' + db_name
            except sp.CalledProcessError:
                print('Could not create blastDB from the reference input file.\nExiting..\n')
                sys.exit()
    else:
        for fe in file_extensions:
            db_files.append(glob.glob(reference + '.*' + fe))
        if not [] in db_files:
            print("Reference given as blastDB basename. Starting analysis\n")
            ref_blast_db = reference
        else:
            print(
                'BLAST database for the given reference does not exist.\n'
                'Trying to construct a BLAST database from {}.'
                .format(reference)
            )
            try:
                db_command = 'makeblastdb -in {} -out {}/{} -dbtype nucl'.format(reference, out_data, db_name)
                sp.run(db_command, shell=True, check=True, stderr=sp.PIPE)
                ref_blast_db = out_data + '/' + db_name
            except sp.CalledProcessError:
                print(
                    'Was not able to construct the BLASTdb for {}\n'
                    'Please check you input.\n'
                    .format(reference)
                )
                sys.exit()

    # check type of query: single FASTA-file, Folder with FASTA-files or txt document with paths to FASTA files
    # TODO: Write summary output to file, per species or per miRNA or both?
    # TODO: Cleanup folders with no hits
    #variable for saving summary output
    hit_dict = {}
    if (
        os.path.isfile(query)
        and query.split('.')[-1] == 'fa'
    ):
        print('Query given as single FASTA-file. Starting anaysis\n')
        query_species = query.split('/')[-1].split('.')[0]
        taxon_out = '{}/{}'.format(output, query_species)
        hits = ncortho(mirnas, models, taxon_out, msl, cpu, query, cm_cutoff, ref_blast_db)
        hit_dict[query_species] = hits
    elif (
        os.path.isfile(query)
        and query.split('.')[-1] == 'txt'
    ):
        print('Found a textfile as the query input.\n'
              'Searching for paths to FASTA-files in {}\n'.format(query))
        with open(query, 'r') as fasta_list:
            query_list = fasta_list.read().splitlines()
            isfasta = []
            for query_path in query_list:
                if (
                    os.path.isfile(query_path)
                    and query_path.split('.')[-1] == 'fa'
                ):
                    isfasta.append(query_path)
                else:
                    print('{} is not a valid FASTA file.\n'
                          'Skipping..'.format(query_path))
                    continue
            print('Found {} paths to valid FASTA files in the query input file\n'.format(len(isfasta)))

            for fasta_path in isfasta:
                query_species = fasta_path.split('/')[-1].split('.')[0]
                taxon_out = '{}/{}'.format(output, query_species)
                print('### Starting search for miRNAs in\n'
                      '### {}\n'.format(query_species))
                hits = ncortho(mirnas, models, taxon_out, msl, cpu, fasta_path, cm_cutoff, ref_blast_db)
                hit_dict[query_species] = hits

    elif os.path.isdir(query):
        print('Query genomes given as directory.\n'
              'Searching for FASTA files in {}'.format(query))
        dir_files = glob.glob(query + '/*')
        isfasta = [file for file in dir_files if file.split('.')[-1] == 'fa']
        print('Found {} valid FASTA files in the query input directory\n'.format(len(isfasta)))
        for fasta_path in isfasta:
            query_species = fasta_path.split('/')[-1].split('.')[0]
            taxon_out = '{}/{}'.format(output, query_species)
            print('### Starting search for miRNAs in\n'
                  '### {}\n'.format(query_species))
            hits = ncortho(mirnas, models, taxon_out, msl, cpu, fasta_path, cm_cutoff, ref_blast_db)
            hit_dict[query_species] = hits
    else:
        print('No valid query found. Exiting..')
        sys.exit()

    # write output
    summary_out = '{}/{}_results.fa'.format(output, mirnas.split('/')[-1].split('.')[0])
    print('# Writing summary output at {}\n'.format(summary_out))
    with open(summary_out, 'w') as of:
        of.write('# Taxon\tmiRNA\tStart in query genome\tEnd in query genome\tstrand\n')
        for taxon in hit_dict:
            mir_dict = hit_dict[taxon]
            for mirna in mir_dict:
                seq = list(mir_dict[mirna])[0]
                info = list(list(mir_dict[mirna])[1])
                info[0] = info[0].split('_')[0]
                info = '\t'.join(info[:-1])
                header = '>{}\t{}'.format(taxon, info)
                of.write(header + '\n')
                of.write(seq + '\n')


    print('### ncOrtho has finished!\n')
    # cleanup
    if cleanup:
        print('# Cleaning up..\n')
        for taxon in hit_dict:
            taxon_out = '{}/{}'.format(output, taxon)
            clean_cmd = 'rm -r {}'.format(taxon_out)
            sp.run(clean_cmd, shell=True)

if __name__ == "__main__":
    main()
