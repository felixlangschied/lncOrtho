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
        '-r', '--reference', metavar='<.fa>', type=str,
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
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

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

    # TODO: include checks for validity of arguments
    # os.getcwd()
    # os.chdir(path)
    # os.path.exists(path)
    # os.path.isfile(path)
    # os.path.isdir(path)

    mirnas = args.ncrna
    models = args.models
    output = args.output
    query = args.query
    reference = args.reference
    cm_cutoff = args.cutoff
    # Not in use yet
    msl = args.msl
    # blast_cutoff = args.blastc

    ##################################
    # Checking Validity of Arguments #
    ##################################

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
            'Reference given as FASTA file, testing if BlastDB exists in {}'
                .format(out_data)
        )
        if not os.path.isdir(out_data):
            mkdir_cmd = 'mkdir {}'.format(out_data)
            sp.call(mkdir_cmd, shell=True)
        for fe in file_extensions:
            db_files.append(glob.glob(out_data + '/' + db_name + '.*' + fe))
        if not [] in db_files:
            print("BLAST database for the reference species found.")
        else:
            print(
                'BLAST database for the reference species does not exist.\n'
                'Constructing BLAST database.'
            )
            # At least one of the BLAST db files is not existent and has to be
            # created.
            db_command = 'makeblastdb -in {} -dbtype nucl -out {}/{}'.format(reference, out_data, fname)
            sp.call(db_command, shell=True)
    else:
        for fe in file_extensions:
            db_files.append(glob.glob(out_data + '/' + fname + '.*' + fe))
        if not [] in db_files:
            print("Reference given as blastDB. Starting analysis")
        else:
            print(
                'BLAST database for the given reference does not exist.\n'
                'Trying to construct a BLAST database from {}.'
                .format(reference)
            )
            try:
                db_command = 'makeblastdb -in {} -out {}/{} -dbtype nucl'.format(reference, out_data, db_name)
                sp.call(db_command, shell=True)
            except:
                print(
                    'Was not able to construct the BLASTdb for {}.\n'
                    'Please check you input.'
                    .format(reference)
                )
                sys.exit()
    ref_blast_db = out_data + '/' + db_name

    ##############################################################################
    # TODO: implement multiple query files
    # check how many query files are given
    # if (
    #         os.path.isfile(query)
    #     and query.endswith
    # )



    hits = ncortho(mirnas, models, output, msl, cpu, query, cm_cutoff, ref_blast_db)

    print(hits)


if __name__ == "__main__":
    main()
