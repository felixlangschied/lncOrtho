#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate pairwise orthologs with OMA standalone for the core set construction of ncOrtho

Requires a working installation of OMA standalone available here:
https://omabrowser.org/standalone/

USAGE:



Created on Wed Nov 25 13:44:22 2020
@author: felixl
"""
import os
import glob
import sys
import subprocess as sp
import argparse
import multiprocessing as mp
from run_oma import run_oma


def oma4core(reference, isfasta, output, jobs=4, cleanup=True):
    for core in isfasta:
        print('# Starting pre-processing for {}'.format(core))
        sys.stdout.flush()
        result_file = run_oma(reference, core, output, cpu=jobs)
        #result_file = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/GCF_000737145_1_ASM73714v1_protein-GCA_000015425_1_ASM1542v1_protein.txt'

        #process_out = result_file.split('.')[0] + '_mod.txt'
        process_out = '{}/Refvs{}.txt'.format(output, core.split('/')[-1].split('.')[0])
        with open(result_file, 'r') as file, open(process_out, 'w') as outfile:
            example = True
            for line in file:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    out_list = line[2:]
                    if ' ' in out_list[0] or ' ' in out_list[1]:
                        if example:
                            print('# Space detected in ProteinID. Trying to remove description')
                            print('# Printing Example:')
                            print('# Your input IDs look like this:')
                            print('\t'.join(out_list[0:2]))
                            # Process ids
                            out_list = [entry.split(' ')[0] for entry in out_list]
                            print('# Now they look like this:')
                            print('\t'.join(out_list[0:2]))
                            print('# Please double check this processing!')
                            example = False
                        else:
                            out_list = [entry.split(' ')[0] for entry in out_list]
                    outstring = '\t'.join(out_list)
                    outfile.write('{}\n'.format(outstring))
        if cleanup:
            cmd = 'rm {}'.format(result_file)
            sp.run(cmd, shell=True)
        print('# OMA run and post-processing done. Output saved at:\n'
              '{}'.format(process_out))


# Allow boolean argument parsing
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        prog='python oma4core.py',
        description='Calculate pairwise orthologs with OMA standalone for the core set construction of ncOrtho',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-r', '--reference', metavar='<.fa>', type=str,
        help='Path to the proteome of the reference species in FASTA format',
    )
    parser.add_argument(
        '-c', '--core', metavar='<PATH>', type=str,
        help="Path to the proteome of a core species in FASTA format, "
             "a directory with the core species' proteomes, "
             "or a .txt file containing paths to the core species proteomes in each line"
    )
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='Output path'
    )
    parser.add_argument(
        '-x', '--cleanup', type=str2bool, metavar='True/False', nargs='?', const=True, default=True,
        help="Delete un-post-processed output file of OMA pairwise (Default=True)"
    )
    parser.add_argument(
        '-t', '--threads', type=int, metavar='<int>', nargs='?', const=4, default=4,
        help="Delete un-post-processed output file of OMA pairwise (Default=True)"
    )

    # Parse arguments
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if args.threads > available_cpu:
        print(
            '# Error: The provided number of CPU cores is higher than the '
            'number available on this system. Exiting...'
        )
        sys.exit(1)
    else:
        jobs = args.threads

    reference = args.reference
    core = args.core
    output = args.output
    cleanup = args.cleanup

    # Argument checks
    # check type of query: single FASTA-file, Folder with FASTA-files or txt document with paths to FASTA files
    isfasta = []
    if (
            os.path.isfile(core)
            and core.split('.')[-1] in ('fa', 'faa')
    ):
        print('# Core proteoms given as single FASTA-file.\n')
        isfasta.append(core)
    elif (
            os.path.isfile(core)
            and core.split('.')[-1] == 'txt'
    ):
        print('# Found a textfile as the core proteome input. Searching for paths to FASTA-files in\n'
              '{}'.format(core))
        with open(core, 'r') as fasta_list:
            query_list = fasta_list.read().splitlines()
            for query_path in query_list:
                if (
                        os.path.isfile(query_path)
                        and query_path.split('.')[-1] in ('fa', 'faa')
                ):
                    isfasta.append(query_path)
                else:
                    print('{} is not a valid FASTA file.\n'
                          'Skipping..'.format(query_path))
                    continue
            print('# Found {} paths to valid FASTA files in the core proteome'
                  ' input file\n'.format(len(isfasta)))
    elif os.path.isdir(core):
        print('# Core proteoms given as directory. Searching for FASTA files in\n'
              '{}'.format(core))
        dir_files = glob.glob(core + '/*')
        isfasta = [
            file for file in dir_files if file.split('.')[-1] in ('fa', 'faa')
        ]
        print('# Found {} valid FASTA files in the query input directory\n'
              .format(len(isfasta)))
    else:
        print('# No valid protein sequences in FASTA format found. Exiting..')
        sys.exit()

    # Check if output folder exists, create it otherwise
    if not os.path.isdir(output):
        print('# Creating output folder')
        cmd = 'mkdir {}'.format(output)
        try:
            sp.run(cmd, shell=True, check=True, stderr=sp.PIPE)
        except sp.CalledProcessError:
            print('# Could not create output folder at:\n'
                  '{}\n')

    # BODY
    oma4core(reference, isfasta, output, jobs=4, cleanup=True)

if __name__ == "__main__":
    main()


# # # YOUR INPUT HERE
# ref_proteome = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/ref_proteome/GCF_000737145.1_ASM73714v1_protein.faa'
# core_proteome = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/proteomes/GCA_000015425_1_ASM1542v1_protein.fa'
# output = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/pairwise'
#
# main(ref_proteome, core_proteome, output)