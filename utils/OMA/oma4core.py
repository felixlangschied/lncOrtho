#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 13:44:22 2020

@author: felixl
"""
import os
import glob
import sys
import subprocess as sp
from run_oma import run_oma


def main(reference, core, output, jobs=4, cleanup=True):
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
    for core in isfasta:
        print('# Starting pre-processing for {}'.format(core))
        sys.stdout.flush()
        result_file = run_oma(reference, core, output, cpu=jobs)
        #result_file = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/GCF_000737145_1_ASM73714v1_protein-GCA_000015425_1_ASM1542v1_protein.txt'


        process_out = result_file.split('.')[0] + '_mod.txt'
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
            cmd = 'mv {} {}'.format(process_out, result_file)
            sp.run(cmd, shell=True)
            print('# OMA run and post-processing done. Output saved at:\n'
                  '{}'.format(result_file))
        else:
            print('# OMA run and post-processing done. Output saved at:\n'
                  '{}'.format(process_out))



# # YOUR INPUT HERE
ref_proteome = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/ref_proteome/GCF_000737145.1_ASM73714v1_protein.faa'
core_proteome = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/proteomes/GCA_000015425_1_ASM1542v1_protein.fa'
output = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/pairwise'

main(ref_proteome, core_proteome, output)