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
import pandas as pd


def main(reference, core, output, jobs=4):
    # check type of query: single FASTA-file, Folder with FASTA-files or txt document with paths to FASTA files
    isfasta = []
    if (
            os.path.isfile(core)
            and core.split('.')[-1] in ('fa', 'faa')
    ):
        print('Core proteoms given as single FASTA-file. Starting anaysis\n')
        isfasta.append(core)
    elif (
            os.path.isfile(core)
            and core.split('.')[-1] == 'txt'
    ):
        print('Found a textfile as the core proteome input.\n'
              'Searching for paths to FASTA-files in {}\n'.format(core))
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
            print('Found {} paths to valid FASTA files in the core proteome'
                  ' input file\n'.format(len(isfasta)))
    elif os.path.isdir(core):
        print('Core proteoms given as directory.\n'
              'Searching for FASTA files in {}'.format(core))
        dir_files = glob.glob(core + '/*')
        isfasta = [
            file for file in dir_files if file.split('.')[-1] in ('fa', 'faa')
        ]
        print('Found {} valid FASTA files in the query input directory\n'
              .format(len(isfasta)))
    else:
        print('No valid protein sequences in FASTA format found. Exiting..')
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
        # TODO: DEVELOPMENT ONLY
        # result_file = run_oma(reference, core, output, cpu=jobs)
        result_file = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/' \
                      'GCF_000737145_1_ASM73714v1_protein-GCA_000015425_1_ASM1542v1_protein.txt'
        headers = ['No1', 'No2', 'id1', 'id2', 'type', 'group']
        df = pd.read_csv(result_file, comment='#', names=headers, sep='\t')
        new_head = ['id1', 'id2', 'type', 'group']
        out = df.filter(items=new_head)
        print(out.head())


# # YOUR INPUT HERE
ref_proteome = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/ref_proteome/GCF_000737145.1_ASM73714v1_protein.faa'
core_proteome = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/proteomes/GCA_000015425_1_ASM1542v1_protein.fa'
output = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma'

main(ref_proteome, core_proteome, output)