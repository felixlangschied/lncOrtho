#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 09:36:03 2020

@author: felixl

Pre-process files and run OMA standalone.
Requires a working installation of OMA standalone available here:
https://omabrowser.org/standalone/

USAGE:
    outfile_name = run_oma(
        input1.faa, input2.faa, output, [parameters], [out_name], [cpu]
        )
Returns the path to the output
(path to pairwise orthologs when run in 'pairwise_only' mode)

Will overwrite files during preprocessing
if the output and input directory are the same

OPTIONAL:
Parameters for the OMA search can be supplied with
parameters=/path/to/OMA/parameters.drw

Per default, the script can be run with parameters='pairwise_only' to only
save the pairwise orthologs and to delete all temporary files (e.g. for ncOrtho)

You can can name the OMA output directory using the parameter name='yourNameHere'

Number of jobs to start with cpu=int (Default=4)

Compatibility tested for OMA 2.4.1

"""
import os
import glob
import sys
import subprocess as sp


def process_input(fasta_list, output):
    for path in fasta_list:
        new_name = ''
        # test if filename contains '.'
        fname = path.split('/')[-1]
        if fname.count('.') > 1:
            print('# Removing dots from filename: {}'.format(fname))
            new_name = '_'.join(fname.split('.')[0:-1])
            if fname.split('.')[-1] != 'fa':
                print('# {} file extension not supported by OMA.'
                      ' Renaming it to fa'.format(fname.split('.')[-1]))
            new_name = new_name + '.fa'
        elif fname.split('.')[-1] != 'fa':
            print('# {} file extension not supported by OMA.'
                  ' Renaming it to fa'.format(fname.split('.')[-1]))
            new_name = fname.split('.')[0] + '.fa'

        path_stem = '/'.join(path.split('/')[0:-1])
        if path_stem == output and new_name:
            os.rename(path, '{}/{}'.format(path_stem, new_name))
        elif path_stem != output and not new_name:
            cmd = 'cp {} {}/{}'.format(path, output, fname)
            try:
                sp.run(cmd, shell=True, check=True, stderr=sp.PIPE)
            except sp.CalledProcessError:
                print('# Could not copy to {}'.format(output))
        elif path_stem != output and new_name:
            cmd = 'cp {} {}/{}'.format(path, output, new_name)
            try:
                sp.run(cmd, shell=True, check=True, stderr=sp.PIPE)
            except sp.CalledProcessError:
                print('# Could not copy to {}'.format(output))


# run OMA standalone between two proteomes given in FASTA format
# per default only the parwise orthologs will be calculated and
# the output will be processed to look like the downloadable content from
# OMA browser
def run_oma(
        query1, query2, output, parameters='pairwise_only', out_name='OMA_out',
        cpu=4):
    if (
            os.path.isfile(query1) and os.path.isfile(query2)
            and query1.split('.')[-1] in ('fa', 'faa')
            and query2.split('.')[-1] in ('fa', 'faa')
    ):
        fastas = [query1, query2]
        # Create DB directory
        DB = '{}/{}/DB'.format(output, out_name)
        tmp_out = '{}/{}'.format(output, out_name)
        cmd = 'mkdir -p {}'.format(DB)
        sp.run(cmd, shell=True)
        # preprocess input and save in DB directory
        process_input(fastas, DB)
        # move modified oma parameters file to tmp_out
        if parameters == 'pairwise_only':
            curr_path = os.path.dirname(os.path.realpath(__file__))
            paramfile = '{}/nco_parameters.drw'.format(curr_path)
            sp.run('cp {} {}/parameters.drw'.format(paramfile, tmp_out), shell=True)
            # run OMA
            os.chdir(tmp_out)
            try:
                print('# Starting OMA standalone run')
                #sp.run('/home/felixl/applications/OMA.2.3.1/bin/oma {} -n {}'.format(paramfile, cpu), shell=True, check=True)
                sp.run('oma -n {}'.format(cpu), shell=True, check=True)
            except sp.CalledProcessError:
                print('# OMA is returning an error:')
                sys.exit()

            # sort output and remove temporary files
            outfile = glob.glob('{}/Output/PairwiseOrthologs/*'.format(tmp_out))
            if len(outfile) > 1:
                print('# Multiple OMA pairwise files detected at:'
                      '{}/Output/PairwiseOrthologs\n'
                      'Exiting..'.format(tmp_out))
                sys.exit()
            else:
                outfile = outfile[0]
            cmd = 'mv {} {}'.format(outfile, output)
            sp.run(cmd, shell=True)
            cmd = 'rm -r {}'.format(tmp_out)
            sp.run(cmd, shell=True)
            new_out = '{}/{}'.format(output, outfile.split('/')[-1])
            return new_out

        elif os.path.isfile(parameters):
            # run OMA
            sp.run('cp {} {}/parameters.drw'.format(parameters, tmp_out), shell=True)
            os.chdir(tmp_out)
            try:
                print('# Starting OMA standalone run')
                sp.run('oma -n {}'.format(cpu), shell=True, check=True)
            except sp.CalledProcessError:
                print('# OMA is returning an error:')
                sys.exit()
            outfile = glob.glob('{}/Output/PairwiseOrthologs/*'.format(tmp_out))
            return outfile
        else:
            print('# No valid path to OMA parameters file supplied')
            sys.exit()
    else:
        print('# Invalid input. Please enter two valid FASTA files')
    # print('# Finished')



# # YOUR INPUT HERE
# ref_proteome = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/ref_proteome/GCF_000737145.1_ASM73714v1_protein.faa'
# core_proteome = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma/proteomes/GCA_000015425_1_ASM1542v1_protein.fa'
# output = '/share/project/felixl/ncOrtho/data/aci_ref_core/oma'
#
# # RUN PROGRAM
# run_oma(ref_proteome, core_proteome, output)