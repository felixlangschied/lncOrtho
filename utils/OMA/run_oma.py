"""
Created on Tue Nov 24 09:36:03 2020

@author: felixl

Run OMA standalone

Will overwrite files during preprocessing
if the output and input directory are the same

Parameters for the OMA search can be supplied with
parameters=/path/to/OMA/parameters.drw

Per default, the script can be run with parameters='pairwise_only' to only
save the pairwise orthologs and to delete all temporary files (e.g. for ncOrtho)
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
            print('# Done')
        elif path_stem != output and not new_name:
            cmd = 'cp {} {}/{}'.format(path, output, fname)
            try:
                sp.run(cmd, shell=True, check=True, stderr=sp.PIPE)
                print('# Done')
            except sp.CalledProcessError:
                print('# Could not copy to {}'.format(output))
        elif path_stem != output and new_name:
            cmd = 'cp {} {}/{}'.format(path, output, new_name)
            try:
                sp.run(cmd, shell=True, check=True, stderr=sp.PIPE)
                print('# Done')
            except sp.CalledProcessError:
                print('# Could not copy to {}'.format(output))


# run OMA standalone pairwise orthologs only
# between two proteomes given in FASTA format
def run_oma(
        query1, query2, output, parameters='pairwise_only', out_name='OMA_out'
):
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
            cmd = ('cp {}/nco_parameters.drw {}/parameters.drw'
                   .format(curr_path, tmp_out))
            sp.run(cmd, shell=True)
            # run OMA
            os.chdir(tmp_out)
            # sp.run('oma', shell=True)
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

        elif os.path.isfile(parameters):
            cmd = 'cp {} {}/parameters.drw'.format(parameters, tmp_out)
            sp.run(cmd, shell=True)
            # run OMA
            os.chdir(tmp_out)
            sp.run('oma', shell=True)
        else:
            print('# No valid path to OMA parameters file supplied')
            sys.exit()
        print('# Finished')
    else:
        print('# Invalid input. Please enter two valid FASTA files')