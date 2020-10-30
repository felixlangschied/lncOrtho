# Construct and calibrate a covariance model for a given ncRNA core set
# alignment in Stockholm format.
# The calibration step is computationally very expensive and should be
# performed on multiple CPU cores.
# TODO: Include license, author etc.

import multiprocessing as mp
import os
import subprocess as sp
import sys


class CmConstructor(object):
    
    def __init__(self, alignment, outpath, name, cpu):
        self.alignment = alignment
        self.outpath = outpath
        self.name = name
        self.cpu = cpu
        self.model = '{0}/{1}.cm'.format(outpath, name)
    
    def construct(self):
        print('# Constructing covariance model for {}.'.format(self.name))
        cmbuild = 'cmbuild'
        construct_command = (
            '{0} -n {1} -o {2}/{1}_cmbuild.log {2}/{1}.cm {3}'
            .format(cmbuild, self.name, self.outpath, self.alignment)
        )
        sp.call(construct_command, shell=True)
        print(
            '# Writing log file to {}/{}_cmbuild.log'
            .format(self.outpath, self.name)
        )
        print('# Finished covariance model construction.')
    
    def calibrate(self):
        print('# Calibrating covariance model for {}.'.format(self.name))
        cmcalibrate = 'cmcalibrate'
        # cmcalibrate = (
        #     '/home/andreas/Applications/infernal-1.1.2-linux-intel-gcc/'
        #     'binaries/cmcalibrate'
        # )
        calibrate_command = (
            '{0} --cpu {1} {2}'.format(cmcalibrate, self.cpu, self.model)
        )
        sp.call(calibrate_command, shell=True)
        print('# Finished covariance model calibration.')


def CreateCm(alignment,output,cpu):

    
    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if cpu > available_cpu:
        print(
            '# Error: The provided number of CPU cores is higher than the '
            'number available on this system. Exiting...'
        )
        sys.exit(1)

    # Check if alignment file exists.
    if not os.path.isfile(alignment):
        print(
            '# Error: The provided path to the alignment file appears to be '
            'invalid. Exiting...'
        )
        sys.exit(1)

    name = alignment.split('/')[-1].split('.')[0]

    # Check if the output folder exists.
    if not os.path.isdir(output):
        mkdir_cmd = 'mkdir {}'.format(output)
        sp.call(mkdir_cmd, shell=True)

    # Initiate covariance model construction and calibration.
    cmc = CmConstructor(alignment, output, name, cpu)
    # Construct the model.
    cmc.construct()
    # Calibrate the model.
    cmc.calibrate()

