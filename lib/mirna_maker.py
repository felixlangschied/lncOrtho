import os
import subprocess as sp


# Central class of microRNA objects
class Mirna(object):
    def __init__(self, name, chromosome, start, end, strand):
        # miRNA identifier
        self.name = name
        # chromosome that the miRNA is located on
        # downloads from miRBase contain "chr" prefix which leads to inconsistency
        self.chromosome = chromosome
        # start position of the pre-miRNA
        self.start = int(start)
        # end position of the pre-miRNA
        self.end = int(end)
        # sense (+) or anti-sense (-) strand
        self.strand = strand

    def loadSeq(self, seq, seq_type):
        if seq_type == 'pre':
            self.pre = seq
        elif seq_type == 'mat':
            self.mat = seq


# TODO: include both mature strands, 5p and 3p, aka mature and star


# mirna_maker: Parses the miRNA data input file and returns a
#              dictionary of Mirna objects.
# Arguments:
# mirpath: path to file with microRNA data
# cmpath: path to covariance models
# output: path for writing temporary files
def mirna_maker(mirpath, cmpath, output, msl):
    mmdict = {}  # will be the return object

    with open(mirpath) as mirna_file:
        mirna_data = [
            line.strip().split() for line in mirna_file
            if not line.startswith('#')
        ]

    for mirna in mirna_data:
        mirid = mirna[0]
        # Check if the output folder exists, otherwise create it.
        if not os.path.isdir('{}/{}'.format(output, mirid)):
            try:
                mkdir = 'mkdir -p {}/{}'.format(output, mirid)
                sp.run(mkdir, shell=True, check=True, stderr=sp.PIPE)
            except sp.CalledProcessError:
                print(
                    '# Cannot create output folder for {}.'
                    'Skipping to next miRNA.'
                )
                continue

        # check if 'chr' prefix from miRBase exists and delete it
        if 'chr' in mirna[1]:
            mirna[1] = mirna[1].replace('chr', '')

        # Obtain the reference bit score for each miRNA by applying it
        # to its own covariance model.
        print('# Calculating reference bit score for {}.'.format(mirid))
        seq = mirna[5]
        query = '{0}/{1}/{1}.fa'.format(output, mirid)
        model = '{0}/{1}.cm'.format(cmpath, mirid)

        # Check if the covariance model even exists, otherwise skip to
        # the next miRNA.
        if not os.path.isfile(model):
            print('# No covariance model found for {}.\n'
                  'Skipping miRNA..'.format(mirid))
            continue

        # Create a temporary FASTA file with the miRNA sequence as
        # query for external search tool cmsearch to calculate
        # reference bit score.
        with open(query, 'w') as tmpfile:
            tmpfile.write('>{0}\n{1}'.format(mirid, seq))
        cms_output = '{0}/{1}/cmsearch_{1}_tmp.out'.format(output, mirid)
        cms_log = '{0}/{1}/cmsearch_{1}.log'.format(output, mirid)
        cms_command = (
            'cmsearch -E 0.01 --noali -o {3} --tblout {0} {1} {2}'
            .format(cms_output, model, query, cms_log)
        )
        sp.call(cms_command, shell=True)
        with open(cms_output) as cmsfile:
            hits = [
                line.strip().split() for line in cmsfile
                if not line.startswith('#')
            ]
            if hits:
                top_score = float(hits[0][14])
            # In case of any issues occuring in the calculation of the bit
            # score, no specific threshold can be determined. The value will
            # set to zero, which technically turns off the filter.
            else:
                print(
                    '# Warning: Self bit score not applicable, '
                    'setting threshold to 0.'
                )
                top_score = 0.0

        tmp_mir = Mirna(*mirna[0:5])
        tmp_mir.bit = top_score
        # add sequences
        tmp_mir.loadSeq(mirna[5], 'pre')
        if (
                len(mirna) > 6
                and mirna[6] != '.'
        ):
            tmp_mir.loadSeq(mirna[6], 'mat')
        else:
            tmp_mir.loadSeq(None, 'mat')

        # Create output.
        mmdict[mirna[0]] = tmp_mir
        # Remove temporary files.
        for rmv_file in [cms_output, cms_log, query]:
            sp.call('rm {}'.format(rmv_file), shell=True)

        #print('Reference Bit-Score determined as: {}'.format(top_score))
    return mmdict
