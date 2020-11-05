c_gtfs = '/share/project/felixl/ncOrtho/data/mouse_ref_core/core_annotation'
#c_gtfs = '/home/felixl/Desktop/tmp/gft_input.txt'
map = '/share/project/felixl/ncOrtho/data/oma/oma-ensembl.txt'
r_gtf = '/share/project/felixl/ncOrtho/data/mouse_ref_core/reference_data/Mus_musculus.GRCm38.101.gtf'

import os
import sys

from oma_lib.parse_coreGTF import coreGTF_parser


def refGTF_parser(r_gtf, map):
    if (
        os.path.isfile(r_gtf)
        and r_gtf.split('.')[-1] == 'gtf'
    ):
        with open(gtf, 'r') as file:
            head = [next(file) for x in range(c)]
            ensembl_id = head[c - 1].split('\t')[-1].split(';')[0].split(' ')[1].replace('"', '')

    else:
        print('No valid GTF file (<.gtf>) found at\n'
              '{}'.format(r_gtf))
        sys.exit()




# Body

core_oma = coreGTF_parser(c_gtfs, map)




if core_oma:
    print('Continuing to build files of OMA pairwise orthologs between\n'
          'the reference species and each core species represented in the OMA database')
