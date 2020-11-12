"""
This tool will write the 1:1 oma orthologs between a reference species and a set of core species to respective files.
Since the OMA groups file contains only 1:1 orthologs, extracting the OMAid-prefix with the findOMAprefix_fromEnsemblGTF
function and then using the genome pair view service on: https://omabrowser.org/oma/genomePW/ might be all you need.
"""

import os
import glob
import re

from utils.oma_lib.ensembl2oma import findOMAprefix_fromEnsemblGTF
from utils.oma_lib.ensembl2oma import map_ensembl2OMA
from utils.oma_lib.ensembl2oma import find_core_ortholog
from utils.oma_lib.ensembl2oma import map_OMA2ensembl


def oma_parser(r_gtf, c_gtf, ens2oma_map, oma_groups, output):
    # parse gtfs
    # load gtf paths into list
    core_paths = []
    if os.path.isdir(c_gtf):
        input_paths = glob.glob(c_gtf + '/*')
        for path in input_paths:
            if path.split('.')[-1] == 'gtf':
                core_paths.append(path)
            else:
                print('{} is not a valid annotation file (<.gtf>)'
                      .format(path))
                continue
    elif os.path.isfile(c_gtf):
        print('File found as core species annotation input. Searching for paths to GTF files..')
        with open(c_gtf, 'r') as fl:
            for line in fl:
                line = line.strip()
                if (
                        os.path.isfile(line)
                        and line.split('.')[-1] == 'gtf'
                ):
                    core_paths.append(line)
                else:
                    print('{} is not a valid annotation file (<.gtf>)'
                          .format(line))
                    continue
    print('Found {} core species annotation files'.format(len(core_paths)))

    if (
            os.path.isfile(r_gtf)
            and r_gtf.split('.')[-1] == 'gtf'
    ):
        print('Found reference species annotation file.')
        ref_path = r_gtf
    else:
        print('No valid GTF file (<.gtf>) found at\n'
              '{}'.format(r_gtf))

    # find ID prefix of taxa in the OMA database
    core_oma_ids = []
    no_oma_taxon = []
    for core in core_paths:
        core_oma = findOMAprefix_fromEnsemblGTF(core, ens2oma_map)
        if core_oma is None:
            no_oma_taxon.append(core)
        else:
            core_oma_ids.append(core_oma)
    ref_oma = findOMAprefix_fromEnsemblGTF(ref_path, ens2oma_map)

    if no_oma_taxon:
        print('No OMA id found for:'
              '{}\n'.format('\n'.join(no_oma_taxon)))
        print('Starting to map reference Proteins to OMA orthologs of the remaining core species')
    else:
        print('\nStarting to extract all ensembl gene ids from the reference species')

    # Extract ensembl gene ids from the reference GTF
    ensembl_ids = []
    with open(r_gtf, 'r') as r_file:
        for line in r_file:
            line = line.strip()
            if (
                    not line.startswith('#')
                    and line.split('\t')[2] == 'CDS'
            ):
                ensembl_id = line.split('\t')[8].split(';')[0].split(' ')[1].replace('"', '')
                ensembl_ids.append(ensembl_id)

    # convert list to set to get rid of duplicates
    ensembl_ids = set(ensembl_ids)
    print('\nFound {} unique ensembl gene ids (Type: CDS) in {}'
          .format(len(ensembl_ids), r_gtf.split('/')[-1]))
    print('\nStarting to map ensembl gene ids of the reference species to the respective OMA ids')

    # map ensembl ids of the reference species to the OMA ids
    map_dict, no_oma = map_ensembl2OMA(ensembl_ids, ens2oma_map)
    if no_oma:
        print('Could not find oma ids for {} ensembl ids\n'.format(len(no_oma)))
    del ensembl_ids


    print('Starting to map reference OMA ids to the OMA groups\n')


    # load dictionary with reference OMAid as key and all other OMAids
    # of the orthologous group as the values
    ref_oma2rest_oma = {}
    with open(oma_groups, 'r') as file:
        for line in file:
            line = line.strip()
            if ref_oma in line:
                hit = re.search(ref_oma + '\d+', line)
                hit = hit.group()
                line_list = line.split('\t')[1:]
                line_list.remove(hit)
                ref_oma2rest_oma[hit] = line_list

    print('Found orthologous groups for {} reference OMA ids\n'.format(len(ref_oma2rest_oma)))


    for taxon_suff in core_oma_ids:
        print('Starting to extract orthologs from the OMA groups for:\n'
              '## {} ##'.format(taxon_suff))
        core_oma2ref_ens, no_hit = find_core_ortholog(map_dict, ref_oma2rest_oma, taxon_suff)

        print('Found {} OMA ids with the {} taxon suffix.\n'.format(len(core_oma2ref_ens), taxon_suff))
        print('Found no hits for these OMAids:\n'
              '{}'.format(no_hit))

        # map OMA ids of the core species back to the ensembl IDs
        core_omas = core_oma2ref_ens.keys()
        core_map = map_OMA2ensembl(core_omas, ens2oma_map)

        final_dict = {}
        for oma_id in core_oma2ref_ens:
            key = core_oma2ref_ens[oma_id]
            value = core_map[oma_id]
            final_dict[key] = value

        #write output
        outpath = '{}/{}.tsv'.format(output, taxon_suff)
        with open(outpath, 'w') as outfile:
            for key, value in final_dict.items():
                outfile.write('{}\t{}\n'.format(key, value))

        #print(final_dict)
        print('Wrote {} pairwise orthologs between {} and the reference species to {}'
              .format(len(final_dict), taxon_suff, outpath))




# testing

#c_gtf = '/share/project/felixl/ncOrtho/data/mouse_ref_core/core_annotation'
c_gtf = '/home/felixl/Desktop/tmp/gft_input.txt'
ens2oma_map = '/share/project/felixl/ncOrtho/data/oma/oma-ensembl.txt'
oma_groups = '/share/project/felixl/ncOrtho/data/oma/oma-groups.txt'
r_gtf = '/share/project/felixl/ncOrtho/data/mouse_ref_core/reference_data/Mus_musculus.GRCm38.101.gtf'
output = '/home/felixl/PycharmProjects/ncOrtho/utils'

oma_parser(r_gtf, c_gtf, ens2oma_map, oma_groups, output)


