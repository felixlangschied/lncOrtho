#c_gtf = '/share/project/felixl/ncOrtho/data/mouse_ref_core/core_annotation'
c_gtf = '/home/felixl/Desktop/tmp/gft_input.txt'
ens2oma_map = '/share/project/felixl/ncOrtho/data/oma/oma-ensembl.txt'
oma_groups = '/share/project/felixl/ncOrtho/data/oma/oma-groups.txt'
r_gtf = '/share/project/felixl/ncOrtho/data/mouse_ref_core/reference_data/Mus_musculus.GRCm38.101.gtf'

import os
import glob
import multiprocessing as mp
import numpy as np

from oma_lib.ensembl2oma import findOMAprefix_fromEnsemblGTF
from oma_lib.ensembl2oma import map_ensembl2OMA

def oma_parser(r_gtf, c_gtf, ens2oma_map, oma_groups):
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


    core_oma_ids = []
    no_oma_id = []
    for core in core_paths:
        core_oma = findOMAprefix_fromEnsemblGTF(core, ens2oma_map)
        if core_oma == None:
            no_oma_id.append(core)
        else:
            core_oma_ids.append(core_oma)
    ref_oma = findOMAprefix_fromEnsemblGTF(ref_path, ens2oma_map)

    # Body
    if no_oma_id:
        print('No OMA id found for:'
              '{}\n'.format('\n'.join(no_oma_id)))
        print('Starting to map reference Proteins to OMA orthologs of the remaining core species')
    else:
        print('\nStarting to extract all ensembl gene ids from the reference species')

    # Extract ensembl gene ids from reference GTF
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
    print('\nFound {} unique ensembl gene ids (Type: CDS) in {}'.format(len(ensembl_ids), r_gtf.split('/')[-1]))
    print('\nStarting to map ensembl gene ids to the OMA ids. This might take a while')

    #map list of ensembl gene ids in the reference species to the OMA id

    # maybe implement multiprocessing
    # output = mp.Queue()
    # num_workers = 3
    # chunks = np.array_split(list(ensembl_ids), num_workers)
    # processes = [mp.Process(target=map_ensembl2OMA, args=(chunks[x], ens2oma_map, output)) for x in range(num_workers)]
    # for p in processes:
    #     p.start()
    # for p in processes:
    #     p.join()
    # results = [output.get() for p in processes]
    # print(results)

    map_dict = map_ensembl2OMA(ensembl_ids, ens2oma_map)

    output = []
    with open('ensembl2oma.txt','w') as file:
        for key, value in map_dict.items():
            line = key + '\t' + value
            output.append(line)
            file.write('\n'.join(output))





    #print(list(map_dict.keys())[0:10])
    #print('Did not find orthologs for {} ensembl ids'.format(len(not_found)))






    print(core_oma_ids)
    print(no_oma_id)
    print(ref_oma)


oma_parser(r_gtf, c_gtf, ens2oma_map, oma_groups)


#
# print('this is ref_oma')
# print(ref_oma)
# print('this is core_ome')
# print(core_oma)
# for value in core_oma:
#     print(value)
#
#
#
# if core_oma:
#     print('Continuing to build files of OMA pairwise orthologs between\n'
#           'the reference species and each core species represented in the OMA database')
