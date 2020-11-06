import multiprocessing as mp


def findOMAprefix_fromEnsemblGTF(gtf, mapfile, repeats=5):
    oma_id = ''
    stop = 0
    tested_ensembl = []

    with open(gtf, 'r') as file:
        while (
                not oma_id
                and stop <= repeats
        ):
            for line in file:
                line = line.strip()
                if (
                    not line.startswith('#')
                    and line.split('\t')[2] == 'CDS'
                ):
                    ensembl_id = line.split('\t')[8].split(';')[0].split(' ')[1].replace('"','')
                    if not ensembl_id in tested_ensembl:
                        tested_ensembl.append(ensembl_id)
                        break
            with open(mapfile, 'r') as mf:
                for line in mf:
                    if ensembl_id in line:
                        oma_id = line.split('\t')[0][0:5]
            stop += 1

    if oma_id:
        if stop == 1:
            print('Needed to test {} ensembl gene id to find a match in the OMA database for {}'
                  .format(stop, gtf.split('/')[-1]))
            return oma_id
        else:
            print('Needed to test {} ensembl gene ids to find a match in the OMA database for {}'
                  .format(stop, gtf.split('/')[-1]))
            return oma_id
    else:
        print('\n## WARNING##')
        print('ensembl2oma parser searched for {} different ensembl gene IDs in the OMA database for:\n'
              '## {} ##\n'
              'You can increase the number different gene IDs to look for with the (repeats=) flag.\n'
              'However, the input species might not be represented in the OMA database.\n'
              'In this case, please run the OMA standalone software, available at:\n'
              'https://omabrowser.org/standalone/#downloads\n'
              .format(repeats, gtf.split('/')[-1]))

#gtf = '/share/project/felixl/ncOrtho/data/mouse_ref_core/reference_data/Mus_musculus.GRCm38.101.gtf'
#mapfile = '/share/project/felixl/ncOrtho/data/oma/oma-ensembl.txt'
#ensembl2oma(gtf, mapfile)

def map_ensembl2OMA(ensembl_ids, ens2oma_map):
    #output = mp.Queue()
    # Map ensembl gene id to oma id
    out_dict = {}
    not_found = []
    c = 0
    myrange = list(range(0, 110, 1))

    with open(ens2oma_map, 'r') as m_file:
        map_data = m_file.readlines()
        perc_done = c / len(ensembl_ids)
        milestone = myrange.pop(0)
        for e_id in ensembl_ids:
            perc_done = c / len(ensembl_ids)
            # print(e_id)
            for line in map_data:
                if e_id in line:
                    oma_id = line.split('\t')[0]
                    out_dict[e_id] = oma_id
            c += 1
            if perc_done >= milestone:
                print('Progress: {}%'.format(perc_done))
                milestone = myrange.pop(0)
    return(out_dict)
    #output.put(out_dict)