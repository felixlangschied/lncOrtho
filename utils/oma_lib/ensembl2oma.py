import re

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
                    ensembl_id = line.split('\t')[8].split(';')[0].split(' ')[1].replace('"', '')
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


def map_ensembl2OMA(ensembl_ids, ens2oma_map):
    # Map ensembl gene id to oma id
    out_dict = {}
    map_dict = {}
    no_oma = []

    #load mapping file as dictionary
    with open(ens2oma_map, 'r') as m_file:
        for line in m_file:
            tmp_oma = []
            if not line.startswith('#'):
                (oma, ens) = line.strip().split('\t')
                ens = ens.split('.')[0]
                # for ensembl_ids with multiple oma_ids, load oma ids as list
                if ens not in map_dict:
                    map_dict[ens] = oma
                else:
                    if type(map_dict[ens]) == str:
                        tmp_oma = map_dict[ens].split(None)
                        tmp_oma.append(oma)
                    elif type(map_dict[ens]) == list:
                        tmp_oma = map_dict[ens]
                        tmp_oma.append(oma)
                    map_dict[ens] = tmp_oma

    # look vor ensembl ids in the mapping dictionary
    for e_id in ensembl_ids:
        try:
            oma_id = map_dict[e_id]
            out_dict[e_id] = oma_id
        except KeyError:
            no_oma.append(e_id)


    # return output
    return out_dict, no_oma


#####################################################################

# find core taxon ortholog

def find_core_ortholog(map_dict, ref_oma2rest_oma, taxon_suff):
    out_dict = {}

    for ens_id in map_dict:
        oma_ids = map_dict[ens_id]
        if type(oma_ids) == str:
            oma_iterator = oma_ids.split(None)
        elif type(oma_ids) == list:
            oma_iterator = oma_ids

        for oma_id in oma_iterator:
            tmp_list = []
            no_hits = []
            try:
                searchstring = str(ref_oma2rest_oma[oma_id])
            except KeyError:
                continue
            hit = re.search(taxon_suff + '\d+', searchstring)
            if hit:
                core_ortholog = hit.group()
                if core_ortholog not in out_dict:
                    out_dict[core_ortholog] = ens_id
                else:
                    print('lala')
                    tmp_list = out_dict[core_ortholog]
                    tmp_list.append(ens_id)
                    out_dict[core_ortholog] = tmp_list
                    print(out_dict[core_ortholog])
            else:
                no_hits.append(oma_id)
    return out_dict, no_hits


##########################################################################

def map_OMA2ensembl(oma_ids, ens2oma_map):
    map_dict = {}
    # load mapping file as dictionary
    with open(ens2oma_map, 'r') as m_file:
        for line in m_file:
            if not line.startswith('#'):
                (oma, ens) = line.strip().split('\t')
                ens = ens.split('.')[0]
                map_dict[oma] = ens
    return map_dict


#########################################################################

