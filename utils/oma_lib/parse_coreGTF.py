import os
import glob
import sys


# c_gtf = path to directory with the annotation of the core spezies/ path to list with gtf annotations
def getOMAid(c_gtf, r_gtf, mapfile, repeats=20):
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
        print('Found reference species annotation file. Adding it to the search.')
        core_paths.append(r_gtf)
    else:
        print('No valid GTF file (<.gtf>) found at\n'
              '{}'.format(r_gtf))

    # open GTFs and load the ensemblIDs in a list
    # omaIDs = set()
    # foundfiles = set()
    ensembl_ids = []
    out_dict = {}
    c = 6
    stop = 0
    while (
            len(out_dict) < len(core_paths)
            and stop <= repeats
    ):
        for gtf in core_paths:
            if gtf not in out_dict.keys():
                with open(gtf, 'r') as file:
                    head = [next(file) for x in range(c)]
                    ensembl_id = head[c - 1].split('\t')[-1].split(';')[0].split(' ')[1].replace('"', '')
                    #ensembl_ids.append(ensembl_id)
                # TODO: Don't go into the mapfile multiple times with the same id (ensembl_ids are multiple times in the GTF)
                # map ensembl-ids to OMA-ids
                with open(mapfile, 'r') as mf:
                    for line in mf:
                        if ensembl_id in line:
                            oma_id = line.split('\t')[0][0:5]
                            out_dict[gtf] = oma_id
                            # omaIDs.add(oma_id)
                            # foundfiles.add(gtf)
        if len(out_dict) < len(core_paths):
            #print('Found only {} of {} omaIDs'.format(len(out_dict), len(core_paths)))
            #print('Checking if the next protein in the core species gtf file is in an OMA group')
            c += 1
            stop += 1
        print(stop)

    # check if reference species is part of the OMA database
    if r_gtf in out_dict:
        print('Reference species is part of the OMA database. The OMA id is:\n'
              '{}'.format(out_dict[r_gtf]))
        #return out_dict.pop(r_gtf)

    else:
        print('Reference species is not part of the OMA database for:\n'
              '{}'
              'Please run OMA standalone to identify the pairwise orthologs\n'
              'of all proteins between the reference species and all core species, respectively.\n'
              'OMA standalone is available at: https://omabrowser.org/standalone/#downloads\n'
              'Exiting..')
        sys.exit()

    # give feedback if all core species are in the oma database
    missing = []
    for gtf in core_paths:
        if gtf not in out_dict:
            missing.append(gtf)

    if not missing:
        print('\nFound the omaID of all core species')
    else:
        print('\nOnly these {} omaIDs were found until time-out:'.format(len(out_dict)))
        print('\n{} of your core species is not in the OMA database.\n'
              'Please run OMA standalone to identify the pairwise orthologs\n'
              'of all proteins of the missing core species and your reference species.\n'
              'OMA standalone is available at: https://omabrowser.org/standalone/#downloads\n'
              .format(len(core_paths) - len(out_dict)))
        print('You are missing pairwise orthologs for these files:')
        for miss in missing:
            print(miss)
    return out_dict.pop(r_gtf), out_dict.values()
