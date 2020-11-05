import os
import glob


# c_gtf = path to directory with the annotation of the core spezies/ path to list with gtf annotations
def coreGTF_parser(c_gtf, map):
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
        with open(c_gtf, 'r') as list:
            for line in list:
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

    # open GTFs and load the ensemblIDs in a list
    omaIDs = set()
    foundfiles = set()
    c = 6
    stop = 0
    while (
            len(omaIDs) < len(core_paths)
            and stop <= 5
    ):
        for gtf in core_paths:
            if not gtf in foundfiles:
                with open(gtf, 'r') as file:
                    head = [next(file) for x in range(c)]
                    ensembl_id = head[c - 1].split('\t')[-1].split(';')[0].split(' ')[1].replace('"', '')

                # map ensembl-ids to OMA-ids
                with open(map, 'r') as mapfile:
                    for line in mapfile:
                        if ensembl_id in line:
                            oma_id = line.split('\t')[0][0:5]
                            omaIDs.add(oma_id)
                            foundfiles.add(gtf)

        if len(omaIDs) < len(core_paths):
            print('Found only {} of {} omaIDs'.format(len(omaIDs), len(core_paths)))
            print('Checking if the next protein in the core species gtf file is in an OMA group')
            c += 1
            stop += 1

    # give feedback if all core species are in the oma database
    missing = set(core_paths).difference(foundfiles)
    if len(omaIDs) == len(core_paths):
        print('\nFound the omaID of all core species')
        print(omaIDs)
    else:
        print('\nOnly these {} omaIDs were found until time-out:'.format(len(omaIDs)))
        print(omaIDs)
        print('\n{} of your core species is not in the OMA database.\n'
              'Please run OMA standalone to identify the pairwise orthologs\n'
              'of all proteins of the missing core species and your reference species.\n'
              'OMA standalone is available at: https://omabrowser.org/standalone/#downloads\n'
              .format(len(core_paths) - len(omaIDs)))
        print('You are missing pairwise orthologs for these files:')
        for miss in missing:
            print(miss)
    return(omaIDs)