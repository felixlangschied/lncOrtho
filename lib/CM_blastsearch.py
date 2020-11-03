import os
import subprocess as sp
import glob
import sys

# blast_search: Perform a reverse BLAST search in the reference genome for a
# candidate.
# s: cmsearch result
# out: output of the main() script
# o: output for the blastsearch
# c: number of threads
# db: path to reference genome (as fasta or as blastDB)
def blast_search(s, db, out, o, c):

    # test if reference DB exists or has to be created
    file_extensions = ['.nhr', '.nin', '.nsq']
    out_data = out + '/data'
    fname = db.split('/')[-1]
    db_name = fname.replace('.fa', '')
    db_files = []
    if (
            os.path.isfile(db)
            and fname.split('.')[-1] == 'fa'
    ):
        print(
            'Reference given as FASTA file, testing if BlastDB exists in {}'
                .format(out_data)
        )
        if not os.path.isdir(out_data):
            mkdir_cmd = 'mkdir {}'.format(out_data)
            sp.call(mkdir_cmd, shell=True)
        for fe in file_extensions:
            db_files.append(glob.glob(out_data + '/' + db_name + '.*' + fe))
        if not [] in db_files:
            print("BLAST database for the reference species found. Starting re-blast")
        else:
            print(
                'BLAST database for the reference species does not exist.\n'
                'Constructing BLAST database.'
            )
            # At least one of the BLAST db files is not existent and has to be
            # created.
            db_command = 'makeblastdb -in {} -dbtype nucl -out {}/{}'.format(db, out_data, fname)
            sp.call(db_command, shell=True)
    else:
        for fe in file_extensions:
            db_files.append(glob.glob(out_data + '/' + fname + '.*' + fe))
        if not [] in db_files:
            print("Reference given as blastDB. Starting re-blast")
        else:
            print(
                'BLAST database for the given reference does not exist.\n'
                'Trying to construct a BLAST database from {}.'
                .format(db)
            )
            try:
                db_command = 'makeblastdb -in {} -out {}/{} -dbtype nucl'.format(db, out_data, db_name)
                sp.call(db_command, shell=True)
            except:
                print(
                    'Was not able to construct the BLASTdb for {}.\n'
                    'Please check you input.'
                    .format(db)
                )
                sys.exit()

    # print('db = {}'.format(db))
    # print('s = {}'.format(s))
    # print('o = {}'.format(o))
    # print('c = {}'.format(c))
    # print('out = {}'.format(out))
    # print('out_data = {}'.format(out_data))
    # print('{}/{}'.format(out_data, fname))
    os.chdir(out_data)
    blast_command = (
        'blastn -task blastn -db {0}/{1} -query {2} '
        '-out {3} -num_threads {4} -outfmt 6'.format(out_data, db_name, s, o, c)
    )
    sp.call(blast_command, shell=True)
