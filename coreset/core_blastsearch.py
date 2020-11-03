import os
import subprocess as sp
import glob
import sys

# for python versions <= 3, the blast commands need an additional keyword:
# encoding='utf8'

#find python version
version = sys.version.split()[0]
gen_version = int(version.split('.')[0])

# Perform reciprocal BLAST search and construct Stockholm alignment
def blast_search(mirna, r_path, o_path, c):
    core_set = {}

    # create blastDB for the reference score, to calculate bit score threshold
    # Check if BLAST database already exists, otherwise create it.
    # Database files are ".nhr", ".nin", ".nsq".
    file_extensions = ['.nhr', '.nin', '.nsq']
    out_data = o_path + '/data'
    if not os.path.isdir(out_data):
        mkdir_cmd = 'mkdir {}'.format(out_data)
        sp.call(mkdir_cmd, shell=True)

    db_files = []
    fname = r_path.split('/')[-1]
    db_name = fname.replace('.fa', '')
    for fe in file_extensions:
        db_files.append(glob.glob(out_data + '/' + db_name + '.*' + fe))

    if not [] in db_files:
        print("BLAST database for the reference species found")
    else:
        print(
            'BLAST database for the reference species does not exist.\n'
            'Constructing BLAST database.'
        )
        # At least one of the BLAST db files is not existent and has to be
        # created.

        db_command = (
            'makeblastdb -dbtype nucl -in {0} -out {1}/{2}'
                .format(r_path, out_data, db_name)
        )
        sp.call(db_command, shell=True)


    print("##### Starting to identify miRNAs in candidate regions #####")

    #for mirna in mirnas:
    mirid = mirna[0]
    # Ensure that output folder exists and change to this folder.
    out_folder = '{}/{}'.format(o_path, mirid)
    if not os.path.isdir(out_folder):
        mkdir_cmd = 'mkdir {}'.format(out_folder)
        sp.call(mkdir_cmd, shell=True)

    print('### {} ###'.format(mirid))
    # Coordinates of the ncRNA
    mchr = mirna[1]
    mstart = int(mirna[2])
    mend = int(mirna[3])

    # Start of reference bit score computation.
    # Convert RNA sequence into DNA sequence.
    preseq = mirna[5].replace('U', 'T')
    # The miRNA precursors can show a low level of complexity, hence it is
    # required to deactivate the dust filter for the BLAST search.
    bit_check = (
        'blastn -num_threads {0} -dust yes -task megablast -db {1}/{2} '
        '-outfmt \"6 bitscore\"'.format(c, out_data, db_name)
    )
    if gen_version == 2:
        ref_bit_cmd = sp.Popen(
            bit_check, shell=True, stdin=sp.PIPE,
            stdout=sp.PIPE, stderr=sp.PIPE
        )
    elif gen_version == 3:
        ref_bit_cmd = sp.Popen(
            bit_check, shell=True, stdin=sp.PIPE,
            stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
        )

    ref_results, err = ref_bit_cmd.communicate(preseq)
    # print(ref_results)
    ref_bit_score = float(ref_results.split('\n')[0].split('\t')[0])
    print("Reference bit score is: {}.".format(ref_bit_score))
    print("Threshold set at: {}".format(ref_bit_score * 0.5))
    # End of reference bit score computation.

    # find output of main() in coreset.py
    fasta = '{0}/{1}/{1}.fa'.format(o_path, mirid)
    #print(fasta)

    # change directory before BLAST
    #os.chdir(out_folder)
    # BLAST
    if os.path.isfile(fasta):
        print('Candidate FASTA file found for {}.'.format(mirna[0]))
        print('Checking BLAST database.')
        # Check if BLAST database already exists, otherwise create it.
        # Database files are ".nhr", ".nin", ".nsq".
        file_extensions = ['.nhr', '.nin', '.nsq']
        # file_extensions = ['.nhr']
        for fe in file_extensions:
            checkpath = '{}{}'.format(fasta, fe)
            if not os.path.isfile(checkpath):
                print(
                    'BLAST database does not exist.\n'
                    'Constructing BLAST database.'
                )
                # At least one of the BLAST db files is not existent and has to be
                # created.
                db_command = (
                    'makeblastdb -dbtype nucl -in {0} -out {0}'
                        .format(fasta)
                )
                sp.call(db_command, shell=True)
                break
        else:
            print('BLAST database already exists for {}.'.format(mirid))

        blastn_cmd = (
            'blastn -num_threads {0} -task blastn -db {1} -outfmt \"6 '
            'sseqid evalue bitscore sseq\"'.format(c, fasta)
        )
        if gen_version == 2:
            blastn = sp.Popen(
                blastn_cmd, shell=True, stdin=sp.PIPE,
                stdout=sp.PIPE, stderr=sp.STDOUT
            )
        elif gen_version == 3:
            blastn = sp.Popen(
                blastn_cmd, shell=True, stdin=sp.PIPE,
                stdout=sp.PIPE, stderr=sp.STDOUT, encoding='utf8'
            )
        results, err = blastn.communicate(preseq)

        print("\nBLAST results:")
        print(results)
        if err != None:
            print(err)

        ##### Collect best hit for each core set species if it is within the accepted bit score range
        core_dict = {}
        result_list = results.split('\n')
        if result_list:
            for hit in result_list:
                if hit:
                    hit_data = hit.split()
                    if (
                            not hit_data[0] in core_dict
                            and float(hit_data[2]) >= 0.5 * ref_bit_score
                    ):
                        core_dict[hit_data[0]] = hit

        #print(core_dict.keys())
        print('##### Starting Re-BLAST #####')

        ##### Re-BLAST #####
        print('Found BLAST-Hits above threshold for:')
        print(', '.join(core_dict.keys()))

        for species in core_dict.copy():
            print('### ' + species + ' ###')
            # Make sure to eliminate gaps
            candidate_seq = core_dict[species].split()[3].replace('-', '')
            reblastn_cmd = (
                'blastn -num_threads {0} -task blastn -db {1}/{2} -outfmt \"6'
                ' sseqid sstart send evalue bitscore\"'
                    .format(c, out_data, db_name)
            )
            if gen_version == 2:
                reblastn = sp.Popen(
                    reblastn_cmd, shell=True, stdin=sp.PIPE,
                    stdout=sp.PIPE, stderr=sp.PIPE
                )
            elif gen_version == 3:
                reblastn = sp.Popen(
                    reblastn_cmd, shell=True, stdin=sp.PIPE,
                    stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf8'
                )

            reresults, reerr = reblastn.communicate(candidate_seq)
            print('Reverse search results:')
            print(reresults.split('\n')[0])
            if reerr:
                print(reerr)

            ##### Check if reverse hit overlaps with reference miRNA
            if reresults:
                first_hit = reresults.split('\n')[0].split()
                rchr = first_hit[0]
                rstart = int(first_hit[1])
                rend = int(first_hit[2])

                if rchr == mchr:
                    print('Same chromosome.')
                    if (
                            (rstart <= mstart and mstart <= rend)
                            or (rstart <= mend and mend <= rend)
                    ):
                        print('Reciprocity fulfilled.')
                elif (
                        (mstart <= rstart and rstart <= mend)
                        or (mstart <= rend and rend <= mend)
                ):
                    print('Reciprocity fulfilled.')
                else:
                    del core_dict[species]
                    print(
                        'No reverse hit for {}. Reciprocity unfulfilled.'
                            .format(mirid)
                    )
            else:
                print('FASTA file not found for {}.'.format(mirna[0]))
                continue

        notreciproc = []
        # write output
        if core_dict:
            corefile = '{}/{}_core.fa'.format(out_folder, mirid)
            with open(corefile, 'w') as outfile:
                outfile.write('>{}\n{}\n'.format(mirid, preseq))
                print('\n##### Constructing core set for {}: #####'.format(mirid))
                for accepted in core_dict:
                    outfile.write(
                        '>{}\n{}\n'
                            .format(accepted, core_dict[accepted].split('\t')[3])
                            .replace('-', '')
                    )
            alignment = '{}/{}.aln'.format(out_folder, mirid)
            stockholm = '{}/{}.sto'.format(out_folder, mirid)
            outtree = '{}/{}.dnd'.format(out_folder, mirid)
            t_coffee = 't_coffee'

            print('Building T-Coffee alignment.')
            tc_cmd_1 = (
                '{} -quiet -multi_core={} -special_mode=rcoffee -in {} '
                '-output=clustalw_aln -outfile={} -newtree={} -remove_template_file=1'
                .format(t_coffee, c, corefile, alignment, outtree)
            )
            sp.call(tc_cmd_1, shell=True)

            print('Adding secondary structure to Stockholm format.')
            tc_cmd_2 = (
                '{} -other_pg seq_reformat -in {} -action +add_alifold -output '
                'stockholm_aln -out {}'
                .format(t_coffee, alignment, stockholm)
            )
            sp.call(tc_cmd_2, shell=True)
        else:
            print('\nNo reverse hit found for any of the core species.\n'
                  'Skipping alignment..'
                  )
            return(mirid)
