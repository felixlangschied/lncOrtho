import os
import subprocess as sp
import sys

from lib.blastparser_felix import BlastParser
from lib.genparser_felix import GenomeParser
from lib.cmsearch_parser import cmsearch_parser
from lib.mirna_maker import mirna_maker

# write_output: Write a FASTA file containing the accepted orthologs.
# Arguments:
# a: dictionary of accepted hits
# o: path for output
def write_output(a, o):
    with open(o, 'w') as outfile:
        for hit in a:
            outfile.write('>{0}\n{1}\n'.format(hit, a[hit]))

# blast_search: Perform a reverse BLAST search in the reference genome for a
# candidate.
# s: cmsearch result
# o: output for the blastsearch
# c: number of threads
# ref_blast_db: blastDB of the reference genome (tested and created in ncortho_main.py)
def blast_search(s, ref_blast_db, o, c):

    blast_command = (
        'blastn -task blastn -db {0} -query {1} '
        '-out {2} -num_threads {3} -outfmt 6'.format(ref_blast_db, s, o, c)
    )
    sp.call(blast_command, shell=True)


def ncortho(mirnas, models, output, msl, cpu, query, cm_cutoff, ref_blast_db):
    # Create miRNA objects from the list of input miRNAs.
    try:
        mirna_dict = mirna_maker(mirnas, models, output, msl)
    except TypeError:
        print('Could not parse the provided miRNA file.'
              'Please check the format of your file at:\n'
              '{}\n'
              'Exiting..'.format(mirnas))
        sys.exit()


    # Create Outdict:
    mirout_dict = {}

    # Identify ortholog candidates.
    for mir_data in mirna_dict:
        #print(mir_data)
        mirna = mirna_dict[mir_data]
        mirna_id = mirna.name
        outdir = '{}/{}'.format(output, mirna_id)
        # Create output folder, if not existent.
        if not os.path.isdir(outdir):
            try:
                sp.run('mkdir {}'.format(outdir), shell=True, check=True, stderr=sp.PIPE)
            except sp.CalledProcessError:
                print('Could not create output folder at\n'
                      '{}'.format(outdir))
        print('\n# Running covariance model search for {}.'.format(mirna_id))
        cms_output = '{0}/cmsearch_{1}.out'.format(outdir, mirna_id)
        # Calculate the bit score cutoff.
        cut_off = mirna.bit * cm_cutoff
        # Calculate the length cutoff.
        len_cut = len(mirna.pre) * msl
        # Perform covariance model search.
        # Report and inclusion thresholds set according to cutoff.
        cms_command = (
            'cmsearch -T {5} --incT {5} --cpu {0} --noali '
            '--tblout {1} {2}/{3}.cm {4}'
                .format(cpu, cms_output, models, mirna_id, query, cut_off)
        )
        sp.call(cms_command, shell=True)
        cm_results = cmsearch_parser(cms_output, cut_off, len_cut, mirna_id)

        print('These are the cm_results')
        print(cm_results)

        # Extract sequences for candidate hits (if any were found).
        if not cm_results:
            print('# No hits found for {}.\n'.format(mirna_id))
            continue
        else:
            gp = GenomeParser(query, cm_results.values())
            candidates = gp.extract_sequences()
            #print(candidates)
            nr_candidates = len(candidates)
            if nr_candidates == 1:
                print(
                    '\n# Covariance model search successful, found 1 '
                    'ortholog candidate above the threshold.\n'
                )
            else:
                print(
                    '\n# Covariance model search successful, found {} '
                    'ortholog candidates above the threshold.\n'
                        .format(nr_candidates)
                )
            print('# Evaluating candidates.\n')

        # Perform reverse BLAST test to verify candidates, stored in
        # a list (accepted_hits).
        accepted_hits = {}
        for candidate in candidates:
            # print('This is a candidate:')
            # print(candidate)
            sequence = candidates[candidate]
            temp_fasta = '{0}/{1}.fa'.format(outdir, candidate)
            # TODO: change BlastParser to take query directly from command-line
            #       to avoid creating temporary files
            with open(temp_fasta, 'w') as tempfile:
                tempfile.write('>{0}\n{1}'.format(candidate, sequence))
            blast_output = '{0}/blast_{1}.out'.format(outdir, candidate)

            blast_search(temp_fasta, ref_blast_db, blast_output, cpu)
            bp = BlastParser(mirna, blast_output, msl)
            if bp.parse_blast_output():
                accepted_hits[candidate] = sequence

        # Write output file if at least one candidate got accepted.
        if accepted_hits:
            nr_orthologs = len(accepted_hits)
            if nr_orthologs == 1:
                print('# ncOrtho found 1 verified ortholog.\n')
            else:
                print(
                    '# ncOrtho found {} verified orthologs.\n'
                        .format(nr_orthologs)
                )
            print('# Writing output of accepted candidates.\n')
            outpath = '{0}/{1}_orthologs.fa'.format(outdir, mirna_id)
            write_output(accepted_hits, outpath)
            print('# Finished writing output.\n')

            print('These are the accepted_hits')
            print(accepted_hits)

            out_dict = {}
            for key in accepted_hits:
                out_dict[key] = (accepted_hits[key], cm_results[key])
            print('this is the out_dict')
            print(out_dict)
            print('Trying to update')
            mirout_dict.update(out_dict)
            print(mirout_dict)
            return(mirout_dict)
        else:
            print(
                '# None of the candidates for {} could be verified.\n'
                    .format(mirna_id)
            )
            print('# No hits found for {}.\n'.format(mirna_id))
        print('# Finished ortholog search for {}.'.format(mirna_id))