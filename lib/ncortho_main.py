import os
import subprocess as sp
import sys

from lib.blastparser_felix import BlastParser
from lib.genparser_felix import GenomeParser
from lib.cmsearch_parser import cmsearch_parser
from lib.mirna_maker import mirna_maker


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

        # check if cm-search was already run, otherwise read from file
        cmparse_results = '{0}/CMresults_{1}.txt'.format(outdir, mirna_id)
        if not os.path.isfile(cmparse_results):
            # Perform covariance model search.
            # Report and inclusion thresholds set according to cutoff.
            cms_command = (
                'cmsearch -T {5} --incT {5} --cpu {0} --noali '
                '--tblout {1} {2}/{3}.cm {4}'
                .format(cpu, cms_output, models, mirna_id, query, cut_off)
            )
            sp.call(cms_command, shell=True)
            cm_results = cmsearch_parser(cms_output, cut_off, len_cut, mirna_id)

            with open(cmparse_results, 'w') as file:
                file.write(str(cm_results))
        else:
            with open(cmparse_results, 'r') as file:
                cm_results = eval(file.read())

        # Extract sequences for candidate hits (if any were found).
        if not cm_results:
            print('# No hits found for {}.\n'.format(mirna_id))
            continue
        else:
            gp = GenomeParser(query, cm_results.values())
            candidates = gp.extract_sequences()
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
        # a dict (accepted_hits).
        accepted_hits = {}
        for candidate in candidates:
            sequence = candidates[candidate]
            blast_cmd = (
                'blastn -task blastn -db {0} -num_threads {1} '
                '-outfmt 6 -query <(echo -e \">{2}\\n{3}\")'.format(ref_blast_db, cpu, candidate, sequence)
            )

            p = sp.Popen(blast_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, executable='/bin/bash')
            blast_output, err = p.communicate()
            if not err:
                blast_output = blast_output.decode('utf-8')
            else:
                print(err)
                exit()

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

            out_dict = {}
            for key in accepted_hits:
                out_dict[key] = (accepted_hits[key], cm_results[key])
            mirout_dict.update(out_dict)

        else:
            print(
                '# None of the candidates for {} could be verified.\n'
                .format(mirna_id)
            )
            print('# No hits found for {}.\n'.format(mirna_id))
        print('# Finished ortholog search for {}.'.format(mirna_id))
    return mirout_dict
