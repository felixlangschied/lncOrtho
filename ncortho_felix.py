'''
# Re-implementation of ncOrtho in Python
# TODO: include license, author information etc
'''

# Modules import

# Python
import argparse
import multiprocessing as mp
import os
import subprocess as sp
import sys

# Internal ncOrtho modules
from andreas.blastparser import BlastParser
from andreas.genparser import GenomeParser
from cmsearch_parser import CmsearchParser

###############################################################################


# Central class of microRNA objects
class Mirna(object):
    def __init__(self, name, chromosome, start, end, strand, pre, bit):
        # miRNA identifier
        self.name = name
        # chromosome that the miRNA is located on
        self.chromosome = chromosome
        # start position of the pre-miRNA
        self.start = int(start)
        # end position of the pre-miRNA
        self.end = int(end)
        # sense (+) or anti-sense (-) strand
        self.strand = strand
        # nucleotide sequence of the pre-miRNA
        self.pre = pre
        # nucleotide sequence of the mature miRNA
        #self.mature = mature
        # reference bit score that miRNA receives by its own
        # covariance model
        self.bit = bit
# TODO: include both mature strands, 5p and 3p, aka mature and star


# mirna_maker: Parses the miRNA data input file and returns a
#              dictionary of Mirna objects.
# Arguments:
# mirpath: path to file with microRNA data
# cmpath: path to covariance models
# output: path for writing temporary files
def mirna_maker(mirpath, cmpath, output, msl):
    
    mmdict = {} # will be the return object
    
    with open(mirpath) as mirna_file:
        mirna_data = [
            line.strip().split() for line in mirna_file
            if not line.startswith('#')
        ]

    for mirna in mirna_data:
        mirid = mirna[0]
        # Check if the output folder exists, otherwise create it.
        if not os.path.isdir('{}/{}'.format(output, mirid)):
            try:
                mkdir = 'mkdir {}/{}'.format(output, mirid)
                sp.call(mkdir, shell=True)
            except:
                print(
                    '# Cannot create output folder for {}.'
                    'Skipping to next miRNA.'
                )
                continue

        # Obtain the reference bit score for each miRNA by applying it
        # to its own covariance model.
        print('# Calculating reference bit score for {}.'.format(mirid))
        seq = mirna[5]
        query = '{0}/{1}/{1}.fa'.format(output, mirid)
        model = '{0}/{1}.cm'.format(cmpath, mirid)

        # Check if the covariance model even exists, otherwise skip to
        # the next miRNA.
        if not os.path.isfile(model):
            print('# No covariance model found for {}.'.format(mirid))
            continue
        
        # Create a temporary FASTA file with the miRNA sequence as
        # query for external search tool cmsearch to calculate
        # reference bit score.
        with open(query, 'w') as tmpfile:
            tmpfile.write('>{0}\n{1}'.format(mirid, seq))
        cms_output = '{0}/{1}/cmsearch_{1}_tmp.out'.format(output, mirid)
        cms_log = '{0}/{1}/cmsearch_{1}.log'.format(output, mirid)
        cms_command = (
            'cmsearch -E 0.01 --noali -o {3} --tblout {0} {1} {2}'
            .format(cms_output, model, query, cms_log)
        )
        sp.call(cms_command, shell=True)
        with open(cms_output) as cmsfile:
            hits = [
                line.strip().split() for line in cmsfile
                if not line.startswith('#')
            ]
            if hits:
                top_score = float(hits[0][14])
            # In case of any issues occuring in the calculation of the bit
            # score, no specific threshold can be determined. The value will
            # set to zero, which technically turns off the filter.
            else:
                print(
                    '# Warning: Self bit score not applicable, '
                    'setting threshold to 0.'
                )
                top_score = 0.0

        mirna.append(top_score)

        # Remove temporary files.
        for rmv_file in [cms_output, cms_log, query]:
            sp.call('rm {}'.format(rmv_file), shell=True)

        # Create output.
        mmdict[mirna[0]] = Mirna(*mirna)

        print('Reference Bit-Score determined as: {}'.format(top_score))

    return mmdict


# blast_search: Perform a reverse BLAST search in the reference genome for a
# candidate.
# s: cmsearch result
# r: reference genome
# o: output name
# c: number of threads
def blast_search(s, r, o, c):
    # Check if BLAST database already exists, otherwise create it.
    # Database files are ".nhr", ".nin", ".nsq".
    file_extensions = ['.nhr', '.nin', '.nsq']
    for fe in file_extensions:
        checkpath = '{}{}'.format(r, fe)
        if not os.path.isfile(checkpath):
        # At least one of the BLAST db files is not existent and has to be
        # created.
            db_command = 'makeblastdb -in {} -dbtype nucl'.format(r)
            sp.call(db_command, shell=True)
            break
    blast_command = (
       'blastn -task blastn -db {0} -query {1} '
       '-out {2} -num_threads {3} -outfmt 6'.format(r, s, o, c)
    )
    sp.call(blast_command, shell=True)


# write_output: Write a FASTA file containing the accepted orthologs.
# Arguments:
# a: dictionary of accepted hits
# o: path for output
def write_output(a, o):
    with open(o, 'w') as outfile:
        for hit in a:
            outfile.write('>{0}\n{1}\n'.format(hit, a[hit]))


# Main function
def main():

    # Print header
    print('\n'+'#'*57)
    print('###'+' '*51+'###')
    print('###   ncOrtho - ortholog search for non-coding RNAs   ###')
    print('###'+' '*51+'###')
    print('#'*57+'\n')
    
    # Parse command-line arguments
    # Define global variables
    parser = argparse.ArgumentParser(
        prog='python ncortho.py', description='ncRNA orthology prediction tool'
    )
    # cpu, use maximum number of available cpus unless specified otherwise
    parser.add_argument(
        '-c', '--cpu', metavar='int', type=int,
        help='number of cpu cores ncOrtho should use', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
    )
    # covariance models folder
    parser.add_argument(
        '-m', '--models', metavar='<path>', type=str,
        help='path to your covariance models'
    )
    # mirna data
    parser.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str,
        help='path to your reference micrornas'
    )
    # output folder
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='path for the output folder'
    )
    # query genome
    parser.add_argument(
        '-q', '--query', metavar='<.fa>', type=str,
        help='path to query genome'
    )
    # reference genome
    parser.add_argument(
        '-r', '--reference', metavar='<.fa>', type=str,
        help='path to reference genome'
    )
    # bit score cutoff for cmsearch hits
    parser.add_argument(
        '-t', '--cutoff', metavar='float', type=float,
        help='cmsearch bit score cutoff', nargs='?', const=0.6, default=0.6
    )
    # length filter to prevent short hits
    parser.add_argument(
        '-l', '--msl', metavar='float', type=float,
        help='hit length filter', nargs='?', const=0.9, default=0.9
    )
    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
    
    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if args.cpu > available_cpu:
        print(
            '# Error: The provided number of CPU cores is higher than the '
            'number available on this system. Exiting...'
        )
        sys.exit(1)
    else:
        cpu = args.cpu

    # TODO: include checks for validity of arguments
    # os.getcwd()
    # os.chdir(path)
    # os.path.exists(path)
    # os.path.isfile(path)
    # os.path.isdir(path)

    mirnas = args.ncrna
    models = args.models
    output = args.output
    query = args.query
    reference = args.reference
    cm_cutoff = args.cutoff
    # Not in use yet
    msl = args.msl
    #blast_cutoff = args.blastc
       
    # Create miRNA objects from the list of input miRNAs.
    mirna_dict = mirna_maker(mirnas, models, output, msl)

    # Identify ortholog candidates.
    for mir_data in mirna_dict:
        print(mir_data)
        mirna = mirna_dict[mir_data]
        mirna_id = mirna.name
        outdir = '{}/{}'.format(output, mirna_id)
        # Create output folder, if not existent.
        if not os.path.isdir(outdir):
            sp.call('mkdir {}'.format(outdir), shell=True)
        print('\n# Running covariance model search for {}.'.format(mirna_id))
        cms_output = '{0}/cmsearch_{1}.out'.format(outdir, mirna_id)
        # Calculate the bit score cutoff.
        cut_off = mirna.bit*cm_cutoff
        # Calculate the length cutoff.
        len_cut = len(mirna.pre)*msl
        # Perform covariance model search.
        # Report and inclusion thresholds set according to cutoff.
        cms_command = (
            'cmsearch -T {5} --incT {5} --cpu {0} --noali '
            '--tblout {1} {2}/{3}.cm {4}'
            .format(cpu, cms_output, models, mirna_id, query, cut_off)
        )
        sp.call(cms_command, shell=True)
        cm_results = CmsearchParser(cms_output, cut_off, len_cut, mirna_id)

        #test
        #print(cm_results)
        #print(cm_results)
        
        # Extract sequences for candidate hits (if any were found).
        if not cm_results:
            print('# No hits found for {}.\n'.format(mirna_id))
            continue
        else:
            gp = GenomeParser(query, cm_results.values())
            candidates = gp.extract_sequences()
            print(candidates)
            nr_candidates = len(candidates)
            if nr_candidates == 1:
                print(
                    '\n# Covariance model search successful, found 1 '
                    'ortholog candidate.\n'
                )
            else:
                print(
                    '\n# Covariance model search successful, found {} '
                    'ortholog candidates.\n'
                    .format(nr_candidates)
                )
            print('# Evaluating candidates.\n')
        
        # Perform reverse BLAST test to verify candidates, stored in
        # a list (accepted_hits).
        accepted_hits = {}
    
        for candidate in candidates:
            print(candidate)
            sequence = candidates[candidate]
            temp_fasta = '{0}/{1}.fa'.format(outdir, candidate)
            # TODO: change BlastParser to take query directly from command-line
            #       to avoid creating temporary files
            with open(temp_fasta, 'w') as tempfile:
                tempfile.write('>{0}\n{1}'.format(candidate, sequence))
            blast_output = '{0}/blast_{1}.out'.format(outdir, candidate)
            blast_search(temp_fasta, reference, blast_output, cpu)
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
        else:
            print(
                '# None of the candidates for {} could be verified.\n'
                .format(mirna_id)
            )
            print('# No hits found for {}.\n'.format(mirna_id))
        print('# Finished ortholog search for {}.'.format(mirna_id))

if __name__ == "__main__":
    main()
