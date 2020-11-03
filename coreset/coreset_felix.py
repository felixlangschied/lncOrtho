# Create a core set of orthologs
# Find the corresponding syntenic regions in reference and core species
# Search for core orthologs by reciprocal BLAST search
# Create Stockholm structural alignment

""" This is a reprogramming of the script written by Andreas"""

import argparse
import multiprocessing as mp
import sys
import glob
import os
import subprocess as sp
import pyfaidx
from core_blastsearch import blast_search
from createcm_felix import create_cm


# Parse a GTF file to store the coordinates for each protein-coding gene in a
# dictionary
def gtf_parser(gtf):
    chr_dict = {}
    chromo = ''

    with open(gtf, 'r') as gtf_file:
        for line in gtf_file:
            if (
                    not line.startswith('#')
                    and line.split()[2] == 'gene'
                    and line.split('gene_biotype ')[1].split('\"')[1] == 'protein_coding'
            ):
                linedata = line.strip().split('\t')

                gene_id = linedata[-1].split('\"')[1]
                contig = linedata[0]
                start = int(linedata[3])
                end = int(linedata[4])
                strand = linedata[6]

                if contig != chromo:
                    i = 1
                    chromo = contig
                chr_dict[gene_id] = (contig, i)

                try:
                    chr_dict[contig][i] = (gene_id, start, end, strand)
                except:
                    chr_dict[contig] = {i: (gene_id, start, end, strand)}
                i += 1
    return chr_dict


###############################################################################
# Try to find the ortholog for a given reference gene in a core set species
def ortho_search(r_gene, ortho_dict):
    orthologs = {}
    for core_taxon in ortho_dict.keys():
        try:
            ortholog = ortho_dict[core_taxon][r_gene]
            orthologs[core_taxon] = ortholog
            print(
                '{0} is the ortholog for {1} in {2}.'
                    .format(ortholog, r_gene, core_taxon)
            )
        except:
            print(
                'No ortholog found for {0} in {1}.'
                    .format(r_gene, core_taxon)
            )
    return orthologs


def main():
    ortho_dict = {}
    mirna_dict = {}
    neighbor_dict = {}

    # Print header
    print('\n' + '#' * 43)
    print('###' + ' ' * 37 + '###')
    print('###   ncOrtho - core set construction   ###')
    print('###' + ' ' * 37 + '###')
    print('#' * 43 + '\n')

    # required arguments
    # input miRNAs, reference gtf, reference genome, core gtf, core genome, OMA orthologs, output, (cpu), (mip)
    # python {} -n {} -r {} -g {} -c {} -q {} -p {} -o {} (-t {}) (-m {})
    # Parse command-line arguments
    # Define global variables
    parser = argparse.ArgumentParser(
        prog='python coreset_felix.py', description='core set construction'
    )
    # mirna data
    parser.add_argument(
        '-n', '--ncrna', metavar='<path>', type=str,
        help='path to your reference micrornas'
    )
    # reference gtf
    parser.add_argument(
        '-r', '--reference', metavar='<.gtf>', type=str,
        help='Path to reference GTF file'
    )
    # reference genome
    parser.add_argument(
        '-g', '--genome', metavar='<.fa>', type=str,
        help='Path to the reference genome in FASTA format'
    )
    # core gtf
    parser.add_argument(
        '-c', '--core', metavar='<path>', type=str,
        help='Path to the GTF files of the core genomes'
    )
    # core genomes
    parser.add_argument(
        '-q', '--query', metavar='<path>', type=str,
        help='Path to the core genomes in FASTA format'
    )
    # pairwise orthologs folder
    parser.add_argument(
        '-p', '--pairwise', metavar='<path>', type=str,
        help='path to pairwise orthologs'
    )
    # output folder
    parser.add_argument(
        '-o', '--output', metavar='<path>', type=str,
        help='path for the output folder'
    )
    # cpu, use maximum number of available cpus if not specified otherwise
    parser.add_argument(
        '-t', '--threads', metavar='int', type=int,
        help='number of CPU cores to use', nargs='?',
        const=mp.cpu_count(), default=mp.cpu_count()
    )
    # Maximum gene insertions
    parser.add_argument(
        '-m', '--mgi', metavar='int', type=int,
        help='maximum number of gene insertions', nargs='?',
        const=3, default=3
    )

    ###############################################################################

    # Show help when no arguments are added.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    # Check if computer provides the desired number of cores.
    available_cpu = mp.cpu_count()
    if args.threads > available_cpu:
        print(
            '# Error: The provided number of CPU cores is higher than the '
            'number available on this system. Exiting...'
        )
        sys.exit(1)
    else:
        cpu = args.threads

    # TODO: include checks for validity of arguments
    # os.getcwd()
    # os.chdir(path)
    # os.path.exists(path)
    # os.path.isfile(path)
    # os.path.isdir(path)

    mirna_path = args.ncrna
    output = args.output
    query = args.query
    ref_gtf_path = args.reference
    ref_genome = args.genome
    oma_paths = glob.glob('{}/*'.format(args.pairwise))
    core_gtf_paths = args.core
    core_fa_paths = args.query
    mgi = args.mgi

    # Create output folder if it not already exists
    if not os.path.isdir(output):
        mkdir_cmd = 'mkdir {}'.format(output)
        sp.call(mkdir_cmd, shell=True)

    ###############################################################################

    # Parse the pairwise orthologs
    for oma_path in oma_paths:
        taxon = oma_path.split('/')[-1]
        with open(oma_path, 'r') as oma_file:
            oma_lines = oma_file.readlines()
            orthologs = {
                ref: core for (ref, core) in [
                    (line.split()[0], line.split()[1])
                    for line in oma_lines
                ]
            }
            ortho_dict[taxon] = orthologs
            # {Taxon: {orth_ref1: orth_core1, orth_ref2: orth_core2...}}

    # Read in the miRNA data
    with open(mirna_path) as mirfile:
        mirnas = [
            line.split() for line in mirfile.readlines()
            if not line.startswith('#')
        ]
    ref_dict = gtf_parser(ref_gtf_path)

    # Determine the position of each miRNA and its neighboring gene(s)
    for mirna in mirnas:
        mirid = mirna[0]
        print('### {0} ###'.format(mirid))
        # Check if output folder exists or create it otherwise
        if not os.path.isdir('{}/{}'.format(output, mirid)):
            sp.call('mkdir {}/{}'.format(output, mirid), shell=True)
        # Workaround for differing naming conventions in miRBase and Ensembl
        if 'chr' in mirna[1]:
            chromo = mirna[1].split('chr')[1]
        else:
            chromo = mirna[1]
        ###
        start = int(mirna[2])
        end = int(mirna[3])
        strand = mirna[4]

        # find left neighbor or check if located inside gene
        # chr_dict[contig][i] = (geneid, start, end, strand)
        ###############################################################################
        # case 1): there is no protein-coding gene on the same contig as the miRNA,
        # so there can be no neighbors (should only occur in highly fragmented
        # assemblies)
        if not chromo in ref_dict.keys():
            print(
                'There are no protein-coding genes on contig {0}. '
                'Synteny around {1} cannot be established.'
                    .format(chromo, mirid)
            )
            continue

        # case 2): miRNA is located left of the first gene and hence has no left
        # neighbor, the first gene is therefore by default the right neighbor
        if end < int(ref_dict[chromo][1][1]):  # ref_dict[chromo][coordinate(starts at 1)][index into values]
            print(
                'There is no left neighbor of {0}, since it is located at the '
                'start of contig {1}.'.format(mirid, chromo)
            )
            print(
                '{0} is the right neighbor of {1}.'
                    .format(ref_dict[chromo][1][0], mirid)
            )
            continue

        # case 3): miRNA is located right to the last gene, so the last gene is the
        # left neighbor and there cannot be a right neighbor
        elif start > int(ref_dict[chromo][len(ref_dict[chromo])][2]):
            print(
                '{0} is the left neighbor of {1}.'
                    .format(ref_dict[chromo][len(ref_dict[chromo])][0], mirid)
            )
            print(
                'There is no right neighbor of {0}, since it is located at the'
                ' end of contig {1}.'.format(mirid, chromo)
            )
            continue

        # case 4): miRNA is located either between two genes or overlapping with (an
        # intron of) a gene, either on the same or the opposite strand
        ###############################################################################
        else:
            solved = False
            for i, gene in enumerate(ref_dict[chromo]):
                gene_data = ref_dict[chromo][gene]  # gene_data = (geneID,start,end,strand)
                # case 4.1): miRNA inside gene
                if (
                        start >= gene_data[1]
                        and end <= gene_data[2]
                        and strand == gene_data[3]
                ):
                    solved = True
                    print(
                        '{0} is located inside the gene {1}.'
                            .format(mirid, gene_data[0])
                    )
                    ortho_hits = ortho_search(gene_data[0], ortho_dict)
                    for core_tax in ortho_hits:
                        try:
                            neighbor_dict[core_tax][mirid] = (
                                ('inside', ortho_hits[core_tax])
                            )
                        except:
                            neighbor_dict[core_tax] = (
                                {mirid: ('inside', ortho_hits[core_tax])}
                            )
                    break
                # case 4.2): miRNA opposite of gene
                elif (
                        start >= gene_data[1]
                        and end <= gene_data[2]
                        and strand != gene_data[3]
                ):
                    solved = True
                    print(
                        '{0} is located opposite of the gene {1}.'
                            .format(mirid, gene_data[0])
                    )
                    ortho_hits = ortho_search(gene_data[0], ortho_dict)
                    for core_tax in ortho_hits:
                        try:
                            neighbor_dict[core_tax][mirid] = (
                                ('opposite', ortho_hits[core_tax])
                            )
                        except:
                            neighbor_dict[core_tax] = (
                                {mirid: ('opposite', ortho_hits[core_tax])}
                            )
                    break
                # case 4.3): miRNA between genes
                elif (
                        int(ref_dict[chromo][gene][2]) < start
                        and ref_dict[chromo][gene + 1][1] > end
                ):
                    solved = True
                    ###############################################################################
                    print(
                        '{1} is the left neighbor of {2}.'
                        .format(gene, ref_dict[chromo][gene][0], mirid)
                    )
                    print(
                        '{1} is the right neighbor of {2}.'
                        .format(gene, ref_dict[chromo][gene + 1][0], mirid)
                    )
                    left_hits = ortho_search(gene_data[0], ortho_dict)
                    right_hits = (
                        ortho_search(ref_dict[chromo][gene + 1][0], ortho_dict)
                    )
                    # save only the hits where both genes have orthologs in a species
                    if left_hits:
                        for taxon in left_hits:
                            if taxon in right_hits:
                                try:
                                    neighbor_dict[taxon][mirid] = (
                                        (
                                            'in-between',
                                            [left_hits[taxon],
                                             right_hits[taxon]]
                                        )
                                    )
                                except:
                                    neighbor_dict[taxon] = (
                                        {mirid: (
                                            'in-between',
                                            [left_hits[taxon],
                                             right_hits[taxon]]
                                        )}
                                    )
                    break
        if not solved:
            print('Unable to resolve synteny for {}.'.format(mirid))

    print('\n#########################')
    print('### Synteny Analysis ###')
    print('#########################')

    if not neighbor_dict:
        print(
            "Not enough orthologs found in the core species for synteny analysis.\n"
            "Exiting..."
        )
        sys.exit()

    # Search for the coordinates of the orthologs and extract the sequences
    for taxon in neighbor_dict:  # neighbor_dict = {'Core_taxon': {mirnaID: ('category', [geneIDleft, geneIDright])}
        print('\nStarting synteny analysis for {}'.format(taxon))
        gtf_path = '{0}/{1}.gtf'.format(core_gtf_paths, taxon)
        fasta_path = glob.glob('{0}/{1}*.fa'.format(core_fa_paths, taxon))
        if len(fasta_path) != 1:
            print('Unable to identify genome file for {}'.format(taxon))
            continue
        genome = pyfaidx.Fasta(fasta_path[0])  # genome = {fasta_geneId: sequence}

        # parse annotation
        # print('Trying to parse GTF file for {}.'.format(taxon))
        try:
            core_gtf_dict = gtf_parser(gtf_path)

            for mirna in neighbor_dict[taxon]:
                # print(mirna)

                style = neighbor_dict[taxon][mirna][0]
                if style == 'inside' or style == 'opposite':
                    ###############################################################################
                    try:
                        ortho_data = (
                            core_gtf_dict[neighbor_dict[taxon][mirna][1]]
                            # ortho_data = (contig_id,coordinate on contig)
                        )
                        positions = list(
                            core_gtf_dict[ortho_data[0]][ortho_data[1]][1:4]  # positions = [start, end, strand]
                        )
                        coordinates = [ortho_data[0]] + positions  # coordinates = [contig_id, start, end, strand]
                        seq = (
                            genome[coordinates[0]]
                            [coordinates[1] - 1:coordinates[2]].seq
                        )
                        try:
                            mirna_dict[mirna][taxon] = seq
                        except:
                            mirna_dict[mirna] = {taxon: seq}
                    except:
                        print('{} not found in GTF file.'.format(mirna[1]))

                elif style == 'in-between':
                    left_data = (
                        core_gtf_dict[neighbor_dict[taxon][mirna][1][0]]
                        # left_data = (contig-id, coordinate on contig)
                    )
                    right_data = (
                        core_gtf_dict[neighbor_dict[taxon][mirna][1][1]]
                        # right_data = (contig-id, coordinate on contig)
                    )

                    # Test to see if the two orthologs are themselves neighbors where their
                    # distance cannot be larger than the selected mgi value. This accounts
                    # for insertions in the core species.
                    # TODO: Apply mgi also to the reference species to account for insertions
                    # in the reference.
                    if (
                            left_data[0] == right_data[0]
                            and abs(left_data[1] - right_data[1]) <= mgi
                    ):
                        # Determine which sequence to include for the synteny-based ortholog search
                        # depending on the order of orthologs. The order of the orthologs in the core
                        # species might be inverted compared to that in the reference species.
                        ###############################################################################
                        # gtf_parser = gene_id: (contig_id,coordinate on contig)
                        # gtf_parser = contig_id: (coordinate: (gene_id,start,end,strand))
                        # left_data = (contig-id, coordinate on contig)

                        if left_data[1] < right_data[1]:
                            contig = left_data[0]
                            seq_start = (
                                core_gtf_dict[contig][left_data[1]][2]
                            )
                            seq_end = (
                                core_gtf_dict[right_data[0]][right_data[1]][1]
                            )
                            seq = genome[contig][seq_start - 1:seq_end].seq
                            try:
                                mirna_dict[mirna][taxon] = seq
                            except:
                                mirna_dict[mirna] = {taxon: seq}

                        elif right_data[1] < left_data[1]:
                            contig = left_data[0]
                            seq_end = (
                                core_gtf_dict[left_data[0]][left_data[1]][2]
                            )
                            seq_start = (
                                core_gtf_dict[right_data[0]][right_data[1]][1]
                            )
                            seq = genome[contig][seq_start - 1:seq_end].seq
                            try:
                                mirna_dict[mirna][taxon] = seq
                            except:
                                mirna_dict[mirna] = {taxon: seq}
                        print('Synteny fulfilled for {}.'.format(mirna))
                    else:
                        print(
                            'No shared synteny for {} in {}.'
                            .format(mirna, taxon)
                        )
        except:
            print('No GTF file found for {}'.format(taxon))
            continue



    def write_fasta():
        for mirna in mirna_dict:
            with open('{0}/{1}/{1}.fa'.format(output, mirna), 'w') as outfile:
                for core_taxon in mirna_dict[mirna]:
                    outfile.write(
                        '>{0}\n{1}\n'
                        .format(core_taxon, mirna_dict[mirna][core_taxon])
                    )

    # write fasta format output of candidate regions of each miRNA in each core species
    write_fasta()

    print('\n#########################')
    print('### Reciprocal BLAST ###')
    print('#########################\n')

    noHits = []
    for mirna in mirnas:
        mirid = mirna[0]
        if os.path.isfile('{0}/{1}/{1}.fa'.format(output, mirid)):
            # calculate and write results of reciprocal BLAST search
            blast_search(mirna, ref_genome, output, cpu)
        else:
            print("No Hits found for {}, skipping...".format(mirid))
            noHits.append(mirid)
        if not os.path.isfile('{0}/{1}/{1}_core.fa'.format(output, mirid)):
            noHits.append(mirid)

    if noHits:
        print('Saving names of  ncRNAs with no hits at {}/noHits.txt'.format(output))
        with open('{}/noHits.txt'.format(output), 'w') as outfile:
            outfile.write('\n'.join(noHits))
            outfile.write('\n')
        for fail in noHits:
            rm_cmd = 'rm -r {0}/{1}'.format(output,fail)
            sp.call(rm_cmd, shell=True)

    # Cleanup extra T-coffee files
    run_path = os.getcwd()
    trash_paths = glob.glob('{}/*.template_list'.format(run_path))
    for trash in trash_paths:
        sp.call('rm {}'.format(trash), shell=True)


    print('\n#########################')
    print('### Creating CMs ###')
    print('#########################\n')

    notreciproc = []
    cm_output = output + '/' + 'core_models'
    for mirna in mirna_dict:
        align_path = '{0}/{1}/{1}.sto'.format(output, mirna)
        if mirna in noHits:
            continue
        elif os.path.isfile('{}/{}.cm'.format(cm_output, mirna)):
            print(
                'CM of {} already exists in {}.\n'
                'Skipping..\n'.format(mirna, cm_output)
            )
            continue
        elif os.path.isfile(align_path):
            with open(align_path, 'r') as infile:
                create_cm(infile.name, cm_output, cpu)
        else:
            print('No alignment file found for {}.\n'
                  'Added to list of ncRNAs that were not the best reciprocal hit'.format(mirna))
            notreciproc.append(mirna)

    if notreciproc:
        print('Saving names of ncRNAs that were not the best reciprocal hit at {}/notReciproc.txt'.format(output))
        with open('{}/notReciproc.txt'.format(output),'w') as outfile:
            outfile.write('\n'.join(notreciproc))
            outfile.write('\n')
        for fail in notreciproc:
            rm_cmd = 'rm -r {0}/{1}'.format(output, fail)
            sp.call(rm_cmd, shell=True)

if __name__ == '__main__':
    main()
