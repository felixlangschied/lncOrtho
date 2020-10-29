import argparse
from Bio import SeqIO
"""
Created on Thu 8.10.20
based on ncortho-python/utils/mirna_input.py

@author: andreas, edited by felix

Create an miRNA input file for ncOrtho (from miRBase input data)
TODO: include both mature sequences (-5p and -3p) if available
(implementation of andreas' script doesn't work yet)

Line 31 is specific to miRBase's gff3 format of v22.1 in giving chromosome numbers as (chrX)
"""


def fasta_parser(fa_path):
    with open(fa_path,'r') as fa_file:
        seq = SeqIO.parse(fa_file,'fasta')
        seqdict = SeqIO.to_dict(seq)
        return(seqdict)


def gff_parser(gff_path,pre_dict):
    with open(gff_path, 'r') as gff_file:
        tmp_dict = dict()
        for line in gff_file:
            if not line.startswith('#'):
                linedata = line.strip().split('\t')
                if linedata[2] == 'miRNA_primary_transcript':
                    chromo = linedata[0].split('chr')[1]
                    start = linedata[3]
                    end = linedata[4]
                    strand = linedata[6]
                    mirna_id = linedata[-1].split('ID=')[1].split(';')[0]
                    name = linedata[-1].split('Name=')[1].split(';')[0]
                    tmp_dict[mirna_id] = [name, chromo, start, end, strand, str(pre_dict[name].seq)]
        return tmp_dict

def write_output(out_path,out_dict):
    with open(out_path,'w') as outfile:
        header = '\t'.join(['# miRBase_id','Chromosome','Start','End','Strand','Sequence']) + '\n'
        outfile.write(header)
        for mirna in out_dict:
            lineout = '\t'.join(out_dict[mirna]) + '\n'
            outfile.write(lineout)

def main():
    #define global variables
    parser = argparse.ArgumentParser(prog='parse_mirna', description='tool to produce ncRNA input data files for ncOrtho')
    #gff3 file from miRBase
    parser.add_argument('-g', '--gff', metavar='<gff3>', type=str, help='Path to gff file')
    #miRNA sequences
    parser.add_argument('-s', '--sequences', metavar='<fa>', type=str, help='Path to (pre-)miRNA sequences')
    #output
    parser.add_argument('-o', '--output', metavar='<path>', type=str, help='Path for the output')

    args = parser.parse_args()
    gff = args.gff
    seqs = args.sequences
    out = args.output

    seq_dict = fasta_parser(seqs)
    out_dict = gff_parser(gff,seq_dict)

    write_output(out,out_dict)

if __name__ == "__main__":
    main()


#python parse_mirna.py -s '/share/project/felixl/ncOrtho/data/mouse_ref_core/mmu_miRNAs/mirbase_v22.1_mmu.fa' -g '/share/project/felixl/ncOrtho/data/mouse_ref_core/test_mirnaSet/test_mmu.gff3' -o '/share/project/felixl/ncOrtho/data/mouse_ref_core/test_mirnaSet/test_mirnas.tsv'



