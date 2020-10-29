'''
ncOrtho submodule
TODO: include license, author details etc

### TODO: calculate length of overlap to see if it is relevant/significant
### also the bit score should be compared
### add a check for the presence of the mature miRNA
'''

# BlastParser object performs reverse BLAST search and reports if a candidate
# fulfills the reverse best hit criterion
class BlastParser(object):
    # init parameters:
    # mirna: Mirna object that holds the location for the reference miRNA
    # blastpath: path to the output file of the reverse BLAST search
    # msl: ncOrtho minimum sequence length threshold
    def __init__(self, mirna, blastpath, msl):
        self.start = mirna.start
        self.end = mirna.end
        self.chromosome = mirna.chromosome
        self.strand = mirna.strand
        del mirna
        self.blastpath = blastpath
        self.msl = msl
        #self.blasthits = []
        #self.top_score = 100
        #self.top_score = blasthits[0][6]

    def parse_blast_output(self,):
        with open(self.blastpath) as blastfile:
            blasthits = [line.strip().split() for line in blastfile]
        # If no BLAST hit was found, the search failed by default.
        if not blasthits:
            return False
        # Otherwise, check if the best hit and the reference miRNA overlap.
        else:
            # Gather coordinates of best BLAST hit.
            tophit = blasthits[0]
            sseqid = tophit[1]
            sstart = int(tophit[8])
            send = int(tophit[9])
            del blasthits
            # Sequences must be on the same contig, otherwise overlap can be
            # ruled out instantaneously
            if not sseqid == self.chromosome:
                return False
            # Contigs match, so overlap is possible.
            else:
                # first within second
                if (
                    (sstart <= self.start and self.start <= send)
                    or (sstart <= self.end and self.end <= send)
                ):
                    return True
                # second within first
                elif (
                    (self.start <= sstart and sstart <= self.end)
                    or (self.start <= send and send <= self.end)
                ):
                    return True
                # No overlap
                else:
                    return False
