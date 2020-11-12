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
        #print('I want it to be on this contig:')
        #print(self.chromosome)

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
            #print('This is the tophit')
            #print(tophit)

            sseqid = tophit[1]
            #print('This is the contig')
            #print(sseqid)

            sstart = int(tophit[8])
            #print('This is the hit start')
            #print(sstart)

            send = int(tophit[9])
            #print('This is the hit end')
            #print(send)

            #print('This is the reference start')
            #print(self.start)

            #print('This is the reference end')
            #print(self.end)

            del blasthits
            # Sequences must be on the same contig, otherwise overlap can be
            # ruled out instantaneously
            if not sseqid == self.chromosome:
                #print('Not on the same contig')
                return False

            # Contigs match, so overlap is possible.
            else:
                # first within second
                if (
                    (sstart <= self.start and self.start <= send)
                    or (sstart <= self.end and self.end <= send)
                ):
                    #print('case 1')
                    return True

                # second within first
                elif (
                    (self.start <= sstart and sstart <= self.end)
                    or (self.start <= send and send <= self.end)
                ):
                    #print('case 2')
                    return True
                # No overlap
                else:
                    #print('case 3')
                    return False
