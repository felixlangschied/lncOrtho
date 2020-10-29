# cmsearch_parser: Parse the output of cmsearch while eliminating
#                  duplicates and filtering entries according to the
#                  defined cutoff.
# Arguments:
# cms: path to cmsearch output
# cmc: cutoff to decide which candidate hits should be included for the
#      reverse BLAST search
# lc: length cutoff
# mirid: name/id of the microRNA
def CmsearchParser(cms, cmc, lc, mirid):
    # Output
    hits_dict = {}
    # Required for finding duplicates, stores hits per chromosome
    chromo_dict = {}
    cut_off = cmc

    with open(cms) as cmsfile:
        # Collect only the hits which satisfy the bit score cutoff.
        hits = [
            line.strip().split() for line in cmsfile
            if not line.startswith('#')
            and float(line.strip().split()[14]) >= cut_off
            and abs(int(line.split()[7])-int(line.strip().split()[8])) >= lc
        ]

        # Add the hits to the return dictionary.
        if hits:
            for candidate_nr, hit in enumerate(hits, 1):
                data = (
                    '{0}_c{1}'.format(mirid, candidate_nr),
                    hit[0], hit[7], hit[8], hit[9], hit[14]
                )

                hits_dict[data[0]] = data
                # Store the hits that satisfy the bit score cutoff to filter
                # duplicates.
                try:
                    chromo_dict[data[1]].append(data)
                except:
                    chromo_dict[data[1]] = [data]

    # Loop over the candidate hits to eliminate duplicates.
    for chromo in chromo_dict:
                nrhits = len(chromo_dict[chromo])
                if nrhits > 1:
                    for hitnr in range(nrhits):
                        start = int(chromo_dict[chromo][hitnr][2])
                        stop = int(chromo_dict[chromo][hitnr][3])
                        strand = chromo_dict[chromo][hitnr][4]
                        score = float(chromo_dict[chromo][hitnr][5])
                        for chitnr in range(hitnr+1, nrhits):
                            if strand != chromo_dict[chromo][chitnr][4]:
                                cstart = int(chromo_dict[chromo][chitnr][2])
                                cstop = int(chromo_dict[chromo][chitnr][3])
                                cscore = float(chromo_dict[chromo][chitnr][5])
                                # Test if the two hits from opposite strands
                                # overlap, which means one of them is
                                # (probably) a false positive.
                                # Out of two conflicting hits, the one with the
                                # highest cmsearch bit score is retained.
                                if (
                                    start in range(cstart, cstop+1)
                                    or stop in range(cstart, cstop+1)
                                    or cstart in range(start, stop+1)
                                    or cstop in range(start, stop+1)
                                ):
                                    if score > cscore:
                                        try:
                                            del hits_dict[chromo_dict[chromo]
                                                [chitnr][0]]
                                        except:
                                            pass
                                    else:
                                        try:
                                            del hits_dict[chromo_dict[chromo]
                                                [hitnr][0]]
                                        except:
                                            pass
    return hits_dict
