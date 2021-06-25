#python3

import logging
import os
import math

def LocationToGene(scaffold, pos, sortedGenes):
    """
    Description:
        Given scaffold, pos, and dict of scaffold to list of sorted genes,
        returns the locusId and the fraction through the gene it is 
        in (as a list of 2 elements) [locusId, fraction]
        
        If the location overlaps multiple genes or no genes, 
        returns locusId = "", f = "".
        
        Each gene should be a hash that contains begin, end, strand, and locusId
        
        This code does not support orfs that wrap around the origin, and it 
        may not give correct results if there are complicated overlaps between ORFs. 
        In particular, it only checks the adjacent ORF on either side to see if there 
        is an overlap.

        If the strand is the "-" strand, we return 1 - f

    """
    if scaffold == "pastEnd":
        return ["",""]
    if scaffold in sortedGenes:
        genelist = sortedGenes[scaffold]
        if genelist == None:
            return ["",""]
    else:
        return ["",""]

    # binary search
    # at all times, either the true index is between lo and hi,
    # or there is no hit
    # We search the middle of current search loc
    nGenes = len(genelist)
    lo = 0
    hi = nGenes - 1
    for nRound in range(100000):
        mid = int((lo+hi)/2)
        iBegin = int(genelist[mid]['begin'])
        iEnd = int(genelist[mid]['end'])
        if pos < iBegin:
            if mid == lo:
                return ["",""]
            hi = mid - 1
        elif pos > iEnd:
            if mid == hi:
                return ["",""]
            lo = mid + 1
        else:
            # Does the previous or next gene also overlap this position?
            if (mid > 0) and (int(genelist[mid-1]['begin']) <= pos) and \
                    pos <= int(genelist[mid-1]['end']):
                        return ["",""]
            if (mid < nGenes - 1) and (int(genelist[mid+1]['begin']) <= pos) \
                and (pos <= int(genelist[mid+1]['end'])):
                        return ["",""]
            if iBegin == iEnd:
                f = 0
            else:
                f = (pos - iBegin)/(iEnd - iBegin)
            strand = genelist[mid]['strand']
            # insertions near N terminus of gene should have f near 0
            # regardless of strand
            if strand == "-":
                f = 1.0 - f
            return [genelist[mid]['locusId'],f]
    raise Exception("Unreachable gene in scf: {}, pos: {}".format(
                       scaffold, pos))


def CheckGeneLocations(sortedGenes):
    """
    CheckGeneLocations Func Def

    Given a dict of scaffold to a list of sorted genes, check to see that
    for almost all genes, the center of the gene does not overlap other genes
    and is identified as the hit for that gene
    Each gene should be a hash that contains begin, end, strand, and locusId
    Writes a short report to STDERR and reports a reference to the list of 
    locusIds that fail (with overlap-adjacent cases excluded)

    Args:
        sortedGenes: (d)
            scaffoldId -> gene_list (list<gene>)
            gene: (d)
                begin:
                end:
                locusId:

    Returns:
        fail: (list<locus_id_str>)
            locus_id_str: (str) LocusId
    """
    ok = []
    wrap = []
    overlap = []
    fail = []
    for k in sortedGenes.keys():
        scaffoldId, geneList = k, sortedGenes[k]
        for i in range(len(geneList)):
            gene = geneList[i]
            if gene['end'] < gene['begin']:
                wrap.append(gene['locusId'])
            else:
                #int in perl returns 'integer portion' of number
                pos = math.floor((int(gene['begin']) + int(gene['end']) + 1 )/2)
                locusHit, f = LocationToGene(scaffoldId, pos, sortedGenes)
                if locusHit == "":
                    # does this position overlap adjacent genes?
                    if (i>0 and (int(geneList[i-1]['end']) >= pos)):
                        overlap.append(gene['locusId'])
                    elif (i< len(geneList) - 1 and 
                          int(geneList[i+1]['begin']) <= pos):
                        overlap.append(gene['locusId'])
                    else:
                        fail.append(gene['locusId'])
                elif locusHit == gene['locusId']:
                    ok.append(gene['locusId'])
                else:
                    fail.append(gene['locusId'])
    logging.info("Check gene locations: success "
            "{} wrap {} overlap-adjacent {} failure {}\n".format(
                len(ok), len(wrap), len(overlap), len(fail)))
    if len(fail) > 0:
        logging.info("Failures: " + " ".join(fail) + '\n')
    return fail
                        




                
