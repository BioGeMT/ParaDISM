def map_alncoords_to_seqcoords(aln):
    """
    Returns a numpy array in which entry [i,j] is the 
    coordinate of character aln[i,j] in sequence i. 
    Essentially removes gaps from the alignment.  
    Gaps before the start of the sequence get assigned a coordinate -1.  
    """
    curr_coord = [-1]*len(aln)
    sequence_coords = np.zeros((len(aln), aln.get_alignment_length()), dtype='int')
    for ci in range(aln.get_alignment_length()):
        for ri in range(len(aln)):
            if aln[ri, ci] != '-': curr_coord[ri] += 1
        sequence_coords[:,ci] = curr_coord
    return sequence_coords

def map_seqcoords_to_alncoords(aln):
    """
    Almost a reverse of map_alncoords_to_seqcoords. 
    Returns a list of lists LL, such that in L = LL[i], L[j] is the column number in the alignment aln
    corresponding to j-th residue in sequence i.
    Works in base zero.  
    """
    aln_coords = [[] for _ in range(len(aln))]
    for ci in range(aln.get_alignment_length()):
        for ri in range(len(aln)):
            if aln[ri, ci] != '-': 
                aln_coords[ri].append(ci)
    return aln_coords

if __name__ == '__main__':
    from Bio import AlignIO
    path_to_test_alignment = ''
    test_aln = AlignIO.read(path_to_test_alignment, 'fasta')
    map_alncoords_to_seqcoords(test_aln)
    map_seqcoords_to_alncoords(test_aln)
