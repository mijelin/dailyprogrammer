import sys

verbose = False

def read_input_seqs(input_fname):
    """Returns a list of input DNA sequences read from filename INPUT_FNAME."""
    f = open(input_fname, 'r')
    input_seqs = [ line[:-1] for line in f ] # Strip \n

    if verbose:
        print "> Reading input sequences from file '%s':" % input_fname
        print ">", input_seqs

    # Sort by string length
    input_seqs.sort(key = lambda s : len(s))

    return input_seqs

def remove_substrings(input_seqs):
    """Remove subsequences that are wholly contained in another subsequence."""
    numSeqs = len(input_seqs)
    input_list = input_seqs[:]

    for i in range(numSeqs):
        for j in range(i + 1, numSeqs):
            if input_seqs[i] in input_seqs[j]:
                if verbose:
                    print "> Removing '%s' which is contained in '%s'." \
                          % (input_seqs[i], input_seqs[j])
                input_list.remove(input_seqs[i])

    return input_list

def find_scs(input_seqs):
    """Returns the shortest common supersequence composed of the subsequences
    in INPUT_SEQS.
    """
    # Sort by string length, descending
    input_seqs.reverse()

    # Start with longest subsequence
    scs = input_seqs[0]
    input_seqs.remove(scs)

    while len(input_seqs) > 0:
        longest_overlap = ""
        next_to_merge = ""

        for seq in input_seqs:

            if verbose:
                print "\n> Current shortest common supersequence:", scs
                print "> Comparing with:", seq

            if verbose:
                print "\n> Prepend"

            # Try prepending seq to scs
            prepend_overlap = find_overlap(seq, scs, True)

            if verbose:
                print "> Prepend '%s' to '%s': %s" % (seq, scs, prepend_overlap)

            if len(prepend_overlap) > len(longest_overlap):
                longest_overlap = prepend_overlap
                next_to_merge = seq

            if verbose:
                print "\n> Append"

            # Try appending seq to scs
            append_overlap = find_overlap(scs, seq, False)

            if verbose:
                print "> Append '%s' to '%s': %s" % (seq, scs, append_overlap)

            if len(append_overlap) > len(longest_overlap):
                longest_overlap = append_overlap
                next_to_merge = seq

        if next_to_merge == "": next_to_merge = seq
        scs = merge_seqs(scs, next_to_merge, longest_overlap)
        input_seqs.remove(next_to_merge)

    return scs

def find_overlap(seqA, seqB, prepend):
    """Finds the overlapping substring between SEQA and SEQB.

    If PREPEND = True, tries to prepend seqA to seqB. Otherwise, tries to
    append seqB to seqA.
    """
    if prepend:
        seqB = seqB[:len(seqA)]
    else:
        seqA = seqA[-len(seqB):]

    if verbose:
        print "%s\n%s\n" % (seqB, seqA)

    if prepend:
        if seqA == "" or seqA == seqB: return seqA
    else:
        if seqB == "" or seqB == seqA: return seqB

    return find_overlap(seqA[1:], seqB[:-1], prepend)

def merge_seqs(seqA, seqB, common):
    """Merge subsequences SEQA and SEQB together, where len(seqA) >= len(seqB),
    and the two strings have COMMON as their longest common subsequence.
    """
    if verbose:
        print "\n> Merging '%s' and '%s', which have '%s' in common." \
              % (seqA, seqB, common)

    if common == "":
        return seqA + seqB

    overlap_start_idx = seqA.find(common)
    return seqA[:overlap_start_idx] + seqB + seqA[overlap_start_idx + len(common):]

fname = sys.argv[1]
input_seqs = remove_substrings(read_input_seqs(fname))
print(find_scs(input_seqs))
