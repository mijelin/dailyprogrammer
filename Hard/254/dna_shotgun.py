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

# Not necessary because find_scs() covers this.
# def remove_substrings(input_seqs):
#     """Remove subsequences that are wholly contained in another subsequence."""
#     numSeqs = len(input_seqs)
#     input_list = input_seqs[:]

#     for i in range(numSeqs):
#         for j in range(i + 1, numSeqs):
#             if input_seqs[i] in input_seqs[j]:
#                 if verbose:
#                     print "> Removing '%s' which is contained in '%s'." \
#                           % (input_seqs[i], input_seqs[j])
#                 input_list.remove(input_seqs[i])

#     return input_list

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

            for i in range(len(seq)):
                for j in range(len(scs)):

                    if seq[i] != scs[j]:
                        continue

                    overlap = ""
                    length = 0

                    while seq[i + length] == scs[j + length]:
                        overlap += seq[i + length]
                        length += 1

                        if (i + length) >= len(seq) or \
                           (j + length) >= len(scs): break

                    if verbose:
                        print "> Overlap is:", overlap

                    if len(overlap) > len(longest_overlap):
                        longest_overlap = overlap
                        next_to_merge = seq

                        if verbose:
                            print "> Longest overlap is now:", longest_overlap

        scs = merge_seqs(scs, next_to_merge, longest_overlap)
        input_seqs.remove(next_to_merge)

    return scs

def merge_seqs(seqA, seqB, common):
    """Merge subsequences SEQA and SEQB together, where len(seqA) >= len(seqB),
    and the two strings have COMMON as their longest common subsequence.
    """
    if verbose:
        print "\n> Merging '%s' and '%s', which have '%s' in common." \
              % (seqA, seqB, common)

    overlap_start_idx = seqA.find(common)
    return seqA[:overlap_start_idx] + seqB + seqA[overlap_start_idx + len(common):]

fname = sys.argv[1]
# input_seqs = remove_substrings(read_input_seqs(fname))
input_seqs = read_input_seqs(fname)
print(find_scs(input_seqs))
