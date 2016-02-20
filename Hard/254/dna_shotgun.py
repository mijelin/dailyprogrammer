import sys

verbose = True

input_seqs = []
lcs_table = {}

def read_input_seqs(input_fname):
    """Returns a list of input DNA sequences read from filename INPUT_FNAME."""
    f = open(input_fname, 'r')
    input_seqs = [ line[:-1] for line in f ] # Strip \n

    if verbose:
        print "> Reading input sequences from file '%s'." % input_fname
        print "> ", input_seqs
        print ""

    return input_seqs

def populate_lcs():
    """Populates a table of longest overlapping subsequences
    between two strings.
    """
    table = {}

    input_seqs.sort(key = lambda s: len(s))
    input_seqs.reverse()

    if verbose:
        print "> Sorted input list by length, descending:"
        print "> ", input_seqs
        print ""

    for i in range(len(input_seqs)):
        seqA = input_seqs[i]
        table[seqA] = {}

        for j in range(i, len(input_seqs)):
            seqB = input_seqs[j]
            table[seqA][seqB] = find_overlap(seqA, seqB)

    return table

def find_overlap(a, b):
    """Returns the longest overlapping substring between strings A and B."""
    if len(a) > len(b):
        longer = a
        shorter = b
    else:
        longer = b
        shorter = a

    longest_overlap = ""
    for i in range(len(shorter)):
        for j in range(len(longer)):

            if shorter[i] != longer[j]:
                continue

            overlap = ""
            length = 0

            while shorter[i + length] == longer[j + length]:
                overlap += shorter[i + length]
                length += 1
                if i + length >= len(shorter) or j + length >= len(longer):
                    break

            if len(overlap) > len(longest_overlap): longest_overlap = overlap

    return longest_overlap

def merge_seqs(seqA, seqB):
    if len(seqA) < len(seqB):
        a = seqB
        b = seqA
    else:
        a = seqA
        b = seqB

    common = lcs_table[a][b]
    overlap_start_idx = a.find(common)
    return a[:overlap_start_idx] + b + a[overlap_start_idx + len(common):]

# input_seqs = read_input_seqs(sys.argv[1])
# lcs_table = populate_lcs()
