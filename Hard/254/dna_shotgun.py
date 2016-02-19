import sys

verbose = True

input_seqs = []

def read_sequences():
    input_fname = sys.argv[1]
    f = open(input_fname, 'r')
    input_seqs = [ line[:-1] for line in f ] # Strip \n

    if verbose:
        print "> Reading input sequences from file '%s'." % input_fname
        print "> ", input_seqs

read_sequences()
