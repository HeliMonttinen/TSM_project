"""
Unifies a given set of sequences so that
all sequences have 'T' instead of 'U'.

Run the script:

    python3 unify_sequence.py input_dir output_dir


    input_dir: a directory on which the .fas-formatted
               input sequences are located
    output_dir: an output directory for the resulting
                dirctory

Author: Heli Mönttinen  (ORCID: 0000-0003-2461-0690)

"""

import os
import sys

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Goes through all the sequences and replaces 'u' with
    't'.
    """

    from common import (split_multiple_fasta_file,
                        filepaths_in_dir)

    directory = sys.argv[1]
    output_directory = sys.argv[2]

    for filepath in filepaths_in_dir(directory):

        filename = filepath.split('/')[-1]

        for fasta in split_multiple_fasta_file(filepath):

            fasta_split = fasta.split('\n')
            for line in fasta_split:
                if '>' in line:
                    identifier = line
                else:
                    seq_lines = line

                    seq_lines = seq_lines.replace('u', 't')
                    seq_lines = seq_lines.replace('U', 'T')

            with open(output_directory + '/' + filename, 'a+') as f:
                f.write(identifier + '\n')
                f.write(seq_lines + '\n')

if __name__ == "__main__":
    main()
