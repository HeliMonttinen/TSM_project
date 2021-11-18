"""
A tool to calculate proportions z-test for
loop sequences and loop lengths. Also calculates
frequencies of background loop sequences and the ratio
to their reverse complements.
"""
import os
import sys
from Bio.Seq import Seq
from decimal import Decimal

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Takes a file containing loop sequences of interest
    and tests their significance to the full set of
    background loop sequences using propotions z-test.

    Argument
    ========
    input_freq: A file which contains counts for loop
                 sequences in the full data set.
                 The file has to be a .csv formatted in which
                 the sequence is in the column three and
                 counts for each loop sequence in the column four
    output : A file identifie for output files.
    immeadiate_loop: A file containing loop sequences of
                     interest.
                     The file has to be .csv formatted in which
                     loop sequences are located in the first column
                     and counts in the fourth column.

    Output
    ======
    Three output files 1) Loop length, 2) Loop motif,
    3) Ratio for loop sequence and its' reverse complement
       among the loop sequence data

    """

    from RNA_loops import identify_seq_fit_to_motif
    from mystatistics import two_proportions_z_test

    input_freq = sys.argv[1]
    output = sys.argv[2]
    immeadiate_loop = sys.argv[3]

    motifs_of_interest = ['GNRA', 'TNCG', 'CTTG', 'GNNA',
                          'GNAB', 'CNNG', 'UUUG', 'TYNC',
                          'CGNA', 'CAAG']

    motifs_of_interest_background = {}
    motifs_of_interest_sample = {}
    loop_sequence_sample = {}
    all_loops_in_sample = 0
    all_loops_background = 0
    loop_sequence_counts = {}
    loop_length = {}
    loop_length_back = {}

    with open(immeadiate_loop, 'r') as f:

        for line in f:
            line_splitted = line.rstrip().split('\t')
            sequence = line_splitted[0]
            count = int(line_splitted[1])

            all_loops_in_sample += count

            if sequence not in loop_sequence_sample:
                loop_sequence_sample[sequence] = count
            else:
                loop_sequence_sample[sequence] += count

            if len(sequence) in loop_length:
                loop_length[len(sequence)] += count
            else:
                loop_length[len(sequence)] = count

            for motif in motifs_of_interest:

                if identify_seq_fit_to_motif(sequence,
                                             motif) is True:

                    if motif in motifs_of_interest_sample:
                        motifs_of_interest_sample[motif] += count
                    else:
                        motifs_of_interest_sample[motif] = count

    with open(input_freq, 'r') as f:
        for line in f:

            line_splitted = line.split('\t')

            all_loops_background += int(line_splitted[3])

            if len(line_splitted[2]) in loop_length_back:
                loop_length_back[len(line_splitted[2])] +=\
                        int(line_splitted[3])
            else:
                loop_length_back[len(line_splitted[2])] =\
                        int(line_splitted[3])

            if line_splitted[2].replace('U', 'T')\
                    not in loop_sequence_counts:
                loop_sequence_counts[
                        line_splitted[2].replace('U', 'T')] =\
                                int(line_splitted[3])
            else:
                loop_sequence_counts[
                        line_splitted[2].replace('U', 'T')] +=\
                                int(line_splitted[3])

            for motif in motifs_of_interest:

                if identify_seq_fit_to_motif(
                        line_splitted[2].replace('U', 'T'),
                        motif) is True:

                    if motif in motifs_of_interest_background:
                        motifs_of_interest_background[motif] +=\
                                int(line_splitted[3])
                    else:
                        motifs_of_interest_background[motif] =\
                                int(line_splitted[3])

    for item in loop_length:

        motif_n_sample = loop_length[item]
        motif_n_back = loop_length_back[item]

        prop_sample = Decimal(motif_n_sample/all_loops_in_sample)
        prop_background = Decimal(motif_n_back/all_loops_background)
        result = two_proportions_z_test(
                motif_n_sample,
                motif_n_back,
                all_loops_in_sample,
                all_loops_background)

        with open(output + '_loop_length_immed_CS.txt', 'a+') as f:

            f.write(str(item) + '\t' + str(motif_n_sample) + '\t' +
                    str(motif_n_back) + '\t' + str(prop_sample) +
                    '\t' + str(prop_background) + '\t' + str(result[0])
                    + '\t' + str(result[1]) + '\n')

    for item in motifs_of_interest_sample:

        motif_n_sample = motifs_of_interest_sample[item]
        motif_n_back = motifs_of_interest_background[item]

        prop_sample = Decimal(motif_n_sample/all_loops_in_sample)
        prop_background = Decimal(motif_n_back/all_loops_background)
        result = two_proportions_z_test(
                motif_n_sample,
                motif_n_back,
                all_loops_in_sample,
                all_loops_background)

        with open(output + '_motif_immed_CS.txt', 'a+') as f:

            f.write(item + '\t' + str(motif_n_sample) + '\t'
                    + str(motif_n_back) + '\t' + str(prop_sample)
                    + '\t' + str(prop_background) + '\t' +
                    str(result[0]) + '\t' + str(result[1]) + '\n')

    loop_sequence_sample = sorted(loop_sequence_counts.items(),
                                  key=lambda x: x[1], reverse=True)

    for item, val in loop_sequence_sample:

        seq = Seq(item)

        reverse_seq = str(seq.reverse_complement())

        if reverse_seq not in loop_sequence_counts:
            reverse_seq_count = str(0)
            prop = str(0)
        else:
            reverse_seq_count = loop_sequence_counts[reverse_seq]
            prop = str(Decimal(int(val)/int(
                loop_sequence_counts[reverse_seq])))

        with open('reverse_complement_counts_170821.txt', 'a+') as f:

            f.write(item + '\t' + str(val) + '\t' + str(reverse_seq) +
                    '\t' + str(reverse_seq_count) + '\t' + prop + '\n')


if __name__ == "__main__":
    main()
