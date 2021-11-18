"""
These scripts are intended for extracting information
about RNA structures.

Author: Heli Mönttinen (ORCID: 0000-0003-2461-0690)

"""

from Bio.Seq import Seq
import common
from itertools import islice
import re


def _join_dot_file_lines(dot_file):
    """
    Joins a dot parenthesis file (dot parenthesis part of a file)
    into one string.

    dot_file:  A filepath to dot parenthesis

    Returns
    dot_parenthesis_string: a dot parenthesis part of
    a file as one string.
    """

    dot_parenthesis_string = ""
    start_options = ['.', '(', '[', ')', ']', '{', '}', '<', '>']

    with open(dot_file, 'r') as data_file:
        for line in islice(data_file, 1, None):
            if line.startswith(tuple(start_options)):
                dot_parenthesis_string += line.rstrip().split()[0]

    return dot_parenthesis_string


def identify_seq_fit_to_motif(sequence, motif):
    """
    Identifies if a given sequence fits to a given motif.

    Arguments
    ===========
    :Sequences: sequence of interest
    :motif: A motif to check

    Returns:
    True if a given motif is found from the given sequence.
    """

    iupac = {'N': {'U': "", 'G': "", 'T': "", 'A': "",
                   'C': "", 'R': "", 'Y': "", 'S': "",
                   'W': "", 'M': "", 'B': "", 'K': "",
                   'D': "", 'H': "", 'V': ""},
             'R': {'G': "", 'A': "", 'R': ""},
             'B': {'B': "", 'T': "", 'U': "", 'G': "",
                   'C': "", 'S': "", 'K': ""},
             'Y': {'C': "", 'T': "", 'U': "", 'Y': ""},
             'S': {'G': "", 'C': "", 'S': ""},
             'W': {'A': "", 'T': "", 'U': "", 'W': ""},
             'K': {'G': "", 'T': "", 'U': ""}}

    if len(sequence) != len(motif):
        return

    for i in range(len(sequence)):

        if sequence[i] == motif[i]:
            continue
        elif motif[i] in iupac and sequence[i] in iupac[motif[i]]:
            continue
        else:
            return
    return True


def _get_sequence_dot_parenthesis(dot_file):
    """
    Joins a sequence part of a dot parenthesis file into one sequence string.

    dot_file: a filepath to the dot parenthesis file

    Reuturns

    sequence: sequence as one string
    """

    sequence = ""

    exclude = ['(', ')', '.', '#', '[', ']', '{', '}']

    with open(dot_file, 'r') as data_file:
        for line in data_file:
            if not common.contains_any(line, exclude):
                sequence += line.rstrip()

    return sequence


def make_reverse_complement(seqs):
    """
    Takes a sequences (a list or string).
    and returns a reverse complement sequence as
    a string or a list dependening on the input.

    seqs: sequences as a string or a list of strings

    output: a reverse complement as a string or a list
    depending on the input

    if an input is not a string or a list it raises an error.
    """

    if isinstance(seqs, str):
        rna_seq = Seq(seqs)

        return rna_seq.reverse_complement()

    elif isinstance(seqs, list):

        reverse_rna_seqs = []

        for seq in seqs:
            rna_seq = Seq(seq, generic_rna)
            reverse_rna_seqs.append(rna_seq.reverse_complement())

        return reverse_rna_seqs

    else:
        raise TypeError(
                'Your Sequence cannot be recognized as a string or a list')


def loop_part_reverse_complement(full_seq, indexes):
    """
    Makes a reverse complement for a loop part of the RNA sequence.
    The indexes are from the full length sequence.

    Requires

    seqs: A sequence in which the loop is located.
    stem_loop_indexes: A indexes (4) for a loop and stems (start of stem,
                        start of a loop, end of the loop, end of the stem)

    Returns

    comp_reverse_loop: The input sequence in which the loop part is
                        reverse complement.

    """

    loop = full_seq[indexes[1]:indexes[2]]
    loop_seq = Seq(loop, generic_rna)
    reverse_loop = loop_seq.reverse_complement()

    return reverse_loop


class RNA_loops:
    """
    Functions for identifying and analysing loop sequences and
    secondary structures.
    """

    def __init__(self, dot_parenthesis, sequence, mode=None):
        """
        Instances for RNA Loops class.
        RNA_id: refers to identifier of a RNA file downloaded
                from the RNASTRAND 2.
        root_dir: A root_directory to look at.
        loop_sequences: Loop sequences as a list
        loop_indexes: Loop  indexes as a list
        """
        self.dot_parenthesis = dot_parenthesis
        self.sequence = sequence
        self.loop_sequences = []
        self.reverse_loop_sequences = []
        self.loop_indexes = []
        self.mode = mode

    def identify_loops(self):
        """
        Identify loops from dot-parenthesis files (for a given
        RNA_identifier). Identifies also in maximum three nucleotides
        around the loop if they belong to a surronding stem.
        The stem nucleotides a needed because they may anchor the
        loop in the sequence alignment.

        Requires
        self (self.dot_parenthesis)

        Result
        loop_indexes: a list of four indexes:
                        1) first nucleotide of the surrounding stem
                           (max -3 from start of the loop),
                        2) start of the loop,
                        3) end of the loop and
                        4) the last nucleotide of the surrounding stem
                           (max +3 from the end of the loop)
        """

        def _identify_stem_anchor(loop_start, loop_end):
            """
            Explores sequence to both directions from the
            given start and end index of the loop until
            a char indicating a neighboring loop is reached.

            Requires
                loop_start: The staring index of a loop
                loop_end: The last index of a loop

            Returns
                stem_first: The first nucleotide
                            of the left stem
                stem_last: The last nucleotide of the
                            latter stem

            """

            stem_first = loop_start
            stem_last = loop_end

            if self.mode is None:
                loop_start_char = '('
                loop_end_char = ')'
            elif self.mode == 'Stockholm':
                loop_start_char = '<'
                loop_end_char = '>'

            for index in range(1, len(self.dot_parenthesis)-1):

                if index == 1:
                    if ((self.dot_parenthesis[
                        loop_start-index] == loop_start_char)
                        and (self.dot_parenthesis[
                            loop_end-1+index] == loop_end_char)):

                        if loop_start-index > 0 and\
                                ((loop_end-1+index) <
                                 len(self.dot_parenthesis)-1):
                            stem_first = loop_start-index
                            stem_last = loop_end+index

                elif ((self.dot_parenthesis[
                    loop_start-index] == loop_start_char or
                            (self.dot_parenthesis[
                                loop_start-index] == '.')) and
                            (self.dot_parenthesis[
                                loop_end-1+index] ==
                                loop_end_char or
                                self.dot_parenthesis[loop_end-1+index]
                                == '.')):

                    if loop_start-index > 0 and\
                            ((loop_end-1+index) < len(self.dot_parenthesis)-1):

                        stem_first = loop_start-index
                        stem_last = loop_end+index

                    else:
                        break
                else:
                    break

            for i in range(stem_last):
                if self.dot_parenthesis[stem_last] == '.' and\
                        self.dot_parenthesis[stem_first] == '.':
                    stem_first = stem_first + 1
                    stem_last = stem_last - 1
                else:
                    break

            return stem_first, stem_last+1

        pattern = re.compile(r"(?<=\()[.-]+(?=\))")

        if self.mode == 'Stockholm':
            pattern = re.compile(r"(?<=\<)[._]+(?=\>)")

        for loop in re.finditer(pattern, self.dot_parenthesis):
            indexes = loop.span()

            loop_start = indexes[0]
            loop_end = indexes[1]

            stem_first, stem_last = _identify_stem_anchor(loop_start, loop_end)

            self.loop_indexes.append(
                    [stem_first, loop_start, loop_end, stem_last])

    def sequence_for_loops(self):
        """
        Creates a list of sequences for identied loop patterns with
        a small part of stem. Also creates a list of sequences in which
        a loop part is a reverse complement

        Requires:
        self (self.identify_loops(), and self.loop_indexes

        Saves resulting sequences of loops into self.loop_sequences list
        and self.reverse_loop_sequences.
        """

        self.identify_loops()

        for index in self.loop_indexes:
            loop_sequence = self.sequence[index[0]:index[3]]
            reverse_loop_sequence = make_reverse_complement(
                    self.sequence[index[0]:index[3]])
            reverse_str = str(
                    reverse_loop_sequence)  # Seq() object to str
            self.reverse_loop_sequences.append(reverse_str)
            self.loop_sequences.append(loop_sequence)
