"""
This package includes tools for analysing a RNA structure.

Author: Heli Mönttinen (ORCId: 0000-0003-2461-0690)

"""
from RNA_loops import (_join_dot_file_lines,
                       RNA_loops,
                       make_reverse_complement)
from tree_parse import mismatches_between_sequences
from RNA_alignments import (unusual_iupac_char,
                            indexes_in_original_seq,
                            indexes_in_alignment)


def identify_alignment_area(structure_file, align_dict,
                            index1, index2):
    """
    On the basis of given indexes the secondary structure
    is identified (dot-parenthesis).

    Arguments
    ==========
    :structure_file: A path to RNA fot-parenthesis file
    :alignment_dict: A dictionary which contains the aligned
                     sequence of the RNA
    :index1: The starting index of the region of interest
    :index2: The ending index of the region of interest

    Returns
    =======
    :ts_region: The extracted sequence of the region of
                interest
    :aligned_struct: The secondary structure in the region
                     if interest. The structure is aligned,
                     contains gaps.
    """

    structure = _join_dot_file_lines(structure_file)

    with open(structure_file, 'r') as f:
        for line in f:
            if '>' in line:
                if '#_' in line:
                    structure_id = line.lstrip('>').rstrip().split('_')[0]
                else:
                    structure_id = line.lstrip('>').rstrip()

    alignment_seq = align_dict[structure_id]

    aligned_struct = ""

    b = 0
    for a in range(len(alignment_seq)):
        if alignment_seq[a] == '-':

            aligned_struct = aligned_struct + '-'
        elif b < len(structure):
            aligned_struct = aligned_struct + structure[b]
            b += 1

    ts_region = aligned_struct[index1:index2]

    return ts_region, aligned_struct


def region_overlaps_with_loop(aligned_structure,
                              alignment_seq,
                              template_switch):
    """

    The function identifies if a given
    region overlaps with a RNA hairpin loop structure.

    Arguments
    ===========
    :aligned_structure: A RNA secondary structure in
                        a dot-parenthesis format. The struture
                        should be alinged (gaps added to
                        follow the sequence alignment).
    :alignment_seq: The aligned RNA sequence
    :template_switch: A set of indexes covering all the
                     residues in the region of interest.

    Returns:
    =======

    :loop_struct: The entire RNA hairpin loop structure
                  (which overlaps with the region of interest)
                  As a dot-parenthesis format.
    :temp_struct: That part of the RNA hairpin structure
                  which covers the region of interest
    :loop_part_seq: The sequence in the loop
    :loop_index_hit: The indexes for the hairpin loop as
                     a list [stem1 starting index, loop part
                     starting index, loop ending index,
                     stem2 ending index]
    :loop_overlap: All indexes as a set which overlap with
                   a loop region and the region of interest.
    :stem1_overlap: All stem1 indexes overlapping with
                    the region of interest as a set
    :stem2_overlap: All stem2 residue indexes overlapping with
                    the region of interest as a set
    :stem1_ind: All stem1 residue indexes as a set
    :stem2_ind: All stem2 residue indexes as a set
    :loop_ind: All loop residue indexes as a set
    :adj_nucleotide1: The adjacent nucleotide to the loop
                      in the stem1
    :adj_nucleotide2: The adjacent nucleotide to the loop
                      in the stem2.

    """

    RNA_seq_orig = alignment_seq.replace('-', '')
    orig_struct = aligned_structure.replace('-', '')

    RNA_sequence = RNA_loops(orig_struct,
                             RNA_seq_orig)

    RNA_sequence.sequence_for_loops()
    loop_inds = RNA_sequence.loop_indexes
    loop_indexes = list()

    for loop_x in loop_inds:
        loop_ind = indexes_in_alignment(loop_x,
                                        alignment_seq)
        loop_indexes.append(loop_ind)

    loop_sequence = RNA_sequence.loop_sequences

    loop_sequence = []
    for loop_i in loop_indexes:
        seq_l = alignment_seq[loop_i[0]:loop_i[-1]]
        loop_sequence.append(seq_l)

    loop_part_seq = ""
    loop_overlap = set()
    stem1_overlap = set()
    stem2_overlap = set()
    loop_ind = set()
    stem2_ind = set()
    stem1_ind = set()
    loop_index_hit = None
    temp_struct = None
    loop_struct = None
    num = 0
    adj_nucleotide1 = ""
    adj_nucleotide2 = ""

    for index in loop_indexes:

        hairpin_intersection = set(list(range(index[0], index[3])))
        set_hairpin_intersection =\
            hairpin_intersection.intersection(template_switch)
        if len(set_hairpin_intersection) < (len(template_switch)/2):
            num += 1

            continue

        loop_index = set(list(range(index[1], index[2])))
        set_loop = loop_index.intersection(template_switch)
        loop_ind1 = index[1]-index[0]
        loop_ind2 = index[2]-index[0]

        if len(set_loop) > 0:

            loop_overlap = set_loop
            loop_index_hit = index

            loop_ind = set(list(range(index[1], index[2])))
            loop_part_seq = loop_sequence[num][loop_ind1:loop_ind2]
            adj_number1 = loop_ind1
            adj_number2 = loop_ind2-1
            for i in range(len(loop_sequence[num])):
                adj_number1 -= 1
                adj_nucleotide1 = loop_sequence[num][adj_number1]
                if adj_nucleotide1 != '-':
                    break

            for i in range(len(loop_sequence[num])):
                adj_number2 += 1
                adj_nucleotide2 = loop_sequence[num][adj_number2]
                if adj_nucleotide2 != '-':
                    break

        stem_1_indexes = set(list(range(index[0], index[1])))

        set_stem1 = stem_1_indexes.intersection(template_switch)
        stem1_ind = stem_1_indexes
        if len(set_stem1) > 0:

            stem1_overlap = set_stem1

            if loop_index_hit is None:
                loop_index_hit = index
                loop_part_seq = loop_sequence[num][loop_ind1:loop_ind2]
                adj_number1 = loop_ind1
                adj_number2 = loop_ind2-1
                for i in range(len(loop_sequence[num])):
                    adj_number1 -= 1
                    adj_nucleotide1 = loop_sequence[num][adj_number1]
                    if adj_nucleotide1 != '-':
                        break

                for i in range(len(loop_sequence[num])):
                    adj_number2 += 1
                    adj_nucleotide2 = loop_sequence[num][adj_number2]
                    if adj_nucleotide2 != '-':
                        break

        stem_2_indexes = set(list(range(index[2], index[3])))

        set_stem2 = stem_2_indexes.intersection(template_switch)
        stem2_ind = stem_2_indexes
        if len(set_stem2) > 0:

            stem2_overlap = set_stem2
            if loop_index_hit is None:
                loop_index_hit = index
                loop_part_seq = loop_sequence[num][loop_ind1:loop_ind2]
                adj_number1 = loop_ind1
                adj_number2 = loop_ind2-1

                for i in range(len(loop_sequence[num])):
                    adj_number1 -= 1
                    adj_nucleotide1 = loop_sequence[num][adj_number1]
                    if adj_nucleotide1 != '-':
                        break

                for i in range(len(loop_sequence[num])):
                    adj_number2 += 1
                    adj_nucleotide2 = loop_sequence[num][adj_number2]
                    if adj_nucleotide2 != '-':
                        break

        num += 1
        if len(set_loop) > 0:
            break

    if loop_index_hit is not None:
        minimum = min(template_switch)
        maximum = max(template_switch)
        temp_struct = aligned_structure[minimum:maximum]
        loop_struct = aligned_structure[loop_index_hit[0]: loop_index_hit[3]]

    return (loop_struct,
            temp_struct,
            loop_part_seq,
            loop_index_hit,
            loop_overlap,
            stem1_overlap,
            stem2_overlap,
            stem1_ind,
            stem2_ind,
            loop_ind,
            adj_nucleotide1,
            adj_nucleotide2)


def identify_if_stem_is_extended(
        loop_index_hit_qry,
        loop_index_hit_ref,
        ref_overlaps,
        qry_overlaps,
        ref_ts_start,
        ref_ts_end,
        qry_ts_start,
        qry_ts_end):
    """
    Compares the structures of reference and query
    sequences in the region of interest. Prints out
    the structure change in a string format.

    Arguments
    =========
    :loop_index_hit_qry: A list of loop indexes in
                         the child seq [stem1_start, loop_start,
                         loop_end, stem2_end]
    :loop_index_hit_ref: A list of loop indeces in the parental
                         seq (structure is similar to a
                         child seq)
    :ref_overlaps: A list, which contains three sets of indexes
                   for the overlapping region between the hairpin
                   and the region of interet.
                   The first includes the indeces for the
                   overlapping region in the stem1

    :qry_overlaps: A list containing three sets of indexes.
    :template_switch_start: The starting index for the
                            region of interest
    :temp_struct_ref: The structure in the region of
                      interest
    :temp_struct_qry: The structure in the region of interest

    Returns
    =======
    :output: The structural changes as a string

    """

    output = []

    qry_loop = set(range(int(
        loop_index_hit_qry[1]),
        int(loop_index_hit_qry[2])))
    ref_loop = set(range(
        int(loop_index_hit_ref[1]),
        int(loop_index_hit_ref[2])))

    qry_stem1 = set(range(int(
        loop_index_hit_qry[0]),
        int(loop_index_hit_qry[1])))
    qry_stem2 = set(range(int(
        loop_index_hit_qry[2]),
        int(loop_index_hit_qry[3])))

    if min(qry_loop) < min(ref_loop) or\
            max(qry_loop) > max(ref_loop):
        return output

    loops_overlap = qry_loop.intersection(ref_loop)

    if (len(qry_loop) < len(ref_loop)) and\
        (len(loops_overlap) == len(qry_loop)) and len(qry_loop) > 0 and\
        ((ref_ts_start in qry_stem1 and qry_ts_end in qry_stem2) or
         (qry_ts_start in qry_stem1 and
          ref_ts_end in qry_stem2)) and\
        ((qry_ts_start in ref_loop and ref_ts_end in ref_loop) or
         (qry_ts_end in ref_loop and ref_ts_start in ref_loop)):

        output.append('A stem is extended toward a loop.')

    elif loop_index_hit_qry[0] is not None and\
            loop_index_hit_ref[0] is not None:

        if (int(loop_index_hit_qry[0]) < int(loop_index_hit_ref[0])) and\
            (len(loops_overlap) == len(qry_loop)) and\
            len(qry_loop) > 0 and\
            (qry_ts_start < loop_index_hit_ref[0] or
             ref_ts_end > loop_index_hit_ref[3]) and\
            ((qry_ts_start in qry_stem1 and ref_ts_end in qry_stem2 and
              ref_ts_start not in qry_stem1 and
              qry_ts_end not in qry_stem2)):
            output.append('A stem is extended at the base.')

    elif (int(loop_index_hit_qry[2]) >
            int(loop_index_hit_ref[2])) and\
            (len(qry_overlaps[2]) > len(ref_overlaps[2])) and\
            (qry_ts_end > loop_index_hit_ref[2] or
             ref_ts_start < loop_index_hit_ref[0]) and\
            (len(loops_overlap) == len(qry_loop)) and\
            len(qry_loop) > 0 and\
            ((ref_ts_start in qry_stem1 and
              qry_ts_end in qry_stem2 and ref_ts_end not in qry_stem2 and
                qry_ts_start not in qry_stem1)):

        output.append('A stem is extended at the base.')

    return output


def identify_if_structure_has_changed(
        loop_index_hit_qry,
        loop_index_hit_ref,
        ref_overlaps,
        qry_overlaps,
        template_switch_start,
        temp_struct_ref,
        temp_struct_qry,
        ref_ts_start,
        ref_ts_end,
        qry_ts_start,
        qry_ts_end):
    """
    Compares the structures of reference and query
    sequences in the region of interest. Prints out
    the structure change in a string format.

    Arguments
    =========
    :loop_index_hit_qry: A list of loop indexes in
                         the child seq [stem1_start, loop_start,
                         loop_end, stem2_end]
    :loop_index_hit_ref: A list of loop indeces in the parental
                         seq (structure is similar to a
                         child seq)
    :ref_overlaps: A list, which contains three sets of indexes
                   for the overlapping region between the hairpin
                   and the region of interet.
                   The first includes the indeces for the
                   overlapping region in the stem1

    :qry_overlaps: A list containing three sets of indexes.
    :template_switch_start: The starting index for the
                            region of interest
    :temp_struct_ref: The structure in the region of
                      interest
    :temp_struct_qry: The structure in the region of interest

    Returns
    =======
    :output: The structural changes as a string

    """

    output = []
    position = 0
    new_index = []

    if loop_index_hit_qry is None or\
            loop_index_hit_ref is None:
        output.append('A structure does not match a loop.')
        output = set(output)
        return output

    qry_stem1 = set(range(int(
        loop_index_hit_qry[0]),
        int(loop_index_hit_qry[1])))

    qry_stem2 = set(range(int(
        loop_index_hit_qry[2]),
        int(loop_index_hit_qry[3])))

    qry_loop = set(range(int(
        loop_index_hit_qry[1]),
        int(loop_index_hit_qry[2])))

    ref_loop = set(range(
        int(loop_index_hit_ref[1]),
        int(loop_index_hit_ref[2])))

    for ind in loop_index_hit_qry:
        new_ind = ind - template_switch_start
        new_index.append(new_ind)

    new_ref_overlaps = []

    for set1 in ref_overlaps:
        new_set = set()

        if set1 != 0:
            for ind in set1:
                new_ind = ind - template_switch_start
                new_set.add(new_ind)
            new_ref_overlaps.append(new_set)

        else:
            new_set.add(None)
            new_ref_overlaps.append(new_set)

    new_qry_overlaps = []

    for set1 in qry_overlaps:
        new_set = set()

        if set1 != 0:
            for ind in set1:
                new_ind = ind - template_switch_start
                new_set.add(new_ind)
            new_qry_overlaps.append(new_set)
        else:
            new_set.add(None)
            new_qry_overlaps.append(new_set)

    loops_overlap = qry_loop.intersection(ref_loop)

    if len(loops_overlap) < 3:
        output.append('A major change')

    else:

        for a, b in zip(temp_struct_ref, temp_struct_qry):

            if a != '-' and b != '-' and a != b:

                if position in new_qry_overlaps[0] and\
                        position in new_ref_overlaps[0] and\
                        temp_struct_ref[position] == '.' and\
                        temp_struct_qry[position] == '(':

                    output.append('A mismatch in stem1 fixed')

                elif position in new_qry_overlaps[2] and\
                        position in new_ref_overlaps[2] and\
                        temp_struct_ref[position] == '.' and\
                        temp_struct_qry[position] == ')':

                    output.append('A mismatch in stem2 fixed.')

                elif position in new_qry_overlaps[0] and\
                        position in new_ref_overlaps[0] and\
                        temp_struct_ref[position] == '(' and\
                        temp_struct_qry[position] == '.':

                    output.append('A new mismatch in stem1.')

                elif position in new_qry_overlaps[2] and\
                        position in new_ref_overlaps[2] and\
                        temp_struct_ref[position] == ')' and\
                        temp_struct_qry[position] == '.':

                    output.append('A new mismatch in stem2.')

                elif ((position in new_qry_overlaps[1] and
                       position in new_ref_overlaps[0]) or
                      (position in new_qry_overlaps[1] and
                       position in new_ref_overlaps[2])) and\
                        len(qry_loop) > len(ref_loop):

                    output.append('A loop has extended.')

                elif ((position in new_qry_overlaps[0] and
                       position in new_ref_overlaps[1]) or
                      (position in new_qry_overlaps[2] and
                       position in new_ref_overlaps[1])) and\
                        len(qry_loop) < len(ref_loop):

                    output.append('A loop has shorten.')

                elif (((position in new_qry_overlaps[0] and
                        position in new_ref_overlaps[1]) or
                       (position in new_qry_overlaps[2] and
                        position in new_ref_overlaps[1])) or
                        ((position in new_qry_overlaps[1] and
                          position in new_ref_overlaps[0]) or
                         (position in new_qry_overlaps[1] and
                          position in new_ref_overlaps[2]))) and\
                        len(qry_loop) == len(ref_loop):
                    output.append('A loop position has shifted.')

            elif a != b and\
                ((ref_ts_start in qry_stem1 and
                  qry_ts_end in qry_stem2 and
                  ref_ts_end not in qry_stem2 and
                  qry_ts_start not in qry_stem1) or
                 (ref_ts_end in qry_stem2 and
                  qry_ts_start in qry_stem1 and
                  qry_ts_end not in qry_stem2 and
                  ref_ts_start not in qry_stem1)):
                chr_set = set(list(temp_struct_qry.replace('-', '')))
                if position in new_qry_overlaps[0] and\
                        position in new_ref_overlaps[0] and\
                        temp_struct_ref[position] == '-' and\
                        temp_struct_qry[position] == '(' and\
                        len(chr_set) == 1:

                    output.append('A stem is extended')

                elif position in new_qry_overlaps[2] and\
                        position in new_ref_overlaps[2] and\
                        temp_struct_ref[position] == '-' and\
                        temp_struct_qry[position] == ')' and\
                        len(chr_set) == 1:

                    output.append('A stem is extended')

            position += 1
    output = set(output)

    return output


def identify_immediate_CS(orig_struct_ref,
                          orig_struct_qry,
                          ref_aligned_seq,
                          qry_aligned_seq,
                          regions_dict,
                          child1=None,
                          child2=None):
    """
    This script identifies the instantaneous CS (i.e. CS
    observed in the very same sequence with the initial
    substitution on the opposite stem).
    The script identifies immeadiate CS only from those loops
    that are in the same position in both sequences,

    Arguments
    =========

    :orig_struct_ref: aligned reference (parental) structure
    :orig_struct_qry: aligned query (descendant) structure
    :ref_aligned_seq: aligned reference (parental) sequence
    :qry_aligned_seq: aligned query (descendant) sequence
    :ancestor_seq: The sequence of the parent's ancestor
    :sister_seq: The sister sequence of the query sequence

    Returns:
    ========

    :count: Count for instantaneus CSs between parent and child
    :RNA_loop: Counts for instantaneous CS events per loop sequence
               (dictionary)
    :adj_nucleotide: The adjacent nucleotides to the loop
    :rev_comp_loop: The count for reverse complement loops
    :mutations: A number of mutations involved in one
                immediate CM affected event. A dictionary.
    :lengths: A dictionary containing lengths for the
              perfect Watson-Crick base-pairing region in which
              the instantaneous CM region is located.
    :mutations_type: A dictionary which contains counts for
                     different base-pairs mutations in stem.
                     The key is a tuple. The first one indicates the
                     original nucleotide pair and the second is the
                     mutated one.

    """

    qry_aligned_seq = qry_aligned_seq.upper()
    ref_aligned_seq = ref_aligned_seq.upper()
    RNA_seq_orig_qry = qry_aligned_seq.replace('-', '')
    RNA_seq_orig_ref = ref_aligned_seq.replace('-', '')

    RNA_sequence_qry = RNA_loops(orig_struct_qry,
                                 RNA_seq_orig_qry)

    RNA_sequence_ref = RNA_loops(orig_struct_ref,
                                 RNA_seq_orig_ref)

    RNA_sequence_qry.identify_loops()
    RNA_sequence_ref.identify_loops()

    loop_inds_qry = RNA_sequence_qry.loop_indexes
    loop_inds_ref = RNA_sequence_ref.loop_indexes

    loop_ind_qry_dict = dict()

    for x in range(len(loop_inds_qry)):

        try:
            loop_ind_qry = indexes_in_alignment(loop_inds_qry[x],
                                                qry_aligned_seq)

            loop_ind_qry_dict[
                str(loop_ind_qry[1]) + ',' + str(loop_ind_qry[2])] =\
                loop_ind_qry

        except:
            continue
    loop_ind_ref_dict = dict()

    for x in range(len(loop_inds_ref)):
        try:
            loop_ind_ref = indexes_in_alignment(
                    loop_inds_ref[x],
                    ref_aligned_seq)
            loop_ind_ref_dict[
                str(loop_ind_ref[1]) + ',' + str(loop_ind_ref[2])] =\
                loop_ind_ref
        except:
            continue

    nucl = {'A': '', 'T': '', 'C': '', 'G': '', 'U': ''}

    count = {}
    adj_nucleotide = dict()
    RNA_loop = dict()
    rev_comp_loop = dict()
    length = dict()
    mutations = dict()
    mutations_type = dict()
    loop_inds = dict()
    rounds = 0

    for loop_ind in loop_ind_qry_dict:
        if loop_ind in loop_ind_ref_dict and\
            len(loop_ind_qry_dict[loop_ind]) == 4 and\
                len(loop_ind_ref_dict[loop_ind]) == 4:

            flag_qry = False

            for reg in regions_dict:
                reg1 = int(reg.split('_')[0])
                reg2 = int(reg.split('_')[1])

                int_first = int(loop_ind_qry_dict[loop_ind][0])
                int_last = int(loop_ind_qry_dict[loop_ind][1])

                if len(set(range(int(int_first), int(int_last))).intersection(
                    set(range(reg1, reg2)))) ==\
                        len(range(int(int_first), int(int_last))):

                    flag_qry = True

                if flag_qry is False:
                    continue
        else:
            continue

        loop_qry_ind_orig = indexes_in_original_seq(
                loop_ind_qry_dict[loop_ind],
                qry_aligned_seq)

        RNA_seq_qry_stem1 = RNA_seq_orig_qry[
                loop_qry_ind_orig[0]:loop_qry_ind_orig[1]]

        loop_ref_ind_orig = indexes_in_original_seq(
                loop_ind_ref_dict[loop_ind],
                ref_aligned_seq)

        RNA_seq_ref_stem1 = RNA_seq_orig_ref[
                loop_ref_ind_orig[0]:loop_ref_ind_orig[1]]
        RNA_seq_qry_stem2 = RNA_seq_orig_qry[
                loop_qry_ind_orig[2]:loop_qry_ind_orig[3]]
        RNA_seq_ref_stem2 = RNA_seq_orig_ref[
                loop_ref_ind_orig[2]:loop_ref_ind_orig[3]]

        if RNA_seq_qry_stem1 != RNA_seq_ref_stem1 and\
                RNA_seq_qry_stem2 != RNA_seq_ref_stem2:
            rev_RNA_seq_qry_stem1 = RNA_seq_qry_stem1[::-1]
            rev_RNA_seq_ref_stem1 = RNA_seq_ref_stem1[::-1]
            RNA_seq_qry_loop = RNA_seq_orig_qry[
                loop_qry_ind_orig[1]:loop_qry_ind_orig[2]]
            RNA_seq_ref_loop = RNA_seq_orig_ref[
                loop_ref_ind_orig[1]:loop_ref_ind_orig[2]]
            minimum = min(set([
                len(rev_RNA_seq_qry_stem1), len(RNA_seq_qry_stem2),
                len(rev_RNA_seq_ref_stem1), len(RNA_seq_ref_stem2)]))

            if child1 is not None:
                RNA_seq_child1_stem1 = child1[
                    loop_qry_ind_orig[0]:loop_qry_ind_orig[1]]
                RNA_seq_child2_stem1 = child2[
                    loop_qry_ind_orig[0]:loop_qry_ind_orig[1]]
                RNA_seq_child1_stem2 = child1[
                    loop_qry_ind_orig[2]:loop_qry_ind_orig[3]]
                RNA_seq_child2_stem2 = child2[
                    loop_qry_ind_orig[2]:loop_qry_ind_orig[3]]
                ch1_stem1_rev = make_reverse_complement(RNA_seq_child1_stem1)
                ch2_stem1_rev = make_reverse_complement(RNA_seq_child2_stem1)

            found = False
            mismatch_count = 0
            length_c = 0

            for ind in range(minimum):

                rev_qry_stem2_nucl = make_reverse_complement(
                    rev_RNA_seq_qry_stem1[ind])
                rev_ref_stem2_nucl = make_reverse_complement(
                    rev_RNA_seq_ref_stem1[ind])

                if rev_RNA_seq_qry_stem1[ind] !=\
                        rev_RNA_seq_ref_stem1[ind] and\
                        (rev_RNA_seq_qry_stem1[ind] in nucl) and\
                        (rev_RNA_seq_ref_stem1[ind] in nucl):

                    rev_qry_stem2_full = make_reverse_complement(
                        RNA_seq_qry_stem1)[0:ind]
                    rev_ref_stem2_full = make_reverse_complement(
                        RNA_seq_ref_stem1)[0:ind]
                    if RNA_seq_qry_stem2[ind] == rev_qry_stem2_nucl and\
                        RNA_seq_ref_stem2[ind] == rev_ref_stem2_nucl and\
                            rev_qry_stem2_full == RNA_seq_qry_stem2[0:ind] and\
                            rev_ref_stem2_full == RNA_seq_ref_stem2[0:ind]:
                        if child1 is not None:

                            if RNA_seq_child1_stem2[0:ind] !=\
                                RNA_seq_qry_stem2[0:ind] and\
                                RNA_seq_child2_stem2[0:ind] !=\
                                RNA_seq_qry_stem2[0:ind] and\
                                    ch1_stem1_rev[0:ind] !=\
                                    RNA_seq_qry_stem2[0:ind] and\
                                    ch2_stem1_rev[0:ind] !=\
                                    RNA_seq_qry_stem2[0:ind]:
                                break

                        if found is False:
                            nuclpairs = str(rev_RNA_seq_ref_stem1[ind]) +\
                                        str(RNA_seq_ref_stem2[ind]) + ',' +\
                                        str(rev_RNA_seq_qry_stem1[ind]) +\
                                        str(RNA_seq_qry_stem2[ind])
                            if rounds not in mutations_type:
                                mutations_type[rounds] = dict()
                            if nuclpairs in mutations_type[rounds]:
                                mutations_type[rounds][nuclpairs] += 1
                                loop_inds[rounds] = str(loop_ind)
                            else:
                                mutations_type[rounds][nuclpairs] = 1
                                loop_inds[rounds] = str(loop_ind)

                            mismatch_count += 1
                            found = True
                            adj_nucl1 = rev_RNA_seq_qry_stem1[0]
                            adj_nucl2 = RNA_seq_qry_stem2[0]
                            rev_qry_loop = make_reverse_complement(
                                    RNA_seq_qry_loop)
                            if rev_qry_loop == RNA_seq_ref_loop and\
                                    rev_qry_loop != RNA_seq_qry_loop:

                                if rounds not in rev_comp_loop:
                                    rev_comp_loop[rounds] = dict()

                                if RNA_seq_qry_loop in rev_comp_loop[rounds]:
                                    rev_comp_loop[rounds][
                                            RNA_seq_qry_loop] += 1
                                else:
                                    rev_comp_loop[rounds][RNA_seq_qry_loop] = 1
                            count[rounds] = 1
                            if rounds not in adj_nucleotide:
                                adj_nucleotide[rounds] = dict()
                                RNA_loop[rounds] = dict()

                            if adj_nucl1 + adj_nucl2 in adj_nucleotide[rounds]:
                                adj_nucleotide[rounds][
                                        adj_nucl1 + adj_nucl2] += 1
                            else:
                                adj_nucleotide[rounds][
                                        adj_nucl1 + adj_nucl2] = 1
                            if RNA_seq_ref_loop in RNA_loop[rounds]:
                                RNA_loop[rounds][RNA_seq_ref_loop] += 1
                            else:
                                RNA_loop[rounds][RNA_seq_ref_loop] = 1

                            for x in range(1, minimum):
                                rev_qry_stem2_full = make_reverse_complement(
                                        RNA_seq_qry_stem1)[:x]
                                rev_ref_stem2_full = make_reverse_complement(
                                        RNA_seq_ref_stem1)[:x]
                                if rev_qry_stem2_full == RNA_seq_qry_stem2[:x]\
                                    and\
                                        rev_ref_stem2_full ==\
                                        RNA_seq_ref_stem2[:x]:
                                    length_c = len(rev_qry_stem2_full)

                                elif found is True:
                                    mismatch_count += 1
                                    nuclpairs =\
                                        str(rev_RNA_seq_ref_stem1[ind]) +\
                                        str(RNA_seq_ref_stem2[ind]) + ',' +\
                                        str(rev_RNA_seq_qry_stem1[ind]) +\
                                        str(RNA_seq_qry_stem2[ind])

                                    if nuclpairs in mutations_type[rounds]:
                                        mutations_type[rounds][nuclpairs] += 1
                                    else:
                                        mutations_type[rounds][nuclpairs] = 1

                            else:
                                break

                if length_c > 0 and length_c == ind:

                    break

            if rounds not in mutations:
                mutations[rounds] = dict()

            if mismatch_count in mutations[rounds] and found is True:
                mutations[rounds][mismatch_count] += 1
            elif found is True:
                mutations[rounds][mismatch_count] = 1

            if rounds not in length:
                length[rounds] = dict()

            if length_c in length[rounds] and found is True:
                length[rounds][length_c] += 1
                loop_inds[rounds] = loop_inds[rounds] + '_' +\
                    RNA_seq_qry_stem2[:(length_c)]
                loop_inds[rounds] = loop_inds[rounds] + '_' +\
                    RNA_seq_ref_stem2[:(length_c)]
            elif found is True:
                length[rounds][length_c] = 1
                loop_inds[rounds] = loop_inds[rounds] + '_' +\
                    RNA_seq_qry_stem2[:(length_c)]
                loop_inds[rounds] = loop_inds[rounds] + '_' +\
                    RNA_seq_ref_stem2[:(length_c)]

        rounds += 1

    return (count, RNA_loop, adj_nucleotide, rev_comp_loop,
            mutations, length, mutations_type, loop_inds)


def identify_compensating_mutations(ref_ts_start,
                                    ref_ts_end,
                                    qry_align_seq,
                                    ref_align_seq,
                                    loop_index_hit_qry,
                                    loop_index_hit_ref,
                                    ref_overlaps,
                                    qry_overlaps,
                                    temp_struct_ref,
                                    temp_struct_qry,
                                    struct_output):
    """
    This script aims to identify the compensating mutations
    associated with template switches.

    Arguments
    ==========
    :ref_ts_start: A TSM source start index in the aligned
                   sequence
    :ref_ts_end: A TSM source index in the aligned sequence
    :qry_align_seq: The child sequence aligned
    :ref_align_seq: The parental sequence aligned
    :loop_index_hit_qry: A list of RNA hairpin indeces in
                         a child seq
                        [stem1 start index, loop start index,
                         loop end index, stem2 end index)
    :loop_index_hit_ref: A list of RNA hairpin indeces in
                        a parental seq
    :ref_overlaps: A list containing the three sets of indexes
                   indicating (in reference sequence)
                   1. TSM overlaps with stem1
                   2. TSM overlaps with a loop and
                   3. TSM overlaps with stem2
    :qry_overlaps: A similar list as above but for the
                   child sequence
    :temp_struct_ref: The structure in the tsm region of the
                      parental sequence
    :temp_struct_qry: The same as above but in the child
                      sequence
    :struct_output: The output from the 'identify_structural_change'
                    function.

    Returns:

    :diff_set: The number of  explained by a TSM on the stem1
    :first_nucleotide: The first nucleotide of the loop (on left)
    :adj_loop_nucl_mutated: True if the nucleotide next to stem is
                            mutated due to TSM
    :diff_set: The number of  explained by a TSM on the stem2
    :mutation_types: The nucleotide pair in the parental and
                     child sequences that has mutated.
                     It is presented in the form parent,child

    """

    ts_or_min = ref_ts_start
    ts_or_max = ref_ts_end

    accepted_mutations = {"GT,GC": "", "GT,AT": "", "TG,CG": "",
                          "TG,TA": "", "TT,TA": "", "TT,AT": "", "AA,TA": "",
                          "AA,AT": "", "GG,CG": "", "GG,GC": "", "CC,CG": "",
                          "CC,GC": "", "AC,AT": "", "AC,GC": "", "CA,CG": "",
                          "CA,TA": "", "GA,GC": "", "GA,TA": "", "AG,CG": "",
                          "AG,AT": "", "CT,CG": "", "CT,AT": "", "TC,GC": "",
                          "TC,TA": ""}

    mutation_types = {}
    ts_orig = set(range(ts_or_min, ts_or_max))

    ts_orig_overlap_stem2 = set(range(
        loop_index_hit_ref[2],
        loop_index_hit_ref[3])).intersection(set(ts_orig))

    ts_orig_overlap_stem1 = set(range(
        loop_index_hit_ref[0],
        loop_index_hit_ref[1])).intersection(set(ts_orig))

    loop_start = loop_index_hit_qry[1]
    loop_end = loop_index_hit_qry[2]
    first_nucleotide = False
    adj_loop_nucl_mutated = False
    stem1_seq_qry = ""
    stem1_seq_ref = ""
    stem2_seq_qry = ""
    stem2_seq_ref = ""

    corresponding_res = set()
    (qry_align_seq,
     ref_align_seq) = unusual_iupac_char(qry_align_seq,
                                         ref_align_seq)

    if len(qry_overlaps[0]) != 0 and ref_ts_start > min(qry_overlaps[0]):

        qry_overlaps_min = min(qry_overlaps[0])
        qry_overlaps_max = max(qry_overlaps[0])

        stem1_seq_qry = qry_align_seq[qry_overlaps_min:(qry_overlaps_max+1)]
        stem1_seq_ref = ref_align_seq[qry_overlaps_min:(qry_overlaps_max+1)]
        if qry_overlaps_max+1 == loop_start:
            first_nucleotide = True
            qry_start_seq = qry_align_seq[loop_start:].replace('-', '')
            ref_start_seq = ref_align_seq[loop_start:].replace('-', '')

            if qry_start_seq[0] != ref_start_seq[0]:
                adj_loop_nucl_mutated = True

        stem1_loop_indexes = []
        stem2_loop_indexes = []

        loop_start_or = indexes_in_original_seq([loop_start],
                                                qry_align_seq)
        for index in qry_overlaps[0]:

            if qry_align_seq[index] != '-':

                ind_or = indexes_in_original_seq([index],
                                                 qry_align_seq)

                fixed_index = abs(ind_or[0] - loop_start_or[0])
                stem1_loop_indexes.append(fixed_index)
        stem1_loop_indexes.sort()

        loop_end_or = indexes_in_original_seq([loop_end],
                                              ref_align_seq)
        for index in ts_orig_overlap_stem2:

            if ref_align_seq[index] != '-':
                ind_or = indexes_in_original_seq([index],
                                                 ref_align_seq)
                fixed_index = abs(ind_or[0] - (loop_end_or[0]-1))
                stem2_loop_indexes.append(fixed_index)
        stem2_loop_indexes.sort()

        corresponding_res = set(
            stem1_loop_indexes).intersection(stem2_loop_indexes)

    elif len(qry_overlaps[2]) != 0 and (max(qry_overlaps[2]) > ref_ts_end):

        qry_overlaps_min = min(qry_overlaps[2])
        qry_overlaps_max = max(qry_overlaps[2])

        stem2_seq_qry = qry_align_seq[qry_overlaps_min:(qry_overlaps_max+1)]
        stem2_seq_ref = ref_align_seq[qry_overlaps_min:(qry_overlaps_max+1)]

        if qry_overlaps_min == loop_end:
            first_nucleotide = True
            qry_seq_end = qry_align_seq[:(loop_end)].replace('-', '')
            ref_seq_end = ref_align_seq[:(loop_end)].replace('-', '')

            if qry_seq_end[-1] != ref_seq_end[-1]:
                adj_loop_nucl_mutated = True

        stem1_loop_indexes = []
        stem2_loop_indexes = []
        loop_end_or = indexes_in_original_seq([loop_end],
                                              qry_align_seq)

        for index in qry_overlaps[2]:

            if qry_align_seq[index] != '-':
                ind_or = indexes_in_original_seq([index],
                                                 qry_align_seq)

                fixed_index = abs(ind_or[0] - (loop_end_or[0]-1))
                stem2_loop_indexes.append(fixed_index)
        stem2_loop_indexes.sort()

        loop_start_or = indexes_in_original_seq([loop_start],
                                                ref_align_seq)

        for index in ts_orig_overlap_stem1:

            if ref_align_seq[index] != '-':
                ind_or = indexes_in_original_seq([index],
                                                 ref_align_seq)
                fixed_index = abs(ind_or[0] - (loop_start_or[0]))
                stem1_loop_indexes.append(fixed_index)
        stem1_loop_indexes.sort()

        corresponding_res = set(
            stem1_loop_indexes).intersection(stem2_loop_indexes)

    if len(qry_overlaps[0]) > 0 and len(ref_overlaps[0]) > 0 and\
            stem1_seq_qry != stem1_seq_ref and\
            stem1_loop_indexes == stem2_loop_indexes and\
            len(corresponding_res) == len(stem1_loop_indexes) and\
            ('A mismatch in stem1 fixed' in struct_output or
             'No change in structure' in struct_output):

        (stem1_seq_qry,
         stem1_seq_ref) = unusual_iupac_char(stem1_seq_qry,
                                             stem1_seq_ref)

        diff_set = 0
        counter_a = -1
        counter_b = -1
        st_index = 0
        if max(qry_overlaps[0]) < loop_index_hit_qry[1]-1:
            st_index = (loop_index_hit_qry[1]-1) - max(qry_overlaps[0])
        ref_stem2 = ref_align_seq[
            (loop_index_hit_qry[2] + st_index):
                loop_index_hit_qry[3]].replace('-', '')
        qry_stem2 = qry_align_seq[
            (loop_index_hit_qry[2] + st_index):
                loop_index_hit_qry[3]].replace('-', '')
        rev_qry_stem1 = stem1_seq_qry[::-1]
        rev_ref_stem1 = stem1_seq_ref[::-1]
        rev_qry_stem1_no_gaps = rev_qry_stem1.replace('-', '')
        rev_ref_stem1_no_gaps = rev_ref_stem1.replace('-', '')

        for a, b in zip(rev_qry_stem1_no_gaps,
                        rev_ref_stem1_no_gaps):
            if a != '-':
                counter_a += 1
            if b != '-':
                counter_b += 1

            if a != b and a != '-' and b != '-':

                try:
                    nucls = rev_ref_stem1_no_gaps[counter_b] +\
                            ref_stem2[counter_b] + ',' +\
                            rev_qry_stem1_no_gaps[counter_a] +\
                            qry_stem2[counter_b]
                except:
                    continue

                if nucls in accepted_mutations:
                    diff_set += 1
                else:
                    continue

                if nucls in mutation_types and '-' not in nucls:
                    mutation_types[nucls] += 1
                elif '-' not in nucls:
                    mutation_types[nucls] = 1

        return (diff_set, first_nucleotide,
                adj_loop_nucl_mutated, 0, mutation_types)

    elif len(qry_overlaps[0]) > 0 and len(ref_overlaps[0]) > 0 and\
            stem1_seq_qry != stem1_seq_ref and\
            len(stem2_loop_indexes) > 0 and\
            len(stem1_loop_indexes) > 0 and\
            len(corresponding_res) != len(stem1_loop_indexes):

        (stem1_seq_qry,
         stem1_seq_ref) = unusual_iupac_char(stem1_seq_qry,
                                             stem1_seq_ref)

        diff_set = 0
        for a, b in zip(stem1_seq_qry,
                        stem1_seq_ref):
            if a != b:
                diff_set += 1

        return (0, first_nucleotide, adj_loop_nucl_mutated,
                diff_set, mutation_types)

    elif len(qry_overlaps[2]) > 0 and len(ref_overlaps[2]) > 0 and\
            stem2_seq_qry != stem2_seq_ref and\
            stem1_loop_indexes == stem2_loop_indexes and\
            len(corresponding_res) == len(stem2_loop_indexes) and\
            ('A mismatch in stem2 fixed' in struct_output or
             'No change in structure' in struct_output):

        (stem2_seq_qry,
         stem2_seq_ref) = unusual_iupac_char(stem2_seq_qry,
                                             stem2_seq_ref)

        diff_set = 0

        counter_a = -1
        counter_b = -1
        st_index = 0

        if min(qry_overlaps[2]) > loop_index_hit_qry[2]:
            st_index = min(qry_overlaps[2]) - loop_index_hit_qry[2]
        ref_stem1 = (ref_align_seq[loop_index_hit_qry[0]: (
            loop_index_hit_qry[1]-st_index)][::-1]).replace('-', '')
        qry_stem1 = (qry_align_seq[loop_index_hit_qry[0]: (
            loop_index_hit_qry[1]-st_index)][::-1]).replace('-', '')

        ref_stem2 = stem2_seq_ref.replace('-', '')
        qry_stem2 = stem2_seq_qry.replace('-', '')

        for a, b in zip(
                qry_stem2,
                ref_stem2):

            if a != '-':
                counter_a += 1
            if b != '-':
                counter_b += 1
            if a != b and a != '-' and b != '-':

                try:
                    nucls = ref_stem1[counter_b] +\
                            ref_stem2[counter_b] + ',' +\
                            qry_stem1[counter_a] +\
                            qry_stem2[counter_a]
                except:
                    continue
                if nucls in accepted_mutations:
                    diff_set += 1
                else:
                    continue

                if nucls in mutation_types and '-' not in nucls:
                    mutation_types[nucls] += 1
                elif '-' not in nucls:
                    mutation_types[nucls] = 1

        return (diff_set, first_nucleotide,
                adj_loop_nucl_mutated, 0, mutation_types)

    elif len(qry_overlaps[2]) > 0 and len(ref_overlaps[2]) > 0 and\
        stem2_seq_qry != stem2_seq_ref and\
        len(stem1_loop_indexes) > 0 and\
        len(stem2_loop_indexes) > 0 and\
            len(corresponding_res) != len(stem2_loop_indexes):

        (stem2_seq_qry,
         stem2_seq_ref) = unusual_iupac_char(stem2_seq_qry,
                                             stem2_seq_ref)

        diff_set = 0
        for a, b in zip(stem2_seq_qry,
                        stem2_seq_ref):
            if a != b:
                diff_set += 1

        return (0, first_nucleotide, adj_loop_nucl_mutated,
                diff_set, mutation_types)

    return 0, first_nucleotide, adj_loop_nucl_mutated, 0, mutation_types


def new_loop(ref_ts_start,
             ref_ts_end,
             loop_index_hit_qry,
             loop_index_hit_ref,
             ref_overlaps,
             qry_overlaps,
             qry_align_seq,
             ts_distance,
             qry_ts_start,
             qry_ts_end):
    """
    This function identifies if a template switch has caused a new loop.

    Arguments
    ----------

    :ref_ts_start: the starting index of the region of interest
    :ref_ts_end: the ending index of the region of interest
    :loop_index_hit_qry:  A list of loop indexes
                          [stem_start, loop_start, loop_end, stem_end]
    :loop_index_hit_ref: A list of loop indexes
    :ref_overlaps: A list containing the three sets of indexes
                        indicating (in reference sequence)
                        1. TSM overlaps with stem1
                        2. TSM overlaps with a loop and
                        3. TSM overlaps with stem2
    :qry_overlaps: A similar list as above but for the
                   child sequence
    :qry_align_seq: aligned query sequence
    :ts_distance: Distance between the template switch source and target
                  regions.
    :qry_ts_start: the first index of a TSM target site
    :qry_ts_end: the last index of a TSM target site

    Output
    ------
    True if a tsm has caused a new loop.
    """

    if loop_index_hit_ref is not None:
        loop_overlapping = set(range(
            loop_index_hit_qry[1], loop_index_hit_qry[2])).intersection(
                set(range(loop_index_hit_ref[1], loop_index_hit_ref[2])))
        if len(loop_overlapping) > 0:

            return False

    if ts_distance > 10:

        return False

    if len(qry_overlaps[2]) >= 4 and\
            qry_ts_start in set(range(
                loop_index_hit_qry[1], (loop_index_hit_qry[3]+1))) and\
            ref_ts_end in set(range(
                loop_index_hit_qry[0], (loop_index_hit_qry[2]+1))) and\
            qry_ts_end >= loop_index_hit_qry[2] and\
            ref_ts_start <= (loop_index_hit_qry[1]):

        return True

    elif len(qry_overlaps[0]) >= 4 and\
            ref_ts_start in set(range(
                loop_index_hit_qry[1], (loop_index_hit_qry[3]+1))) and\
            qry_ts_end in set(range(
                loop_index_hit_qry[0], (loop_index_hit_qry[2]+1))) and\
            qry_ts_start <= loop_index_hit_qry[1] and\
            (ref_ts_end) >= loop_index_hit_qry[2]:

        return True

    else:
        return False


def loop_seq_to_stem(ref_overlaps,
                     qry_overlaps,
                     qry_align_seq,
                     ref_align_seq,
                     ref_ts_start,
                     ref_ts_end,
                     loop_indexes_hit_ref,
                     loop_indexes_hit_qry):
    """
    Identifies TSM-cases between a loop sequence source and
    stem sequence target. The fucntion requires that position of
    the loop has to be the same between the parent and
    child sequences.

    Arguments:
    ----------
    ref_overlaps: TSM-target overlap in a reference sequence.
                  Ref_overlaps is a list with three sets of sequence
                  indexes 1. a TS overlap in stem1, 2. a TS-overlap
                  in loop. a TS-overlap in stem2.
    qry_overlaps: Similar to ref_overlaps, but for a query structure.
    qry_align_seq: An aligned query sequence
    ref_align: An aligned reference sequence
    ref_ts_start: A TSM sequence source region
    ref_ts_end: A TSM sequence target region
    :loop_index_hit_qry: A list of RNA hairpin indeces in
                         a child seq
                        [stem1 start index, loop start index,
                         loop end index, stem2 end index)
    :loop_index_hit_ref: A list of RNA hairpin indeces in
                        a parental seq

    Returns
    -------

    Number of mismatches that a loop to sequence TSM has caused.

    """

    if loop_indexes_hit_ref is None or loop_indexes_hit_qry is None:
        return 0
    if loop_indexes_hit_ref[1] == loop_indexes_hit_qry[1] and\
            loop_indexes_hit_ref[2] == loop_indexes_hit_qry[2]:

        ts_orig = set(range(ref_ts_start, ref_ts_end))

        if len(ts_orig.intersection(set(range(
                loop_indexes_hit_ref[1], loop_indexes_hit_ref[2])))) >= 3:
            if len(qry_overlaps[2]) >= 3:
                minimum = min(qry_overlaps[2])
                maximum = max(qry_overlaps[2])+1
                mismatches, length = mismatches_between_sequences(
                    qry_align_seq[minimum:maximum],
                    ref_align_seq[minimum:maximum])
                return mismatches
            elif len(qry_overlaps[0]) >= 3:
                minimum = min(qry_overlaps[0])
                maximum = max(qry_overlaps[0])+1
                mismatches, length = mismatches_between_sequences(
                    qry_align_seq[minimum:maximum],
                    ref_align_seq[minimum:maximum])
                return mismatches
    return 0


def inverted_stem(ref_overlaps,
                  qry_overlaps,
                  qry_align_seq,
                  ref_align_seq,
                  ref_ts_start,
                  ref_ts_end,
                  loop_indexes_hit_ref,
                  loop_indexes_hit_qry):
    """
    Identifies if the TSM source and target sequences are within the
    same stem.

    Arguments
    =========

    ref_overlaps: TSM-target overlap in a reference sequence.
                  Ref_overlaps is a list with three sets of sequence
                  indexes 1. a TS overlap in stem1, 2. a TS-overlap
                  in loop. a TS-overlap in stem2.
    qry_overlaps: Similar to ref_overlaps, but for a query structure.
    qry_align_seq: An aligned query sequence
    ref_align: An aligned reference sequence
    ref_ts_start: A TSM sequence source region
    ref_ts_end: A TSM sequence target region
    :loop_index_hit_ref: A list of RNA hairpin indeces in
                         a child seq
                        [stem1 start index, loop start index,
                         loop end index, stem2 end index)
    :loop_index_hit_qry: A list of RNA hairpin indeces in
                         a parental seq

    Returns
    =======

    Number of mismatches between query and reference sequences that a
    TSM has caused.

    """

    if loop_indexes_hit_ref is None or loop_indexes_hit_qry is None:
        return 0

    if loop_indexes_hit_ref[1] == loop_indexes_hit_qry[1] and\
            loop_indexes_hit_ref[2] == loop_indexes_hit_qry[2]:

        ts_orig = set(range(ref_ts_start, ref_ts_end))

        if len(qry_overlaps[0]) >= 6:
            if len(ts_orig.intersection(set(range(
                   loop_indexes_hit_ref[0],
                   loop_indexes_hit_ref[1])))) == len(ts_orig):
                minimum = min(qry_overlaps[0])
                maximum = max(qry_overlaps[0])+1

                if len(ts_orig.intersection(set(range(
                    loop_indexes_hit_qry[0],
                    loop_indexes_hit_qry[1])))) == len(ts_orig) and\
                        len(qry_overlaps[1]) == 0 and\
                        len(qry_overlaps[2]) == 0:
                    mismatches, length = mismatches_between_sequences(
                            qry_align_seq[minimum:maximum],
                            ref_align_seq[minimum:maximum])
                    return mismatches

        if len(qry_overlaps[2]) >= 6:
            if len(ts_orig.intersection(set(range(
                loop_indexes_hit_ref[2],
                    loop_indexes_hit_ref[3])))) == len(ts_orig):
                if len(ts_orig.intersection(set(range(
                       loop_indexes_hit_qry[2],
                       loop_indexes_hit_qry[3])))) == len(ts_orig) and\
                               len(qry_overlaps[0]) == 0 and\
                               len(qry_overlaps[1]) == 0:
                    minimum = min(qry_overlaps[2])
                    maximum = max(qry_overlaps[2])

                    mismatches, length = mismatches_between_sequences(
                        qry_align_seq[minimum:maximum],
                        ref_align_seq[minimum:maximum])

                    return mismatches
    return 0


def inverted_loop(ref_overlaps,
                  qry_overlaps,
                  ref_ts_start,
                  ref_ts_end,
                  ref_stem1_indexes,
                  ref_loop_indexes,
                  ref_stem2_indexes,
                  qry_stem1_indexes,
                  qry_loop_indexes,
                  qry_stem2_indexes,
                  qry_seq_aligned,
                  ref_seq_aligned):
    """
    identifies if template switch contains an inverted loop sequence.

    Arguments
    =========

    :ref_overlaps: A list of three sets covering ts_region overlapping
                   stem1, loop and stem2.
    :qry_overlaps: A list of three sets covering ts_region overlapping
                   stem1, loop and stem2.
    :ref_ts_start: A TSM sequence source region
    :ref_ts_end: A TSM sequence target region
    :ref_stem1_indexes: Alignment indexes belonging to stem1 in
                        the parental sequence
    :ref_loop_indexes:  Alignement loop indexes in the parental sequence
    :ref_stem2_indexes: Alignment indexes belonging to stem2 in
                        the parental sequence
    :qry_stem1_indexes: Alignment indexes belonging to stem1 in
                        the child sequence
    :qry_loop_indexes: Alignment indexes belonging to the loop in
                        the child sequence
    :qry_stem2_indexes: Alignment indexes belonging to stem2 in
                        the child sequence
    :qry_seq_aligned: An aligned query sequence
    :ref_seq_aligned: An aligned reference sequence


    Returns
    ---------

    inverted_loop_status: the length of an inverted loop
    palindromic: >0 if a loop inversion is palindromic
    symmetric: True, if a loop inversion is symmetric


    """

    ts_or_min = ref_ts_start
    ts_or_max = ref_ts_end

    ts_orig = set(range(ts_or_min, ts_or_max))
    symmetric = False

    if len(qry_overlaps[1]) > 0 and len(ref_overlaps[1]) > 0:
        qry_ts_min = min(qry_loop_indexes)
        qry_ts_max = max(qry_loop_indexes)+1
        ref_ts_min = min(ref_loop_indexes)
        ref_ts_max = max(ref_loop_indexes)+1
        loop_ref = ref_seq_aligned[ref_ts_min:ref_ts_max]
        loop_qry = qry_seq_aligned[qry_ts_min:qry_ts_max]
    else:
        return 0, 0, symmetric

    if len(ref_stem1_indexes.intersection(qry_stem1_indexes)) > 0 and\
            len(ref_stem2_indexes.intersection(qry_stem2_indexes)) > 0 and\
            len(loop_ref.replace('-', '')) ==\
            len(loop_qry.replace('-', '')) and\
            len(qry_overlaps[1]) != 0:

        qry_ts_min = min(qry_overlaps[1])
        qry_ts_max = max(qry_overlaps[1])+1

        loop_ref = ref_seq_aligned[qry_ts_min:qry_ts_max]
        loop_qry = qry_seq_aligned[qry_ts_min:qry_ts_max]
        loop_ref, loop_qry = unusual_iupac_char(loop_ref, loop_qry)

        ts_orig = set(range(ts_or_min, ts_or_max))

        if len(qry_overlaps[0]) > 0 and\
                len(qry_overlaps[0]) >= len(qry_overlaps[2]) and\
                len(ts_orig.intersection(ref_stem2_indexes)) > 0:
            overlap1_qry_min = min(qry_overlaps[0])
            overlap1_qry_max = max(qry_overlaps[0])+1
            overlap2_ref_min = min(ts_orig.intersection(ref_stem2_indexes))
            overlap2_ref_max = max(ts_orig.intersection(ref_stem2_indexes))+1
            ref_stem2_seq = ref_seq_aligned[
                    overlap2_ref_min:overlap2_ref_max].replace('-', '')
            qry_stem1_seq = qry_seq_aligned[
                    overlap1_qry_min:overlap1_qry_max].replace('-', '')
            if len(ref_stem2_seq) == len(qry_stem1_seq):
                symmetric = True

            elif len(qry_overlaps[2]) > 0 and\
                    len(qry_overlaps[2]) >= len(qry_overlaps[0]) and\
                    len(ts_orig.intersection(ref_stem1_indexes)) > 0:
                overlap2_qry_min = min(qry_overlaps[2])
                overlap2_qry_max = max(qry_overlaps[2])+1
                overlap1_ref_min = min(ts_orig.intersection(
                    ref_stem1_indexes))
                overlap1_ref_max = max(ts_orig.intersection(
                    ref_stem1_indexes))+1
                ref_stem1_seq = ref_seq_aligned[
                        overlap1_ref_min:overlap1_ref_max].replace('-', '')
                qry_stem2_seq = qry_seq_aligned[
                        overlap2_qry_min:overlap2_qry_max].replace('-', '')
                if len(ref_stem1_seq) == len(qry_stem2_seq):
                    symmetric = True

            ts_ref_overlap_stem2 = ref_stem2_indexes.intersection(set(ts_orig))
            ts_ref_overlap_loop = ref_loop_indexes.intersection(set(ts_orig))
            ts_ref_overlap_stem1 = ref_stem1_indexes.intersection(set(ts_orig))

            if len(qry_overlaps[1]) > 0 and len(qry_overlaps[0]) > 0 and\
                    len(ts_ref_overlap_loop) > 0 and\
                    len(ts_ref_overlap_stem2) > 0:

                ind_set = set()
                for ind in qry_overlaps[1]:
                    if qry_seq_aligned[ind] != '-':
                        ind_set.add(ind)

                    if loop_ref == loop_qry:

                        return len(ind_set), len(loop_qry.replace(
                            '-', '')), symmetric
                    else:
                        return len(ind_set), 0, symmetric

                if len(qry_overlaps[1]) > 0 and len(qry_overlaps[2]) > 0 and\
                        len(ts_ref_overlap_loop) > 0 and\
                        len(ts_ref_overlap_stem1) > 0:

                    ind_set = set()
                    for ind in qry_overlaps[1]:
                        if qry_seq_aligned[ind] != '-':
                            ind_set.add(ind)

                    if loop_ref == loop_qry:
                        return (len(ind_set),
                                len(loop_qry.replace('-', '')),
                                symmetric)

                    else:
                        return len(ind_set), 0, symmetric

    return 0, 0, symmetric
