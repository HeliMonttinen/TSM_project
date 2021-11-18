"""
A tool to analyse observed TSMs and to save
them into a database.

Author: Heli Mönttinen (ORCID: 0000-0003-2461-0690)

"""

import copy
import os
import sys
from Bio import SeqIO
from collections import (defaultdict,
                         OrderedDict)
import multiprocessing as mp
from mongoengine import (connect,
                         disconnect)


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def TSM_analysis(i, text, filename, align_dict, tree,
                 structure_dir, template_switch_clusters,
                 tupl):
    """
    A function that analyses TSMs and links them to
    the RNA secondary structure.
    This function is run via a python multiprocessing
    async function.

    Requires:
    :i: An index required by the multiprocessing async
        parallizing function.
    :text: Information on taxonomy and RNA family
    :filename: The name of the ancestral tree
    :align_dict: Alignment dictionary
    :tree: The ancestral tree read py the ete3 tree tool.
    :structure_dir; A fullpath to the directory containing
     .db structure files
    :template_switch_clusters: All TSM events within a tree.
    Outpit of the map_template_switch_in_same_region function as
     a list format.
    :tupl: (Position, dictionary containing information on high quality
            regions

    output:

    A saved document about a TSM in a given database

    """

    from tree_parse import mismatches_between_sequences

    from RNA_alignments import unusual_iupac_char

    from structure import (identify_alignment_area,
                           identify_if_stem_is_extended,
                           identify_if_structure_has_changed,
                           region_overlaps_with_loop,
                           identify_compensating_mutations,
                           new_loop,
                           inverted_loop,
                           loop_seq_to_stem,
                           inverted_stem)

    results_dict = {}
    position = tupl[0]
    regions_dict = tupl[1]

    results_dict["identifier"] = filename.split('/')[-1].rstrip('_fpa')\
        + '_' + str(position)
    results_dict["cluster"] = filename.split('/')[-1].rstrip('_fpa')
    results_dict["position"] = position
    results_dict["ts_distance"] = 1000
    results_dict["ts_quality"] = False
    results_dict["ts_length"] = 0
    results_dict["ts_mismatch"] = 0
    results_dict["structural_quality"] = False
    results_dict["CM_present"] = False
    results_dict["one_CM"] = False
    results_dict["multiple_CM"] = False
    results_dict["CM_and_loop_change"] = False
    results_dict["CM_and_inverted_loop"] = False
    results_dict["inverted_loop"] = False
    results_dict["new_loop"] = False
    results_dict["asymmetric_TSM_between_stems"] = False
    results_dict["TSM_loop_to_stem"] = False
    results_dict["TSM_within_one_stem"] = False
    results_dict["loop_to_seq"] = 0
    results_dict["adj_nucl1"] = ""
    results_dict["adj_nucl2"] = ""
    results_dict["adj_nucl1_ref"] = ""
    results_dict["adj_nucl2_ref"] = ""
    results_dict["tsm_within_stem"] = 0
    results_dict["loop_seq"] = ""
    results_dict["loop_seq_ref"] = ""
    results_dict["major_change"] = False
    results_dict["no_change"] = False
    results_dict["loop_extension"] = False
    results_dict["mismatch_fixed"] = False
    results_dict["leaf"] = False
    results_dict["new_mismatch"] = False
    results_dict["new_loop"] = False
    results_dict["inverted_loop"] = False
    results_dict["compensating_mutations"] = 0
    results_dict["alignment_file"] = filename.rstrip('_fpa') + '_pagan.fas'
    results_dict["tree_file"] = filename.rstrip('_fpa') + '_pagan.anctree'
    results_dict["family"] = ""
    results_dict["taxonomy"] = ""
    results_dict["loop_extension"] = False
    results_dict["stem_base_extension"] = False
    results_dict["stem_extended"] = False

    for line in text.split('\n'):

        if 'Family' in line:

            family = line.split('Family: ')[1].rstrip()
            results_dict["family"] = family

        elif 'Taxonomy' in line:
            taxonomy = line.split('Taxonomy: ')[1].rstrip()
            results_dict["taxonomy"] = taxonomy

    nodes = {}
    ref_node = template_switch_clusters[position][2][0]

    nodes[ref_node] = ""
    query_node = template_switch_clusters[position][2][1]

    nodes[query_node] = ""
    if not (tree&ref_node).is_root():
        ancestor_node = (tree&template_switch_clusters[position][2][0]).up.name
    else:
        print(filename, position, "reference is a root.")
        return i, None, None, None, None, None, None, None

    # Parsing the TSM cluster input: indexes of TSM target
    # and source in the alignment
    qry_ts_start = int(template_switch_clusters[position][2][4])
    qry_ts_end = int(template_switch_clusters[position][2][5])
    ref_ts_start = int(template_switch_clusters[position][2][7])
    ref_ts_end = int(template_switch_clusters[position][2][6])

    index_info = [qry_ts_start, qry_ts_end, ref_ts_start, ref_ts_end]

    align_ts_qry = str(align_dict[query_node])[
            int(qry_ts_start):int(qry_ts_end)]
    align_ts_ref = str(align_dict[ref_node])[
            int(qry_ts_start):int(qry_ts_end)]
    align_orig_ref = str(align_dict[ref_node])[
            int(ref_ts_start):int(ref_ts_end)]
    align_orig_qry = str(align_dict[query_node])[
            int(ref_ts_start):int(ref_ts_end)]
    qry_child_structure_node1 = None
    qry_child_structure_node2 = None

    (align_ts_qry_c,
     align_ts_ref_c) = unusual_iupac_char(align_ts_qry,
                                          align_ts_ref)

    (align_orig_qry_c,
     align_orig_ref_c) = unusual_iupac_char(align_orig_qry,
                                            align_orig_ref)

    if ref_ts_start < qry_ts_start:

        if ref_ts_end < qry_ts_start:

            if align_orig_qry_c != align_orig_ref_c:

                return i, None, None, None, None, None, None, None
        else:
            diff = qry_ts_start - ref_ts_start

            if align_orig_qry_c[0:diff] != align_orig_ref_c[0:diff]:
                return i, None, None, None, None, None, None, None

    elif qry_ts_start < ref_ts_start:
        if qry_ts_end < ref_ts_start:

            if align_orig_qry_c != align_orig_ref_c:
                return i, None, None, None, None, None, None, None
        else:
            diff = ref_ts_end - qry_ts_end
            if align_orig_qry_c[diff:] != align_orig_ref_c[diff:]:
                return i, None, None, None, None, None, None, None

    mismatches, length = mismatches_between_sequences(align_ts_ref_c,
                                                      align_ts_qry_c,
                                                      strict_mode=True)

    length = len(align_ts_qry.replace('-', ''))

    if length < 6 or mismatches < 1:
        return i, None, None, None, None, None, None, None

    if qry_ts_end < ref_ts_start:

        dist_seq = str(align_dict[query_node])[qry_ts_end:ref_ts_start]
        dist_seq_w_gaps = dist_seq.replace('-', '')
        ts_distance = len(dist_seq_w_gaps)

    elif ref_ts_end < qry_ts_start:

        dist_seq = str(align_dict[query_node])[ref_ts_end:qry_ts_start]
        dist_seq_w_gaps = dist_seq.replace('-', '')
        ts_distance = len(dist_seq_w_gaps)

    else:
        ts_distance = 0

    results_dict["ts_distance"] = ts_distance

    # If the query is not a leaf node, collect the information on
    # child nodes and the location in the alignment
    qry_child_structure_node1 = None
    qry_child_structure_node2 = None
    results_dict["structureFileChild1"] = None
    results_dict["structureFileChild2"] = None

    if not (tree&query_node).is_leaf():
        qry_child_structure_node1 = (
                tree&template_switch_clusters[position][2][1]).children[0].name
        qry_child_structure_node2 = (
                tree&template_switch_clusters[position][2][1]).children[1].name
        nodes[qry_child_structure_node1] = ""
        nodes[qry_child_structure_node2] = ""
        children_structs = []
        child_struct1 = structure_dir +\
            filename.split('/')[-1].rstrip('_fpa')\
            + '_' + qry_child_structure_node1 + '.db'
        child_struct2 = structure_dir +\
            filename.split('/')[-1].rstrip('_fpa')\
            + '_' + qry_child_structure_node2 + '.db'
        results_dict["structureFileChild1"] = child_struct1
        results_dict["structureFileChild2"] = child_struct2
        children_structs.append(qry_child_structure_node1)
        children_structs.append(qry_child_structure_node2)

    regions = regions_dict

    flag_ref = False
    flag_qry = False
    ch_flag = False
    for reg in regions:
        reg1 = int(reg.split('_')[0])
        reg2 = int(reg.split('_')[1])

        if len(set(range(qry_ts_start, qry_ts_end)).intersection(
            set(range(reg1, reg2)))) ==\
                len(range(qry_ts_start, qry_ts_end)):
            flag_qry = True

        if len(set(range(ref_ts_start, ref_ts_end)).intersection(
            set(range(reg1, reg2)))) ==\
                len(range(ref_ts_start, ref_ts_end)):
            flag_ref = True

        if not (tree&query_node).is_leaf() and flag_qry is True:
            for ch in children_structs:
                (align_ch_orig,
                 align_qrych_orig) =\
                    unusual_iupac_char(
                        align_dict[ch][
                            int(ref_ts_start):int(ref_ts_end)],
                        align_dict[query_node][
                            int(ref_ts_start):int(ref_ts_end)])

                (align_ch_ts,
                 align_qrych_ts) =\
                    unusual_iupac_char(
                        align_dict[ch][
                            int(qry_ts_start):int(qry_ts_end)],
                        align_dict[query_node][
                            int(qry_ts_start):int(qry_ts_end)])

                if align_ch_orig == align_qrych_orig and\
                        align_ch_ts == align_qrych_ts:
                    ch_flag = True

    quality_info = False
    if (tree&query_node).is_leaf() and flag_qry is True and\
            flag_ref is True:
        results_dict["ts_quality"] = True
        results_dict["structural_quality"] = True
        quality_info = True

    elif ch_flag is True and flag_ref is True:
        results_dict["ts_quality"] = ch_flag
        results_dict["structural_quality"] = ch_flag
        quality_info = True

    # Collect the information to the sequence quality check

    ref_tree = tree&template_switch_clusters[position][2][0]

    child_node1 = ref_tree.children[0].name
    child_node2 = ref_tree.children[1].name

    if child_node1 != template_switch_clusters[position][2][1]:

        sister_node = child_node1
    else:
        sister_node = child_node2

    nodes[sister_node] = ""

    if (tree&query_node).is_leaf():
        results_dict["leaf"] = True
    else:
        results_dict["leaf"] = False

    results_dict["sister_node"] = sister_node
    results_dict["child_node1"] = qry_child_structure_node1
    results_dict["child_node2"] = qry_child_structure_node2
    results_dict["ancestor_node"] = ancestor_node

    # Check that structural information exists for all nodes

    structure_file_ref = structure_dir +\
        filename.split('/')[-1].rstrip('_fpa') + '_' + ref_node + '.db'
    if not os.path.isfile(structure_dir +
                          filename.split('/')[-1].rstrip('_fpa') +
                          '_' + ref_node + '.db'):
        return i, None, None, None, None, None, None, None
    structure_file_qry = structure_dir +\
        filename.split('/')[-1].rstrip('_fpa') + '_' + query_node + '.db'
    if not os.path.isfile(structure_dir +
                          filename.split('/')[-1].rstrip('_fpa') +
                          '_' + query_node + '.db'):
        return i, None, None, None, None, None, None, None
    structure_file_sister = structure_dir +\
        filename.split('/')[-1].rstrip('_fpa') + '_' + sister_node + '.db'
    if not os.path.isfile(structure_dir +
                          filename.split('/')[-1].rstrip('_fpa') +
                          '_' + sister_node + '.db'):
        return i, None, None, None, None, None, None, None

    results_dict["structureRef"] = structure_file_ref
    results_dict["structureQry"] = structure_file_qry
    results_dict["structureSister"] = structure_file_sister

    # Structural quality check for reference sequence and
    # for the query sequence if query is not a leaf.

    minimum = min(set([qry_ts_start, qry_ts_end,
                       ref_ts_end, ref_ts_start]))
    maximum = max(set([qry_ts_start, qry_ts_end,
                       ref_ts_end, ref_ts_start]))
    # identify aligned structure
    ts_region_qry, aligned_struct_qry = identify_alignment_area(
            structure_file_qry,
            align_dict,
            minimum,
            maximum)

    ts_region_ref, aligned_struct_ref = identify_alignment_area(
            structure_file_ref,
            align_dict,
            minimum,
            maximum)

    # Identify if the TSM region overlaps with a hairpin loop.
    # Identification is done separately for both the query and reference.

    template_switch = set(range(template_switch_clusters[position][0],
                          template_switch_clusters[position][1]))

    (align_qry_seq,
     align_ref_seq) = unusual_iupac_char(
             str(align_dict[query_node]),
             str(align_dict[ref_node]))

    adj_nucl1 = "None"
    adj_nucl2 = "None"
    adj_nucl1_ref = "None"
    adj_nucl2_ref = "None"

    (loop_struct,
     temp_struct,
     loop_seq_qry,
     loop_indexes,
     loop_overlap,
     stem1_overlap,
     stem2_overlap,
     qry_stem1,
     qry_stem2,
     qry_loop,
     adj_nucl1,
     adj_nucl2) = region_overlaps_with_loop(aligned_struct_qry,
                                            align_qry_seq,
                                            template_switch)
    (loop_struct_ref,
     temp_struct_ref,
     loop_seq_ref,
     loop_indexes_ref,
     loop_overlap_ref,
     stem1_overlap_ref,
     stem2_overlap_ref,
     ref_stem1,
     ref_stem2,
     ref_loop,
     adj_nucl1_ref,
     adj_nucl2_ref) = region_overlaps_with_loop(aligned_struct_ref,
                                                align_ref_seq,
                                                template_switch)

    if loop_indexes is not None:
        loop_seq = loop_seq_qry.replace('-', '')
        if len(loop_seq) == 0:
            loop_seq = 'None'
        if len(adj_nucl1) == 0:
            adj_nucl1 = 'None'
        if len(adj_nucl2) == 0:
            adj_nucl2 = 'None'

        results_dict["loop_seq"] = loop_seq
        results_dict["adj_nucl1"] = adj_nucl1
        results_dict["adj_nucl2"] = adj_nucl2

    if loop_indexes_ref is not None:
        loop_seq_ref = loop_seq_ref.replace('-', '')
        if len(loop_seq_ref) == 0:
            loop_seq_ref = 'None'
        if len(adj_nucl1_ref) == 0:
            adj_nucl1_ref = 'None'
        if len(adj_nucl2_ref) == 0:
            adj_nucl2_ref = 'None'

        results_dict["loop_seq_ref"] = loop_seq_ref
        results_dict["adj_nucl1_ref"] = adj_nucl1_ref
        results_dict["adj_nucl2_ref"] = adj_nucl2_ref

    ts_ref_struct = aligned_struct_ref[qry_ts_start:qry_ts_end]
    ts_qry_struct = aligned_struct_qry[qry_ts_start:qry_ts_end]

    inverted_loop_status = 0
    palindromic = 0
    symmetric = False

    # Identify if a loop sequence is inverted in RNA hairpin
    # overlapping with a identified TSM

    (inverted_loop_status,
     palindromic,
     symmetric) = inverted_loop(
            [stem1_overlap_ref, loop_overlap_ref, stem2_overlap_ref],
            [stem1_overlap, loop_overlap, stem2_overlap],
            ref_ts_start,
            ref_ts_end,
            ref_stem1,
            ref_loop,
            ref_stem2,
            qry_stem1,
            qry_loop,
            qry_stem2,
            align_qry_seq,
            align_ref_seq)

    # Identify if TSM has swapped a loop sequence or part of it
    # to the stem.

    loop_seq_mismatch = loop_seq_to_stem(
            [stem1_overlap_ref, loop_overlap_ref, stem2_overlap_ref],
            [stem1_overlap, loop_overlap, stem2_overlap],
            align_qry_seq,
            align_ref_seq,
            ref_ts_start,
            ref_ts_end,
            loop_indexes_ref,
            loop_indexes)

    # identify a TMS has inverted a part of a stem

    within_stem_tsms = inverted_stem(
            [stem1_overlap_ref, loop_overlap_ref, stem2_overlap_ref],
            [stem1_overlap, loop_overlap, stem2_overlap],
            align_qry_seq,
            align_ref_seq,
            ref_ts_start,
            ref_ts_end,
            loop_indexes_ref,
            loop_indexes)

    if inverted_loop_status == len(loop_seq_qry.replace('-', '')) and\
        inverted_loop_status > 0 and palindromic == 0 and\
            len(loop_seq_ref.replace('-', '')) ==\
            len(loop_seq_qry.replace('-', '')):
        full_loop_inversion = True
    else:
        full_loop_inversion = False

    # Identify if the structure has changed in the TSM-linking area.
    if ts_ref_struct == ts_qry_struct:
        structural_change = 'No change in structure'
        results_dict["no_change"] = True

    elif loop_indexes is not None and temp_struct_ref is not None:

        output = identify_if_structure_has_changed(
            loop_indexes,
            loop_indexes_ref,
            [stem1_overlap_ref, loop_overlap_ref, stem2_overlap_ref],
            [stem1_overlap, loop_overlap, stem2_overlap],
            int(template_switch_clusters[position][0]),
            temp_struct_ref,
            temp_struct,
            ref_ts_start,
            ref_ts_end,
            qry_ts_start,
            qry_ts_end)

        for out in output:

            if 'new' in out:
                results_dict["new_mismatch"] = True
            elif 'fixed' in out:
                results_dict['mismatch_fixed'] = True
            if 'stem is' in out:
                results_dict["stem_extended"] = True
            if 'major' in out and len(output) == 1:
                results_dict["major_change"] = True

        structural_change = ', '.join(output)
        if structural_change == 0:
            structural_change = 'A major change'

    elif temp_struct != temp_struct_ref:
        structural_change = 'A major change'
        results_dict["major_change"] = True

    else:
        structural_change = 'major change'
        results_dict["major_change"] = True

    extension_output = []
    if loop_indexes is not None and loop_indexes_ref is not None:

        extension_output = identify_if_stem_is_extended(
                loop_indexes,
                loop_indexes_ref,
                [stem1_overlap_ref, loop_overlap_ref, stem2_overlap_ref],
                [stem1_overlap, loop_overlap, stem2_overlap],
                ref_ts_start,
                ref_ts_end,
                qry_ts_start,
                qry_ts_end)

    for ext in extension_output:
        if "A stem is extended toward a loop" in ext:
            results_dict["loop_extension"] = True

        elif "A stem is extended at the base." in ext:
            results_dict["stem_base_extension"] = True
    # Identify compensating substitutions linking to a TSM

    compensating = 0
    unsymmetric_subs = 0
    mutation_types = {}
    if (loop_indexes is not None) and (loop_indexes_ref is not None) and\
            loop_indexes[1] == loop_indexes_ref[1] and\
            loop_indexes[2] == loop_indexes_ref[2]:

        (compensating,
         first_nucleotide,
         adj_loop_nucleotide,
         unsymmetric_subs,
         mutation_types) = identify_compensating_mutations(
                ref_ts_start,
                ref_ts_end,
                align_qry_seq,
                align_ref_seq,
                loop_indexes,
                loop_indexes_ref,
                [stem1_overlap_ref, loop_overlap_ref, stem2_overlap_ref],
                [stem1_overlap, loop_overlap, stem2_overlap],
                temp_struct_ref,
                temp_struct,
                structural_change)

        results_dict["compensating_mutations"] = compensating

    else:
        results_dict["compensating_mutations"] = compensating

    # Identify if a TSM is linked to a new loop in a query sequence
    new_loop_status = False

    if loop_indexes is not None:

        new_loop_status = new_loop(
            ref_ts_start,
            ref_ts_end,
            loop_indexes,
            loop_indexes_ref,
            [stem1_overlap_ref,
             loop_overlap_ref,
             stem2_overlap_ref],
            [stem1_overlap,
             loop_overlap,
             stem2_overlap],
            align_qry_seq,
            ts_distance,
            qry_ts_start,
            qry_ts_end)

    results_dict["query_node"] = query_node
    results_dict["ref_node"] = ref_node
    results_dict["qry_ts_start"] = qry_ts_start
    results_dict["qry_ts_end"] = qry_ts_end
    results_dict["ref_ts_start"] = ref_ts_start
    results_dict["ref_ts_end"] = ref_ts_end
    results_dict["ts_mismatch"] = mismatches
    results_dict["ts_length"] = length

    # Collect information from the TSM-linked changes in the RNA hairpin
    if compensating > 0:
        results_dict["CM_present"] = True

        if (compensating == 1 and
            ((inverted_loop_status == 0) or
             (inverted_loop_status > 0 and
                 palindromic > 0))):
            results_dict["one_CM"] = True

        if (compensating > 1 and
            ((inverted_loop_status == 0) or
                (inverted_loop_status > 0 and
                    palindromic > 0))):

            results_dict["multiple_CM"] = True

        if inverted_loop_status > 0 and\
            palindromic == 0 and\
                full_loop_inversion is False and\
                results_dict["multiple_CM"] is False:

            results_dict["CM_and_loop_change"] = True

        if full_loop_inversion is True and symmetric is True:
            results_dict["CM_and_inverted_loop"] = True

    if compensating == 0 and\
        full_loop_inversion is True and\
            palindromic == 0 and\
            symmetric is True:
        results_dict["inverted_loop"] = True

    if new_loop_status is True and\
            length >= 8 and\
            mismatches >= 2:
        results_dict["new_loop"] = True

    if loop_seq_mismatch > 0 and unsymmetric_subs == 0 and\
            mismatches >= 2 and\
            results_dict["stem_extended"] is False and\
            results_dict["CM_present"] is False:
        results_dict["TSM_loop_to_stem"] = True

    if unsymmetric_subs > 0 and mismatches >= 2 and\
            results_dict["TSM_loop_to_stem"] is False and\
            results_dict["CM_present"] is False and\
            results_dict["stem_extended"] is False:
        results_dict["asymmetric_TSM_between_stems"] = True

    if within_stem_tsms > 0 and mismatches >= 2 and\
            results_dict["asymmetric_TSM_between_stems"] is False and\
            results_dict["CM_present"] is False:
        results_dict["TSM_within_one_stem"] = True

    return (i, results_dict, quality_info, index_info,
            mutation_types, align_ts_qry_c, align_orig_ref_c)


def main():
    """
    This script identifies, if the TSMs found by FPA tool
    are linked to RNA hairpins and if they have caused:

    1) CSs
    2) partial or full loop inversions
    3) Asymmetric TSM from the loop to stem
       within the same hairpin
    4) Asymmetric TSMs within the same stem
    5) Insertion in stem half

    PLEASE NOTE! The script expects a directory that
    contains alignment file and tree and fpa results files.
    They should be named as identifier_pagan.fas, identifieri_pagan.anctree and
    identifier_fpa, respectively. The alignment file
    fasta headers have to match the nodes of the tree files.

    Arguments
    =========

    fpa_dir: A directory that contains alignment file and
             directory and fpa results file. They should be
             named as identifier.fas, identifier.anctree and
             identifier_fpa, respectively. The alignment file
             fasta headers have to match the nodes of the tree
             files.

    structure_dir: A directory that contains dot-parenthesis
                   structure files for all nodes in the tree.
                   Structure files should be named as
                   identfier_nodename.db. The node name has to
                   match alignment file fasta headers and
                   node names in the tree.

    filename: The fpa file that contains output from
              run_fpa.py

    outpufile: A file where the output is written.

    Output
    =======

    A TSM document in a given database
    """

    from tree_parse import (get_node_db,
                            read_tree,
                            make_family_tax_template_switch_title,
                            map_template_switch_in_same_region)

    from results_database import add_information

    structure_dir = sys.argv[1]
    fpa_directory = sys.argv[2]
    mutation_output = sys.argv[3]
    database = sys.argv[4]
    cpu = sys.argv[5]
    regions_file = sys.argv[6]

    connect(database, "default")

    regions_dict = defaultdict(dict)
    mutations_dict = {}

    with open(regions_file, 'r') as f:

        for line in f:
            line_splitted = line.rstrip().split('\t')
            identifier = line_splitted[0].split('_')[0]
            parent = line_splitted[1].rstrip()
            regions = line_splitted[6].split(';')

            for reg in regions:
                if len(reg) > 0:

                    reg_0 = int(reg.lstrip('[').split(', ')[0])
                    reg_1 = int(reg.rstrip(']').split(', ')[1])
                    if identifier not in regions_dict:
                        regions_dict[identifier] = dict()
                    if parent not in regions_dict[identifier]:
                        regions_dict[identifier][parent] = dict()
                    regions_dict[identifier][parent][
                            str(reg_0) + '_' + str(reg_1)] = ""

    for subdir, dirs, files in os.walk(fpa_directory):
        for file in files:
            filename = fpa_directory + os.sep + file

            if filename.endswith('_fpa'):

                identifier = file.rstrip('_fpa')
                used_indexes = set()

                tree = read_tree(filename.rstrip('_fpa') + '_pagan.anctree')

                bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                    filename.rstrip('_fpa') + '_pagan.fas',
                                    "fasta"))

                features_dict = get_node_db(
                        tree, 'Rfam_loops_v2', 'rfam_loops')

                text = make_family_tax_template_switch_title(
                        tree,
                        'Rfam_family',
                        'rfam',
                        features_dict)

                align_dict = OrderedDict()

                for seq in bio_align_dict:
                    align_dict[seq] = str(bio_align_dict[seq].seq)

                # Identification of the TSMs in the alignment

                template_switch_clusters = map_template_switch_in_same_region(
                        filename,
                        filename.rstrip('_fpa') + '_pagan.fas',
                        do_not_merge=True)

                # Going through all identified TSMs in the alignment
                position_list = []

                for position in template_switch_clusters:

                    ref_node = template_switch_clusters[
                            position][2][0].rstrip("'").lstrip("'")
                    try:

                        regs = copy.deepcopy(
                                regions_dict[identifier][ref_node])

                    except:
                        continue

                    tupl = [position, regs]
                    position_list.append(tupl)

                if len(position_list) == 0:
                    continue
                if len(position_list) < int(cpu):
                    pool = mp.Pool(len(position_list))
                else:
                    pool = mp.Pool(int(cpu))

                print("start", filename)
                results_objects = [pool.apply_async(
                    TSM_analysis,
                    args=(i, text, filename,
                          align_dict, tree,
                          structure_dir,
                          template_switch_clusters,
                          tupl)) for i, tupl in enumerate(position_list)]

                results1 = [r.get()[1] for r in results_objects]
                try:
                    results2 = [r.get()[2] for r in results_objects]
                except:
                    continue

                try:
                    results3 = [r.get()[3] for r in results_objects]

                except:
                    continue

                results4 = [r.get()[4] for r in results_objects]
                results5 = [r.get()[5] for r in results_objects]
                results6 = [r.get()[6] for r in results_objects]

                for a in range(len(results1)):

                    if results2[a] is True and\
                        (str(results3[a]) + ',' +
                         results5[a].replace('-', ''))\
                            not in used_indexes:
                        add_information(results1[a])
                        used_indexes.add(
                            str(results3[a]) + ',' +
                            results5[a].replace('-', ''))
                        used_indexes.add(
                            str([results3[a][2], results3[a][3],
                                results3[a][0], results3[a][1]]) + ','
                            + results6[a].replace('-', ''))

                        if results4[a] is not None and\
                                results1[a]["ts_distance"] < 9:
                            for mutation in results4[a]:
                                if mutation in mutations_dict:
                                    mutations_dict[mutation] +=\
                                            results4[a][mutation]
                                else:
                                    mutations_dict[mutation] =\
                                            results4[a][mutation]

                    elif results1[a] is False:
                        add_information(results1[a])

                pool.close()
                pool.join()

        with open(mutation_output, 'w') as f:
            for a in mutations_dict:

                f.write(a + '\t' + str(mutations_dict[a]) + '\n')

        disconnect("default")


if __name__ == "__main__":
    main()
