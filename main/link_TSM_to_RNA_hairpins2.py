"""
A tool to analyse observed TSMs and to save
them into a database and a .csv file.

Author: Heli Mönttinen (ORCID: 0000-0003-2461-0690)

"""

import os
import sys
from Bio import SeqIO
from collections import OrderedDict
import multiprocessing as mp
from mongoengine import (connect,
                         disconnect)


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def TSM_analysis(i, filename, align_dict, tree,
                 structure_dir, template_switch_clusters,
                 position):
    """
    A function that analyses TSMs and links them to
    the RNA secondary structure.
    This function is run via a python multiprocessing
    async function.

    Requires:
    :i: An index required by the multiprocessing async
        parallizing function.
    :filename: The name of the ancestral tree
    :align_dict: Alignment dictionary
    :tree: The ancestral tree read py the ete3 tree tool.
    :structure_dir; A fullpath to the directory containing 
     .db structure files
    :template_switch_clusters: All TSM events within a tree.
    Outpit of the
     map_template_switch_in_same_region function as
     a list format.
    :position: An index for the TSM events.

    output:
    :result_txt: One line for the .csv results file. The columns
    are separated by a tab.

    Line format (python index in the beginning of each column)
    [0]filename + [1]tsm_event_identifier + [2]ref_node +
    [3]query_node+ [4]ancestor_node +
    [5]qry_ts_start + [6]qry_ts_end + [7]ref_ts_start +
    [8]ref_ts_end + [9]mismatches+
    [10]length + [11]ts_distance + [12]sister_node +
    [13]query_child_node1 +
    [14]query_child_node2 + [15]query_is_leaf(true/false) +
    [16]sequence_quality(true/false) + [17]ancestral_tree_file
    + [18]alignment_file + [19]ref_struct_file +
    [20]query_struct_file + [21]sister_node_file +
    [22]child_node_file1 + [23]child_node_file2 +
    [24]structural_quality(true/false) + [25]qry_loopseq +
    [26]qry_closing_nucl_stem1 + [27]qry_closing_nucl_stem2 +
    [28]ref_loopseq + [29]ref_closing_nucl_stem1 +
    [30]ref_closing_nucl_stem2 + [31]structure_change +
    [32]number of compensating mutations + [33]contains_CM(true/false)+
    [34]one_CM(true_false + [35]multiple_CMs +
    [36]CM and partial loop inversion (true/false) +
    [37]CM and full loop inversion + [38]full loop inversion only +
    [39]new_loop +  [40]asymmetric CM  betweenstems +
    [41]loop_to_stem_TSM + [42]TSM_within_a_stem

    """

    from tree_parse import (read_tree,
                            map_template_switch_in_same_region,
                            reference_sequence_quality,
                            compare_predicted_sequence,
                            reference_sequence_quality,
                            mismatches_between_sequences)

    from RNA_alignments import unusual_iupac_char

    from structure import (identify_alignment_area,
                           identify_if_structure_has_changed,
                           region_overlaps_with_loop,
                           identify_compensating_mutations,
                           new_loop,
                           inverted_loop,
                           structure_quality,
                           loop_seq_to_stem,
                           inverted_stem)


    results_dict = {}

    results_dict["identifier"] = filename.split('/')[-1].rstrip('_fpa') + '_' + str(position)
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

    nodes = {}
    ref_node = template_switch_clusters[position][2][0]
    result_text = filename.rstrip('_fpa').split('/')[-1]
    result_text = result_text + '\t' + result_text +'_' + str(position)
    result_text = result_text + '\t' +ref_node

    nodes[ref_node] = ""
    query_node = template_switch_clusters[position][2][1]
    result_text = result_text + '\t' + query_node

    nodes[query_node] = ""
    if not (tree&ref_node).is_root():
        ancestor_node = (tree&template_switch_clusters[position][2][0]).up.name
        result_text = result_text + '\t' + ancestor_node
    else:
        print(filename, position, "reference is a root.")
        return i, None, None

    #Parsing the TSM cluster input: indexes of TSM target and source in the alignment
    qry_ts_start = int(template_switch_clusters[position][2][4])
    qry_ts_end = int(template_switch_clusters[position][2][5])
    ref_ts_start = int(template_switch_clusters[position][2][7])
    ref_ts_end = int(template_switch_clusters[position][2][6])

    result_text = result_text + '\t' + str(qry_ts_start) +'\t' + str(qry_ts_end) +\
            '\t' + str(ref_ts_start) + '\t' + str(ref_ts_end) 

    align_ts_qry = str(align_dict[query_node])[int(qry_ts_start):int(qry_ts_end)]
    align_ts_ref = str(align_dict[ref_node])[int(qry_ts_start):int(qry_ts_end)]
    align_orig_ref = str(align_dict[ref_node])[int(ref_ts_start):int(ref_ts_end)] 
    align_orig_qry = str(align_dict[query_node])[int(ref_ts_start):int(ref_ts_end)]
    align_orig_anc = str(align_dict[ancestor_node])[int(ref_ts_start):int(ref_ts_end)]
    align_ts_anc = str(align_dict[ancestor_node])[int(qry_ts_start):int(qry_ts_end)]
    qry_child_structure_node1 = None
    qry_child_structure_node2 = None

    (align_ts_qry_c,
     align_ts_ref_c) = unusual_iupac_char(align_ts_qry,
                                          align_ts_ref)

    mismatches, length = mismatches_between_sequences(align_ts_qry_c,
                                                      align_ts_ref_c)

    length = len(align_ts_qry.replace('-',''))

    result_text = result_text + '\t' + str(mismatches) + '\t' + str(length) 

    relaxed_status=None

    if length < 6 and mismatches < 1:
        return i, None, None

    if qry_ts_end < ref_ts_start:

        dist_seq = str(align_dict[query_node])[qry_ts_end:ref_ts_start]
        dist_seq_w_gaps = dist_seq.replace('-','')
        ts_distance = len(dist_seq_w_gaps)


    elif ref_ts_end < qry_ts_start:

        dist_seq = str(align_dict[query_node])[ref_ts_end:qry_ts_start]
        dist_seq_w_gaps = dist_seq.replace('-','')
        ts_distance = len(dist_seq_w_gaps)

    else:
        ts_distance = 0

    results_dict["ts_distance"] = ts_distance

    result_text = result_text + '\t' + str(ts_distance)


    #If the query is not a leaf node, collect the information on
    # the child nodes and the location in the alignment
    qry_child_structure_node1 = None
    qry_child_structure_node2 = None

    if not (tree&query_node).is_leaf():
        qry_child_structure_node1 = (tree&template_switch_clusters[position][2][1]).children[0].name
        qry_child_structure_node2 = (tree&template_switch_clusters[position][2][1]).children[1].name
        nodes[qry_child_structure_node1] = ""
        nodes[qry_child_structure_node2] = ""
        children_lists = []
        children_structs = []
        child_leaves = False
        if not (tree&qry_child_structure_node1).is_leaf():
            children_list = []
            child1_seq = str(align_dict[qry_child_structure_node1])[int(qry_ts_start):int(qry_ts_end)]
            children_list.append(child1_seq)
            children_lists.append(children_list)
            children_structs.append([qry_child_structure_node1])
        if not (tree&qry_child_structure_node2).is_leaf():
            children_list = []
            child2_seq = str(align_dict[qry_child_structure_node2])[int(qry_ts_start):int(qry_ts_end)]
            children_list.append(child2_seq)
            children_lists.append(children_list)
            children_structs.append([qry_child_structure_node2])
        if len(children_lists) == 0:
            child_leaves = True
            children_list = []
            child1_seq = str(align_dict[qry_child_structure_node1])[int(qry_ts_start):int(qry_ts_end)]
            child2_seq = str(align_dict[qry_child_structure_node2])[int(qry_ts_start):int(qry_ts_end)]
            children_list.append(child1_seq)
            children_list.append(child2_seq)
            children_lists.append(children_list)
            children_structs.append([qry_child_structure_node1, qry_child_structure_node2])


    #Collect the information to the sequence quality check

    ref_tree = tree&template_switch_clusters[position][2][0]

    child_node1 = ref_tree.children[0].name
    child_node2 = ref_tree.children[1].name

    if child_node1 != template_switch_clusters[position][2][1]:

        sister_node = child_node1
    else:
        sister_node = child_node2

    nodes[sister_node] = ""

    align_ts_sister = str(align_dict[sister_node])[int(qry_ts_start):int(qry_ts_end)]
    align_orig_sister = str(align_dict[sister_node])[int(ref_ts_start):int(ref_ts_end)]

    result_text = result_text + '\t' + str(sister_node) + '\t' +\
            str(qry_child_structure_node1) + '\t' + str(qry_child_structure_node2)

    if (tree&query_node).is_leaf():
        result_text = result_text + '\t' + 'true'
        results_dict["leaf"] = True
    else:
        result_text = result_text + '\t' + 'false'
        results_dict["leaf"] = False

    flag = False

    results_dict["sister_node"] = sister_node
    results_dict["child_node1"] = qry_child_structure_node1
    results_dict["child_node2"] = qry_child_structure_node2
    results_dict["ancestor_node"] = ancestor_node


    #Quality check of the reference sequence prediction
    if reference_sequence_quality(align_orig_ref,
                                  align_ts_ref,
                                  align_orig_sister,
                                  align_ts_sister,
                                  align_orig_qry,
                                  align_orig_anc,
                                  align_ts_anc,
                                  relaxed=relaxed_status) is True:
        flag = True
    
    #Quality check of the query sequence prediction, if
    # query is not a leaf
    if (tree&query_node).is_leaf() is False and flag is True:
        for children_list in children_lists:
            if compare_predicted_sequence(align_ts_qry,
                    set(children_list),
                    relaxed=relaxed_status,
                    iupac=True) is True:
                flag = True
                break
            else:
                flag = False
    if flag is True:
        result_text = result_text + '\t' + 'true'
        results_dict["ts_quality"] = True
    else:
        result_text = result_text + '\t' + 'false'
        results_dict["ts_quality"] = False

    #Check that structural information exists for all nodes
    struct_list = []
    structure_file_ref = structure_dir + filename.split('/')[-1].rstrip('_fpa') + '_' + ref_node + '.db'
    if not os.path.isfile(structure_dir + filename.split('/')[-1].rstrip('_fpa') + '_' + ref_node + '.db'):
        return i, None, None
    structure_file_qry = structure_dir + filename.split('/')[-1].rstrip('_fpa') + '_' + query_node + '.db'
    if not os.path.isfile(structure_dir + filename.split('/')[-1].rstrip('_fpa') + '_' + query_node + '.db'):
        return i, None, None
    structure_file_sister = structure_dir + filename.split('/')[-1].rstrip('_fpa') + '_' + sister_node + '.db'
    if not os.path.isfile(structure_dir + filename.split('/')[-1].rstrip('_fpa') + '_' + sister_node + '.db'):
        return i, None, None

    results_dict["structureRef"] = structure_file_ref
    results_dict["structureQry"] = structure_file_qry
    results_dict["structureSister"] = structure_file_sister
    results_dict["structureFileChild1"] = None
    results_dict["structureFileChild2"] = None

    result_text = result_text + '\t' + filename.split('/')[-1].rstrip('_fpa') + '_' + '_pagan.anctree'
    result_text = result_text + '\t' + filename.split('/')[-1].rstrip('_fpa') + '_' + '_pagan.fas'
    result_text = result_text + '\t' + filename.split('/')[-1].rstrip('_fpa') + '_' + ref_node + '.db'
    result_text = result_text + '\t' + filename.split('/')[-1].rstrip('_fpa') + '_' + query_node + '.db'
    result_text = result_text  + '\t' + filename.split('/')[-1].rstrip('_fpa') + '_' + sister_node + '.db'
    #Structural quality check for reference sequence and
    # for the query sequence if query is not a leaf.
    try:
        if (tree&query_node).is_leaf() is False:
            struct_flag = False
            found = False
            for child_list in children_structs:
                child2 = None
                double_child = None
                child1 = structure_dir + filename.split('/')[-1].rstrip('_fpa') + '_' + child_list[0]  + '.db'
                if found is False:
                    result_text = result_text + '\t' + filename.split('/')[-1].rstrip('_fpa') + '_' + child_list[0] + '.db'
                    results_dict["structureFileChild1"] = child1
                if len(child_list)==2 and found is False:
                    double_child = True
                    child2 = structure_dir + filename.split('/')[-1].rstrip('_fpa') + '_' + child_list[1]  + '.db'
                    results_dict["structureFileChild2"] = child2
                    result_text = result_text + '\t' + child2
                    found = True

                elif found is False:
                    result_text = result_text + '\t' + "None"
                    found = True

                if structure_quality(qry_ts_start,
                                     qry_ts_end,
                                     align_dict,
                                     structure_file_ref,
                                     structure_file_sister,
                                     structure_file_qry,
                                     qry_child_structure_file1=child1,
                                     qry_child_structure_file2=child2,
                                     double_child=double_child,
                                     leaf_mode=None,
                                     relaxed=True) is True:
                    struct_flag = True

            if struct_flag is False:                         

                result_text = result_text + '\t' + 'false'

            else:
                results_dict["structural_quality"] = True
                result_text = result_text + '\t' + 'true'

        else:

            if structure_quality(qry_ts_start,
                                 qry_ts_end,
                                 align_dict,
                                 structure_file_ref,
                                 structure_file_sister,
                                 structure_file_qry,
                                 leaf_mode=True,
                                 relaxed=True) is False:
                results_dict["structural_quality"] = False
                result_text = result_text + '\t' + "None" + '\t' + "None" + '\t' +'false'

            else:
                results_dict["structural_quality"] = True
                result_text = result_text + '\t' + "None" + "\t" + "None" +'\t' + 'true'

    except:

        results_dict["structural_quality"] = False
        print(filename, position, "a problem with structures.")
        return i, None, None

    minimum = min(set([qry_ts_start, qry_ts_end, ref_ts_end, ref_ts_start]))
    maximum = max(set([qry_ts_start, qry_ts_end, ref_ts_end, ref_ts_start]))
    #identify aligned structure
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

    #Identify if the TSM region overlaps with a hairpin loop.
    #Identification is done separately for both the query and reference.

    template_switch = set(range(template_switch_clusters[position][0],
                          template_switch_clusters[position][1]))

    align_qry_seq = str(align_dict[query_node])
    align_ref_seq = str(align_dict[ref_node])

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

        result_text = result_text + '\t' + loop_seq
        result_text = result_text + '\t' + adj_nucl1 + '\t' + adj_nucl2
        results_dict["loop_seq"] = loop_seq
        results_dict["adj_nucl1"] = adj_nucl1
        results_dict["adj_nucl2"] = adj_nucl2
    else:
        result_text = result_text + '\t' + 'None' + '\t' + 'None' + '\t' + 'None'

    if loop_indexes_ref is not None:
        loop_seq_ref = loop_seq_ref.replace('-', '')
        if len(loop_seq_ref) == 0:
            loop_seq_ref = 'None'
        if len(adj_nucl1_ref) == 0:
            adj_nucl1_ref = 'None'
        if len(adj_nucl2_ref) == 0:
            adj_nucl2_ref = 'None'

        result_text = result_text + '\t' + loop_seq_ref
        result_text = result_text + '\t' + adj_nucl1_ref + '\t' + adj_nucl2_ref
        results_dict["loop_seq_ref"] = loop_seq_ref
        results_dict["adj_nucl1_ref"] = adj_nucl1_ref
        results_dict["adj_nucl2_ref"] = adj_nucl2_ref


    else:
        result_text = result_text + '\t' + 'None' + '\t' + 'None' + '\t' + 'None'

    ts_ref_struct = aligned_struct_ref[qry_ts_start:qry_ts_end]
    ts_qry_struct = aligned_struct_qry[qry_ts_start:qry_ts_end]

    inverted_loop_status = 0
    palindromic = 0
    symmetric = False
    mutation = 0
    outside_overlap = 0


    #Identify if a loop sequence is inverted in RNA hairpin
    #overlapping with a identified TSM

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


    #Identify if TSM has swapped a loop sequence or part of it
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

    #identify a TMS has inverted a part of a stem

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
            inverted_loop_status >0 and palindromic == 0 and\
            len(loop_seq_ref.replace('-','')) == len(loop_seq_qry.replace('-','')):
    
                full_loop_inversion = True
    else:
        full_loop_inversion = False


    #Identify if the structure has changed in the TSM-linking area.
    if ts_ref_struct == ts_qry_struct:
        structural_change = 'No change in structure'
        result_text = result_text + '\t' + 'no_change'

    elif loop_indexes is not None and temp_struct_ref is not None:
        output = identify_if_structure_has_changed(
                loop_indexes,
                loop_indexes_ref,
                [stem1_overlap_ref, loop_overlap_ref, stem2_overlap_ref],
                [stem1_overlap, loop_overlap, stem2_overlap],
                int(template_switch_clusters[position][0]),
                temp_struct_ref,
                temp_struct)

        for out in output:

            if 'extended' in out:
                results_dict["loop_extension"] = True
            elif 'new' in out:
                results_dict["new_mismatch"] = True
            elif 'fixed' in out:
                results_dict['mismatch_fixed'] = True
            elif 'major' in out:
                results_dict["major_change"] = True

        structural_change = ', '.join(output)
        if len(structural_change) != 0:
            result_text = result_text + '\t' + structural_change
        else:
            structural_change = 'A major change'
            result_text = result_text + '\t' + 'A major change'

    elif temp_struct != temp_struct_ref:
        structural_change = 'A major change'
        results_dict["major_change"] = True
        result_text = result_text + '\t' + structural_change

    else:
        structural_change = 'major change'
        results_dict["major_change"] = True
        result_text = result_text + '\t' + structural_change


    #Identify compensating substitutions linking to a TSM

    compensating = 0
    unsymmetric_subs = 0
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

        result_text = result_text + '\t' + str(compensating)
        results_dict["compensating_mutations"] = compensating

    else:
        result_text = result_text + '\t' + str(compensating)
        results_dict["compensating_mutations"] = compensating

    #Identify if a TSM is linked to a new loop in a query sequence
    new_loop_status = False
    if loop_indexes is not None:

        new_loop_status = new_loop(ref_ts_start,
                                ref_ts_end,
                                loop_indexes,
                                [stem1_overlap_ref,
                                 loop_overlap_ref,
                                 stem2_overlap_ref],
                                [stem1_overlap,
                                 loop_overlap,
                                 stem2_overlap],
                                align_qry_seq,
                                ts_distance)


    results_dict["query_node"] = query_node
    results_dict["ref_node"] = ref_node
    results_dict["qry_ts_start"] = qry_ts_start
    results_dict["qry_ts_end"] = qry_ts_end
    results_dict["ref_ts_start"] = ref_ts_start
    results_dict["ref_ts_end"] = ref_ts_end
    results_dict["ts_mismatch"] = mismatches
    results_dict["ts_length"] = length 


    #Collect information from the TSM-linked changes in the RNA hairpin
    if compensating > 0:
        result_text = result_text + '\t' + 'True'
        results_dict["CM_present"] = True

        if (compensating == 1 and\
                inverted_loop_status == 0) or\
                (inverted_loop_status > 0 and\
                 palindromic > 0):
                    results_dict["one_CM"] = True
                    result_text = result_text + '\t' + 'True'

        else:
            result_text = result_text + '\t' +'False'

        if (compensating > 1 and\
                ((inverted_loop_status == 0) or\
                  (inverted_loop_status > 0 and\
                  palindromic > 0))):

                    results_dict["multiple_CM"] = True
                    result_text = result_text + '\t' + 'True'
        else:
            result_text = result_text + '\t' +'False'

        if inverted_loop_status > 0 and\
                palindromic == 0 and\
                full_loop_inversion is False and\
                results_dict["multiple_CM"] is False:
            
                    result_text = result_text + '\t' + 'True'
                    results_dict["CM_and_loop_change"] = True
            
        else:
            result_text = result_text + '\t' + 'False'

        if full_loop_inversion is True and symmetric is True:
            results_dict["CM_and_inverted_loop"] = True
            result_text = result_text + '\t' + 'True'
        else:
            result_text = result_text + '\t' + 'False'

    else:
        result_text = result_text + '\t' + 'False' + '\t' + 'False'+\
                '\t' +'False' +'\t'+ 'False' + '\t' + 'False'


    if compensating == 0 and\
        full_loop_inversion is True and\
            palindromic == 0 and\
            symmetric is True:
                results_dict["inverted_loop"] = True
                result_text = result_text + '\t' + 'True'

    else:
        result_text = result_text + '\t' + 'False'


    if new_loop_status is True and\
            length >= 8 and\
            mismatches >= 2:
                results_dict["new_loop"] = True
                result_text = result_text + '\t' + 'True'
    else:
        result_text = result_text + '\t' + 'False'

    if unsymmetric_subs > 0 and mismatches >= 2:
        results_dict["asymmetric_TSM_between_stems"] = True
        result_text = result_text + '\t' + 'True'
    else:
        result_text = result_text + '\t' + 'False'

    if loop_seq_mismatch > 0 and unsymmetric_subs == 0 and\
            mismatches >= 2:
        results_dict["TSM_loop_to_stem"] = True
        result_text = result_text + '\t' + 'True'
    else:
        result_text = result_text + '\t' + 'False'

    if within_stem_tsms > 0 and mismatches >= 2:
        results_dict["TSM_within_one_stem"] = True
        result_text = result_text + '\t' + 'True'
    else:
        result_text = result_text + '\t' + 'False'
    print(results_dict)

    return i, result_text, results_dict


def main():
    """
    This script identifies, if the TSMs found by FPA tool
    are linked to RNA hairpins and if they have caused:

    1) CSs
    2) partial or full loop inversions
    3) Asymmetric TSMs between the loop stems within
       the same hairpin
    4) Asymmetric TSM from the loop to stem
       within the same hairpin
    5) Asymmetric TSMs within the same stem

    PLEASE NOTE! The script expects a directory that 
    contains alignment file and tree and fpa results files.
    They should be named as identifier.fas, identifier.anctree and
    identifier_fpa, respectively. The alignment file
    fasta headers have to match the nodes of the tree files.

    python link_TSM_to_RNA_hairpins.py structure_dir/ full_path/fpa_results_file outputfile

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

    Writes the result into a .csv file

    The columns in the .csv file are separated by tabs, and
    follows the format (columns are separated with || ):

    reference_name ||query_name || result_text || tsm target in
    the alignment, the index of the last nucleotide ||
    tsm target in the alignment, the index of the last nucleotide ||
    tsm source in the alignment, the index of the last nucleotide ||
    tsm source in the alignment, the index of the last nucleotide ||
    mismatches caused by a single TSM event || length of a TSM ||
    query node is leaf (True/False) || sequences passed quality control
    (true/false) || structures passed quality control (true/false) ||
    information on the structural change || the query loop sequence,
    if TSM overlaps a hairpin in the query sequence || the reference
    loop sequence, if TSM overlaps a hairpin in the reference sequence ||
    structural change || number of compensating mutations in the
    TSM target region || TSM caused one CS (true/false) || TSM caused
    multiple CS in a stem (true/false) || TSM caused CS and partial
    loop inversion  (true/false) || TSM caused CS and a loop inversion ||
    TSM caused a loop inversion without a CS (true/false) ||
    TSM asymmetrically between stems (true/false) || TSM between loop
    and stem || TSM caused an inversion within a stem.

    """

    from tree_parse import (read_tree,
                            map_template_switch_in_same_region)

    from results_database import add_information

    structure_dir = sys.argv[1]
    fpa_directory = sys.argv[2]
    output_file = sys.argv[3]
    database = sys.argv[4]
    cpu = sys.argv[5]

    connect(database, "default")


    for subdir, dirs, files in os.walk(fpa_directory):
        for file in files:
            filename = fpa_directory + os.sep + file

            if filename.endswith('_fpa'):

                tree = read_tree(filename.rstrip('_fpa') + '_pagan.anctree')

                bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                    filename.rstrip('_fpa') + '_pagan.fas',
                                                            "fasta"))
                align_dict = OrderedDict()

                for seq in bio_align_dict:
                    align_dict[seq] = str(bio_align_dict[seq].seq)

                #Identification of the TSMs in the alignment

                template_switch_clusters = map_template_switch_in_same_region(
                        filename,
                        filename.rstrip('_fpa') + '_pagan.fas',
                        do_not_merge=True)

                #Going through all identified TSMs in the alignment
                position_list = []
                for position in template_switch_clusters:

                    position_list.append(position)

                sets = set(position_list)
                
                pool = mp.Pool(int(cpu))

                results_objects = [pool.apply_async(TSM_analysis, args=(i, filename, align_dict,
                                                                        tree,
                                                                        structure_dir,
                                                                        template_switch_clusters,
                                                                        position)) for i, position in enumerate(sets)]

                results = [r.get()[1] for r in results_objects]
                
                for a in results:

                    if a is not None:

                        with open(output_file, 'a+') as f:
                            f.write(a + '\n')

                results = [r.get()[2] for r in results_objects]
                for a in results:

                    if a is not None:

                        add_information(a)

                pool.close()
                pool.join()

        disconnect("default")

if __name__ == "__main__":
    main()
