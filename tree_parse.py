"""
Tools for parsing phylogenetic trees.

Author: Heli Mönttinen (ORCID: 0000-0003-2461-0690) 

"""

import re

from Bio import SeqIO
from collections import (defaultdict,
                         OrderedDict)
from decimal import Decimal
from ete3 import Tree, PhyloTree
import numpy

from  pymongo import MongoClient
import random

from RNA_alignments import (indexes_in_original_seq,
                            indexes_in_alignment,
                            unusual_iupac_char)
from fpa_tools import fpa_parsing
import statistics


def read_tree(tree_file):
    """
    Read a tree from a file,

    Requires

    tree_file: A full path to the tree file

    Returns

    a loaded tree

    """

    tree = Tree(tree_file, format=1)

    return tree


def go_through_branches(dictomy_tree):
    """
    This tool goes systematically through tree branches
    and samples systematically the tree by returning
    pairs of branches.

    Requires:

    tree: a tree from which polytomy is resolved

    Yields:
        Node name, identifier of sequence1 and identifier of sequence 2
    """
    node_number = 0

    for node in dictomy_tree.traverse("postorder"):

        if node.is_leaf() is not True:
            sample_seq1 = node.children[0].name
            sample_seq2 = node.children[1].name

            yield node.name, sample_seq1, sample_seq2

        node_number += 1


def write_tree(tree, outfile):
    """
    Writes a tree into a file.
    """

    tree.write(outfile=outfile)


def location_of_template_switch(fpa_file, alignment_file):
    """
    The location of the template switch in the alignment.
    
    Parameters
    -------------
    :fpa_output: The output file of fpa
    :alignment_file: Alignment as a dictionary format.

    Yields
    --------
    Indexes as tuple
    (The identifier od the parent,
     the identifier of the child,
     the sequence of the parent,
     the sequence of the child,
     the starting index of the tsm source region in the alignment,
     the ending index of the tsm source region in the alignment,
     the starting index of the tsm target region in the alignment,
     the ending index of the tsm target region in the alignment)

    """

    fpa_dictionary = fpa_parsing(fpa_file)
    align_dict = SeqIO.to_dict(SeqIO.parse(alignment_file, "fasta"))

    for hit in fpa_dictionary:

        for case in fpa_dictionary[hit]:

            if int(case["sw_end"]) > 0:

                template_switch_indexes = [int(case["start"]),
                                           int(case["end"]),
                                           int(case["sw_start"]),
                                           int(case["sw_end"])-1]
            else:
                template_switch_indexes = [int(case["start"]),
                                           int(case["end"]),
                                           int(case["sw_start"]),
                                           int(case["sw_end"])]

            alignment_indexes2 = indexes_in_alignment(
                    template_switch_indexes[0:2], align_dict[case["query"]].seq)

            alignment_indexes = indexes_in_alignment(
                    template_switch_indexes[2:], align_dict[case["ref"]].seq)

            if len(alignment_indexes2) > 1 and len(alignment_indexes) > 1:
                yield (case["ref"],
                       case["query"],
                       case["seq_ref"],
                       case["seq_qry"],
                       alignment_indexes2[0],
                       alignment_indexes2[1],
                       alignment_indexes[0],
                       alignment_indexes[1])


def map_template_switch_in_same_region(fpa_file,
                                        alignment_file,
                                        do_not_merge=None):
    """
    Identifies if there are several template switch in the same
    region of the alignment. Returns a regions dictionary with
    region_numbers as keys. The values are tuple-formatted containing
    reference, query, sequence, ind1, ind2, ind3 and ind4.

    Requires
    --------
    :fpa_file: a fpa_file output with sequence pair information
    :alignment_file: an alignment file

    Returns
    -------

    :regions: Regions defaultdict list. Each key corresponds
              to one template switch region. A list contains
              information on the node pairs that has a template
              switch in the same region. The first index of the
              list tells the starting index of the region, and
              second the ending index of the region.

    """

    regions = defaultdict(list)

    template_switch_region = 0
    used = []

    for ref, query, seq_qry, seq_ref, start, end, sw_start, sw_end in\
            location_of_template_switch(fpa_file, alignment_file):
                
        template_switch = set(list(range(start, end+1)))

        found = False

        if do_not_merge is True:
            minimum = int(start)
            maximum = int(end)

            regions[template_switch_region].append(minimum)
            regions[template_switch_region].append(maximum)
            regions[template_switch_region].append(
                    (ref, query, seq_qry, seq_ref, start, end, sw_start, sw_end))

            template_switch_region += 1

        if do_not_merge is None:
            for region in regions:
                region_range = list(
                        range(regions[region][0], regions[region][1]+1))

                if len(template_switch.intersection(region_range)) > 0:

                    minimum = min(set(region_range).union(template_switch))
                    maximum = max(set(region_range).union(template_switch))

                    regions[region][0] = minimum
                    regions[region][1] = maximum

                    regions[region].append(
                            (ref, query, seq_qry, seq_ref, start, end, sw_start, sw_end))
                    found = True
                    break

            if found is False:
                regions[template_switch_region].append(start)
                regions[template_switch_region].append(end)
                regions[template_switch_region].append(
                        (ref, query, seq_qry, seq_ref, start, end, sw_start, sw_end))

                template_switch_region += 1        

    return regions



def fix_indexes_based_on_position(position_list,
                                  index,
                                  is_list=False):
    """
    Takes a position list (designed especially for the subtree
    visualization) and substracts from given indexes the number of
    those positions that are smaller than the given index.

    :position_list: list of indexes that are removed form the original
                    alignment
    :index_list: Indexes that should be adjusted to the altered
                 alignment length.
    """

    smallerThan = lambda x,y: [i for i in x if i<y]

    if is_list is True:
        fixed_index = []
        for index1 in index:

            smaller = smallerThan(position_list, index1)

            new_index = index1 - len(smaller)

            fixed_index.append(new_index)
        return fixed_index

    else:
        smaller = smallerThan(position_list, index)
        fixed_index = index - len(smaller)

        return fixed_index


def subtree_for_visualization(tree,
                              ancestor_node,
                              alignment_dict):
    """
    This script is designed for large trees that are impossible
    to visualize. Extracts a subtree and creates a new alignment
    dictionary for a subset of sequences.

    :subtree: a sub tree of interest
    :alignment_dict: information of alignment dictonary
    :tree_leaves: optional parameter that identifies, how many leaves
                  a tree has to have that this script is run. default: 200

    Returns:
        new_align_dict: A dictionary of aligned sequences from which
                        all the positions
                        removed that are gaps in all the aligned sequences.
    """

    new_align_dict = OrderedDict()


    subtree = tree&ancestor_node
    names = set()
    for node in subtree.traverse():
        names.add(node.name)

    for node in alignment_dict:
        if node in names:
            new_align_dict[node] = alignment_dict[node]
            last = node

    positions = []
    for position in range(len(new_align_dict[last])):

        gap = True

        for seq in new_align_dict:

            if new_align_dict[seq][position] != '-':
                gap = False
                continue

        if gap is True:
            positions.append(position)

    positions.reverse()

    for pos in positions:

        for seq in new_align_dict:

             new_align_dict[seq] = new_align_dict[seq][0 : pos : ] + new_align_dict[seq][pos + 1 : :]

    return new_align_dict, positions


def mismatches_between_sequences(seq1,
                                 seq2,
                                 strict_mode=None):
    """
    Returns the number of mismatches
    and the length of alignment between two sequences (does not
    take into account those positions which both have a gap).


    Arguments
    ==========

    Seq1: an aligned sequence 1
    Seq2: an aligned sequence 2
    strict_mode: if true, gaps are not calculated as mismatches
                 (default: None)

    Reuturns
    =========
    :mismatches: Number of  mismatches
    :length: Length of the TSM

    """

    mismatches = 0
    length = 0

    if strict_mode is True:
        for a, b in zip(seq1, seq2):

            if a != b and a != '-' and b != '-':

                mismatches += 1

            if a != '-' or b != '-':

                length += 1

        return mismatches, length


    for a, b in zip(seq1, seq2):

        if a != b:

            mismatches += 1

        if a != '-' or b != '-':

            length += 1


    return mismatches, length


def sequence_quality(seq1, seq2, relaxed=None, iupac=None):
    """
    The aim of this function is to confirm sequence
    quality of a predicted sequence.
    By default the mismatches between these
    two sequences should be less than 10 %.
    In the relaxed mode the mismatches has to be less than
    30%.  The script should return True, the sequences
    fulfil the requirement.

    :seq1: The first sequence to compare (string)
    :seq2: The second sequence to compare (string)
    :relaxed: If True, a relaxed mode is applied.
    :iupac: if True, the uncertain iupac chars in the same
            position of the alignment are interpreted as
            the same base between seq1 and seq2, if the the groups
            the possible bases overlap. 

    Returns :
        Boolean
    """

    if iupac != None:
        seq1, seq2 = unusual_iupac_char(seq1, seq2)
    
    mismatches, length = mismatches_between_sequences(seq1, seq2)

    if length == 0:
        return False

    if relaxed is True:
        if Decimal(mismatches)/Decimal(length) <= 0.3:
            return True

    if Decimal(mismatches)/Decimal(length) <= 0.1:
        return True

    return False


def compare_predicted_sequence(seq,
                               children_seqs,
                               ancestor_seq=None,
                               relaxed=None,
                               iupac=None):
    """
    The aim is to confirm the quality of a predicted seq.
    Compares a sequence to the ancestral sequence
    and its possible child sequences. The sequence has to be always
    similar to the ancestral sequence, one of the child
    sequences.

    seq: The predicted sequence (aligned) which quality is wanted
         to check
    ancestor_seq: The ancestor sequence (aligned) of the seq
    children seqs: A set or a list of child sequences

    Returns:
    True if the sequence is of high quality
    """

    if ancestor_seq:
        if sequence_quality(seq, ancestor_seq, relaxed=relaxed, iupac=iupac) is False:

            return False

    for child in children_seqs:
        if sequence_quality(seq, child, relaxed=relaxed, iupac=iupac) is False:
            return False

    return True



def reference_sequence_quality(reference_orig,
                               reference_ts,
                               sister_orig,
                               sister_ts,
                               qry_orig,
                               ancestor_orig,
                               ancestor_ts,
                               relaxed=None):

    """
    Ensures that ..
    1) Ts-region of reference sequence is similar to sister
       sequence and to the ancestral sequence

    2) Source ts-region has to be similar to the sister
       sequence or ts_sequence and ancestral sequence.

    Arguments
    ----------
    
    reference_orig: The origin sequence of the ts in the reference
    reference_ts: The template switch sequence in the reference
    sister_orig: The origin sequence in the sister sequence
    sister_ts: The template switch in the sister sequence
    ancestor_orig: The origin sequence in the ancestor of
    the reference ancestor_ts: The template switch sequence
    in the ancestor of reference sequence

    Output
    --------

    True if all the requirements are fulfilled. False, if some of them fails.

    """

    if compare_predicted_sequence(reference_ts,
                                  set([sister_ts]),
                                  ancestor_seq=ancestor_ts,
                                  relaxed=relaxed,
                                  iupac=True) is False:
        return False


    if compare_predicted_sequence(reference_orig,
                                  set([sister_orig]),
                                  ancestor_seq=ancestor_orig,
                                  relaxed=relaxed,
                                  iupac=True) is False:


        if compare_predicted_sequence(reference_orig,
                                      set([qry_orig]),
                                      ancestor_seq=ancestor_orig,
                                      relaxed=relaxed,
                                      iupac=True) is False:

            return False

    return True


def events_per_substitution(database,
                            collection,
                            level=1,
                            alternative_query={}):
    """
    Counts TSM events per tree branch in the given dataset.
     

    :database: in which all the events are reported.
    :collection: In which all the events are reported.
    :level: level from which the template switch and branch distances
            are calculated. Level 1 means terminal leaf nodes.
            Default=1.
    
    """

    mongo_client = MongoClient()
    db = mongo_client[database]

    collection = db[collection]

    query_list = []

    for criteria in ["multiple_CM", "CM_and_loop_change",
            "CM_and_inverted_loop", "inverted_loop", "new_loop" ]:

        query_dict = {}
        query_dict["ts_quality"] = True
        query_dict["structural_quality"] = True
        query_dict[criteria] = True
        query_dict["ts_mismatch"] = {"$gt": 1 }
        query_list.append(query_dict)

    query = {}
    tar_nodes_list = []

    if len(alternative_query) != 0:
        for field, values in alternative_query.items():
            query[field] = values
    
    tsm_collection = defaultdict(dict)
    tsm_substitution = defaultdict(dict) 
    align_dict = {}
    trees_dict = {}

    for query in query_list:
        cursor = collection.find(query)

        count = 0
        for document in cursor:

            tree_file = document["tree_file"].split('//')
            tree_file = '/'.join(tree_file)
            case = ""
            for x in query:
                if x in ["multiple_CM", "CM_and_loop_change",
            "CM_and_inverted_loop", "inverted_loop", "new_loop" ]:
                    case = x

            query_node = document["query_node"]
            trees_dict[tree_file] = ""


            if case not in tsm_collection:
                tsm_collection[case] = defaultdict(dict)
                tsm_substitution[case] = defaultdict(dict)

                tsm_substitution[case][tree_file][query_node] =\
                        document["ts_mismatch"]

                tsm_collection[case][tree_file][query_node] = 1


            elif tree_file not in tsm_collection[case]:

                tsm_substitution[case][tree_file][query_node] =\
                                document["ts_mismatch"]
                tsm_collection[case][tree_file][query_node] = 1


            elif query_node not in tsm_collection[case][tree_file]: 
                tsm_substitution[case][tree_file][query_node] =\
                                document["ts_mismatch"]
                tsm_collection[case][tree_file][query_node] = 1


            else:
                tsm_substitution[case][tree_file][query_node] +=\
                        document["ts_mismatch"]
                                
                tsm_collection[case][tree_file][query_node] += 1


    tsm = 0
    tsm_subs = 0
    branch_length = 0
    branch_count = 0
    root_count = 0
    leaf_count = 0

    levels = defaultdict(set)
    for lev in range(0,level):
        levels[lev] = set()

    target_nodes = set()
    for treefile in trees_dict:   
        tree = PhyloTree(treefile, format=1) 

        for node in tree.traverse("levelorder"):
            if node.is_leaf():
                leaf_count += 1

                node_name = node.name

                target = node
                levels[0].add(treefile + '=' + target.name)
                for i in range(1,level):
                    if not target.is_root():
                        target = target.up
                        levels[i].add(treefile + '=' + target.name)
                    else:
                        continue

    level_dict = {}

    tsm = defaultdict(dict)
    tsm_subs = defaultdict(dict)
    tar_nodes = defaultdict(dict)

    for case in ["multiple_CM", "CM_and_loop_change",
            "CM_and_inverted_loop", "inverted_loop", "new_loop"]:

        tsm[case][lev] = 0
        tsm_subs[case][lev] = 0
        tar_nodes[case] = defaultdict(list)
        for lev in [0, 1, 2]:
            tsm[case][lev] = 0
            tsm_subs[case][lev] = 0
            tar_nodes[case] = defaultdict(list)


    for lev in levels:
        for tar_node in levels[lev]:
            treex = tar_node.split('=')[0]

            node_name = tar_node.split('=')[-1]
            
            for case in ["multiple_CM", "CM_and_loop_change",
            "CM_and_inverted_loop", "inverted_loop", "new_loop"]:

                if treex in tsm_collection[case]:
                    if node_name in tsm_collection[case][treex]:
                        tsm[case][lev] += tsm_collection[case][treex][node_name]
                        tsm_subs[case][lev] += tsm_substitution[case][treex][node_name]

                        tar_nodes[case][lev].append(tsm_collection[case][treex][node_name])

                    else:
                        tar_nodes[case][lev].append(0)
                else:
                    tar_nodes[case][lev].append(0)


    for case in ["multiple_CM", "CM_and_loop_change",
            "CM_and_inverted_loop", "inverted_loop", "new_loop"]:
        for lev in [0, 1, 2]:
            rate3 = round((Decimal(tsm[case][lev])/Decimal(len(levels[lev]))), 10)
            std = round(numpy.std(tar_nodes[case][lev]), 10)
            mean = round(numpy.mean(tar_nodes[case][lev]),  10)
            print(case, lev, tsm[case][lev], len(levels[lev]),  rate3, len(trees_dict))
            
            yield len(trees_dict), case, lev, tsm[case][lev], tsm_subs[case][lev], rate3, mean, std


def cut_align_seq_position(align_dict,
                           qry_ts_start,
                           qry_ts_end,
                           ref_ts_start,
                           ref_ts_end):
    """
    Cut sequences to match the sequence region of interest. 
    Leaves -+ 20 residues around the region of interest.
    Returns a dictionary containing the cut sequences.

    Arguments:
    ------------
    :align_dict: Alignment file in dictionary format
    :qry_ts_start: The starting index of the TSM target region
    :qry_ts_end: The ending index of the TSM target source
    :ref_ts_start: The starting index of the TSM source region
    :ref_ts_end: The encing index of the TSM source region

    Returns:
    --------

    :minimum: The smallest index within TSM source
              and target regions
    :maximum: The largest index within the TSM source and
              target regions
    :cut_start: The first index of the cut sequence 
    :cut_end: The last index of the cut seqeucne
    :cut_align_dict: A dictionary containing the cut sequences
    """

    minimum = 6000
    maximum = 0


    num_set = set([qry_ts_start, qry_ts_end, ref_ts_start, ref_ts_end])

    if  min(num_set) < minimum:
        minimum = min(num_set)

    if max(num_set) > maximum:
        maximum = max(num_set)

    cut_align_dict = {}
    cut_start = 0
    cut_end = 0

    for align in align_dict:

        new_name = align.replace('_', '.')

        if minimum > 20 and (len(align_dict[align]) - maximum) > 20:
            cut_align_dict[new_name] =\
                    str(align_dict[align][minimum - 20:maximum + 20])
                    
            if cut_end == 0:
                cut_start = minimum - 20
                cut_end = maximum + 20

        elif  minimum <= 20 and (len(align_dict[align]) - maximum) <= 20:
            cut_align_dict[new_name] =\
                    str(align_dict[align][0:-1])

            if cut_end == 0:
                cut_start = 0
                cut_end = len(align_dict[align]) - 1

        elif minimum <= 20:
            cut_align_dict[new_name] =\
                str(align_dict[align][0:maximum + 20])

            if cut_end == 0:
                cut_start = 0
                cut_end = maximum + 20

        elif (len(align_dict[align]) - maximum) <= 20:
            cut_align_dict[new_name] =\
                    str(align_dict[align][minimum-20:-1])

            if cut_end == 0:
                cut_start = minimum -20
                cut_end = len(align_dict[align]) - 1

    return minimum, maximum, cut_start, cut_end, cut_align_dict
