"""
This script is designed for drawing .svg figures from selected
xases in database. Requires the TSM identifier
pointing to the database.

    python3 draw_figs_from_database.py output_dir
                    database collection identifier

:output_dir: The name of a directory, in which the figure is saved
:database: A TSM database name
:collection: A collection name
:identifier: The identifier of the TSM

"""
import os
import sys
from Bio import SeqIO
from ete3 import Tree
from collections import OrderedDict
from pymongo import MongoClient

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """creates figures from ts cases that match the query.
       Change the #query = {} line to search a case of
       interest from a database
    """

    from tree_parse import (fix_indexes_based_on_position,
                            subtree_for_visualization,
                            cut_align_seq_position)

    from structure import identify_alignment_area

    from visualize_svg import (hit,
                               files,
                               params,
                               make_fig)

    output_dir = sys.argv[1]
    database = sys.argv[2]
    collection = sys.argv[3]

    query = {"identifier":"16602_24"}

    client = MongoClient()
    db = client[database]
    collection = db[collection]

    for document in collection.find(query):

        hit_dict = {}

        for key, value in document.items():

            hit_dict[key] = value

        orig_tree = Tree(hit_dict["tree_file"], format=1)

        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                    hit_dict["alignment_file"],
                                    "fasta"))
        orig_align_dict = OrderedDict()

        for seq in bio_align_dict:
            orig_align_dict[seq] = str(bio_align_dict[seq].seq)

        position_0 = min(set([hit_dict["qry_ts_start"],
                              hit_dict["qry_ts_end"],
                              hit_dict["ref_ts_start"],
                              hit_dict["ref_ts_end"]]))

        position_1 = max(set([hit_dict["qry_ts_start"],
                              hit_dict["qry_ts_end"],
                              hit_dict["ref_ts_start"],
                              hit_dict["ref_ts_end"]]))

        ancestor_node = hit_dict["ancestor_node"]
        qry_ts_start = hit_dict["qry_ts_start"]
        qry_ts_end = hit_dict["qry_ts_end"]
        ref_ts_start = hit_dict["ref_ts_start"]
        ref_ts_end = hit_dict["ref_ts_end"]
        query_node = hit_dict["query_node"]
        ref_node = hit_dict["ref_node"]
        sister_node = hit_dict["sister_node"]
        if "child_node1" in hit_dict:

            qry_child_structure_node1 = hit_dict["child_node1"]
            structure_file_qry_child1 = hit_dict["structure_child1"]
        else:
            qry_child_structure_node1 = None
            structure_file_qry_child1 = None

        if "child_node2" in hit_dict:

            qry_child_structure_node2 = hit_dict["child_node2"]
            structure_file_qry_child2 = hit_dict["structure_child2"]

        else:
            qry_child_structure_node2 = None
            structure_file_qry_child2 = None

        structure_file_ref = hit_dict["structure_ref"]
        structure_file_qry = hit_dict["structure_qry"]
        structure_file_sister = hit_dict["structure_sister"]

        if len(orig_align_dict) > 60:

            align_dict, position_list = subtree_for_visualization(
                    orig_tree,
                    ancestor_node,
                    orig_align_dict)
            tree = orig_tree&ancestor_node
            qry_ts_start = fix_indexes_based_on_position(position_list,
                                                         qry_ts_start)
            qry_ts_end = fix_indexes_based_on_position(position_list,
                                                       qry_ts_end)

            ref_ts_start = fix_indexes_based_on_position(position_list,
                                                         ref_ts_start)

            ref_ts_end = fix_indexes_based_on_position(position_list,
                                                       ref_ts_end)

            position_0 = fix_indexes_based_on_position(position_list,
                                                       position_0)

            position_1 = fix_indexes_based_on_position(position_list,
                                                       position_1)

        else:
            tree = orig_tree
            align_dict = orig_align_dict

        try:
            (minimum,
             maximum,
             cut_start,
             cut_end,
             cut_align_dict) = cut_align_seq_position(
                    align_dict,
                    qry_ts_start,
                    qry_ts_end,
                    ref_ts_start,
                    ref_ts_end)

            ts_region_qry, aligned_struct_qry = identify_alignment_area(
                    structure_file_qry,
                    align_dict,
                    cut_start,
                    cut_end)

            ts_region_ref, aligned_struct_ref = identify_alignment_area(
                    structure_file_ref,
                    align_dict,
                    cut_start,
                    cut_end)
        except:
            continue

        outputfile = output_dir + hit_dict["cluster"] + '_'\
            + str(hit_dict["position"]) + '.svg'

        if not (tree&query_node).is_leaf():

            showseq = [ancestor_node, ref_node, query_node,
                       sister_node, qry_child_structure_node1,
                       qry_child_structure_node2]
            showstruct = [ref_node, query_node, sister_node,
                          qry_child_structure_node1,
                          qry_child_structure_node2]

            struct_files = {
                    ref_node: structure_file_ref,
                    query_node: structure_file_qry,
                    sister_node: structure_file_sister,
                    qry_child_structure_node1: structure_file_qry_child1,
                    qry_child_structure_node2: structure_file_qry_child2}

        else:

            showseq = [ancestor_node, ref_node, query_node, sister_node]
            showstruct = [ref_node, query_node, sister_node]

            struct_files = {ref_node: structure_file_ref,
                            query_node: structure_file_qry,
                            sister_node: structure_file_sister}

        if qry_ts_end > len(aligned_struct_ref):
            qry_ts_end = len(aligned_struct_ref)
            ref_ts_end = len(aligned_struct_ref)

        hit_class = hit(cut_start,
                        cut_end,
                        ref_node,
                        query_node,
                        qry_ts_start,
                        qry_ts_end,
                        ref_ts_start,
                        ref_ts_end,
                        aligned_struct_ref,
                        aligned_struct_qry,
                        showseq,
                        showstruct)

        par = params()

        par.lmar = 350

        files_class = files(align_dict,
                            outputfile,
                            struct_files)

        make_fig(files_class,
                 hit_class,
                 par,
                 tree)


if __name__ == "__main__":
    main()
