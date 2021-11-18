"""
This script is for identifing instantaneous CSs in which
the initial mutations and CM are in the same sequence.
Goes automatically through all the .fas files in the directory.
All fasta files have to be in a format identifier.fas or identifier_pagan.fas.
The same directory has to contain also tree files following
format identifier.anctree.

Script requires python3.7

Arguments
==========

tree_dir: A directory in which all tree files are located
structure_dir: A directory, where the predicted secondary structures are
               located
output_name: An output file identidier that is used for specifying the written
             output files.
number_of_cores: The wanted number of the cores

regions_file: A .csv file containing information on high quality regions
              for each internal node sequence

alignment_dir: Directory, where the alignments are located

leaf_mode: True or False. True if instantaneous CMs are collected from
           leaf nodes.

Author: Heli Mönttinen  (ORCID: 0000-0003-2461-0690)

"""
from Bio import SeqIO
import copy
import os
import sys
import multiprocessing as mp

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def count_loops(i,
                tree,
                fasta_dict,
                alignment_file,
                structure_dir,
                tuplei):
    """
    This function is for identifying instantaneous CSs
    between the parent and child nodes. This function is
    designed to be given for a process that runs node comparisons
    within a tree in parallel.

    Arguments
    ==========
    i: A subprocess count needed for the paralellization
    tree: A phylogenetic tree read by the ete3 library tool

    fasta_dict: Alignment dictionary, where fasta headers are
                the keys and they has to correspond the node
                names of the treei

    structure_dir: A path to a directory containing the
                   dot-parenthesis files. The file names
                   have in correct format identifier_nodename.db

    tuplei: a set that includes the parental node name, the child 1
            node name and regions dictionary containig information on
            high quality regions.

    Returns
    =======

    Function returns a tuple containing:

    i: The subprocess count
    count: a count for successful comparisons
    loop_count: a total count for the loops
    adj_nucl_count: count for different closing pairs (dictionary)
    rev_comp: If an instantaneous CS is linked to a loop inversion,
              they are collected into this dictionary
    leaf_count: A count for intantaneous CS containing loops in
                leaf nodes:
    nonleaf_count: A count for instantaneous CS containing loops
                   in the internal nodes
    mismatches: A dictionary counting frequencies, how many
                mutations has occurred in a single instaneous CS event
    lengths: A dictionary for counting cases of different lengths.
             Shows how long perfect Watson-crick base pairing region
             was on the both sides of the loop, in which the
             instantaneous CS was located.

    """
    from RNA_loops import _join_dot_file_lines
    from structure import identify_immediate_CS

    dict_structs = {}

    anc_seqid = tuplei[0]
    sample_seq1 = tuplei[1]
    regions = tuplei[2]

    structure1 = _join_dot_file_lines(
            structure_dir + alignment_file + '_' + anc_seqid + '.db')
    structure2 = _join_dot_file_lines(
            structure_dir + alignment_file + '_' + sample_seq1 + '.db')
    dict_structs[anc_seqid] = structure1
    dict_structs[sample_seq1] = structure2
    child1 = None
    child2 = None

    if not (tree&sample_seq1).is_leaf():
        child1 = fasta_dict[(tree&anc_seqid).children[0].name]
        child2 = fasta_dict[(tree&anc_seqid).children[1].name]

    (count1, loop, adj1, revs, mismatch,
     length, mut_types, used_ind) = identify_immediate_CS(
            dict_structs[anc_seqid],
            dict_structs[sample_seq1],
            fasta_dict[anc_seqid],
            fasta_dict[sample_seq1],
            regions,
            child1=child1,
            child2=child2)

    return (i, count1, loop, adj1, revs, mismatch,
            length, mut_types, used_ind, child1)


def main():
    """
    This function paralellizes the sequence comparisons
    collects the results and writes the output files.

    """

    from collections import (defaultdict,
                             OrderedDict)

    from tree_parse import (read_tree,
                            go_through_branches)

    tree_dir = sys.argv[1]
    structure_dir = sys.argv[2]
    output_name = sys.argv[3]
    number_of_cores = int(sys.argv[4])
    regions_file = sys.argv[5]
    alignment_dir = sys.argv[6]
    leaf_mode = str(sys.argv[7])

    all_count = 0
    leaf_count = 0
    nonleaf_count = 0
    loop_count = {}
    adj_nucl_count = {}
    mismatches = {}
    lengths = {}
    mutation_types = {}

    regions_dict = defaultdict(dict)

    print("done")
    with open(regions_file) as f:

        for line in f:
            line_splitted = line.rstrip().split('\t')
            identifier = line_splitted[0].split('_')[0]
            parent = line_splitted[1].rstrip()
            try:
                regions = line_splitted[6].split(';')
            except:
                continue

            kid1 = line_splitted[2]
            kid2 = line_splitted[3]

            if leaf_mode == "True":
                if '#' in kid1 and '#' in kid2:
                    continue
            else:
                if ('#' not in kid1) and ('#' not in kid2):
                    continue

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
    print("done")

    for subdir, dirs, files3 in os.walk(tree_dir):
        for file in files3:

            if file.endswith('.anctree'):
                used = set()

                identifier = file.split('/')[-1].rstrip('_pagan.anctree')
                print(identifier)
                if identifier not in regions_dict:
                    continue

                if os.path.exists(alignment_dir + identifier + '.fas'):
                    filename = alignment_dir + identifier + '.fas'
                else:
                    filename = alignment_dir + identifier + '_pagan.fas'

                bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                    filename,
                                    "fasta"))

                fasta_dict = OrderedDict()

                for seq in bio_align_dict:
                    fasta_dict[seq] =\
                            str(bio_align_dict[seq].seq).replace('+', '-')

                if len(fasta_dict) > 19:

                    tree = read_tree(tree_dir + identifier + "_pagan.anctree")

                    sets = list()
                    for anc_seqid, sample_seq1, sample_seq2 in\
                            go_through_branches(tree):

                        if anc_seqid not in regions_dict[identifier]:
                            continue

                        regs = copy.deepcopy(
                                regions_dict[identifier][anc_seqid])

                        samp1_leaf = (tree&sample_seq1).is_leaf()
                        samp2_leaf = (tree&sample_seq2).is_leaf()

                        if leaf_mode == "True":

                            if samp1_leaf is True:
                                sets.append((anc_seqid, sample_seq1, regs))
                            if samp2_leaf is True:
                                sets.append((anc_seqid, sample_seq2, regs))
                        else:
                            if samp1_leaf is False:
                                sets.append((anc_seqid, sample_seq1, regs))
                            if samp2_leaf is False:
                                sets.append((anc_seqid, sample_seq2, regs))

                    pool = mp.Pool(number_of_cores)

                    results_objects = [pool.apply_async(
                        count_loops,
                        args=(i, tree, fasta_dict,
                              identifier, structure_dir,
                              tuplei)) for i, tuplei in enumerate(sets)]

                    results1 = [r.get()[1] for r in results_objects]
                    results2 = [r.get()[2] for r in results_objects]
                    results3 = [r.get()[3] for r in results_objects]
                    results5 = [r.get()[5] for r in results_objects]
                    results6 = [r.get()[6] for r in results_objects]
                    results7 = [r.get()[7] for r in results_objects]
                    results8 = [r.get()[8] for r in results_objects]
                    results9 = [r.get()[9] for r in results_objects]

                    for a in range(len(results8)):
                        for rx in results8[a]:
                            if results8[a][rx] in used:
                                continue
                            elif len(results8[a]) == 0:
                                continue
                            else:
                                used.add(results8[a][rx])

                                for x in results2[a][rx]:
                                    if x in loop_count:
                                        loop_count[x] += results2[a][rx][x]
                                    else:
                                        loop_count[x] = results2[a][rx][x]

                                for x in results3[a][rx]:

                                    if x in adj_nucl_count:
                                        adj_nucl_count[x] += results3[a][rx][x]
                                    else:
                                        adj_nucl_count[x] = results3[a][rx][x]

                                if results9[a] is None:

                                    leaf_count += results1[a][rx]
                                else:

                                    nonleaf_count += results1[a][rx]

                                all_count += results1[a][rx]

                                for x in results5[a][rx]:

                                    if x in mismatches:
                                        mismatches[x] += results5[a][rx][x]
                                    else:
                                        mismatches[x] = results5[a][rx][x]

                                for x in results6[a][rx]:
                                    if x in lengths:
                                        lengths[x] += results6[a][rx][x]
                                    else:
                                        lengths[x] = results6[a][rx][x]

                                for x in results7[a][rx]:
                                    if x in mutation_types:
                                        mutation_types[x] += results7[a][rx][x]
                                    else:
                                        mutation_types[x] = results7[a][rx][x]

                    pool.close()
                    pool.join()

    with open("classical_CS_loops_" + output_name + ".txt", 'w') as f:
        for i in loop_count:
            f.write(str(i) + '\t' + str(loop_count[i]) + '\n')
    with open("classical_CS_closing_pair" + output_name + ".txt", 'w') as f:
        for i in adj_nucl_count:
            f.write(str(i) + '\t' + str(adj_nucl_count[i]) + '\n')

    with open("mismatches_" + output_name, 'w') as f:
        f.write("mismatches\n")
        for i in mismatches:
            f.write(str(i) + '\t' + str(mismatches[i]) + '\n')
        f.write("lenghts\n")
        for i in lengths:
            f.write(str(i) + '\t' + str(lengths[i]) + '\n')
        f.write("mutation_types\n")
        for i in mutation_types:
            f.write(str(i) + '\t' + str(mutation_types[i]) + '\n')

    with open("node_info_" + output_name, 'w') as f:

        f.write("all cases: " + str(all_count) + '\n')
        f.write("Classical CS in leafs: " + str(leaf_count) + '\n')
        f.write("Classical CS in non-leafs: " + str(nonleaf_count) + '\n')


if __name__ == "__main__":
    main()
