"""
This script is designed for going through a set of phylogenetic trees
and extracting RNA loop sequences from their interal nodes. In addition,
the script ensures the quality of the predicted RNA hairpin by comparing
the hairpin sequence to the parent and child nodes. In addition, either of 
the child nodes have to have a loop in the corresponding region.

The script is run: 

    python3 extract_information_all_loops.py cluster_dir output_file_loop_count ouput_file_loop_seq_count structure_dir cpu

cluster_dir = A directory, in which alignment_file, tree file and fpa output file
              The alignment and tree files have to also contain the ancestral sequences.
              The tree file has to contain names for the internal nodes too.

output_file_loop_count = An output file to which the loop counts are\
                         are written

ouput_file_loop_seq_count = An output to which the loop sequences and their\
        counts are written

structure_dir = A directory, which contains dot-parenthesis (.dp) files for\
                the all the nodes of the phylogenetic tree.

cpu = number of cpus


Author: Heli Mönttinen (ORCID: 0000-0003-2461-0690)


"""

from Bio import SeqIO
from decimal import Decimal
from collections import defaultdict
from collections import OrderedDict
import numpy as np
import os
import re
from  pymongo import MongoClient
import sys
import multiprocessing as mp

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)



def run_loop_study(i, align_dict, structure_dir, filepath, cluster, tuplei):
    """
    Extracts loops from the reference sequence and identifies their quality by
    comparing the loop sequence to the corresponding sequence regions in the 
    child sequence. The reference loop sequence has to be at least 90% similar to 
    either of the child sequences. In addition, reference loop cannot contain 
    uncertain iupac characters.

    Arguments
    ==========

    :i: an index required for pooling
    :align_dict: A dictionary containing alignments
    :structure_dir: A path to directory containing protein structures
    :filepath: A full path to a tree file (output from pagan software)
               
    :cluster: identifier of the sequence cluster
    :tuplei: A tuple, which contains (identifier for a reference sequence,
                                      identifier for the first child of the reference,
                                      identifier for the second child of the reference,
                                      identifier for the ancestor of the reference)


    Output
    ======

    :i: an index required for pooling
    :accept_pairs: A dictionary containing counts for 
                   the nucleotide pairs in stem region.
                   Only one base-pair long mismatches are allowed.
                   If longer, counting breaks.
    :loop_refs: A dictionary containing 

    """
    from tree_parse import compare_predicted_sequence
    from RNA_alignments import indexes_in_alignment
    from RNA_loops import (RNA_loops,
                           _join_dot_file_lines,
                           _get_sequence_dot_parenthesis,
                           make_reverse_complement)

    ref = tuplei[0]
    child1 = tuplei[1]
    child2 = tuplei[2]
    ancestor = tuplei[3]

    loop_refs = {}
    accept_pairs = {"GU":0, "UG":0, "AU":0, "UA":0,
                    "GC":0, "CG":0}


    if not os.path.isfile(structure_dir + filepath.split('/')[-1].rstrip('_pagan.anctree') + '_' + ref  + '.db'):

        return (i, accept_pairs, "none")

    dot_parenthesis = _join_dot_file_lines(structure_dir + cluster + '_' + ref + ".db")
    rna_sequence = _get_sequence_dot_parenthesis(structure_dir + cluster + '_' + ref + ".db")

    RNA_sequence = RNA_loops(dot_parenthesis, rna_sequence)
    
    RNA_sequence.sequence_for_loops()

    loop_inds_ref = RNA_sequence.loop_indexes

    loop_ind_qry_dict = {}

    for child in [child1, child2]:
        if os.path.exists(structure_dir + cluster + '_' + child + ".db"):

            child_seq = align_dict[child]
            dot_child_parenthesis = _join_dot_file_lines(structure_dir + cluster + '_' + child + ".db")
            rna_child_sequence = _get_sequence_dot_parenthesis(structure_dir + cluster + '_' + child + ".db")
            RNA_sequence_c = RNA_loops(dot_child_parenthesis, rna_child_sequence)
            RNA_sequence_c.sequence_for_loops()
            loop_inds_child = RNA_sequence_c.loop_indexes


            for x in range(len(loop_inds_child)):
                loop_ind_qry = indexes_in_alignment(loop_inds_child[x],
                                                    child_seq)

                loop_ind_qry_dict[str(loop_ind_qry[1]) + ',' + str(loop_ind_qry[2])] = loop_ind_qry
        else:
            return (i, accept_pairs, "none")


    set_loops = list()
    for index, seq in list(zip(loop_inds_ref,
        RNA_sequence.loop_sequences)):

        orig_inds = indexes_in_alignment(index,
                align_dict[ref])

        ind_key = str(orig_inds[1]) + ',' + str(orig_inds[2])
        if ind_key not in loop_ind_qry_dict.keys():

            continue


        flag = False
        seq_ok = True

        seq_set = set(seq)
        for x in seq_set:
            if x not in {'A':"", 'G':"", 'C':"", 'T':"", 'U':""}:
                seq_ok = False
        if seq_ok is False:
                
            continue


        for kid in {child1, child2}:
            child_seq = set([align_dict[kid][orig_inds[0]:orig_inds[3]]])
            ancestor_seq = align_dict[ancestor][orig_inds[0]:orig_inds[3]]
            if compare_predicted_sequence(align_dict[ref][orig_inds[0]:orig_inds[3]],
                                        child_seq,
                                        ancestor_seq=ancestor_seq,
                                        iupac=True) is True:
                flag = True
            
            if flag is False:
                return (i, accept_pairs, "none")


        loop_start = index[1] - index[0]
        loop_end = index[2] - index[0]


        adj_nucl1 = seq[loop_start-1]
        adj_nucl2 = seq[loop_end]
        stem2 = seq[loop_end:-1]
        stem1 = seq[:loop_start][::-1]

        for n1, n2 in list(zip(stem1, stem2)):

            if n1+n2 in accept_pairs:

                accept_pairs[n1+n2] +=1
            else:
                break

        loop_sequence = seq[loop_start:loop_end]

        loop_reverse = make_reverse_complement(loop_sequence)
        loop_part = make_reverse_complement(loop_sequence[0])
        loop_ref = ""
        if loop_sequence == loop_reverse:
            loop_ref = loop_sequence + '/' + adj_nucl1 + '/' + adj_nucl2 + '/' + 'fullrev'
            if loop_ref in loop_refs:
                loop_refs[loop_ref] += 1
            else:
                loop_refs[loop_ref] = 1

        elif loop_part == loop_sequence[-1]:
            loop_ref = loop_sequence + '/' + adj_nucl1 + '/' + adj_nucl2 + '/' + 'part_rev'

            if loop_ref in loop_refs:
                loop_refs[loop_ref] += 1
            else:
                loop_refs[loop_ref] = 1

        else:
            loop_ref = loop_sequence + '/' + adj_nucl1 + '/' + adj_nucl2 + '/' +'non_rev'
            if loop_ref in loop_refs:
                loop_refs[loop_ref] += 1
            else:
                loop_refs[loop_ref] = 1


    return (i, accept_pairs, loop_refs)

        

def main():
    """
    The main script for running the analysis.
    Requires user's input.
    """

    from common import return_file
    from RNA_alignments import indexes_in_alignment
    from RNA_loops import (RNA_loops,
                           _join_dot_file_lines,
                           _get_sequence_dot_parenthesis,
                           make_reverse_complement)
    from tree_parse import (read_tree,
                            compare_predicted_sequence)

    cluster_dir = sys.argv[1]
    output_file_loop_count = sys.argv[2]
    ouput_file_loop_seq_count = sys.argv[3]
    structure_dir = sys.argv[4]
    cpu = int(sys.argv[5])

    all_loops = 0
    loop_leafes = 0
    loop_non_leafes = 0
    loops_palindromes = 0
    loop_leafes_palin = 0
    loop_non_leafes_palin = 0
    loop_leaf_palin_partial = 0
    loop_non_leafes_palin_partial = 0

    reverse_complement_sequences = {}
    non_reverse_complement_sequences = {}
    partial_rev_comp = {}

    path_list = []
    sizes_list = []
    all_sizes = 0

    for filename in return_file(cluster_dir):

        if filename.endswith("_pagan.anctree"):
  
            number = filename.rstrip("_pagan.anctree").split('/')[-1]

            try:
                if os.path.exists(cluster_dir + os.sep + str(number) + "_pagan.fas"):
                    name = cluster_dir + os.sep + str(number) + "_pagan.fas"
                    count = 0
                    with open(name, 'r') as f:
                        for line in f:
                            if '>' in line and '#' not in line:
                                count += 1
                    path_list.append(cluster_dir + os.sep + str(number) + "_pagan.anctree")
                    sizes_list.append(count)
                    all_sizes += count

            except:
                continue


    file_count = 0
    palin_seq = 0
    node_count = 0
    loop_palin_partial = 0

    accept_pairs = {"GU":0, "UG":0, "AU":0, "UA":0,
                    "GC":0, "CG":0}


    for filename in path_list:
        file_count += 1

        filepath = filename

        tree = read_tree(filepath)

        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                filepath.rstrip('_pagan.anctree') + '_pagan.fas',
                "fasta"))

        align_dict = OrderedDict()

        for seq in bio_align_dict:
            align_dict[seq] = str(bio_align_dict[seq].seq)
            cluster =  filepath.split('/')[-1].rstrip('_pagan.anctree')

        del bio_align_dict

        set_loops = list()
        for node in tree.traverse("postorder"):

            if node.is_root() or node.is_leaf():
                continue


            if not node.is_root():
                node_count+= 1
                ancestor = (node.up).name
                child1 = (node.children[0]).name
                child2 = (node.children[1]).name

                set_loops.append((node.name, child1, child2, ancestor))


        if len(set_loops) < cpu:
            pool = mp.Pool(len(set_loops))
        else:
            pool = mp.Pool(cpu)

        results_objects = [pool.apply_async(run_loop_study, args=(i, align_dict, structure_dir, filepath, cluster, tuplei)) for i, tuplei in enumerate(set_loops)]

        all_loops += 1

        results = [r.get()[1] for r in results_objects]

        for r in results:

            for x in r:
                accept_pairs[x] += r[x]


        results = [r.get()[2] for r in results_objects]

        for r in results:

            if not isinstance(r, str) :

                for x in r:
                    if 'full' in x:
                        string = '/'.join(x.split('/')[:-1])
                        if string in reverse_complement_sequences:
                            reverse_complement_sequences[string] += r[x]
                        else:
                            reverse_complement_sequences[string] = r[x]
                        loops_palindromes += 1
                        all_loops += 1

                    elif 'part' in x:
                        string = '/'.join(x.split('/')[:-1])

                        if string in partial_rev_comp:
                            partial_rev_comp[string] += r[x]
                        else:
                            partial_rev_comp[string] = r[x]
                        all_loops += 1

                    elif "none" not in x:
                        string = '/'.join(x.split('/')[:-1])

                        if string in non_reverse_complement_sequences:
                            non_reverse_complement_sequences[string] += r[x]
                        else:
                            non_reverse_complement_sequences[string] = r[x]

                        all_loops += 1

        pool.close()
        pool.join()

    with open(ouput_file_loop_seq_count, 'w') as f:
        
        rev_sorted = sorted(reverse_complement_sequences.items(), key=lambda x:x[1], reverse=True)
        for seq, val in rev_sorted:
            adj1 = seq.split('/')[1]
            adj2 = seq.split('/')[2]
            seq_letters = seq.split('/')[0]
            f.write(adj1 + '\t' + adj2 + '\t' + seq_letters + '\t' + str(val) + '\t' + 'reverse_complement' + '\n')

        partial_sorted = sorted(partial_rev_comp.items(), key=lambda x:x[1], reverse=True)
        for seq, val in partial_sorted:
            adj1 = seq.split('/')[1]
            adj2 = seq.split('/')[2]
            seq_letters = seq.split('/')[0]
            f.write(adj1 + '\t' + adj2 + '\t' + seq_letters + '\t' + str(val) + '\t' + 'partial_reverse_complement' + '\n')

        non_reverse_sorted = sorted(non_reverse_complement_sequences.items(), key=lambda x:x[1], reverse=True)
        for seq, val in non_reverse_sorted:
            adj1 = seq.split('/')[1]
            adj2 = seq.split('/')[2]
            seq_letters = seq.split('/')[0]
            f.write(adj1 + '\t' + adj2 + '\t' + seq_letters + '\t' + str(val) + '\t' + 'non_reverse_complement' + '\n')
   
    with open(output_file_loop_count, 'w+') as f3:
        f3.write('file_count: ' + str(file_count) + '\n')
        f3.write('all_nodes: ' + str(node_count) + '\n')
        f3.write('all_loops: ' + str(all_loops) + '\n')
        f3.write('sequence_count: ' + str(all_sizes) + '\n')
        f3.write('palin_seq: ' +  str(palin_seq) + '\n')
        f3.write('loop_palin_partial: ' +  str(loop_palin_partial) + '\n')

        for c in accept_pairs:
            f3.write(c + '\t' + str(accept_pairs[c]) + '\n')


if __name__ == "__main__":
        main()
