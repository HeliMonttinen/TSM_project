"""
This script is for identifing instantaneous CSs in which
the initial mutations and CM are in the same sequence. 
Goes automatically through all the .fasta files in the directory.
All fasta files have to be in a format identifier.fasta. The
same directory has to contain also tree files following
format identifier.anctree.

Script requires python3.7

The script is run:

python find_classical_CS.py fpa_dir/ structure_dir outputfile number_of_cores

Arguments
==========

alignment_dir: Directory, where the alignments and trees are located.  
structure_dir: A directory, where the predicted secondary structures are
               located
output_name: An output file identidier that is used for specifying the written
             output files.
number_of_cores: The wanted number of the cores


Author: Heli Mönttinen  (ORCID: 0000-0003-2461-0690)

"""


from Bio import SeqIO
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
            node name, the child 2 node name, information if a child 1
            is a leaf node (true/false) and information if a child 2
            is a leaf node (true/false).


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
                leaf nodes
    nonleaf_count: A count for instantaneous CS containing loops
                   in the internal nodes
    mismatches: A dictionary counting frequencies, how many
                mutations has occurred in a single instaneous CS event
    lengths: A dictionary for counting cases of different lengths.
             Shows how long perfect Watson-crick base pairing region
             was on the both sides of the loop, in which the
             instantaneous CS was located.

    """

    from RNA_alignments import unusual_iupac_char
    from RNA_loops import _join_dot_file_lines
    from tree_parse import (read_tree,
                            go_through_branches)
    from structure import identify_immediate_CS


    count = 0
    loop_count = {}
    adj_nucl_count = {}
    dict_structs = {}
    leaf_count = 0
    nonleaf_count = 0
    rev_comp = {}
    mismatches = {}
    lengths = {}

    anc_seqid = tuplei[0]
    sample_seq1 = tuplei[1]
    sample_seq2 = tuplei[2]
    s1_leaf = tuplei[3]
    s2_leaf = tuplei[4]

    sample_list = []

    for seq in [anc_seqid, sample_seq1, sample_seq2]:

        
        anc_up = ((tree&anc_seqid).up).name

        try:
            structure = _join_dot_file_lines(structure_dir  + alignment_file + '_' + seq + '.db')
            dict_structs[seq] = structure
            if seq != anc_seqid:
                sample_list.append(seq)

        except:
            if seq == anc_seqid:
                break
            continue


    child1 = None
    child2 = None
    
    for x in sample_list:

        if x == sample_seq1:
            sis_seq = sample_seq2
        else:
            sis_seq = sample_seq1
        if not (tree&x).is_leaf():
            child1= fasta_dict[(tree&x).children[0].name]
            child2= fasta_dict[(tree&x).children[1].name]

        count1, loop, adj1, revs, mismatch, length = identify_immediate_CS(
                dict_structs[anc_seqid],
                dict_structs[x],
                fasta_dict[anc_seqid],
                fasta_dict[x],
                fasta_dict[anc_up],
                fasta_dict[sis_seq],
                child1=child1,
                child2=child2)

        count += count1
        for loopx in loop:
            if x == sample_seq1 and s1_leaf is True:
                leaf_count += loop[loopx]
            elif x == sample_seq1 and s1_leaf is False:
                nonleaf_count += loop[loopx]
            if x == sample_seq2 and s2_leaf is True:
                leaf_count += loop[loopx]
            elif x == sample_seq2 and s2_leaf is False:
                nonleaf_count += loop[loopx]

            if loopx in loop_count:
                loop_count[loopx] += loop[loopx]
            else:
                loop_count[loopx] = loop[loopx]

        for adjx in adj1:
            if adjx in adj_nucl_count:

                adj_nucl_count[adjx] += adj1[adjx]
            else:
                adj_nucl_count[adjx] = adj1[adjx]

        for r in revs:
            if r in rev_comp:
                rev_comp[r] += revs[r]
            else:
                rev_comp[r] = revs[r]

        for r in mismatch:
            if r in mismatches:
                mismatches[r] += mismatch[r]
            else:
                mismatches[r] = mismatch[r]

        for r in length:
            if r in lengths:
                lengths[r] += length[r]
            else:
                lengths[r] = length[r]



    return (i, count, loop_count, adj_nucl_count, rev_comp, leaf_count, nonleaf_count, mismatches, lengths)


def main():
    """
    This function paralellizes the sequence comparisons
    collects the results and writes the output files.

    """

    from tree_parse import (read_tree,
                            go_through_branches)
                            

    fpa_dir = sys.argv[1]
    structure_dir = sys.argv[2]
    output_name = sys.argv[3]
    number_of_cores = sys.argv[4]


    all_count = 0
    leaf_count = 0
    nonleaf_count = 0
    loop_count = {}  
    adj_nucl_count = {}
    mismatches = {}
    lengths = {}



    alignment_file_list = []

    for subdir, dirs, files3 in os.walk(fpa_dir):
        for file in files3:
            if file.endswith('.fas'):

                alignment_file = file.split('/')[-1].rstrip('_pagan.fas')

                bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                                    filename.rstrip('_fpa') + '.fas',
                                                            "fasta"))
                fasta_dict = OrderedDict()

                for seq in bio_align_dict:
                    fasta_dict[seq] = str(bio_align_dict[seq].seq)
                
                if len(fasta_dict) > 19:
                    
                    tree = read_tree(fpa_dir + alignment_file + ".anctree")

                    sets = list()
                    for anc_seqid, sample_seq1, sample_seq2 in go_through_branches(tree):

                        if (tree&anc_seqid).is_root():
                            continue

                        samp1_leaf = (tree&sample_seq1).is_leaf()
                        samp2_leaf = (tree&sample_seq2).is_leaf()

                        sets.append((anc_seqid, sample_seq1, sample_seq2, samp1_leaf, samp2_leaf))

                    pool = mp.Pool(number_of_cores)   

                    results_objects = [pool.apply_async(
                        count_loops,
                        args=(i, tree, fasta_dict, alignment_file, structure_dir, tuplei)) for i, tuplei in enumerate(sets)]

                    results = [r.get()[2] for r in results_objects]

                    for r in results:
                        for x in r:
                            if x in loop_count:
                                loop_count[x] += r[x]
                            else:
                                loop_count[x] = r[x]
    
                    results = [r.get()[3] for r in results_objects]

                    for r in results:
                        for x in r:
                                
                            if x in adj_nucl_count:
                                adj_nucl_count[x] += r[x]
                            else:
                                adj_nucl_count[x] = r[x]



                    results = [r.get()[5] for r in results_objects]

                    for r in results:
                        leaf_count += r

                    results = [r.get()[6] for r in results_objects]

                    for r in results:
                        nonleaf_count += r

                    results = [r.get()[7] for r in results_objects]

                    for r in results:
                        for x in r:
                            if x in mismatches:
                                mismatches[x] += r[x]
                            else:
                                mismatches[x] = r[x]


                    results = [r.get()[8] for r in results_objects]

                    for r in results:
                        for x in r:
                            if x in lengths:
                                lengths[x] += r[x]
                            else:
                                lengths[x] = r[x]



                    pool.close()
                    pool.join()

    with open("classical_CS_loops_" + output_name + ".txt", 'w') as f:
        for i in loop_count:
            f.write(str(i) + '\t' + str(loop_count[i]) + '\n')
    with open("classical_CS_closing_pair" + output_name + "txt", 'w') as f:
        for i in adj_nucl_count:
            f.write(str(i) + '\t' + str(adj_nucl_count[i]) + '\n')

    with open("mismatches_" + output_name, 'w') as f:
        f.write("mismatches\n")
        for i in mismatches:
            f.write(str(i) + '\t' + str(mismatches[i]) + '\n')
        f.write("lenghts\n")
        for i in lengths:
            f.write(str(i) + '\t' + str(lengths[i]) + '\n')

    with open("node_info_" + output_name, 'w') as f:

        f.write(str("all cases: " + all_count) + '\n')
        f.write(str("Classical CS in leafs: " + leaf_count) + '\n')
        f.write(str("Classical CS in non-leafs: " + nonleaf_count) + '\n')


if __name__ == "__main__":
    main()
