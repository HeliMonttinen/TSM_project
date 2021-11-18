from Bio import SeqIO
from collections import (defaultdict,
                         OrderedDict,
                         Counter)
import os
import sys
import multiprocessing as mp

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def identify_conserved_areas(i,
                             align_dict,
                             info):
    """
     Identifies high quality regions for a parental
    sequence.
    Identifies identical bases between a parent and child
    sequences, and between a parent and an ancestral
    sequence. Identifies also identical characters between
    parent and child dot-pranthesis structures. Uncertain
    IUPAC charaters are unified before the comparison.
    A node sequence residue is considered to be a part of
    a high quality region if it is involved in a 10-residues-long
    window containing >=9/10 identical bases with either of child
    sequences, >=9/10 identical bases with an ancestral
    sequence and >=7/10 identical secondary structure characters
    between the node and either of child sequences.

    Arguments:
    ==========
    :i: index for a Pool worker
    :align_dict: A dictionary containing aligned sequences
    :info: A tuple containing an aligned node (parent) sequence,
           its' both child sequences, its' parent (ancestor),
           secondary structure (dor-parenthesis) for both children
           and the node structure.

    Returns
    =======
    :good_indexes: Indexes of high quality regions in an aligned sequence
    :info[7]: node name
    :info[9]: child1 name
    tset[10]: child2 name
    good_indexes_orig: Indexes for high quality regions in non-aligned
                       node sequence
    parent_len: An effective length for a parental sequence


    """

    from RNA_alignments import (unusual_iupac_char,
                                indexes_in_alignment)

    from structure import identify_alignment_area

    aligned_parent = info[0]
    aligned_seq1 = info[1]
    aligned_seq2 = info[2]
    anc_ancestor = info[3]
    struct_seq1 = info[4]
    struct_seq2 = info[5]
    struct_parent = info[6]

    child1 = []
    child2 = []
    anc_seq = []
    child1_struct = []
    child2_struct = []

    parent_numbers = []

    _, aligned_struct_seq1 = identify_alignment_area(
            struct_seq1,
            align_dict,
            0,
            10)

    _, aligned_struct_seq2 = identify_alignment_area(
            struct_seq2,
            align_dict,
            0,
            10)

    _, aligned_struct_parent = identify_alignment_area(
            struct_parent,
            align_dict,
            0,
            10)

    for base in range(len(aligned_parent)):
        if aligned_parent[base] == '-':
            continue

        seq1_char = aligned_seq1[base]
        seq2_char = aligned_seq2[base]
        parent_char_seq1 = aligned_parent[base]
        parent_char_seq2 = aligned_parent[base]

        parent_char_seq1, seq1_char = unusual_iupac_char(aligned_parent[base],
                                                         seq1_char)

        parent_char_seq2, seq2_char = unusual_iupac_char(aligned_parent[base],
                                                         seq2_char)

        parent_char, anc_char = unusual_iupac_char(aligned_parent[base],
                                                   anc_ancestor[base])

        if seq1_char == parent_char_seq1:

            child1.append(1)
        else:
            child1.append(0)

        if seq2_char == parent_char_seq2:

            child2.append(1)
        else:
            child2.append(0)

        if parent_char == anc_char:
            anc_seq.append(1)

        else:
            anc_seq.append(0)

        if aligned_struct_seq1[base] == aligned_struct_parent[base]:
            child1_struct.append(1)
        else:
            child1_struct.append(0)

        if aligned_struct_seq2[base] == aligned_struct_parent[base]:
            child2_struct.append(1)
        else:
            child2_struct.append(0)

    parent_len = 0
    for number in range((len(child1)-10)+1):
        flag = False
        if number > 10:
            up_num = 10
        else:
            up_num = number

        for i in range(0, up_num):

            count_seq1 = sorted(
                    (Counter(child1[number-i:number+11-i])).items())
            count_seq2 = sorted(
                    (Counter(child2[number-i:number+11-i])).items())
            count_anc = sorted(
                    (Counter(anc_seq[number-i:number+11-i])).items())

            count_seq1_str = sorted(
                    (Counter(child1_struct[number-i:number+11-i])).items())
            count_seq2_str = sorted(
                    (Counter(child2_struct[number-i:number+11-i])).items())

            if count_seq1[0][0] == 1:
                item = (0, 0)
                prev_count = count_seq1[0]
                count_seq1[0] = item
                count_seq1.append(prev_count)

            if count_seq2[0][0] == 1:
                item = (0, 0)
                prev_count = count_seq2[0]
                count_seq2[0] = item
                count_seq2.append(prev_count)

            if count_anc[0][0] == 1:
                item = (0, 0)
                prev_count = count_anc[0]
                count_anc[0] = item
                count_anc.append(prev_count)

            if count_seq1_str[0][0] == 1:
                item = (0, 0)
                prev_count = count_seq1_str[0]
                count_seq1_str[0] = item
                count_seq1_str.append(prev_count)

            if count_seq2_str[0][0] == 1:
                item = (0, 0)
                prev_count = count_seq2_str[0]
                count_seq2_str[0] = item
                count_seq2_str.append(prev_count)

            if int(count_seq1[0][1]) < 2 and int(count_anc[0][1]) < 2 and\
                    int(count_seq1_str[0][1]) < 4 and flag is False:

                flag = True
                parent_numbers.append(1)
                parent_len += 1

            elif int(count_seq2[0][1]) < 2 and int(count_anc[0][1]) < 2 and\
                    int(count_seq2_str[0][1]) < 4 and flag is False:

                flag = True
                parent_numbers.append(1)
                parent_len += 1

        if flag is False:
            parent_numbers.append(0)

    good_indexes = []
    good_indexes_orig = []

    good_ind = []

    good = False
    good_length = 0

    for i in range(len(parent_numbers)):

        if parent_numbers[i] == 1 and\
                good is False:
            good_ind.append(i)
            good = True

        elif parent_numbers[i] == 0 and\
                good is True:

            good_ind.append(i)
            good = False

        elif i == (len(parent_numbers)-1) and\
                good is True:
            good_ind.append(i)
            good = False

        if len(good_ind) == 2:

            orig_alignment = indexes_in_alignment(good_ind, aligned_parent)

            good_length += (orig_alignment[1]-orig_alignment[0])

            good_indexes.append(good_ind)
            good_indexes_orig.append(orig_alignment)
            good_ind = []

    return (i, good_indexes, info[7], info[9],
            info[10], good_indexes_orig, parent_len)


def main():
    """
        Runs a identify_conserved_area function
    for a internal tree nodes. Output is written into
    a .csv file.

    Arguments
    =========
    structure_dir: A directory where dot_parenthesis
                   files are located
                   The structures has to be named as
                   'clusteridentifier_nodename.db'
    tree_directory: A directory, in which tree
                    files are located.
                    Tree files has to be names as
                    'clusteridentifier_xxx.anctree'
    outputfile: A name for output file
    level: A level how deeply a given tree is iterated
    cpu: Number of CPUs for the run
    alignment_dir: A directory in which alignment files are located.

    Output
    ======

    .csv file

    Columns are: identifier_level, a node name, a child node1,
                 a child node2, effective length of a parent,
                 indexes for high quality region in unaligned
                 node sequence and indexes for high quality regions
                 in aligned node sequence.

    """

    from tree_parse import read_tree

    structure_dir = sys.argv[1]
    tree_directory = sys.argv[2]
    outputfile = sys.argv[3]
    level = int(sys.argv[4])
    cpu = int(sys.argv[5])
    alignment_dir = sys.argv[6]

    for subdir, dirs, files in os.walk(tree_directory):

        for file in files:

            if not str(file).endswith('anctree'):
                continue

            filename = tree_directory + os.sep + file
            tree = read_tree(filename)

            identifier = str(file).split('_')[0]

            if not os.path.exists(alignment_dir + identifier + '.fas'):

                print(file, "file missing")
                continue

            bio_align_dict = SeqIO.to_dict(SeqIO.parse(
                alignment_dir + identifier + '.fas',
                "fasta"))

            align_dict = OrderedDict()

            for seq in bio_align_dict:
                align_dict[seq] = str(
                        bio_align_dict[seq].seq).replace('+', '-')

            set_data = []
            used_tar = set()
            levels = defaultdict(set)
            for node in tree.iter_leaves():

                target = node

                for i in range(1, level):

                    if not target.is_root():
                        target = target.up
                        if target.is_root():
                            continue
                    else:
                        continue
                    if target.name in used_tar:
                        continue
                    else:
                        used_tar.add(target.name)
                    seq1 = align_dict[((target.children)[0]).name]
                    seq2 = align_dict[((target.children)[1]).name]
                    parent_seq = align_dict[target.name]
                    ancestor = align_dict[(target.up).name]
                    seq1_struct = structure_dir + identifier +\
                        '_' + ((target.children)[0]).name + '.db'
                    seq2_struct = structure_dir + identifier +\
                        '_' + ((target.children)[1]).name + '.db'
                    parent_struct = structure_dir + identifier +\
                        '_' + target.name + '.db'
                    if (not os.path.exists(seq1_struct)) or\
                        (not os.path.exists(seq2_struct)) or\
                            (not os.path.exists(parent_struct)):
                        print(file, target.name, "node missing")
                        continue

                    anc_seqid = (target.up).name
                    sample_seq1 = ((target.children)[0]).name
                    sample_seq2 = ((target.children)[1]).name

                    to_pool = (parent_seq, seq1, seq2, ancestor,
                               seq1_struct, seq2_struct,
                               parent_struct, target.name, anc_seqid,
                               sample_seq1, sample_seq2)

                    if i not in levels:
                        levels[i] = set()
                    levels[i].add(to_pool)

            for lev in levels:

                set_data = list(levels[lev])

                if len(set_data) < cpu:
                    pool = mp.Pool(len(set_data))

                else:
                    pool = mp.Pool(cpu)

                good_indexes = []
                print(to_pool[7])

                results_objects = [pool.apply_async(
                    identify_conserved_areas, args=(i, align_dict, tuplei))
                            for i, tuplei in enumerate(set_data)]

                results1 = [r.get()[1] for r in results_objects]
                results2 = [r.get()[2] for r in results_objects]
                results3 = [r.get()[3] for r in results_objects]
                results4 = [r.get()[4] for r in results_objects]
                results5 = [r.get()[5] for r in results_objects]
                results6 = [r.get()[6] for r in results_objects]

                for a in range(len(results1)):
                    good_indexes = []
                    good_indexes_orig = []

                    for r in results1[a]:

                        good_indexes.append(r)

                    parent = results2[a]

                    ch1 = results3[a]

                    ch2 = results4[a]
                    par_len = results6[a]

                    for r2 in results5[a]:
                        good_indexes_orig.append(r2)

                    with open(outputfile, 'a') as f:

                        f.write(identifier + '_' + str(lev) + '\t')
                        f.write(parent + '\t')
                        f.write(ch1 + '\t')
                        f.write(ch2 + '\t')

                        f.write(str(par_len) + '\t')

                        for item in good_indexes:

                            f.write(str(item) + ';')

                        f.write('\t')

                        for item in good_indexes_orig:

                            f.write(str(item) + ';')

                        f.write('\n')

                pool.close()
                pool.join()


if __name__ == "__main__":
    main()
