"""
This runs sequence simulations for a set of phylogenetic trees
and alignments. The script goes through each given tree and simulates
two child sequences for each parental sequence.
In addition, TSMs are identified between parental and child sequences
by running FPA. If the TSMs are observed containing more than
two mismatches, secondary structures are predicted.
If the parent and children pass the quality control, the TSM is considered
to occur.
"""

from Bio import SeqIO
from collections import defaultdict
from decimal import Decimal
import numpy as np
import os
import re
import random
import sys
from treetime import GTR

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Runs a script for simulations.

    output_dir: A directory for output files
    filename_inf: A treefile to which the simulations
                  are run
    log_directory: A directory to .iqtree log
                   outputs from which the tree parameters
                   are extracted
    path_to_script: A full path to lambda.pl script of
                    dawg package
    """

    from collections import OrderedDict
    from tree_parse import (read_tree,
                            tree_params,
                            simulate_sequences,
                            define_insertion_and_deletion_rates)

    # connect(database, "default")

    output_dir = sys.argv[1]
    filename_int = sys.argv[2]
    log_directory = sys.argv[3]
    path_to_script = sys.argv[4]

    iupac = defaultdict(list)
    iupac['R'].extend(['A', 'G'])
    iupac['Y'].extend(['C', 'T'])
    iupac['S'].extend(['G', 'C'])
    iupac['W'].extend(['A', 'T'])
    iupac['K'].extend(['G', 'T'])
    iupac['M'].extend(['A', 'C'])
    iupac['B'].extend(['C', 'G', 'T'])
    iupac['D'].extend(['A', 'G', 'T'])
    iupac['H'].extend(['A', 'C', 'T'])
    iupac['V'].extend(['A', 'C', 'G'])
    iupac['N'].extend(['A', 'C', 'G', 'T'])

    distance_dict = defaultdict(dict)
    identifier = filename_int.split('/')[-1].split('_')[0]

    if filename_int.endswith("_pagan.anctree"):

        identifier = filename_int.rstrip("_pagan.anctree").split('/')[-1]

        bio_align_dict = SeqIO.to_dict(SeqIO.parse(
            filename_int.rstrip('_pagan.anctree') + '_pagan.fas',
            "fasta"))

        align_dict = OrderedDict()
        iupac_dict = OrderedDict()

        for seq in bio_align_dict:
            align_dict[seq] = str(bio_align_dict[seq].seq)

            iu_seq = ""
            for a in align_dict[seq]:
                if a in iupac:
                    a = random.choice(iupac[a])
                    iu_seq += a
                else:
                    iu_seq += a

            iupac_dict[seq] = iu_seq

        del bio_align_dict

        tree = read_tree(filename_int)

        file_num = filename_int.rstrip('_pagan.anctree').split('/')[-1]

        (rate_parameters,
         base_frequencies,
         gamma_shape_alpha) = tree_params(
                 log_directory + '/' + file_num + '_pagan_align.fas.log')

        GTR_c = GTR(alphabet='nuc')
        matrix = np.matrix([['-', rate_parameters["AC"],
                            rate_parameters["AG"],
                            rate_parameters["AT"]],
                            [rate_parameters["AC"],
                             '-',
                             rate_parameters["CG"],
                             rate_parameters["CT"]],
                            [rate_parameters["AG"],
                             rate_parameters["CG"],
                             '-',
                             rate_parameters["GT"]],
                            [rate_parameters["AT"],
                             rate_parameters["CT"],
                             rate_parameters["GT"],
                             '-']])

        pi = [base_frequencies["A"],
              base_frequencies["C"],
              base_frequencies["G"],
              base_frequencies["T"]]

        GTR_c = GTR_c.custom(mu=1.0, pi=pi, W=matrix)

        leaf_dict = {}
        for node in tree.traverse("postorder"):
            if node.is_root():
                continue
            if node.is_leaf():
                leaf_dict[node.name] = ""

            distance_dict[node.name] = GTR_c.optimal_t(
                    np.array(list(iupac_dict[(node.up).name])),
                    np.array(list(iupac_dict[node.name])))

            node.dist = distance_dict[node.name]

        with open(filename_int, 'r') as f:
            for line in f:
                tree_text = line.rstrip()
        tree_text = tree_text.replace(')', '&')

        for item in distance_dict:
            match1 = re.match(
                    r'(.+%s\s*).+?(\s*%s.+)' % (item + ':', ','),
                    tree_text, re.DOTALL)
            match2 = re.match(
                    r'(.+%s\s*).+?(\s*%s.+)' % (item + ':', '&'),
                    tree_text, re.DOTALL)
            try:
                match1.groups()
            except:
                print("no groups")
            match2.groups()
            replacement = str(distance_dict[item])
            if match1 is None:
                tree_text = match2.group(1) + replacement + match2.group(2)

            elif len(match1.group(2)) > len(match2.group(2)):
                tree_text = match1.group(1) + replacement + match1.group(2)
            else:
                tree_text = match2.group(1) + replacement + match2.group(2)

        tree_text = tree_text.replace('&', ')')
        with open(output_dir + "all/" + identifier + ".nw", 'w') as f:
            f.write(tree_text)

        (gapModel,
         gapRate,
         SeqLen,
         indelRate) = define_insertion_and_deletion_rates(
                 output_dir + "all/" + identifier + ".nw",
                 filename_int.rstrip('_pagan.anctree') + '_pagan_align.fas',
                 path_to_script)

        choose_root = random.choice(list(leaf_dict))
        ins_rate = Decimal(indelRate)/2
        del_rate = Decimal(indelRate)/2

        simulate_sequences(base_frequencies,
                           rate_parameters,
                           gamma_shape_alpha,
                           SeqLen,
                           tree_text,
                           identifier,
                           iupac_dict[choose_root],
                           SeqLen,
                           output_dir + "all/",
                           gapModel=gapModel,
                           gapRate=gapRate,
                           ins_rate=ins_rate,
                           del_rate=del_rate)

        with open(output_dir + "all/" + identifier + '.fasta', 'r') as f1:

            for line in f1:

                if '>' in line:

                    with open(output_dir + "all/" +
                              identifier + '.fas', 'a+') as f:
                        f.write(line)

                elif '>' not in line:
                    with open(output_dir + "all/" +
                              identifier + '.fas', 'a+') as f:

                        f.write(line)


if __name__ == "__main__":
    main()
