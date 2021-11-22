"""
Tools for working on phylogenetic trees.

Author: Heli Mönttinen (ORCID: 0000-0003-2461-0690)

"""

from Bio import SeqIO
from collections import (defaultdict,
                         OrderedDict)
from decimal import Decimal
from ete3 import Tree, PhyloTree
import numpy

from pymongo import MongoClient
from subprocess import check_output

from RNA_alignments import indexes_in_alignment
from fpa_tools import fpa_parsing


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


def get_node_db(tree, database, collection):
    """
    get features for nodes from a database.

    Parameters
    ----------
    :tree:  tree in a read format
    :database: a database from which data is retrived.

    Returns

    :feature_dict: a dictionary of node features in a tree

    """

    mongo_client = MongoClient()
    db = mongo_client[database]

    collection = db[collection]

    features_dict = defaultdict(dict)

    unknown_species = 1

    for node in tree:

        node_dict = {}

        if node.is_leaf():

            name = node.name

            if '-' in node.name:

                name = rename_identifier(node.name)

            document = collection.find_one({"identifier": name})

            for key, value in document.items():
                node_dict[key] = value

        if 'species' not in node_dict:
            node_dict['species'] = 'unknown species ' + str(unknown_species)
            unknown_species += 1

        node_dict['species'] = node_dict['species'].split('\t')[0]

        features_dict[node.name] = node_dict

    return features_dict


def rename_identifier(node_name):
    """
    Renames a node identifier to match that of in a database.

    :node_name: a node name identifier to rename

    Returns

    :renamed_node: A renamed node
    """

    char_found = False
    renamed_node = ""

    rev = range(len(node_name))

    for i in rev[::-1]:

        if node_name[i] == '-':
            char_found = True

        if char_found is True and node_name[i] == '.':
            renamed_node = '/' + renamed_node
            char_found = False
        elif char_found is True and node_name[i] == '_':
            renamed_node = '/' + renamed_node
            char_found = False

        else:
            renamed_node = node_name[i] + renamed_node

    return renamed_node


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


def make_family_tax_template_switch_title(tree,
                                          database,
                                          collection,
                                          feature_dict):

    """
    Makes a title and legend for a tree.
    """

    family = common_family_description(feature_dict, database, collection)

    feature_tree = add_features_to_node(tree, feature_dict)

    tax1, tax2, tax3 = highest_common_taxonomy(feature_tree)

    if len(tax3) > 0:

        taxonomy = tax1 + ';' + tax2 + ';' + tax3

    else:
        taxonomy = tax1 + ';' + tax2

    text = 'Family: ' + family + '\n' + 'Taxonomy: ' + taxonomy

    return text


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

    if isinstance(fpa_file, dict):
        fpa_dictionary = fpa_file
    else:
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
                    (ref, query, seq_qry, seq_ref,
                     start, end, sw_start, sw_end))

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
                            (ref, query, seq_qry, seq_ref,
                             start, end, sw_start, sw_end))
                    found = True
                    break

            if found is False:
                regions[template_switch_region].append(start)
                regions[template_switch_region].append(end)
                regions[template_switch_region].append(
                        (ref, query, seq_qry, seq_ref,
                         start, end, sw_start, sw_end))

                template_switch_region += 1

    return regions


def common_family_description(feature_dict, database, collection):
    """
    Returns a shared family for a tree and description
    for it.

    :feature_dict1: features_dict from identifier database
    :feature_dict2: features_fict from family database
    """

    family = []
    description = []
    highest_count = {}
    for feature in feature_dict:
        fam_id = feature_dict[feature]["family"]
        family.append(fam_id)
        if fam_id in highest_count and fam_id != "Silva":
            highest_count[fam_id] += 1
        elif fam_id != "Silva":
            highest_count[fam_id] = 1

    mongo_client = MongoClient()
    db = mongo_client[database]

    collection = db[collection]
    family_set = list(set(family))

    for fam in family_set:

        try:

            document = collection.find_one({"family": fam})

            description.append(document["description"])
        except:
            if fam == 'Silva':
                description.append('Silva Ribosomal RNA')

    prev_highest = 0
    highest = ""
    for item in highest_count:
        if item != 'Silva' and highest_count[item] > prev_highest:
            highest = item
            prev_highest = highest_count[item]

    text = ""
    for i in range(len(family_set)):

        if family_set[i] == highest:
            text = family_set[i] + ' - ' + description[i]

        elif family_set[i] == 'Silva' and len(family_set) == 1:

            text = family_set[i] + ' - ' + description[i]

    return text


def add_features_to_node(tree,
                         feature_dict,
                         regions=False,
                         alignment_dict=None):
    """
    Adds features to the node and identifies.

    Parameters:

    :tree: tree read
    :features_dict: features in a dictionary format
    :regions: if regions is True also corresponding secondary
              structure is affed
    """

    for node in tree:

        if node.name in feature_dict:

            for feature, value in feature_dict[node.name].items():

                node.add_feature(feature, value)

    return tree


def highest_common_taxonomy(tree_features):
    """
    Searches a tree and returns the highest taxonomy
    level that is shared by all the nodes.

    :tree_features: a tree with added features

    Returns
    -------

    :taxonomic_level: taxonomic level shared by all the leaves

    """

    upper_taxes = []
    all_lineages = []
    life_domain_count = {}

    for node in tree_features:

        if node.get_ascii(attributes=["tax_lineage"]):
            lineages = node.get_ascii(
                    attributes=[
                        "tax_lineage"]).lstrip().lstrip('--').split(';')

            if lineages[0] == 'cellular organisms':
                del lineages[0]

                cleaned_lineage = []
                unclassified = None
                level = 0

                main_lin = lineages[0].lstrip().rstrip()
                if main_lin in life_domain_count:
                    life_domain_count[main_lin] += 1
                elif 'unclass' not in main_lin:
                    life_domain_count[main_lin] = 1

                for item in lineages:
                    item = item.lstrip().rstrip()
                    if len(item) > 0:
                        cleaned_lineage.append(item)
                    if 'unclassified' in item and level < 2:
                        unclassified = True
                    if 'environmental' in item and level < 2:
                        unclassified = True
                    if 'uncultured' in item and level < 2:
                        unclassified = True
                    if 'metagenome' in item and level < 2:
                        unclassified = True
                    if 'unknown' in item and level < 2:
                        unclassified = True
                    if level == 3:
                        break
                    level += 1

                if unclassified is None:
                    all_lineages.append(cleaned_lineage)

    highest_count = 0
    life_dom = ""
    for item in life_domain_count:
        if life_domain_count[item] > highest_count:
            highest_count = life_domain_count[item]
            life_dom = item

    all_to_study = all_lineages
    for item in all_lineages:
        if item[0] != life_dom:
            all_to_study.remove(item)

    elements_in_all = []

    if len(all_to_study) > 0:
        shared_elements = set.intersection(*map(set, all_to_study))
        for element in all_to_study[0]:
            if element in shared_elements and\
                    len(element) > 0:
                elements_in_all.append(element)

        elements_in_all = list(elements_in_all)

    if len(elements_in_all) > 0:

        first = elements_in_all[0]
        if len(elements_in_all) > 1:
            second = elements_in_all[1]
        else:
            second = "Unknown"
        if len(elements_in_all) > 2:
            third = elements_in_all[2]
        else:
            third = "Unknown"
        return first, second, third

    else:
        if len(life_dom) > 0:
            first = life_dom
        else:
            first = "Unknown"
        second = 'Unknown'

        upper_taxes = []
        for lineages in all_to_study:
            if len(lineages) > 2:
                up_tax = ';'.join(lineages[0:3])
                if up_tax not in upper_taxes:
                    upper_taxes.append(up_tax)
        third = ', '.join(upper_taxes)

        return first, second, third


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

    smallerThan = lambda x, y: [i for i in x if i < y]

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

            new_align_dict[seq] =\
                new_align_dict[seq][0: pos:] + new_align_dict[seq][pos + 1::]

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

            if a != b and b != '-':

                mismatches += 1

            if b != '-':

                length += 1

        return mismatches, length

    for a, b in zip(seq1, seq2):

        if a != b:

            mismatches += 1

        if a != '-' or b != '-':

            length += 1

    return mismatches, length


def tree_params(tree_log):
    """
    Takes a dictionary of sequences and returns frequencies
    for each base as a dictionary. If a base in the sequence
    is uncertain. The script raffles off a base.

    Arguments
    ---------
    :seqs: A sequence dictionary

    Returns:
    --------
    Freqs: Frequencies for base frequecies
    """

    with open(tree_log, 'r') as f:

        rate_parameters = {}
        base_frequencies = {}
        gamma_shape_alpha = ""

        for line in f:

            if 'Rate parameters:' in line:
                line_splitted = line.rstrip().split()
                print(line_splitted)
                rate_parameters["AC"] = Decimal(line_splitted[3])
                rate_parameters["AG"] = Decimal(line_splitted[5])
                rate_parameters["AT"] = Decimal(line_splitted[7])
                rate_parameters["CG"] = Decimal(line_splitted[9])
                rate_parameters["CT"] = Decimal(line_splitted[11])
                rate_parameters["GT"] = Decimal(line_splitted[13])

            elif 'Base frequencies:' in line:
                line_splitted = line.rstrip().split()
                base_frequencies["A"] = Decimal(line_splitted[3])
                base_frequencies["C"] = Decimal(line_splitted[5])
                base_frequencies["G"] = Decimal(line_splitted[7])
                base_frequencies["T"] = Decimal(line_splitted[9])

            elif 'Gamma shape alpha:' in line:
                line_splitted = line.rstrip().split()
                gamma_shape_alpha = Decimal(line_splitted[3])

    return rate_parameters, base_frequencies, gamma_shape_alpha


def define_insertion_and_deletion_rates(tree,
                                        fasta,
                                        path_to_script):
    """
    Runs a dawg script to identify the insertion and
    deletion rates. The script parses the output and
    returns suitable models, and insertion and deletion
    rates.

    Arguments
    ---------
    tree: A phylogenetic tree file
    fasta: A fasta file to study
    path_to_script: A path to the lambda.pl script (dawg package)

    """

    output = check_output(
            ["perl",
             "{path_to_script}lambda.pl".format(
                 path_to_script=path_to_script),
             "{tree_file}".format(tree_file=tree),
             "{fasta_file}".format(fasta_file=fasta)],
            encoding='utf-8')

    output_split = output.split("Most Likely Model\n")[-1]
    output_split = output_split.split()
    GapModel = output_split[2].rstrip('"').lstrip('"')
    if GapModel == "PL":
        GapModel = "POWER-LAW"
    elif GapModel == "NB":
        GapModel = "GEO"
    GapRate = output_split[5].lstrip('{').rstrip(', ')
    SeqLen = output_split[-1].rstrip('}')
    gapRate = output.split("Lambda Estimate is ")[-1].split('.\n')[0]

    return GapModel, GapRate, SeqLen, gapRate


def simulate_sequences(base_frequencies,
                       substitution_rate,
                       gamma_shape_alpha,
                       seq_len,
                       tree,
                       identifier,
                       parental_sequence,
                       output_dir,
                       gapModel=None,
                       gapRate=None,
                       ins_rate=None,
                       del_rate=None,
                       max_ins=None,
                       max_del=None):
    """
    Simulates sequences under a given substitution
    model.

    First, it takes a substitution model under, a subtree (with
    a parent and two_children) and a parental sequence and creates
    a .dawg file. Based on this information the script simulates
    two child sequences.

    Arguments
    ---------
    :Base frequencies: Base frequencies as dictionary
    :Substitution_rates: A dictionary
    :Gamma shape alpha: Gamma shape alpha number
    :seq_len: Alignment length
    :tree: A tree in newick formar
    :identifier: Identifier for a cluster
    :parental_sequence: A root sequence for a simulation
    :outpudir: A directory, in which files are written

    Output
    -------

    Generates a .dawg file and runs dawg

    """

    with open(output_dir + identifier + '.dawg', 'w') as f:

        f.write("Output.markins = 1" + '\n')
        f.write("Output.keepempty = 0" + '\n')
        f.write("Root.Length = " + str(
            len(parental_sequence.replace('-', ''))) + '\n')
        f.write("Root.Code = 1" + '\n')
        f.write("Root.Seq = " + "\"" +
                parental_sequence.replace('-', '') + "\"" + '\n')
        f.write("Sim.Reps = 1\n")
        f.write("Sim.Seed = 42\n")

        f.write('\n')
        f.write("[Indel]" + '\n')
        f.write("Params.Ins = " + str(gapRate) + ", " + seq_len + "\n")
        f.write("Model.Ins = \"" + gapModel + "\"\n")
        f.write("Rate.Ins = " + str(ins_rate) + '\n')
        f.write("Params.Del = " + str(gapRate) + ", " + seq_len + "\n")
        f.write("Model.Del = \"" + gapModel + "\"\n")
        f.write("Rate.Del = " + str(del_rate) + "\n")

        f.write("[[-]]" + '\n')
        f.write("Tree.Tree = " + "\"" + tree + "\"" + '\n')
        f.write("Subst.Model = \"gtr\"\n")
        f.write("Subst.Freqs  = ")
        for a in ["A", "C", "G", "T"]:
            if a != "T":
                f. write(str(round(base_frequencies[a], 3)) + ', ')
            else:
                f. write(str(round(base_frequencies[a], 3)))

        f.write('\n')
        f.write("Subst.params = ")
        for a in ["AC", "AG", "AT", "CG", "CT", "GT"]:
            if a != "GT":
                f.write(str(round(substitution_rate[a], 3)) + ', ')
            else:
                f.write(str(round(substitution_rate[a], 3)))
        f.write('\n')
        f.write("Subst.Shape.Model = \"GAMMA\"\n")
        f.write("Subst.Shape.Params = " + str(gamma_shape_alpha) + '\n')

    p = check_output(
            ["dawg",
             "{output_dir}{identifier}.dawg".format(output_dir=output_dir,
                                                    identifier=identifier),
             "-o",
             "{output_dir}{identifier}.fasta".format(output_dir=output_dir,
                                                     identifier=identifier)])


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
                     "CM_and_inverted_loop", "inverted_loop",
                     "stem_extended"]:

        query_dict = {}
        query_dict["ts_quality"] = True
        query_dict["structural_quality"] = True
        query_dict[criteria] = True
        if criteria != "stem_extended":
            query_dict["ts_mismatch"] = {"$gt": 1}
        query_dict["ts_distance"] = {"$lt": 9}
        query_list.append(query_dict)

    query = {}

    if len(alternative_query) != 0:
        for field, values in alternative_query.items():
            query[field] = values

    tsm_collection = defaultdict(dict)
    tsm_substitution = defaultdict(dict)
    trees_dict = {}

    for query in query_list:
        cursor = collection.find(query)

        for document in cursor:

            tree_file = document["tree_file"].split('//')
            tree_file = '/'.join(tree_file)
            case = ""
            for x in query:
                if x in ["multiple_CM", "CM_and_loop_change",
                         "CM_and_inverted_loop", "inverted_loop",
                         "stem_extended"]:
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
    leaf_count = 0

    levels = defaultdict(set)
    for lev in range(0, level):
        levels[lev] = set()

    for treefile in trees_dict:
        tree = PhyloTree(treefile, format=1)

        for node in tree.traverse("levelorder"):
            if node.is_leaf():
                leaf_count += 1

                node_name = node.name

                target = node
                levels[0].add(treefile + '=' + target.name)
                for i in range(1, level):
                    if not target.is_root():
                        target = target.up
                        levels[i].add(treefile + '=' + target.name)
                    else:
                        continue

    tsm = defaultdict(dict)
    tsm_subs = defaultdict(dict)
    tar_nodes = defaultdict(dict)

    for case in ["multiple_CM", "CM_and_loop_change",
                 "CM_and_inverted_loop", "inverted_loop",
                 "stem_extended"]:

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
                         "CM_and_inverted_loop", "inverted_loop",
                         "stem_extended"]:

                if treex in tsm_collection[case]:
                    if node_name in tsm_collection[case][treex]:
                        tsm[case][lev] +=\
                            tsm_collection[case][treex][node_name]
                        tsm_subs[case][lev] +=\
                            tsm_substitution[case][treex][node_name]

                        tar_nodes[case][lev].append(
                                tsm_collection[case][treex][node_name])

                    else:
                        tar_nodes[case][lev].append(0)
                else:
                    tar_nodes[case][lev].append(0)

    for case in ["multiple_CM", "CM_and_loop_change",
                 "CM_and_inverted_loop", "inverted_loop",
                 "stem_extended"]:
        for lev in [0, 1, 2]:
            rate3 = round(
                    (Decimal(tsm[case][lev])/Decimal(len(levels[lev]))), 10)
            std = round(numpy.std(tar_nodes[case][lev]), 10)
            mean = round(numpy.mean(tar_nodes[case][lev]),  10)

            yield (len(trees_dict), case, lev, tsm[case][lev],
                   tsm_subs[case][lev], rate3, mean, std)


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

    num_set = set(
        [qry_ts_start, qry_ts_end, ref_ts_start, ref_ts_end])

    if min(num_set) < minimum:
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

        elif minimum <= 20 and (len(align_dict[align]) - maximum) <= 20:
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
                cut_start = minimum - 20
                cut_end = len(align_dict[align]) - 1

    return minimum, maximum, cut_start, cut_end, cut_align_dict
