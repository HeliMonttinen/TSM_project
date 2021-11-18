"""
Counts frequencies of different mutation types
in empirical and simulated data, and tests if
a mutation type is significantly (using t-test)
more frequent in empirical data than in simulated data.
"""
from collections import defaultdict
from decimal import Decimal
import os
import pandas as pd
from pymongo import MongoClient
import researchpy as rp
import sys

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Arguments
    =========

    directory: A path to drectory in which tree files
               are located. Should be formatted as
               ".anctree".
    orig_db: A name of the database containing empirical
             data results
    simulation_db: A  name of the database containing
                   simulated data results
    mode: leaf or nonleaf. Leaf means that mutations in
          terminal nodes are studied
    regions_file_orig: A file containing high_quality regions
                       for the original data
    regions_file: A file containing high quality regions
                  for simulated sequences

    Output
    ======
    prints out correction factors and results of t_test
    for each mutation type
    """

    from common import return_file
    from tree_parse import read_tree

    directory = sys.argv[1]
    orig_db = sys.argv[2]
    simulation_db = sys.argv[3]
    mode = sys.argv[4]
    regions_file_orig = sys.argv[5]
    regions_file = sys.argv[6]

    regions_dict_orig = {}
    regions_dict = {}

    with open(regions_file_orig) as f:
        for line in f:
            line_splitted = line.rstrip().split('\t')
            identifier = line_splitted[0].split('_')[0]
            parent = line_splitted[1].rstrip()
            length = line_splitted[4]
            kid1 = line_splitted[2]
            kid2 = line_splitted[3]
            if mode == "leaf":
                if '#' in kid1 and '#' in kid2:
                    continue
            else:
                if ('#' not in kid1) and ('#' not in kid2):
                    continue

            if identifier not in regions_dict_orig:
                regions_dict_orig[identifier] = dict()
            if parent not in regions_dict_orig[identifier]:
                regions_dict_orig[identifier][parent] = dict()
                regions_dict_orig[identifier][parent] = int(length)

    with open(regions_file) as f:

        for line in f:
            line_splitted = line.rstrip().split('\t')
            identifier = line_splitted[0].split('_')[0]
            parent = line_splitted[1].rstrip()

            kid1 = line_splitted[2]
            kid2 = line_splitted[3]
            if mode == "leaf":
                if '#' in kid1 and '#' in kid2:
                    continue
            else:
                if ('#' not in kid1) and ('#' not in kid2):
                    continue

            length = line_splitted[4]
            if identifier not in regions_dict:
                regions_dict[identifier] = dict()
            if parent not in regions_dict[identifier]:
                regions_dict[identifier][parent] = dict()
                regions_dict[identifier][parent] = int(length)

    node_length = 0
    node_length_sim = 0

    for ident in regions_dict:
        for item in regions_dict[ident]:
            length = regions_dict[ident][item]
            node_length += length

            if ident not in regions_dict_orig:
                continue
            elif item not in regions_dict_orig[ident]:
                continue
            length = regions_dict_orig[ident][item]
            node_length_sim += length

    for ident in regions_dict_orig:
        for item in regions_dict_orig[ident]:
            if ident not in regions_dict:
                length = regions_dict_orig[ident][item]
                node_length_sim += length
            elif item not in regions_dict[ident]:
                length = regions_dict_orig[ident][item]
                node_length_sim += length

    factor = Decimal(node_length/node_length_sim)

    categories = ["CM_present", "one_CM", "multiple_CM",
                  "CM_and_loop_change", "CM_and_inverted_loop",
                  "inverted_loop", "asymmetric_TSM_between_stems",
                  "TSM_loop_to_stem", "TSM_within_one_stem",
                  "stem_extended"]

    sim_categories = defaultdict(list)
    orig_categories = defaultdict(list)

    for path in return_file(directory):

        if not path.endswith(".anctree"):

            continue

        tree = read_tree(path)
        node_length = 0
        node_length_sim = 0
        identifier = path.split('/')[-1].split('_')[0]
        if identifier not in regions_dict:
            continue

        for node in tree.traverse("postorder"):

            if node.is_root():
                continue
            parent = (node.up).name
            if parent not in regions_dict[identifier]:
                continue
            if mode == "leaf":
                if not node.is_leaf():
                    continue
            else:
                if node.is_leaf():
                    continue

        for database in [orig_db, simulation_db]:

            if mode == "leaf":
                l_mode = True
            else:
                l_mode = False

            query = {"cluster": identifier,
                     "ts_quality": True,
                     "structural_quality": True,
                     "leaf": l_mode,
                     "ts_distance": {"$lt": 9},
                     "ts_length": {"$gt": 9}}

            client = MongoClient()
            db = client[database]
            collection = db['results__r_n_a_ts']

            categories = {"CM_present": 0, "one_CM": 0,
                          "multiple_CM": 0,
                          "CM_and_loop_change": 0,
                          "CM_and_inverted_loop": 0,
                          "inverted_loop": 0,
                          "asymmetric_TSM_between_stems": 0,
                          "TSM_loop_to_stem": 0,
                          "TSM_within_one_stem": 0,
                          "stem_extended": 0}

            for document in collection.find(query):

                hit_dict = {}

                for key, value in document.items():

                    hit_dict[key] = value

                if hit_dict["ref_node"] not in regions_dict[identifier]:
                    continue

                for category in categories:

                    if hit_dict[category] is True:
                        if category in ["inverted_loop"]:
                            if hit_dict["ts_mismatch"] > 1:
                                categories[category] += 1
                        else:
                            categories[category] += 1

            for category in categories:

                if database == orig_db:

                    prop = Decimal(categories[category])
                    orig_categories[category].append(float(prop))

                elif database == simulation_db:

                    prop = Decimal(categories[category])*factor
                    sim_categories[category].append(float(prop))

    for category in categories:

        print(category)
        orig_data = orig_categories[category]

        sim_data = sim_categories[category]

        d = {'orig_data': orig_data, 'sim_data': sim_data}
        df = pd.DataFrame(data=d)

        des, res = rp.ttest(group1=df['orig_data'],
                            group1_name="orig_data",
                            group2=df['sim_data'],
                            group2_name="sim_data")
        print(des)
        print(res)


if __name__ == "__main__":
    main()
