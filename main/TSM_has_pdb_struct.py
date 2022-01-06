"""
Compares BLAST results and identified CMs to identify, if
a sequence with a CM has a PDB structure 
(with >95% sequence similarity).
"""

from pymongo import MongoClient
import os
import sys
from Bio import SeqIO
from collections import OrderedDict

dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():

    
    blastres_file = sys.argv[1]
    database = sys.argv[2]
    collection = sys.argv[3]
    client = MongoClient()
    db = client[database]
    collection = db[collection]

    id_set = {}
    with open(blastres_file, 'r') as f:

        for line in f:

            identifier = line.split('\t')[0]
            identifier = identifier.replace('/', '_')

            id_set[identifier] = line.split('\t')[1]

        for identifier in id_set:

            query = {"query_node": identifier,
                     "ts_quality": True,
                     "structural_quality": True,
                     "CM_present": True}
            for document in collection.find(query):


                hit_dict = {}

                for key, value in document.items():

                    hit_dict[key] = value


                print(identifier, id_set[identifier], hit_dict["identifier"],\
                        hit_dict["qry_ts_start"], hit_dict["qry_ts_end"])

            query = {"child_node2": identifier,
                     "ts_quality": True,
                     "structural_quality": True,
                     "CM_present": True}

            for document in collection.find(query):


                hit_dict = {}

                for key, value in document.items():

                    hit_dict[key] = value


                print(identifier, id_set[identifier], hit_dict["identifier"],\
                        hit_dict["qry_ts_start"], hit_dict["qry_ts_end"])


            query = {"child_node1": identifier,
                     "ts_quality": True,
                     "structural_quality": True,
                     "CM_present": True}
            for document in collection.find(query):


                hit_dict = {}

                for key, value in document.items():

                    hit_dict[key] = value


                print(identifier, id_set[identifier], hit_dict["identifier"],\
                        hit_dict["qry_ts_start"], hit_dict["qry_ts_end"])
if __name__ == "__main__":
    main()

