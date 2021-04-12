"""
This script is intended for adding information 
about the identified TSMs to a mongo database.
"""

from mongoengine import (Document, StringField,
                         BooleanField, IntField)


class Results_RNA_ts(Document):
    """
    A Document for an information of the identified TSM. 
    """

    identifier = StringField(required=True)
    cluster = StringField(required=True)
    position = IntField()
    ts_distance = IntField()
    ts_quality = BooleanField()
    ts_length = IntField()
    ts_mismatch = IntField()
    structural_quality = BooleanField()
    loop_seq = StringField()
    loop_seq_ref = StringField()
    adj_nucl1 = StringField()
    adj_nucl2 = StringField()
    adj_nucl1_ref = StringField()
    adj_nucl2_ref = StringField()
    CM_present = BooleanField()
    one_CM = BooleanField()
    multiple_CM = BooleanField()
    CM_and_loop_change = BooleanField()
    CM_and_inverted_loop = BooleanField()
    inverted_loop = BooleanField()
    new_loop = BooleanField()
    asymmetric_TSM_between_stems = BooleanField()
    TSM_loop_to_stem = BooleanField()
    TSM_within_one_stem = BooleanField()
    loop_to_seq = IntField()
    tsm_within_stem = IntField()
    major_change = BooleanField()
    no_change = BooleanField()
    loop_extension = BooleanField()
    mismatch_fixed = BooleanField()
    leaf = BooleanField()
    new_mismatch = BooleanField()
    new_loop = BooleanField()
    full_loop_inversion = BooleanField()
    compensating_mutation = IntField()
    reverse_complement = IntField()
    query_node = StringField(required=True)
    ref_node = StringField(required=True)
    sister_node = StringField(required=True)
    ancestor_node = StringField(required=True)
    child_node1 = StringField()
    child_node2 = StringField()
    alignment_file = StringField(required=True)
    qry_ts_start = IntField(required=True)
    qry_ts_end = IntField(required=True)
    ref_ts_start = IntField(required=True)
    ref_ts_end = IntField(required=True)
    structure_ref = StringField(required=True)
    structure_qry = StringField(required=True)
    structure_sister = StringField(required=True)
    structure_child1 = StringField()
    structure_child2 = StringField()
    tree_file = StringField()


    meta = {
        'db_alias': 'default',
        'indexes': [
            'identifier',
            '$identifier',  # text index
            '#identifier',  # hashed index
        ]
    }


def add_information(result_dict):
    """
    Adds information about an identified TSM to mongo database.

    Requires:
    result_dict: A dictionary, which contains all information
    on the identified TSM.

    """

    result_document = Results_RNA_ts(
            identifier=result_dict["identifier"],
            cluster=result_dict["cluster"],
            ts_distance = result_dict["ts_distance"],
            position=result_dict["position"],
            ts_quality=result_dict["ts_quality"],
            ts_length=result_dict["ts_length"],
            ts_mismatch=result_dict["ts_mismatch"],
            loop_seq = result_dict["loop_seq"],
            loop_seq_ref = result_dict["loop_seq_ref"],
            structural_quality=result_dict["structural_quality"],
            CM_present=result_dict["CM_present"],
            one_CM=result_dict["one_CM"],
            multiple_CM=result_dict["multiple_CM"],
            CM_and_loop_change=result_dict["CM_and_loop_change"],
            CM_and_inverted_loop=result_dict["CM_and_inverted_loop"],
            inverted_loop=result_dict["inverted_loop"],
            new_loop=result_dict["new_loop"],
            asymmetric_TSM_between_stems=result_dict["asymmetric_TSM_between_stems"],
            TSM_loop_to_stem=result_dict["TSM_loop_to_stem"],
            TSM_within_one_stem=result_dict["TSM_within_one_stem"],
            loop_to_seq = result_dict["loop_to_seq"],
            adj_nucl1 = result_dict["adj_nucl1"],
            adj_nucl2 = result_dict["adj_nucl2"],
            adj_nucl1_ref = result_dict["adj_nucl1_ref"],
            adj_nucl2_ref = result_dict["adj_nucl2_ref"],
            major_change=result_dict["major_change"],
            no_change=result_dict["no_change"],
            loop_extension=result_dict["loop_extension"],
            mismatch_fixed=result_dict["mismatch_fixed"],
            leaf=result_dict["leaf"],
            new_mismatch=result_dict["new_mismatch"],
            compensating_mutation=result_dict["compensating_mutations"],
            query_node = result_dict["query_node"],
            ref_node = result_dict["ref_node"],
            sister_node = result_dict["sister_node"],
            ancestor_node = result_dict["ancestor_node"],
            child_node1 = result_dict["child_node1"],
            child_node2 = result_dict["child_node2"],
            alignment_file = result_dict["alignment_file"],
            tree_file = result_dict["tree_file"],
            qry_ts_start = result_dict["qry_ts_start"],
            qry_ts_end = result_dict["qry_ts_end"],
            ref_ts_start = result_dict["ref_ts_start"], 
            ref_ts_end = result_dict["ref_ts_end"],
            structure_ref = result_dict["structureRef"],
            structure_qry = result_dict["structureQry"],
            structure_sister = result_dict["structureSister"],
            structure_child1 = result_dict["structureFileChild1"],
            structure_child2 = result_dict["structureFileChild2"])

    result_document.save()


