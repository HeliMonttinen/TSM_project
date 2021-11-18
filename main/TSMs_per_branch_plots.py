"""
A tool for drawing a barplot about two-step
CMs.

Author: Heli Mönttinen (ORCID: 0000-0003-2461-0690)

"""
import os
import sys
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    This script goes through trees in a directory,
    counts tsms/branch and creates a bar graph.

    :input: path to an input file containing information
            on the TSMs
    :output: .svg output file
    """

    from tree_parse import events_per_substitution

    output = sys.argv[1]
    database = sys.argv[2]
    collection = sys.argv[3]

    classes_of_interest = ["multiple_CM",
                           "CM_and_loop_change",
                           "CM_and_inverted_loop",
                           "inverted_loop",
                           "stem_extended"]

    ind = np.arange(len(classes_of_interest))  # the x locations for the groups
    width = 0.27
    datadict = defaultdict(list)

    for class_x in classes_of_interest:

        alternative = {"ts_quality": True,
                       "structural_quality": True,
                       class_x: True,
                       "ts_distance": {"$lt": 9}}

    for (trees_dict, case, level, tsm,
         tsm_subs, rate3, mean, std) in events_per_substitution(
            database,
            collection,
            level=3,
            alternative=alternative):

        datadict[level].append(rate3)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    times = 0

    colors = ['#f1cf07', '#92eb66', '#b688ff']

    rects_list = defaultdict(list)

    for level in datadict:

        rects1 = ax.bar(ind+width*times,
                        datadict[level],
                        width,
                        color=colors[times])

        rects_list[times].extend(rects1)
        times += 1

    ax.set_ylabel('TSMs per branch')
    ax.legend((rects_list[0][0],
               rects_list[1][0],
               rects_list[2][0]),
              ('Terminal node', '-1', '-2'))

    ax.set_xticks(ind+width)

    ax.set_xticklabels([
        'Multiple CMs',
        'CM +\nparallel\nmut. in loop',
        'CM +\nloop\ninversion',
        'Loop\ninversion\nonly',
        'Insertion\nin stem half'])

    plt.savefig(output, bbox_inches='tight')


if __name__ == "__main__":
    main()
