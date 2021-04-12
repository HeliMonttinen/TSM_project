"""
This tool draws a heatmap from CM-linked base-pairs.

The script is run: 
    python3 heatmaps_instantaneous_CMs.py input output

    input: A CSV formatted result file
            (stem-bp-in-parent,stem-bp-in-child  count)
    output: an output file 


Author: Heli Mönttinen  (ORCID: 0000-0003-2461-0690)
"""


import os
import sys
import numpy as np
from collections import OrderedDict
import matplotlib
import matplotlib.pyplot as plt
import math
import seaborn as sns


dirname = os.path.dirname(__file__)
grand_parent_dir = os.path.dirname(dirname)
filen = os.path.join(grand_parent_dir)
sys.path.append(filen)


def main():
    """
    Reads data for heatmaps.
    The value pairs in the input file should be formatted as
    x,y  number

    """

    input_f = sys.argv[1]
    output = sys.argv[2]

    labels = set()
    input_data = {}

    with open(input_f, 'r') as f:

        for line in f:
            line_split = line.rstrip().split()
            count = line_split[1]
            y_val = line_split[0].split(',')[1]
            x_val = line_split[0].split(',')[0]

            labels.add(str(y_val))
            labels.add(str(x_val))

            input_data[line_split[0]] = line_split[1]

    labels = list(labels)
    labels.sort()
    np_list = np.array([])
    res_matrix_list = []
    res_matrix_list_log = []

    count= 0
    for x in labels:
        list_np = np.array([])
        list_np_log = np.array([])

        for y in labels:

            if x + ',' + y in input_data:

                list_np = np.append(list_np, int(input_data[x + ',' + y]))
                list_np_log = np.append(list_np_log, math.log(int(input_data[x + ',' + y])))

            else:
                list_np = np.append(list_np, np.nan)
                list_np_log = np.append(list_np_log, np.nan)

        res_matrix_list.append(list_np)
        res_matrix_list_log.append(list_np_log)
        count += 1

    res_matrix = np.array(res_matrix_list)
    res_matrix_log = np.ma.masked_array(res_matrix_list_log, mask=np.isnan(res_matrix_list_log))

    fig, ax = plt.subplots()
    im = ax.imshow(res_matrix_log)

    ax.set_ylim(sorted(ax.get_xlim(), reverse=True))

    ax.set_xticks(np.arange(res_matrix_log.shape[1]), minor=False)
    ax.set_yticks(np.arange(res_matrix_log.shape[0]), minor=False)
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("log10", rotation=-90, va="bottom")

    for i in range(len(labels)):
        for j in range(len(labels)):

            try:
                text = ax.text(j, i, int(res_matrix[i, j]),
                        ha="center", va="center", color="w")
            except:
                text = ax.text(j, i, str("NaN"),
                        ha="center", va="center", color="w")

    plt.savefig(output, format="svg")
    plt.show()

if __name__ == "__main__":
    main()
                

                

            





