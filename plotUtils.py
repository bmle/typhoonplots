"""
plotUtils.py
Author: bmle

Contains common functions used by the plot-generating scripts in this repository.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches

def calculate_reads(reads, start, end, max_len):
    """
    Takes a pysam.AlignedSegment object (which has read information), coverts the reads to template length,
    and builds an array of read values.
    :param reads: A pysam.AlignedSegment object containing the reads to convert
    :param start: Start position to capture
    :param end: End position capture
    :param max_len: Maximum length of read to capture
    :return: An array of read values
    """

    # Create temporary array to store alignment data
    temp = np.zeros((max_len + 1, end - start + 1))

    # Loop over all reads in file
    for r in reads:

        # Calculates read width based on template length
        r_len = r.template_length

        # Shrinks template length (was present in the original R script)
        # Appears to make the bars in the plot shorter in width
        offset_start = round(r_len / 5)  # default: 0
        offset_end = offset_start * 3  # default: 0

        # Calculate adjusted start and stop positions
        r_start = r.reference_start + offset_start
        r_end = min(r_start + offset_end, end)  # min() prevents reads from extending past the region of interest

        # Adds read positions to temporary array
        for i in range(r_start, r_end + 1):
            temp[r_len][i - start] += 1

    temp[temp > 15] = 15  # Cap all values at 15 so the plot has higher contrast
    return temp


def cmap():
    """
    Generates a color map to be used in the typhoon plots
    :return: A LinearSegmentedColormap that defines a white-to-blue color gradient.
    """

    colordict = {"red":   ((0, 1, 1),
                           (1, 0, 0)),
                 "green": ((0, 1, 1),
                           (1, 0, 0)),
                 "blue":  ((0, 1, 1),
                           (1, 1, 1))}
    return LinearSegmentedColormap("mymap", colordict)


def make_gene_bodies(ax, file ="./genebodies.txt"):
    """
    Annotates genes on the typhoon plot
    :param ax: The Axes object to make the plot on
    :param file: A comma-delimited file specifying gene bodies
    :return: A list of matplotlib.patches.Rectangle representing gene bodies
    """

    # Parses gene body file
    genes = [line.strip() for line in open(file)]
    for gene in genes:

        # Check if line is a comment (starts with #)
        if gene[0] == "#": continue

        # Extract name, start, and end of gene
        temp = gene.split(",")
        name = temp[0]
        left = int(temp[2])
        length = int(temp[3]) - int(temp[2])
        right = left + length

        # Draws gene body
        p = patches.Rectangle((left, 15), length, 12, alpha=0.6, facecolor="0.6", edgecolor="none")
        ax.add_patch(p)

        # Adds text label
        # Problem: text will remain centered on gene body, even if outside of drawing canvas
        plt.text(0.5*(left+right), 21, name, ha="center", va="center")

    return ax
