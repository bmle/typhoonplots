#!/usr/bin/env python

"""
typhoonPlot_multiWindows.py
Author: bmle

Makes typhoon plots from sequence alignment data, based on how many windows are specified.
Based off an R script that also creates typhoon plots.

Adjust the parameters in the execution section, then run the script to produce
one or more files containing typhoon plots at the specified genomic windows.

"""

import argparse
import numpy as np
import pysam
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import math


def main(chr_name, origin, extent, limit, out, file):
	# noinspection PyUnresolvedReferences
	f = pysam.AlignmentFile(file, "rb")  # Load data file

	# Calculate genomic regions to search in, each region having a width of "extent"
	plots = []
	start = origin
	while start < limit:
		end = start + extent
		plots.append((start, end))
		start = end
	num_plots = len(plots)

	max_len = 250	# Maximum length of read to capture (250 should capture everything)

	# Make the color scheme for the plot
	colordict = {"red": ((0, 1, 1),
	                     (1, 0, 0)),
	             "green": ((0, 1, 1),
	                       (1, 0, 0)),
	             "blue": ((0, 1, 1),
	                      (1, 1, 1))}
	my_map = LinearSegmentedColormap("mymap", colordict)

	# =========================================================================

	# Splits plots into multiple files (because of memory issues when generating large files)
	n = 30  # Number of typhoon plots to be included in each file
	plots = [plots[i*n : (i+1)*n] for i in range(math.ceil(len(plots) / n))]

	# Loop over all files
	for x, plot in enumerate(plots, start=1):

		# =========================================================================
		# Calculates read overlaps for each genomic window, using each read's
		# template length
		# =========================================================================
		print("Calculating overlaps for file " + str(x) + " of " + str(len(plots)) + "...")

		data = []  # List to store overlaps

		# Loop over all regions in each file
		for y, region in enumerate(plot, start=1):
			print("\tCalculating overlaps for region " + str((n*(x-1))+y) + " of " + str(num_plots) + "...")

			# Create temporary array to store alignment data
			temp = np.zeros((max_len + 1, region[1] - region[0] + 1))

			# Filter for region of interest
			reg_reads = f.fetch(chr_name, region[0], region[1])

			# Loop over all reads
			for read in reg_reads:

				# Calculates read width based on template length
				r_len = read.template_length

				# Shrinks template length (was present in the original R script)
				# Appears to make the bars in the plot shorter in width
				offset_start = round(r_len/5)   # default: 0
				offset_end = offset_start * 3   # default: 0

				# Calculate proper start and stop positions
				r_start = read.reference_start + offset_start
				r_end = min(r_start + offset_end, region[1])  # min() prevents reads from extending past the region of interest

				# Adds read positions to temporary array
				for z in range(r_start, r_end + 1):
					cur_val = temp[r_len][z - region[0]]
					temp[r_len][z - region[0]] = min(cur_val + 1, 30)   # Cap all values at 30 so plot is displayed correctly

			# Append temp array to the list of all overlaps
			data.append(temp)

		# =========================================================================
		# Makes a "typhoon" plot from alignment data
		# =========================================================================
		print("Generating plots for file " + str(x) + "...")

		# Make the heatmap
		fig = plt.figure(figsize=(15, 6 * len(data)))

		# Add each data set as a subplot
		for y, (datum, region) in enumerate(zip(data, plot), start=1):
			ax = fig.add_subplot(len(data), 1, y)
			ax.imshow(datum, cmap=my_map, aspect="auto", extent=[region[0], region[1], max_len, 0])

			# Flip plot vertically so that higher read lengths are towards the top
			ax.invert_yaxis()

			# Set labels
			ax.set_ylabel("Fragment Length (bp)")
			ax.set_xlabel(" Position from " + chr_name + ":" + str(region[0]))

		# Save the plot
		plt.tight_layout()
		plt.savefig(out + "_" + str(x) + ".png")
		plt.close()
		print("Plot " + str(x) + " of " + str(len(plots)) + " created!")


# === Execution ===============================================================

if __name__ == "__main__":

	p = argparse.ArgumentParser(description="Create multiple typhoon plots from a single read file")
	p.add_argument("chr_name",				help="Name of the chromosome of interest")
	p.add_argument("origin",	type=int,	help="Position to start plotting from")
	p.add_argument("extent",	type=int,	help="How large each window should be")
	p.add_argument("limit",		type=int,	help="Furthest search extent wanted (or max size of chromosome)")
	p.add_argument("out",					help="Filename of outputted plot")
	p.add_argument("file",					help="BAM/SAM read file to make plots out of")
	args = p.parse_args()

	main(vars(args))
