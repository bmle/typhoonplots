#!/usr/bin/env python

"""
typhoonPlot.py
Author: bmle

Makes typhoon plots from sequence alignment data. Based off an R script that also creates typhoon plots.

Produces a single file containing one or more typhoon plots at the specified genomic coordinates.

"""

import argparse
import numpy as np
import pysam
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def main(chr_name, origin, offset, out, files, condensed = False):

	# =========================================================================
	# Calculate read overlaps for each genomic position, using each read's
	# template length
	# =========================================================================
	print("Calculating overlaps...")

	# Make matrix to store overlaps, an entry for each file,
	# and a row for each read length and a column for each nucleosomal position
	start = origin - offset
	end = origin + offset
	data = []

	max_len = 250	# Maximum length of read to capture (250 should capture everything)

	# Iterate over all data files
	for file in files:

		# Load each file and filter for region of interest
		# noinspection PyUnresolvedReferences
		tempf = pysam.AlignmentFile(file, "rb")
		f = tempf.fetch(chr_name, start, end)

		# Create temporary array to store alignment data
		temp = np.zeros((max_len + 1, end - start + 1))

		# Loop over all reads in file
		for r in f:

			# Calculates read width based on template length
			r_len = r.template_length

			# Shrinks template length (was present in the original R script)
			# Appears to make the bars in the plot shorter in width
			offset_start = round(r_len/5)   # default: 0
			offset_end = offset_start * 3   # default: 0

			# Calculate adjusted start and stop positions
			r_start = r.reference_start + offset_start
			r_end = min(r_start + offset_end, end)  # min() prevents reads from extending past the region of interest

			# Adds read positions to temporary array
			for i in range(r_start, r_end + 1):
				# temp[r_len][i - start] += 1
				cur_val = temp[r_len][i - start]
				temp[r_len][i - start] = min(cur_val + 1, 15)   # Cap all values at 15 so plot has higher contrast

		# Append temp array to the list of all overlaps
		data.append(temp)

	print("Finished calculating overlaps!")


	# =========================================================================
	# Make "typhoon" plot(s) from alignment data
	# =========================================================================
	print("Generating plot...")

	# Make the color scheme
	colordict = {"red":   ((0, 1, 1),
	                       (1, 0, 0)),
	             "green": ((0, 1, 1),
	                       (1, 0, 0)),
	             "blue":  ((0, 1, 1),
	                       (1, 1, 1))}
	my_map = LinearSegmentedColormap("mymap", colordict)

	# Make the heatmap
	h = 2 if condensed else 6
	fig = plt.figure(figsize=(15, len(data) * h))
	if condensed: fig.suptitle(out.split("/")[-1], size=16)

	# Add each data set as a subplot
	for i, (datum, f_name) in enumerate(zip(data, files), start=1):
		ax = fig.add_subplot(len(data), 1, i)
		ax.imshow(datum, cmap=my_map, aspect="auto", extent=[origin - offset, origin + offset, max_len, 0])

		# Flip plot vertically so that higher read lengths are towards the top
		ax.invert_yaxis()

		# Set labels
		if not condensed:
			ax.set_title(f_name.split("/")[-1])
			ax.set_ylabel("Fragment Length (bp)")
			ax.set_xlabel(" Position on chromosome " + chr_name + " (center: " + str(origin) + ")")

	# Save the plot
	r = [0, 0, 1, 0.97] if condensed else [0, 0, 1, 1]
	plt.tight_layout(rect=r)
	plt.savefig(out, dpi=300)
	print("Plot created!")


# === Execution ===============================================================

if __name__ == "__main__":

	p = argparse.ArgumentParser(description="Create a typhoon plot of protein locations from read data")
	p.add_argument("chr_name",				help="Name of the chromosome of interest")
	p.add_argument("origin",	type=int,	help="Position of interest (in bp)")
	p.add_argument("offset",	type=int,	help="Distance to capture (in either direction of origin)")
	p.add_argument("out",					help="Filename of outputted plot")
	p.add_argument("files",		nargs="+",	help="List of BAM/SAM read files to make plots out of")
	p.add_argument("condensed",	nargs="?", default=False, type=bool, help="Specify whether the plot should be visually condensed or not")
	args = p.parse_args()

	main(vars(args))
