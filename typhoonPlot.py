#!/usr/bin/env python

"""
typhoonPlot.py
Author: bmle

Makes typhoon plots from sequence alignment data. Based off an R script that also creates typhoon plots.

"""

import argparse
import pysam
import matplotlib.pyplot as plt
import plotUtils

def main(chr_name, origin, offset, out, files, condensed = False):
	"""
	Create a typhoon plot of reads (using template length)
	:param chr_name: Name of the chromosome (as stated in the BAM/SAM file)
	:param origin: Genomic location to start capturing from
	:param offset: How many bases (on either side of `origin`) to capture
	:param out: Name of outputted file
	:param files: BAM/SAM files to analyze
	:param condensed: Whether the plot should be visually condensed or not (optional)
	:return: A file at `out` containing one or more typhoon plots at the specified genomic coordinates
	"""

	print("Calculating overlaps...")

	# Make matrix to store overlaps, an entry for each file,
	# and a row for each read length and a column for each nucleosomal position
	data = []
	start = origin - offset
	end = origin + offset
	max_len = 250	# Maximum length of read to capture (250 should capture everything)

	# Iterate over all data files
	for file in files:

		# Load each file and filter for region of interest
		temp_reads = pysam.AlignmentFile(file, "rb")
		reads = temp_reads.fetch(chr_name, start, end)

		# Convert reads to template lengths, and append to list of all overlaps
		temp_array = plotUtils.calculate_reads(reads, start, end, max_len)
		data.append(temp_array)

	print("Finished calculating overlaps!")


	# =========================================================================
	# Make "typhoon" plot(s) from alignment data
	# =========================================================================
	print("Generating plot...")

	# Make the color scheme for the plot
	my_map = plotUtils.cmap()

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

		# Annotate gene bodies
		ax = plotUtils.make_gene_bodies(ax)

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
