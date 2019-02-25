#!/usr/bin/env python

"""
typhoonPlot.py
Author: bmle

Makes a kernel density plot from sequence alignment data at the specified chromosome, origin, and offset.
Uses a gaussian kernel.

Produces a single file containing one or more kernel density plots for the specified files.

"""

import argparse
import numpy as np
import pysam
import matplotlib.pyplot as plt
from scipy.stats import norm

def main(chr_name, origin, offset, out, files, bandwidth = 30):

	# =========================================================================
	# Calculate midpoints for each read, using each read's template length
	# =========================================================================
	print("Calculating midpoints...")

	# Matrix to store midpoints, a value in "data" for each file
	data = []
	start = origin - offset
	end = origin + offset

	# Iterate over all data files
	for file in files:

		# Load each file and filter for region of interest
		# noinspection PyUnresolvedReferences
		tempf = pysam.AlignmentFile(file, "rb")
		f = tempf.fetch(chr_name, start, end)

		# Loop over all reads in file, calculate midpoints
		temp = []
		for r in f:
			midpoint = r.reference_start + round(r.template_length / 2)
			temp.append(midpoint)
		data.append(temp)

	print("Finished calculating midpoints!")


	# =========================================================================
	# Make kernel plots from alignment data
	# =========================================================================
	print("Generating plot...")

	# Make the main figure
	fig = plt.figure(figsize=(15, 6 * len(data)))
	fig.suptitle("DM936-DM941", size=16, y=0.99)

	x_d = np.arange(start, end, bandwidth)

	# Add each data set as a subplot
	for i, (datum, f_name) in enumerate(zip(data, files), start=1):
		print("Making subplot " + str(i) + " of " + str(len(data)) + "...")
		ax = fig.add_subplot(len(data), 1, i)

		# Make the plot using kernel density estimation
		density = sum(norm(xi).pdf(x_d) for xi in datum)  # gaussian kernel
		ax.fill_between(x_d, density, alpha=0.6)
		ax.axis([start, end, 0, 50])

		# Set labels
		ax.set_ylabel("Read counts")
		ax.minorticks_on()
		ax.grid(which="both")

	# Save the plot
	plt.tight_layout(rect=[0, 0, 1, 0.98])
	plt.savefig(out, dpi=300)
	print("Plot created!")


# === Execution ===============================================================

if __name__ == "__main__":

	p = argparse.ArgumentParser(description="Create a kernel density plot from read data")
	p.add_argument("chr_name",				help="Name of the chromosome of interest")
	p.add_argument("origin",	type=int,	help="Position of interest (in bp)")
	p.add_argument("offset",	type=int,	help="Distance to capture (in either direction of origin)")
	p.add_argument("out",					help="Filename of outputted plot")
	p.add_argument("files",		nargs="+",	help="List of BAM/SAM read files to make plots out of")
	p.add_argument("bandwidth", nargs="?", type=int, default=30, help="Specify the bandwidth of the gaussian kernel (default=30)")
	args = p.parse_args()

	main(vars(args))
