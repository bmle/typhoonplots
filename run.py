#!/usr/bin/env python

"""
run.py
Author: bmle

Runs the scripts in this repository.

"""

import time
import kernelDensityPlot
import typhoonPlot
import typhoonPlot_multiWindows

# BAM/SAM files and genomic coordinates to analyze
data = ["/data/home/vt26/DM936/DM936_chr2_pho5_m1_2018-11-29-13-46.bam",
        "/data/home/vt26/DM937/DM937_chr2_pho5_m1_2018-11-29-14-11.bam",
        "/data/home/vt26/DM938/DM938_chr2_pho5_m1_2018-11-29-14-37.bam",
        "/data/home/vt26/DM939/DM939_chr2_pho5_m1_2018-11-29-15-00.bam",
        "/data/home/vt26/DM940/DM940_chr2_pho5_m1_2018-11-29-15-24.bam",
        "/data/home/vt26/DM941/DM941_chr2_pho5_m1_2018-11-29-15-48.bam",]
chr_name = "tpg|BK006936.2|"  # Name of the chromosome (as stated in the BAM/SAM file)
origin = 277555  # Location to start capturing from
offset = 2750    # How many bases (on either side of `origin`) to capture

# Settings specifically for kernelDensityPlot.main
bandwidth = 30  # Size of kernel to use
out_kde = "../../analyses/newVisuals/936-941_kde.png"

# Settings specifically for typhoonPlot.main
condensed = False   # Whether the plot should be visually condensed or not
out_tp = "../../analyses/newVisuals/936-941.png"

# Settings specifically for typhoonPlot_multiWindows.main
extent = 10000  # How wide each window should be
limit = 500000  # How many bases in total to capture (starting from `origin`)
file = "/data/home/vt26/DM936/DM936_chr2_pho5_m1_2018-11-29-13-46.bam"
out_mw = "../../analyses/newVisuals/mwtest.png"


# =============================================================================
# Run scripts (comment out/in the lines as needed)
# =============================================================================
start_time = time.clock()


# kernelDensityPlot.main(chr_name, origin, offset, out_kde, data, bandwidth)

# typhoonPlot.main(chr_name, origin, offset, out_tp, data, condensed)

typhoonPlot_multiWindows.main(chr_name, origin, extent, limit, out_mw, file)


print("Completed in %s seconds!" % (time.clock() - start_time))
