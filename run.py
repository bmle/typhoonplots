#!/usr/bin/env python

"""
run.py
Author: bmle

Runs the scripts in this repository.

"""

import kernelDensityPlot
import typhoonPlot
import typhoonPlot_multiWindows

# Data to analyze
data = ["/data/home/vt26/DM936/DM936_chr2_pho5_m1_2018-11-29-13-46.bam",
        "/data/home/vt26/DM937/DM937_chr2_pho5_m1_2018-11-29-14-11.bam",
        "/data/home/vt26/DM938/DM938_chr2_pho5_m1_2018-11-29-14-37.bam",
        "/data/home/vt26/DM939/DM939_chr2_pho5_m1_2018-11-29-15-00.bam",
        "/data/home/vt26/DM940/DM940_chr2_pho5_m1_2018-11-29-15-24.bam",
        "/data/home/vt26/DM941/DM941_chr2_pho5_m1_2018-11-29-15-48.bam",]
out = "../../analyses/newVisuals/936-941.png"

# Set up coordinates
chr_name = "tpg|BK006936.2|"
origin = 430250
offset = 2750

# Settings specifically for kernelDensityPlot.main
bandwidth = 30

# Settings specifically for typhoonPlot.main
condensed = True

# Settings specifically for typhoonPlot_multiWindows.main
extent = 10000
limit = 500000
file = "/data/home/vt26/DM936/DM936_chr2_pho5_m1_2018-11-29-13-46.bam"
# out = "../../analyses/newVisuals/multiWindow_test.png"

# =============================================================================
# Run scripts (comment out/in the lines below to run them)
# =============================================================================

# kernelDensityPlot.main(chr_name, origin, offset, out, data, bandwidth)

typhoonPlot.main(chr_name, origin, offset, out, data, condensed)

# typhoonPlot_multiWindows.main(chr_name, origin, extent, limit, out, file)