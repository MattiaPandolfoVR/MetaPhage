#!/usr/bin/env python3
# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
import os, sys
import re

report = sys.argv[1]
run_name = sys.argv[2]
copying = True
with open(report, 'rt') as inf, open("./MetaPhage_" + run_name + "_report.html", 'wt') as outf:
    for line in inf:
        if copying:
            if line.startswith('<div class="mqc-toolbox collapse">'):
                copying = False
            else:
                outf.write(line)
        elif line.startswith('</div>'):
            copying = True
