#!/usr/bin/env python3

import os
import sys

sys.stderr.write(f"\n\n### Checking if you have the following dependencies in your $PATH:\n\n")
Dependencies = ["samtools", "augustus", "bioawk", "awk", "perl", "python3", "bedtools"]
sys.stderr.write(str(Dependencies) + "\n\n\n")

sys.stderr.write(f"### Paths found:\n\n")
for dependency in Dependencies:
    os.system(f"which {dependency}")

sys.stderr.write("\n")
