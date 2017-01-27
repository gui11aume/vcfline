#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

from gzopen import gzopen

KEEP = frozenset([
    "2L", "2R", "3L", "3R", "4", "X", "mitochondrion_genome",
])

def main(f):

  # Print only the chromosomes in KEEP.
  printing = False
  buffer = ''

  for line in f:
    # Extract chromosome from fasta header.
    if line[0] == '>':
      chrom = re.sub(r'\s.*', '', line[1:])
      if chrom in KEEP:
        sys.stdout.write(buffer)
        sys.stdout.write(line)
        buffer = '\n'
        printing = True
      else:
        printing = False
    elif printing:
      # If chromosome in KEEP, print sequence.
      # Remove new line character for convenience.
      sys.stdout.write(line.rstrip())

  # Print the final new line character.
  sys.stdout.write(buffer)

if __name__ == '__main__':
  with gzopen(sys.argv[1]) as f:
    main(f)
