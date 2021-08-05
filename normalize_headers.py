#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

from gzopen import gzopen


KEEP = {
        " 2L:1": "2L",
        " 2R:1": "2R",
        " 3L:1": "3L",
        " 3R:1": "3R",
        " 4:1":  "4",
        " X:1":  "X",
        " Y:1":  "Y",
        " dmel_mitochondrion_genome:1": "mt",
    }


def read_genome(f, tag):
  '''Read a fasta file and return a dictionary whose keys are the
  sequence names and values are the sequences in text format.'''

  genome = dict()

  # Read everything.
  txt = f.read()

  # Separate chromsomes.
  segments = txt.split('>')

  for segment in segments:
    if not segment: continue
    # Extract header.
    (header, seq) = segment.split('\n', 1)
    for (code, chrom) in KEEP.items():
       if code in header:
          header = tag + ':' + chrom
          # Remove new lines characters.
          genome[chrom] = ('>%s' % header, seq.replace('\n', ''))
          continue

  return genome


def main(f, tag):
  # Read genome (and throw away weird chromosomes).
  genome = read_genome(f, tag)

  # Reorder chromosomes.
  for name in ('2L', '2R', '3L', '3R', '4', 'X', 'Y', 'mt'):
    (header, seq) = genome[name]
    print header
    print seq

if __name__ == '__main__':
  tag = sys.argv[2]
  with gzopen(sys.argv[1]) as f:
    main(f, tag)
