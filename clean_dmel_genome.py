#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

from gzopen import gzopen


KEEP = frozenset([
    "2L", "2R", "3L", "3R", "4", "X", "dmel_mitochondrion_genome",
])

def read_genome(f):
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
    (header,seq) = segment.split('\n', 1)
    name = re.sub(r'\s.*', '', header)
    if name not in KEEP: continue
    # Remove new lines characters.
    genome[name] = ('>%s' % header, seq.replace('\n', ''))

  return genome


def main(f):
  # Read genome (and throw away weird chromosomes).
  genome = read_genome(f)

  # Reorder chromosomes.
  for name in ('2L', '2R', '3L', '3R', '4', 'X',
        'dmel_mitochondrion_genome'):
    (header, seq) = genome[name]
    print header
    print seq

if __name__ == '__main__':
  with gzopen(sys.argv[1]) as f:
    main(f)
