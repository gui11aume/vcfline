#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

MINLEN = 26

subst_A6 = {
    '1 2L:1': '>2L_A6\n',
    '2 2R:1': '>2R_A6\n',
    '3 3L:1': '>3L_A6\n',
    '4 3R:1': '>3R_A6\n',
    '5 4:1': '>4_A6\n',
    '6 X:1': '>X_A6\n',
    '7 dmel_mitochondrion_genome:1': '>mt_A6\n',
}

subst_A7 = {
    '1 2L:1': '>2L_A7\n',
    '2 2R:1': '>2R_A7\n',
    '3 3L:1': '>3L_A7\n',
    '4 3R:1': '>3R_A7\n',
    '5 4:1': '>4_A7\n',
    '6 X:1': '>X_A7\n',
    '7 dmel_mitochondrion_genome:1': '>mt_A7\n',
}

def read_genome(f, substdict):
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
    name = substdict[header]
    genome[name] = seq.replace('\n', '')

  return genome


def chop(genome):
  for (chrom,seq) in sorted(genome.items()):
    start = 1
    chrom = chrom.strip()
    for segment in seq.upper().split('GATC'):
      end = start + len(segment) - 1
      if end - start > MINLEN:
        print '%s_%d_%d\n%s' % (chrom,start,end,segment)
      start = end+5

if __name__ == '__main__':
  with open(sys.argv[1]) as f:
    chop(read_genome(f, subst_A7))
  with open(sys.argv[2]) as f:
    chop(read_genome(f, subst_A6))
