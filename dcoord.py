#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

def read_genome(f):
   '''Read a fasta file and return a dictionary whose keys are the
   sequence names and values are the sequences in text format.
   Remove pT2 and weird chromosomes.'''

   genome = dict()
   txt = f.read()
   segments = txt.split('>')
   for segment in segments:
      if not segment: continue
      (header,seq) = segment.split('\n', 1)
      name = re.sub(r'\s.*', '', header)
      # Remove "chrUn_GL456385' etc. and pT2
      if '_' in name or 'pT2' in name: continue
      genome[name] = seq.replace('\n', '') 
   return genome

if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      genome = read_genome(f)
   chrom = sys.argv[2]
   start = int(sys.argv[3])
   end = int(sys.argv[4])
   print genome[chrom][start:end]
