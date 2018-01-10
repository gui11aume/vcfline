#!/usr/bin/env python
# -*- coding:utf-8 -*-

import re
import sys

from gzopen import gzopen

dmel = frozenset([
    '2L', '3L', '2R', '3R', '4', 'X', 'Y',
    'dmel_mitochondrion_genome',
])

def norm_mt(txt):
   # Replace "r'dmel_mitochondrion_genome" by "mt".
   return re.sub(r'dmel_mitochondrion_genome', 'mt', txt)


def process_comment(line):
   txt = ''
   if line.startswith('##contig=<ID=dmel_mitochondrion_genome'):
      txt = '##contig=<ID=mt,length=19524>\n'
   elif line.startswith('##contig'):
      (contig,) = re.search(r'ID=([^,]+)', line).groups()
      # Print only main contigs.
      if contig in dmel: txt = line
   # Keep following header lines.
   elif line.startswith('##'): txt = line
   elif line.startswith('#CHROM'): txt = line
   return txt


def main(f):
  SNPs = 0
  for line in f:
     if line[0] == '#': line = process_comment(line)
     else: SNPs += 1
     sys.stdout.write(norm_mt(line))

  return SNPs

if __name__ == '__main__':
   SNPs = main(gzopen(sys.argv[1]))
   sys.stderr.write('SNPs: %d\n' % SNPs)
