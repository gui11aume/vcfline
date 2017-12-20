#!/usr/bin/env python
# -*- coding:utf-8 -*-

'''
Transforms a vcf file into a dummy vcf that can be used for
calibration purposes. The changes are deterministic: indels
are replaced by their complement sequence and SNPs are
replaced by a nucleotide that is neither the reference nor
the orginal one.
'''

import string
import sys

NT = frozenset('GATC')
TR = string.maketrans('GATC', 'TCGA')

def change(seq):
   return seq.translate(TR)

def main(f):
   for line in f:
      if line[0] == '#':
         sys.stdout.write(line)
         continue
      items = line.split('\t')
      if items[4] != '.':
         # Variation from reference.
         if len(items[4]) + len(items[3]) > 2:
            # Insertion or deletion.
            items[4] = change(items[4])
         else:
            # SNP. In principle this command will always produce
            # the same variant if the file is re-run.
            items[4] = next(iter(NT.difference(items[3:5])))
        sys.stdout.write('\t'.join(items))
      else:
         sys.stdout.write(line)


if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      main(f)
