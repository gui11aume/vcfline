#!/usr/bin/env python
# -*- coding:utf-8 -*-

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
         if len(items[4]) + len(items[3]) > 2:
            # Insertion or deletion.
            items[4] = change(items[4])
         else:
            # SNP
            items[4] = next(iter(NT.difference(items[3:5])))
            sys.stdout.write('\t'.join(items))
      else:
         sys.stdout.write(line)


if __name__ == '__main__':
   with open(sys.argv[1]) as f:
      main(f)
