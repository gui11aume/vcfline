#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from gzopen import gzopen

JohnPoolCode = {
    '2': 'mitochondrion_genome',
    '3': '2L',
    '4': 'X',
    '5': '3L',
    '6': '4',
    '7': '2R',
    '8': '3R',
}

contigCode = {
    '##contig=<ID=2,length=19517>':
        '##contig=<ID=mitochondrion_genome,length=19517>',
    '##contig=<ID=3,length=23011544>': '##contig=<ID=2L,length=23011544>',
    '##contig=<ID=4,length=22422827>': '##contig=<ID=X,length=22422827>',
    '##contig=<ID=5,length=24543557>': '##contig=<ID=3L,length=24543557>',
    '##contig=<ID=6,length=1351857>': '##contig=<ID=4,length=1351857>',
    '##contig=<ID=7,length=21146708>': '##contig=<ID=2R,length=21146708>',
    '##contig=<ID=8,length=27905053>': '##contig=<ID=3R,length=27905053>',
}

def parse_and_print_comment(line):
  if line.startswith('##contig'):
    line = contigCode.get(line, '')
  sys.stdout.write(line)

def main(f):
  for line in f:
    if line[0] == '#':
      parse_and_print_comment(line)
      continue
    items = line.split()
    if items[3] == items[4]:
      continue
    if items[0] in JohnPoolCode:
      items[0] = JohnPoolCode[items[0]]
      print '\t'.join(items)

if __name__ == '__main__':
  try:
    f = gzopen(sys.argv[1])
  except IndexError:
    f = sys.stdin
  main(f)
