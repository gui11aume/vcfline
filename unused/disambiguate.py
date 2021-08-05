#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from collections import defaultdict
from gzopen import gzopen

def collect(f):
   current = ''
   increment = 0
   score = defaultdict(int)
   for line in f:
      if line[0] == '@': continue
      items = line.split()
      for item in items[10:]:
         if item.startswith('AS'):
            AS = item.replace('AS:i:', '')
      if items[0] == current:
         increment += 1
      else:
         increment = 1
         current = items[0]
      score['%s_%d' % (current,increment)] = int(AS)
   return score


if __name__ == '__main__':
   with gzopen(sys.argv[1]) as f:
      scores1 = collect(f)
   with gzopen(sys.argv[2]) as f:
      scores2 = collect(f)
   for key in set(scores1) | set(scores2):
      print scores1[key] - scores2[key]
