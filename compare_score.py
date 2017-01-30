import sys

from collections import defaultdict
from math import log

def get_scores(f):
    Scores = defaultdict(float)
    for line in f:
        items = line.split()
        readnm = items[0]
        for item in items:
            if item.startswith('AS'):
                Scores[readnm] += log(1 + float(item.split(':')[2]))
                break
    return Scores

with open(sys.argv[1]) as f:
    S1 = get_scores(f)
with open(sys.argv[2]) as f:
    S2 = get_scores(f)

indis = 0
dis = 0
for nm in S1:
    if nm not in S2: continue
    if S1[nm] == S2[nm]: indis += 1
    else:                dis += 1

print indis
print dis
