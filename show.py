#!/usr/bin/python

# Parses out the logP column from a 'trace' file.
# Used for quick check of convergence

# Usage: show.py trace1 [trace2 ...]


from matplotlib import pylab as plt
import sys

def series_from_src(src):
  data=[]
  for line in src:
    if line[0] == '#': continue
    fields = line.strip().split('\t')
    val = float(fields[2])
    data.append(val)
  return data

colors = list('rkbgc')
markers = list('ov^<>+')

if len(sys.argv) > 1:
  for i,src in enumerate(sys.argv[1:]):
    data = series_from_src(open(src,'r'))
    if len(data) > 0:
      markeridx = int(float(i) / len(colors)+0.5)
      col =  colors[ i % len(colors)] 
      plt.plot(xrange(len(data)), data, col + markers[markeridx % len(markers)], mec=col, label=src[-30:])
plt.legend()
plt.show()
