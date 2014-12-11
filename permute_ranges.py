#!/usr/bin/env python
# encoding: utf-8

usage = """%prog 
Partition a range vector into two range vectors, for all valid permutations.
usage: %prog [options] indexstr N
\tindexstr\t- a comma-delimited string of values, or an undelimited string of chars
\tN\t\t- number of ranges to put in first output range vector. second one gets the rest

Example: %prog abcdef 2
({a,b,c,d}, {e,f});
({a,b,e,f}, {c,d});
({c,d,e,f}, {a,b});"""

from itertools import permutations
import sys, optparse

def r_ok(rvec, cmp=lambda x,y: x<=y):
  """returns true if the range vector is ordered"""
  return all( [cmp(rvec[i], rvec[i+1]) for i in xrange(len(rvec)-1) ])

def v_ok(varlist, cmp=lambda x,y: x<=y):
  """returns true if the variables appear in (lex) order"""
  return r_ok(varlist, cmp) # same functionality as r_ok, but want them 'separate'

def flat(rvec):
  """flatten the list of pairs into a flat list"""
  ret=[]
  for r in rvec:
    ret.extend(r)
  return ret

def mild_in(item, container):
  """Returns true iff item, or its reverse, is already in container"""
  assert len(item) == 2
  return item in container or tuple(reversed(item)) in container

def permute_vars(varlist, child1len):
  """Generator that yields the ordered permutations of two arg vectors, the first one having child1len args.
       e.g. permute_ranges( ('a', 'b', 'c', 'd', 'e', 'f'), 2)
  """
  assert isinstance(varlist, list)
  assert v_ok(varlist) #input must be ordered for our constraining tricks to work
  assert child1len > 0 and child1len < len(varlist)

  cache = set() #to hold what's already been generated

  for p in permutations(varlist):
    child_vars = (p[:child1len], p[child1len:])
    if all( [ v_ok(r) for r in child_vars ]) and not mild_in(child_vars, cache):
      cache.add(child_vars)
      yield child_vars

def permute_ranges(indices, child1len):
  """Generator that yields the ordered permutations of two range vectors, the first one having child1len ranges
  e.g. permute_ranges( (1,2,4,6,8,9), 2)
       permute_ranges( ('a', 'b', 'c', 'd', 'e', 'f'), 2)
  """
  assert len(indices) % 2 == 0
  assert r_ok(indices) #input must be ordered
  rvec = [ (indices[i], indices[i+1]) for i in xrange(0, len(indices)-1, 2)]
  assert all( [r_ok(r, cmp=lambda x,y: x<y) for r in rvec] ) # single range can't have left and right equal
  assert child1len > 0 and child1len < len(rvec)

  cache = set() #to hold what's already been generated

  for p in permutations(rvec):
    child_rvecs = (p[:child1len], p[child1len:])
    if all( [ r_ok(r) for r in child_rvecs ]) and not mild_in(child_rvecs, cache):
      cache.add(child_rvecs)
      yield [flat(rvec) for rvec in child_rvecs]

def decorate(rvecs):
  """Output range vectors into some desired string format"""
  return ', '.join(['{%s}' % ','.join([str(x) for x in rvec]) for rvec in rvecs])

if __name__ == '__main__':
  parser = optparse.OptionParser(usage=usage)
  parser.add_option("-p", "--prefix", dest="prefix", default="(", help="Prepend this to each yield")
  parser.add_option("-s", "--suffix", dest="suffix", default=");", help="Append this to each yield")
  (options,args) = parser.parse_args()

  index_str = args[0]
  if len(index_str.split(',')) == 1:
    index_str = list(index_str)
  else:
    index_str = index_str.split(',')
  if all([i.isdigit() for i in index_str]):
    index_str = [int(i) for i in index_str]
  
  num_ranges_first_batch = int(args[1])
  statements = []
  for x in permute_ranges(   index_str, num_ranges_first_batch ):
    statements.append( options.prefix + decorate(x) + options.suffix )

  print '\n'.join(statements)
