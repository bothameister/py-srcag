#!/usr/bin/env python
# encoding: utf-8

usage = """%prog [options] < inputgrammar > outputgrammar

(c) Jan Botha

Performs conversion between rule data formats.
Possible input formats (may be mixed in same file):
  cfg:     A -> B C               # Both formats can have 
  mcfg:    A(x,zy) -> B(x,z) C(y) # comments like this,
           A(a,b) -> eps          # although they're stripped from the output
           A(abc) -> B(a,b,c)
           A(abcd) -> B3(*) C1(*) # will generate the valid range permutations (3 ranges x 1 range)
           A(abcd) -> B3(*) C(*)  # and default to dim 1 when it's not given

The first NT may be preceded by up to 3 numbers,
which are passed through to the output unmodified.
Each argument in a terminating mcfg-rule must be one terminal only,
so 'A(ab,c) -> eps' is not allowed.

Terminating rules won't be detected in CFG format.

Output format (suitable for input to pymcfg binary):
           A (((0 0)) ((1 0))) -> B C # for input A(x,y) -> B(x) C(y)"""

import re, sys, optparse
from collections import defaultdict
import permute_ranges

class Grammar:

  def __init__(self, options, instream=None):
    self.__options = options
    self.rules = []
    self.startsym = None

    if not instream is None:
      self.__rules_from_stream(instream)

  def __rules_from_stream(self, instream):
    arrow = self.__options.arrow
    for (pre, to_read) in split_optional_prefix(instream):
      self.add_rules(Rule.read_str(to_read, pre, arrow))
 
  def __repr__(self):
    return '\n'.join( [r.to_pymcfg_fmt() + '\t# ' + r.to_mcfg_fmt(False) for r in self.all_rules()] )

  def as_pymcfg(self):
    return repr(self)

  def as_mcfg(self):
    return '\n'.join( [r.to_mcfg_fmt() for r in self.all_rules()] )

  def add_rules(self, rules):
    """ Add the list of Rules to grammar, creating new dummy rules as necessary
    to ensure unique mapping of rhs to lhs"""
    arrow = self.__options.arrow
    default_prerule = self.__options.default_prerule
    for rule in rules:
      clashset = self.__clash_set(rule.lhs, rule.rhs)
      if len(clashset) == 0:
        if self.startsym is None: self.startsym = rule.lhs
        self.rules.append(rule)
      else:
        if len(clashset) == 1:
          self.__raise_existing_singleton(clashset[0])
        self.__transform_and_add(rule, self.__dummy_name(rule.lhs,rule.rhs))

  def __transform_and_add(self, mainrule, dummy):
    """ Transforms mainrule "v A -> B"
    into "v A -> dummy" and "d dummy -> B", where v the potential prefix str
    and d is the default prefix string. Then adds both to the grammar."""
    arrow = self.__options.arrow
    default_prerule = self.__options.default_prerule

    num_lhs_ranges_needed = max([idx[0] for arg in mainrule.lhs_vars for idx in arg]) + 1
    if len(mainrule.rhs) == 2: 
      # only handle binary rules if the parent has a single range
      assert num_lhs_ranges_needed==1
      dummy_lhs_vars = trivial_lhs_vars(1)
    else:
      dummy_lhs_vars = trivial_lhs_vars(1, num_lhs_ranges_needed)

    # Make new          prefstr A -> dummy
    to_add = Rule(mainrule.prefixstr, arrow, lhs=mainrule.lhs, lhs_vars=dummy_lhs_vars, rhs=[dummy,])
    
    # Make mainrule =   defstr dummy -> B
    mainrule.lhs = dummy
    mainrule.prefixstr = default_prerule
    self.add_rules((to_add, mainrule))

  def remove_rule(self, rule):
    """ Removes rule from grammar """
    self.rules.remove(rule) 

  def __clash_set(self, lhs, rhs):
    """ Returns set of existing rules with which the proposed lhs -> rhs would clash """
    concat = lhs + 'l' + ''.join(rhs) 
    return [r for r in self.rules if (r.lhs==lhs and r.rhs==rhs) or (r.lhs.startswith(concat))]

  def __raise_existing_singleton(self, existing):
    """ Promote existing singleton to make use of a counting symbol name """
    self.remove_rule(existing)
    newlhs = self.__dummy_name_render(existing.lhs, existing.rhs, 1)
    self.__transform_and_add(existing, newlhs)

  def __dummy_name_render(self, lhs, rhs, num):
    return lhs + 'l' + ''.join(rhs) + '%.3d' % num

  def __dummy_name(self, lhs, rhs):
    num = 1 + len(self.__clash_set(lhs,rhs))
    return self.__dummy_name_render(lhs,rhs,num)

  def all_rules(self, plain=False):
    """ Generator """
    assert not self.startsym is None

    for r in self.rules:
      if r.lhs == self.startsym or plain:
        yield r

    if not plain:
      for r in self.rules:
        if r.lhs != self.startsym:
          yield r


  def add_corpus_terminals(self, src):
    """ Adds terminal rules for each unique non-whitespace symbol in src stream """
    arrow = self.__options.arrow
    preterm_sym = self.__options.preterm_sym
    default_prerule = self.__options.default_prerule
    whitespace = ' \t\n\v\r'
    termset = set()
    for line in src:
      lineterms = set([term for term in line.strip().split(' ') 
                          if len(term)>0 and not term in whitespace])
      termset.update(lineterms)
    conflicts = 0
    for c in '#_':
      if c in termset:
        sys.stderr.write("Error: Corpus contains '%s', which has a special meaning in py-mcfg. Please replace it.\n" % c)
        conflicts += 1
  
    if conflicts > 0:
      sys.stderr.write("Some free characters are: " + ' '.join(unused(termset, conflicts)) + "\n")
      sys.exit(1)

    for term in sorted(termset):
      self.add_rules((Rule.from_terminal(preterm_sym, term, default_prerule, arrow),))

def unused(termset,n=3):
  ret=[]
  for start in (ord('a'), ord('A')):
    for i in xrange(26):
      c=chr(start+i) 
      if not c in termset:
        ret.append(c)
      if len(ret) == n:
        return ret
  return ret

class Rule:
  """ All information about a rule / clause
  lhs, rhs1, rhs2 are predicate labels
  Assumes each LHS-arg is one or more vars, each RHS-arg is a single var, and each arg in a terminating rule is a single terminal
  
  lhs_vars has the structure of the rhs and maps each element to its position on the lhs.
  lhs_vars is a list with n symbolmappings if the rhs has n grammar symbols.
  A symbolmapping is a list of m pairs (a,b) if the particular rhs symbol has m arguments.
  A pair (a,b) says that the particular rhs argument (i.e. single var) goes into 
        lhs argument a, position b.

  e.g. A(xy) -> B(x) C(y) => lhs_vars = (((0 0)) ((0 1)))
       A(x,y) -> B(x) C(y) => lhs_vars = (((0 0)) (1 0)))
       A(xy,z) -> B(x) C(y,z) => lhs_vars = (((0 0)) ((0 1) (1 0)))
       A(xyz) -> B(x,z) C(y) => lhs_vars = (((0 0) (0 2)) ((0 1)))

  Except for terminating rules, where lhs_vars stores a flat list of terminal ids
  """
  
  def __init__(self, prefixstr='', arrow = '->', lhs='', lhs_vars=[], rhs=[],terminating=False):
    self.lhs = lhs
    self.lhs_vars = lhs_vars
    self.rhs = rhs
    self.arrow = arrow
    self.prefixstr = prefixstr
    self.terminating = terminating

    assert isinstance(rhs, list)
    assert isinstance(lhs_vars, list)
    assert isinstance(lhs, str)

  ############## Input #####################
  @classmethod
  def from_terminal(cls, preterminal, terminal, prefixstr='', arrow = '->'):
    """ Returns a single Rule instance for a terminating rule """
    return cls(prefixstr, arrow, lhs=preterminal, rhs=[terminal,], 
        lhs_vars=trivial_lhs_vars(1), terminating=True)

  @classmethod
  def from_mcfg_str(cls, string, prefixstr='', arrow = '->'):
    """ Read a Rule from the string
        string is a line like one of these
         A(xy,z) -> B(x) C(y,z)  #comment
         C(a) -> eps
         This method returns a list of rules (possible due to abstraction).
    """
    full_lhs, full_rhs = cls.__split_left_right(string, arrow)
    # LHS
    lhs, lhs_var_str = cls.__get_pred(full_lhs) # C  and  a,b

    if full_rhs == "eps": # terminating rule, assuming single arg in predicate
      assert len(lhs_var_str)==1
      rhs = ''.join(lhs_var_str)
      return (cls.from_terminal(lhs, rhs, prefixstr, arrow),)
    else:
      rhs_parts = [cls.__get_pred(pred_str) for pred_str in full_rhs.split(' ')]
      if all( [ rhs_part[1] == '*' for rhs_part in rhs_parts] ): # A(abc) -> B2(*) C(*) abstr
        return cls.__instantiate_range_abstraction(full_lhs, rhs_parts, prefixstr, arrow)
      else: # no range abstraction
        varmap = cls.__index_var_str(lhs_var_str)
        rhs_symbols = []
        lhs_vars = []
        for sym, var_str in cls.order(rhs_parts)[0]: # A  and  x,z
          rhs_symbols.append(sym)
          arg_map = []
          for var in var_str.split(','):
            assert len(var) == 1
            assert var in varmap
            arg_map.append(varmap[var])
          lhs_vars.append(arg_map)
        rhs = rhs_symbols
        return (cls(prefixstr, arrow, lhs, lhs_vars, rhs), )

  @classmethod
  def order(cls, rhs_parts):
    # get dimensions encoded in symbols and rearrange if necessary so highest dim comes first
    childlens = [get_sym_dim(rhs_part) for rhs_part in rhs_parts]
    if len(rhs_parts) == 2 and (childlens != sorted(childlens, reverse=True)):
      childlens[0], childlens[1] = childlens[1], childlens[0]
      rhs_parts[0], rhs_parts[1] = rhs_parts[1], rhs_parts[0]
    return rhs_parts, childlens

  @classmethod
  def __instantiate_range_abstraction(cls, full_lhs, rhs_parts, prefixstr, arrow):
    """ Create the variable permutations requested through range abstraction.
    rhs_parts is a list of two (rhs_sym, rhs_varlist) pairs
    returns a list of Rules.
    Recall that permute_ranges only returns 'one half' of permutations when dim(rhs1)==dim(rhs2)
    (motivated by making the parser itself more efficient).
        a,c & b,d
        b,d & a,c  << not returned by permute_ranges, therefore not parsed by pymcfg. Keep it that way.
    But 'sometimes' we need to be able to have both kinds when instantiating rules, e.g.
        A(abcd) -> B(a,c) C(b,d)
        A(abcd) -> C(a,c) B(b,d)
    These are not equivalent derivational steps because that's a function of symbol & terminal yield.
    'sometimes' is when the dims are equal and symbols different. This function implements these things too.
    """
    lhs, lhs_var_str = cls.__get_pred(full_lhs) # C  and  a,b
    assert len(rhs_parts) == 2   # range abstraction is only defined for binary rules
    assert len(lhs_var_str.split(',')) == 1 # range abstraction is only defined for rules with continuous lhs

    rhs_parts, childlens = cls.order(rhs_parts)
    assert sum(childlens) == len(lhs_var_str)
    ret=[]
    for rhsvars in permute_ranges.permute_vars(list(lhs_var_str), childlens[0]):
      # Standardly
      orderings = [lambda x:x]
      if rhs_parts[0][0]!=rhs_parts[1][0] and childlens[0] == childlens[1]: # 'sometimes', see above
        #also output an instantiation for the current permutation but with reversed symbols
        orderings.append(lambda x:reversed(x)) 

      for f in orderings:
        rulestr = cls.__make_rule_str(full_lhs, arrow, [r[0] for r in f(rhs_parts)], rhsvars)
        ret.extend(cls.from_mcfg_str(rulestr, prefixstr, arrow))

    return ret
  
  def __make_mcfg_str(self, withprefix=True):
    """ Render a string representation of this rule in mcfg format
    by reconstruct variable names and place them in predicates.
    """
    if self.terminating:
      rulestr = self.__get_mcfg_str_terminating()
    else:
      # variable names to use
      varnames2 = [chr(ord('a')+i) for i in range(len([x for y in self.lhs_vars for x in y]))]
    
      vardict = {} # maps an item in lhs_vars (e.g. (0 0) ) to the variable name used for it
      new_lhs_vars = []
      mx = len(varnames2)
      #for each sensible (id1, id2), see if it is somewhere in lhs_vars, pick a name for it and place in new_lhs_vars
      for id1 in xrange(mx):
        if len(varnames2) == 0: break
        new_lhs_vars.append([])
        for id2 in xrange(mx):
          if len(varnames2) == 0: break
          for i, rhs_arglist in enumerate(self.lhs_vars):
            for idxpair in rhs_arglist:
              if tuple(idxpair) == (id1,id2):
                var = varnames2.pop(0)
                new_lhs_vars[id1].append(var)
                vardict[(id1,id2)] = var
      
      # create rhs's simply by replacing from the vardict
      new_rhs_vars = [ [] for x in self.lhs_vars ]
      for i, rhs_arglist in enumerate(self.lhs_vars):
        for idxpair in rhs_arglist:
          new_rhs_vars[i].append( vardict[tuple(idxpair)])

      rulestr=Rule.__make_rule_str(Rule.__make_pred_str(self.lhs, new_lhs_vars), self.arrow, self.rhs, new_rhs_vars)

    if withprefix: 
      return friendly_join((self.prefixstr, rulestr))
    else:
      return rulestr

  def __get_mcfg_str_terminating(self):
    assert self.terminating
    return ' '.join((Rule.__make_pred_str(self.lhs, (self.rhs,) ), self.arrow, 'eps'))

  @classmethod
  def __make_rule_str(cls, full_lhs, arrow, rhs_syms, rhs_vars):
    """ Renders the supplied params as an mcfg-formatted rule string (excluding prefixstr) """
    new_rhs=[cls.__make_pred_str(sym, varlist) for (sym, varlist) in zip(rhs_syms, rhs_vars)]
    new_full_rhs = ' '.join(new_rhs)
    return ' '.join((full_lhs, arrow, new_full_rhs))

  @staticmethod
  def __make_pred_str(sym, varlist):
    """ Renders a predicate as a string. """
    if all([isinstance(v, str) for v in varlist]):
      # right-hand sides are flat comma separated single vars
      return '%s(%s)' % (sym, ','.join(varlist))
    else:
      # left-hand sides have command-separated multiple vars
      return '%s(%s)' % (sym, ','.join([''.join(mvars) for mvars in varlist]))

  @classmethod
  def from_cfg_str(cls, string, prefixstr='', arrow = '->'):
    """ Read a Rule from the string, which is a straightforward CFG rule, e.g. A -> B C """
    full_lhs, full_rhs = cls.__split_left_right(string, arrow)
    rhs = full_rhs.split(' ')
    lhs_vars = trivial_lhs_vars(len(rhs))
    return (cls(prefixstr, arrow, full_lhs, lhs_vars, rhs),)

  @classmethod
  def read_str(cls, to_read, prefixstr, arrow):
    """ Figures out which format is being input and passes on the call appropriately """
    # detect input format: alphanum(...) is sufficient to indicate MCFG-format
    if re.search(r'\w+\([^)]+\)', to_read) is not None:
      return cls.from_mcfg_str(to_read, prefixstr, arrow)
    else:
      return cls.from_cfg_str(to_read, prefixstr, arrow)

  ############# Output #####################################
  def to_pymcfg_fmt(self, withprefix=True):
    """ Outputs the rule in pymcfg format
    """
    return self.__make_pymcfg_str(withprefix)

  def to_mcfg_fmt(self, withprefix=True):
    """ Outputs the rule in mcfg format
    """
    return self.__make_mcfg_str(withprefix)
  

  ############# Input Helpers ##############################
  
  @staticmethod
  def __get_pred(string):
    """ parses predicate from start of string, 
        e.g. C(y,z)  or A(x)  or B(xy) or C(*)
        The name is any alphanum string """
    pred_matcher = re.compile(r'(\w+)\((.+)\)').match
    matches = pred_matcher(string)
    if matches is None:
      print "Could not parse predicate from string '%s'" % string
      return None
    assert len(matches.groups()) == 2
    return matches.groups() # C  and  a,b
  
  @staticmethod
  def __split_left_right(string, arrow):
    """ Splits rule on arrow into left and right """
    string = string.split('#')[0] # strip comment

    sides = [s.strip() for s in string.split(arrow)]
    assert len(sides) == 2
    return sides

  @staticmethod 
  def __index_var_str(string):
    """ Parses a string of variables into a dictionary denoting the position (i,j) or a var
        i is argument index, j is index within arg i
    """
    varmap = {}
    for (i, arg) in enumerate(string.split(',')):
      for (j, var) in enumerate(list(arg)):
        assert var not in varmap
        varmap[var] = (i, j)
    return varmap

  ############## Output Helpers #######################
  def __myrepr(self, s, before='(', after=')', sep=' '):
    """ Recursively output the contents of s, enclosing each item in before and after.
    Consecutive items at the same depth are separated by sep """
    if isinstance(s, list) or isinstance(s, tuple):
      return before + sep.join( map(self.__myrepr, s)) + after
    elif isinstance(s, int) or isinstance(s, float):
      return str(s)
    elif isinstance(s, str):
      return s
    else:
      sys.stderr.write('Untreated case: val=%s, type=%s\n' % (s, type(s)))
      assert False

  def __make_pymcfg_str(self, withprefix=True):
    main_rule = ' '.join( (self.lhs, self.__myrepr(self.lhs_vars), self.arrow, ' '.join(self.rhs)) )
    if withprefix:
      return friendly_join((self.prefixstr, main_rule))
    else:
      return main_rule
  
  ############## Grammar Transformation #################### 
############# end of Rule def

def get_sym_dim((symstr,varstr)):
  """ Internal function that tries to parse dimensionality from symbol name
      e.g. AB2_123 returns 123. Defaults to 1 if sym doesn't end with num 
  """
  if varstr == '*':
    mo = re.match(r'.+?(\d+)$', symstr) #skip non-greedily until capturing final nums
    if mo:
      return int(mo.group(1))
    return 1
  else:
    return len(varstr.split(','))

def trivial_lhs_vars(n, nvars=1):
  """ Creates lhs_vars appropriate for n RHS symbols with nvars var each, trivial order"""
  assert n==1 or nvars==1 #otherwise it's ambiguous what the 'trivial' case is
  return [ tuple([(v,i) for v in xrange(nvars)]) for i in xrange(n) ]  


def friendly_join(args):
  return ' '.join( [x for x in args if not x is None and x!=''] ) 


def split_optional_prefix(source):
  """ Generator that separates the optional rule prefix (of numerical values) 
      from the actual rule content, taking first alphanum character as boundary
  """
  matcher = re.compile(r"""
      \ *         # allow initial whitespace
      (           # capture rule prefix, i.e.
        (?:(?:[0-9e.-]+\ ){0,3}) # up to 3 reps of general numbers
      )(\w+.*)    # then capture the rest of the str which has to start with some \w s
      """, re.VERBOSE).match

  for line in source:
    line = line.strip()
    if re.match(r'\s*(#|$)', line): continue

    mo = matcher(line)
    if mo is None:
      print "FAILED ON '%s'" % line
    else:
      assert len(mo.groups()) == 2
      pre, rulecontent = [s.strip() for s in mo.groups()]
      yield pre, rulecontent

if __name__ == '__main__':
  parser = optparse.OptionParser(usage=usage)
  parser.add_option("-T", "--new-preterm-char", dest="preterm_sym",
                    help="Use this as preterminal symbol when on is needed", default = 'T')
  parser.add_option("-a", dest="arrow",
                    help="Arrow used in rules", default = '-->')
  parser.add_option("-f", dest="fmt",
                    help="Output format. Possible values are mcfg and pymcfg (default)", default = 'pymcfg')
  parser.add_option("-p", "--rule-prefix-default", dest="default_prerule",
                    help="Rule prefix for new preterminal rules", default = '1 1')
  parser.add_option("-t", "--terminals", dest="terms_from_corpus",
                    help="Generate terminal rules for all unique, non-whitespace tokens in the specified file; uses arg of -T for pre-terminal symbol", default = None)
  (options,args) = parser.parse_args()

  g = Grammar(options, sys.stdin)

  if not options.terms_from_corpus is None:
    g.add_corpus_terminals(open(options.terms_from_corpus,'r'))
  
  if options.fmt.lower() == 'pymcfg':
    print g.as_pymcfg()
  elif options.fmt.lower() == 'mcfg':
    print g.as_mcfg()
