# An Implementation of Pitman-Yor Simple Range Concatenating Adaptor Grammars

# Summary
This repository contains extensions over Mark Johnson's original adaptor
grammar code to implement a constrained version of simple range concatenating
grammars (SRCGs). The contraints are on the number of symbols that can be
handled in an input sequence, and particularities around the shapes of SRCG
rules supported.

# Attribution
Please cite the following papers as appropriate:
* [Adaptor Grammars for Learning Non−Concatenative Morphology](http://aclweb.org/anthology//D/D13/D13-1034.pdf), Jan A. Botha and Phil Blunsom.  In *Proceedings of the 2013 Conference on Empirical Methods in Natural Language Processing*, October, 2013. Association for Computational Linguistics.
* [Adaptor Grammars: A Framework for Specifying Compositional Nonparametric Bayesian Models,
Mark Johnson, Thomas L. Griffiths, and Sharon Goldwater](http://papers.nips.cc/paper/3101-adaptor-grammars-a-framework-for-specifying-compositional-nonparametric-bayesian-models.pdf). In *Advances in Neural Information Processing Systems*, 2007
* [Improving nonparameteric Bayesian inference: Experiments on unsupervised word segmentation with adaptor grammars](http://www.aclweb.org/anthology/N09-1036), Mark Johnson and Sharon Goldwater. In *Proceedings of NAACL-HLT*, 2009, Association for Computational Linguistics.

# Dependencies
Boost Serialization is used to save and load models, but the software can be
built without this functionality by passing including
`NOSERIALIZE=-DNOSERIALIZE` on the build command line.

# Configure and Build
1. Update the paths for `BOOST_LIBDIR` and `BOOST_INCLUDEDIR` in the `Makefile`
   to reflect your local setup.
2. Build with `make`.

# Grammar file format
Grammar files can be specified in a 'friendly' format and converted into
converted into the final format by `preprocess-grammar.py`, or they can be specified directly in the final format.
In both cases, the overall form of a rule is as follows:
```
[w [a [b]]] X --> Y Z
```
where:
* X is a nonterminal and Y and Z are either terminals or nonterminals (see details further down),
* w is the Dirichlet hyper-parameter (i.e., pseudo-count) associated with this rule (a positive real)
* a is the PY "a" (discount) constant associated with X (a positive real less than 1). If a=0, the prior reduces to a CRP. If a=1, X is *unadapted*.
* b is the PY "b" (strength) constant associated with X (a positive real)

Default values for a and b can also be specified via command line arguments. As
an example, when specifying defaults and using hyperparameter inference, a rule
for an adapted non-terminal could be specified simply as
```
X --> Y Z
```
while a rule using an unadapted non-terminal could be written as
```
1 1 Q --> P S
```
In the rest of this section, we omit the prefix parameters and focus on the rule symbols only.

## Friendly format
It is 'friendly' because you can use named variables and have some flexibility.
* CFG rule: `X --> Y Z` (equivalent to `X(ab) --> Y(a) Z(b)`)
* SRCG rule: `X(abc) --> V(a,c) Z(b)`
* Terminal CFG rule: `T(t) --> eps`, where `t` is a literal.

## Final format
Running the previous three rules through `preprocess-grammar.py` should produce the following:
```
X (((0 0)) ((0 1))) --> Y Z
X (((0 0) (0 2)) ((0 1))) --> V Z 
T (((0 0))) --> t 
```
You don't want to write this by hand. Roughly: The variables are replaced with an indexing scheme. (TODO: elaborate).

# Restrictions
The SRCG parser makes the following restrictions:
1. Rules can only have one or two symbols on the right-hand side. (This is a
  restriction the original adaptor grammar code did not have, but introducing
  it made the extensions simpler.)
2. Terminal symbols are only handled in unary CFG rules like the one given
  earlier. Terminal symbols should not appear in non-CFG rules.
3. The number of arguments on a left-hand side non-terminal symbol is **at most 5** (e.g.
  `A(a,b,c,d,e) --> ...`).
4. The total arguments across both right-hand side non-terminal symbols is **at most 9* and the symbol with higher arity should come first (e.g. `A(abcdefghi) --> B(a,c,e,g,i) C(b,d,f,h)`).
5. Concatenation of variables within a LHS argument is limited. In a binary
   rule, all the variables should combine into a single argument (example in
   (4.)) or remain as one variable per argument (example in 3.). Partial
   concatenation of variables is only allowed for unary rules, e.g.
   `A(xy,z) -> B(x,y,z)` and `A(x,yz) -> B(x,y,z)` are allowed,
   but `A(xy,z) -> B(x,z), C(y)` is not.
6. (There may be other subtle restrictions, but you could determine that from
   the code or expect the program to misbehave or warn you if they are violated.)

As for the sampler, it does not support *table label resampling* for discontiguous non-terminals.

## Input length limit
The method used to address chart cells impose some environment-specific contraints on the maximum input length. This maximum is set in `consts.h` and defaults to 21.

Run the `indexing` program to assess your situation. It will alert you to key clashes that would require changing the maximum input length.

# Example Script
(TODO)

Basic invocation: `./py-srcag preprocessed\_grammar < corpus`

