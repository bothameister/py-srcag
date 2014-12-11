// vim: tabstop=8:softtabstop=8:shiftwidth=8:noexpandtab
// py-cky.h
//
// (c) Mark Johnson, 27th January 2006, last modified 22nd September 2009
// Further additions by Jan Botha (01/2013) to support SRCGs

#ifndef PY_CKY_H
#define PY_CKY_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
// #include <ext/hash_map>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <vector>
#include <tr1/unordered_map>
#include <unordered_map>

#ifndef NOSERIALIZE
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/unordered_map.hpp>
//#include <boost/serialization/unordered_set.hpp>
#endif

#include "gammadist.h"
#include "mt19937ar.h"
#include "slice-sampler.h"
#include "sym.h"
#include "xtree.h"
#include "trie.h"
#include "utility.h"
#include "nested_list.h"
#include "consts.h"

extern int debug;
const mcfg::ushint MAX_RANGEVECTOR_DIM_IMPL = get_MAX_RANGEVECTOR_DIM_IMPL; //max implemented

//! Suppose there are n samples occupying m tables.
//! Then the probability that the n+1 sample occupies
//! table 1 <= k <= m is:
//!
//!  P(x_{n+1} = k) = (n_k - a)/(n + b)
//!
//! and the probability that the n+1 sample occupies
//! the new table m+1 is:
//!
//!  P(x_{n+1} = m+1) = (m*a + b)/(n + b)
//!
//! The probability of a configuration in which a 
//! restaurant contains n customers at m tables,
//! with n_k customers at table k is:
//!
//!
//!  a^{-m} G(m+b/a)  G(b)                 G(n_k-a)
//!         -------- ------  \prod_{k=1}^m --------
//!          G(b/a)  G(n+b)                 G(1-a)
//!
//! where G is the Gamma function.

inline float power(float x, float y) { return powf(x, y); }
inline double power(double x, double y) { return pow(x, y); }
// inline long double power(long double x, long double y) { return powl(x, y); }

typedef double F;
typedef symbol S;
typedef std::vector<S> Ss;

typedef std::map<S,F> S_F;
// typedef ext::hash_map<S,F> S_F;

//typedef mcfg::Mapping::value_type L;
typedef mcfg::Mapping L;
typedef std::map<S, L> S_L;

typedef std::pair<S,Ss> SSs;
typedef std::map<SSs,F> SSs_F;
typedef std::map<SSs,L> SSs_L;

//! readline_symbols() reads all of the symbols on the current
//! line into syms
//
std::istream& readline_symbols(std::istream& is, Ss& syms) {
  syms.clear();
  std::string line;
  if (std::getline(is, line)) {
    std::istringstream iss(line);
    std::string s;
    while (iss >> s)
      syms.push_back(s);
  }
  return is;
}  // readline_symbols()


//! A default_value_type{} object is used to read an object from a stream,
//! assigning a default value if the read fails.  Users should not need to
//! construct such objects, but should use the default_value() function instead.
//
template <typename object_type, typename default_type>
struct default_value_type {
  object_type& object;
  const default_type defaultvalue;
  default_value_type(object_type& object, const default_type defaultvalue)
    : object(object), defaultvalue(defaultvalue) { }
};

//! default_value() is used to read an object from a stream, assigning a
//! default value if the read fails.  It returns a default_value_type{}
//! object, which does the actual reading.
//
template <typename object_type, typename default_type>
default_value_type<object_type,default_type>
default_value(object_type& object, const default_type defaultvalue=default_type()) {
  return default_value_type<object_type,default_type>(object, defaultvalue);
}

//! This operator>>() reads default_value_type{} from an input stream.
//
template <typename object_type, typename default_type>
std::istream& operator>> (std::istream& is, 
			  default_value_type<object_type, default_type> dv) {
  if (is) {
    if (is >> dv.object)
      ;
    else {
      is.clear(is.rdstate() & ~std::ios::failbit);  // clear failbit
      dv.object = dv.defaultvalue;
    }
  }
  return is;
}

// inline F random1() { return rand()/(RAND_MAX+1.0); }
inline F random1() { return mt_genrand_res53(); }

//! A pycfg_type is a CKY parser for a py-srcag
//
struct pycfg_type {

  pycfg_type(F default_weight=1, F default_pya=1, F default_pyb=1) 
    : estimate_theta_flag(false), default_weight(default_weight), 
      default_pya(default_pya), default_pyb(default_pyb),
      pya_beta_a(0), pya_beta_b(0), pyb_gamma_s(0), pyb_gamma_c(0), 
      max_rangevector_dim(MAX_RANGEVECTOR_DIM_IMPL) { }

  typedef mcfg::U U;
  typedef std::pair<U,U> UU;

  typedef std::map<S,U> S_U;

  typedef std::map<S,UU> S_UU;

  typedef tr1::unordered_map<S,S_F> S_S_F; //for unary rules (symbols)

  typedef trie<S, S_F> St_S_F;		   //for n-ary rules (symbols)
  typedef St_S_F::const_iterator Stit;     //to index into an n-ary rule (symbols)

  typedef tr1::unordered_map<S,S_L> S_S_L; //for unary rules (var-mapping)
                                                                                   
  typedef trie<S, S_L> St_S_L;             //for n-ary rules (var-mapping)
  typedef St_S_L::const_iterator SLtit;    //to index into an n-ary rule (var-mapping)

  //typedef catcounttree_type tree;
  typedef catcountyldtree_type tree;

  typedef std::set<tree*> sT;

  //typedef trie<S,sT> St_sT;
  typedef trie<Ss, sT> St_sT;
  typedef ::yld_type Sss;

  typedef std::vector<tree*> Ts;

  typedef std::map<S,Ts> S_Ts;

  //! If estimate_theta_flag is true, then we estimate the generator 
  //! rule weights using a Dirichlet prior
  //
  bool estimate_theta_flag;

  //! start is the start symbol of the grammar
  //
  S start;

  //! rhs_parent_weight maps the right-hand sides of rules
  //! to rule parent and rule weight 
  //
  St_S_F rhs_parent_weight;

  //! unarychild_parent_weight maps unary children to a vector
  //! of parent-weight pairs
  //
  S_S_F unarychild_parent_weight;

  //! parent_weight maps parents to the sum of their rule weights
  //
  S_F parent_weight;

  //! default_weight is the default weight for rules with no explicit
  //! weight.  Used when grammar is read in.
  //
  F default_weight;

  //! rule_priorweight is the prior weight of rule
  //
  SSs_F rule_priorweight;

  //! rhs_parent_varmap stores the variable map for a rule parent->rhs,
  //
//  SSs_L rhs_parent_varmap;
  St_S_L rhs_parent_varmap;
  S_S_L unarychild_parent_varmap; //need to keep them separate

  //! parent_priorweight is the prior weight the parent
  //
  S_F parent_priorweight;

  //! terms_pytrees maps terminal strings to their PY trees
  //
  St_sT terms_pytrees;

  //! parent_pyn maps parents to the number of times they have been expanded
  //
  S_U parent_pyn;

  //! parent_pym maps parents to the number of distinct PY tables for parent
  //
  S_U parent_pym;

  F default_pya;   //!< default value for pya
  F default_pyb;   //!< default value for pyb

  F pya_beta_a;    //!< alpha parameter of Beta prior on pya
  F pya_beta_b;    //!< beta parameter of Beta prior on pya

  F pyb_gamma_s;   //!< s parameter of Gamma prior on pyb
  F pyb_gamma_c;   //!< c parameter of Gamma prior on pyb

  S_F parent_pya;  //!< pya value for parent
  S_F parent_pyb;  //!< pyb value for parent

  mcfg::ushint max_rangevector_dim; //max value to parse; determined through reading grammar

#ifndef NOSERIALIZE
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & estimate_theta_flag;
    ar & start;
    ar & rhs_parent_weight;
    ar & unarychild_parent_weight;
    ar & parent_weight;
    ar & default_weight;
    ar & rule_priorweight;
    ar & rhs_parent_varmap;
    ar & unarychild_parent_varmap;
    ar & parent_priorweight;
    ar & terms_pytrees;
    ar & parent_pyn;
    ar & parent_pym;
    ar & default_pya;
    ar & default_pyb;
    ar & pya_beta_a;
    ar & pya_beta_b;
    ar & pyb_gamma_s;
    ar & pyb_gamma_c;
    ar & parent_pya;
    ar & parent_pyb;
    ar & max_rangevector_dim;

//    std::cerr << "---pycfg_type::serialize()---\n";
//    std::cerr << "estimate_theta_flag" << estimate_theta_flag << std::endl << std::endl;
//    std::cerr << "start" << start << "\n\n";
//    std::cerr << "rhs_parent_weight" << rhs_parent_weight << "\n\n";
//    std::cerr << "unarychild_parent_weight" << unarychild_parent_weight << "\n\n";
//    std::cerr << "parent_weight" << parent_weight << "\n\n";
//    std::cerr << "default_weight" << default_weight << "\n\n";
//    std::cerr << "rule_priorweight" << rule_priorweight << "\n\n";
//    std::cerr << "rhs_parent_varmap" << rhs_parent_varmap << "\n\n";
//    std::cerr << "unarychild_parent_varmap" << unarychild_parent_varmap << "\n\n";
//    std::cerr << "parent_priorweight" << parent_priorweight << "\n\n";
//    std::cerr << "terms_pytrees" << terms_pytrees << "\n\n";
//    std::cerr << "parent_pyn" << parent_pyn << "\n\n";
//    std::cerr << "parent_pym" << parent_pym << "\n\n";
//    std::cerr << "default_pya" << default_pya << "\n\n";
//    std::cerr << "default_pyb" << default_pyb << "\n\n";
//    std::cerr << "pya_beta_a" << pya_beta_a << "\n\n";
//    std::cerr << "pya_beta_b" << pya_beta_b << "\n\n";
//    std::cerr << "pyb_gamma_s" << pyb_gamma_s << "\n\n";
//    std::cerr << "pyb_gamma_c" << pyb_gamma_c << "\n\n";
//    std::cerr << "parent_pya" << parent_pya << "\n\n";
//    std::cerr << "parent_pyb" << parent_pyb << "\n\n";
//    std::cerr << "logpcorpus=" << logPcorpus() << "\n\n";
//    std::cerr << "max_rangevector_dim=" << max_rangevector_dim << "\n\n";
  }
#endif

  //! get_pya() returns the value of pya for this parent
  //
  F get_pya(S parent) const { 
    S_F::const_iterator it = parent_pya.find(parent);
    return (it == parent_pya.end()) ? default_pya : it->second;
  }  // pycfg_type::get_pya()

  //! set_pya() sets the value of pya for this parent, returning
  //! the old value for pya
  //
  F set_pya(S parent, F pya) {
    F old_pya = default_pya;
    S_F::iterator it = parent_pya.find(parent);
    if (it != parent_pya.end())
      old_pya = it->second;
    if (pya != default_pya)
      parent_pya[parent] = pya;
    else // pya == default_pya
      if (it != parent_pya.end())
	parent_pya.erase(it);
    return old_pya;
  }  // pycfg_type::set_pya()

  //! get_pyb() returns the value of pyb for this parent
  //
  F get_pyb(S parent) const { 
    S_F::const_iterator it = parent_pyb.find(parent);
    return (it == parent_pyb.end()) ? default_pyb : it->second;
  }  // pycfg_type::get_pyb()

  //! sum_pym() returns the sum of the pym for all parents
  //
  U sum_pym() const {
    U sum = 0;
    cforeach (S_U, it, parent_pym)
      sum += it->second;
    return sum;
  }  // pycfg_type::sum_pym()

  //! terms_pytrees_size() returns the number of trees in terms_pytrees.
  //
  U terms_pytrees_size() const {
    U size = 0;
    terms_pytrees.for_each(terms_pytrees_size_helper(size));
    return size;
  }  // pycfg_type::terms_pytrees_size()

  struct terms_pytrees_size_helper {
    U& size;
    terms_pytrees_size_helper(U& size) : size(size) { }

    template <typename Words, typename TreePtrs>
    void operator() (const Words& words, const TreePtrs& tps) {
      size += tps.size();
    }  // pycfg_type::terms_pytrees_size_helper::operator()    

  };  // pycfg_type::terms_pytrees_size_helper{}

  //! rule_weight() returns the weight of rule parent --> rhs
  //
  template <typename rhs_type>
  F rule_weight(S parent, const rhs_type& rhs) const {
    assert(!rhs.empty());
    if (rhs.size() == 1) {
      S_S_F::const_iterator it = unarychild_parent_weight.find(rhs[0]);
      if (it == unarychild_parent_weight.end())
	return 0;
      else
	return dfind(it->second, parent);
    }
    else {  // rhs.size() > 1
      Stit it = rhs_parent_weight.find(rhs);
      if (it == rhs_parent_weight.end())
	return 0;
      else
	return dfind(it->data, parent);
    }
  }  // pycfg_type::rule_weight()

  //! rule_prob() returns the probability of rule parent --> rhs
  //
  template <typename rhs_type>
  F rule_prob(S parent, const rhs_type& rhs) const {
    assert(!rhs.empty());
    F parentweight = afind(parent_weight, parent);
    F ruleweight = rule_weight(parent, rhs);
    assert(ruleweight > 0);
    assert(parentweight > 0);
    return ruleweight/parentweight;
  }  // pycfg_type::rule_prob()

  //! tree_prob() returns the probability of the tree under the current
  //! model
  //
  F tree_prob(const tree* tp) const {
    if (tp->children.empty()) 
      return 1;
    F pya = get_pya(tp->cat);
    if (pya == 1) { // no cache
      F prob = 1;
      Ss children;
      cforeach(tree::ptrs_type, it, tp->children) {
	children.push_back((*it)->cat);
	prob *= tree_prob(*it);
      }
      prob *= rule_prob(tp->cat, children);
      return prob;
    }
    F pyb = get_pyb(tp->cat);
    U pym = dfind(parent_pym, tp->cat);
    U pyn = dfind(parent_pyn, tp->cat);
    if (tp->count > 0) { // existing node
      assert(tp->count <= pyn);
      assert(pym > 0);
      F prob = (tp->count - pya)/(pyn + pyb);
      assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
      return prob;
    }
    // new node
    F prob = (pym * pya + pyb)/(pyn + pyb);
    assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
    Ss children;
    cforeach(tree::ptrs_type, it, tp->children) {
      children.push_back((*it)->cat);
      prob *= tree_prob(*it);
    }
    prob *= rule_prob(tp->cat, children);
    if (prob < 0)
      std::cerr << "## pycfg_type::tree_prob(" << *tp << ") = " 
		<< prob << std::endl;
    assert(finite(prob)); assert(prob <= 1); assert(prob >= 0);
    // assert(prob > 0); 
    return prob;
  }  // pycfg_type::tree_prob()

  //! incrrule() increments the weight of the rule parent --> rhs,
  //! returning the probability of this rule under the old grammar.
  //
  template <typename rhs_type>
  F incrrule(S parent, const rhs_type& rhs, F weight = 1) {
    assert(!rhs.empty());
    assert(weight >= 0);
    F& parentweight = parent_weight[parent];
    F parentweight0 = parentweight;
    F rhsweight0;
    parentweight += weight;
    if (rhs.size() == 1) {
      F& rhsweight = unarychild_parent_weight[rhs[0]][parent];
      rhsweight0 = rhsweight;
      rhsweight += weight;
    }
    else {  // rhs.size() > 1
      F& rhsweight = rhs_parent_weight[rhs][parent];
      rhsweight0 = rhsweight;
      rhsweight += weight;
    }
    assert(parentweight0 >= 0);
    assert(rhsweight0 >= 0);
    return rhsweight0/parentweight0;
  }  // incrrule()

  //! decrrule() decrements the weight of rule parent --> rhs,
  //! returning the probability of this rule under the new grammar,
  //! and deletes the rule if it has weight 0.
  //
  template <typename rhs_type>
  F decrrule(S parent, const rhs_type& rhs, F weight = 1) {
    assert(weight >= 0);
    assert(!rhs.empty());
    F rhsweight;
    F parentweight = (parent_weight[parent] -= weight);
    assert(parentweight >= 0);
    if (parentweight == 0)
      parent_weight.erase(parent);
    if (rhs.size() == 1) {
      S_F& parent1_weight = unarychild_parent_weight[rhs[0]];
      rhsweight = (parent1_weight[parent] -= weight);
      assert(rhsweight >= 0);
      if (rhsweight == 0) {
	parent1_weight.erase(parent);
	if (parent1_weight.empty())
	  unarychild_parent_weight.erase(rhs[0]);
      }
    }
    else {  // non-unary rule
      S_F& parent1_weight = rhs_parent_weight[rhs];
      rhsweight = (parent1_weight[parent] -= weight);
      if (rhsweight == 0) {
	parent1_weight.erase(parent);
	if (parent1_weight.empty())
	  rhs_parent_weight.erase(rhs);
      }
    }
    return rhsweight/parentweight;
  }  // pycfg_type::decrrule()

  //! incrtree() increments the cache for tp, increments
  //! the rules if the cache count is appropriate, and returns
  //! the probability of this tree under the original model.
  //
  F incrtree(tree* tp, U weight = 1) {
    if (tp->children.empty())  
      return 1;  // terminal node
    assert(weight >= 0);
    F pya = get_pya(tp->cat);    // PY cache statistics
    F pyb = get_pyb(tp->cat);
    if (pya == 1) { // don't table this category
      F prob = 1;
      {
	Ss children;
	cforeach (tree::ptrs_type, it, tp->children)
	  children.push_back((*it)->cat);
	prob *= incrrule(tp->cat, children, estimate_theta_flag*weight);
      }
      cforeach (tree::ptrs_type, it, tp->children)
	prob *= incrtree(*it, weight);
      return prob;
    }
    else if (tp->count > 0) {  // old PY table entry
      U& pyn = parent_pyn[tp->cat];
      F prob = (tp->count - pya)/(pyn + pyb);
      assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
      tp->count += weight;              // increment entry count
      pyn += weight;                    // increment PY count
      return prob;
    } 
    else { // new PY table entry
      {
	const Sss& terms = tp->terminals();
	assert(terms.size() > 0 && terms[0].size()>0);
	bool inserted = terms_pytrees[terms].insert(tp).second;
	assert(inserted);
      }
      U& pym = parent_pym[tp->cat];
      U& pyn = parent_pyn[tp->cat];
      F prob = (pym*pya + pyb)/(pyn + pyb);  // select new table
      assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
      tp->count += weight;              // increment count
      pym += 1;                         // one more PY table entry
      pyn += weight;                    // increment PY count
      {
	Ss children;
	cforeach (tree::ptrs_type, it, tp->children)
	  children.push_back((*it)->cat);
	prob *= incrrule(tp->cat, children, estimate_theta_flag*weight);
      }
      cforeach (tree::ptrs_type, it, tp->children)
	prob *= incrtree(*it, weight);
      return prob;
    }
  }  // pycfg_type::incrtree()

  //! decrtree() decrements the cache for tp, decrements
  //! the rules if the cache count is appropriate, and returns
  //! the probability of this tree under the new model.
  //
  F decrtree(tree* tp, U weight = 1) {
    if (tp->children.empty())  
      return 1;  // terminal node
    F pya = get_pya(tp->cat);    // PY cache statistics
    if (pya == 1) {  // don't table this category
      F prob = 1;
      {
	Ss children;
	cforeach (tree::ptrs_type, it, tp->children)
	  children.push_back((*it)->cat);
	F ruleprob = decrrule(tp->cat, children, estimate_theta_flag*weight);
	assert(ruleprob > 0);
	prob *= ruleprob;
      }
      cforeach (tree::ptrs_type, it, tp->children) 
	prob *= decrtree(*it, weight);
      return prob;
    }
    assert(weight <= tp->count);
    tp->count -= weight;
    assert(afind(parent_pyn, tp->cat) >= weight);
    const U pyn = (parent_pyn[tp->cat] -= weight);
    F pyb = get_pyb(tp->cat);
    if (tp->count > 0) {  // old PY table entry
      assert(pyn > 0);
      F prob = (tp->count - pya)/(pyn + pyb);
      assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
      return prob;
    } 
    else { // tp->count == 0, remove PY table entry
      {
	const Sss& terms = tp->terminals();
	sT& pytrees = terms_pytrees[terms];
	sT::size_type nerased = pytrees.erase(tp);
	assert(nerased == 1);
	if (pytrees.empty()) 
	  terms_pytrees.erase(terms);
      }
      // Bug: when pym or pyn goes to zero and the parent is erased, 
      // and then the reference to pym or pyn becomes a dangling reference
      // U& pym = parent_pym[tp->cat];
      // pym -= 1;                         // reduce cache count
      assert(parent_pym.count(tp->cat) > 0);
      const U pym = --parent_pym[tp->cat];
      if (pym == 0) 
	parent_pym.erase(tp->cat);
      if (pyn == 0)
	parent_pyn.erase(tp->cat);
      F prob = (pym*pya + pyb)/(pyn + pyb);  // select new table
      assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
      {
	Ss children;
	cforeach (tree::ptrs_type, it, tp->children)
	  children.push_back((*it)->cat);
	prob *= decrrule(tp->cat, children, estimate_theta_flag*weight);
      }
      assert(prob > 0);
      cforeach (tree::ptrs_type, it, tp->children)
	prob *= decrtree(*it, weight);
      // assert(prob > 0);
      return prob;
    }
  }  // pycfg_type::decrtree()

  //! read() reads a grammar from an input stream (implements >> )
  //
  std::istream& read(std::istream& is) {
    start = symbol::undefined();
    F weight;
    F pya;
    F pyb;
    S parent;
    L varmap;
    mcfg::ushint max_rvec_dim_sofar = 0; 
    while (is >> default_value(weight, default_weight) 
	      >> default_value(pya, default_pya)
	      >> default_value(pyb, default_pyb)
	      >> parent 
	      >> varmap
	      >> " -->") {
      if (weight<=0)
	weight=default_weight;
      if (start.is_undefined())
	start = parent;
      Ss rhs;
      readline_symbols(is, rhs);
      for (auto it = rhs.begin(); it != rhs.end(); ++it)
	if (it->string_reference()[0] == '#') {
	  if (debug >= 200000) std::cerr << "Ignoring comment in rhs-string'" << rhs << "'\n";
	  rhs.erase(it, rhs.end());
	  break;
	}
      if (debug >= 100000)
	std::cerr << "# " << weight << '\t' << parent << " --> " << rhs << std::endl;
      mcfg::intpair maxdims = mcfg::get_max_pred_dim(varmap);
      max_rvec_dim_sofar = std::max(max_rvec_dim_sofar, std::max(maxdims.first, maxdims.second));
      if (max_rvec_dim_sofar > MAX_RANGEVECTOR_DIM_IMPL) {
	std::cerr << "Error: Implementation can only handle up to " << MAX_RANGEVECTOR_DIM_IMPL 
		  << " arguments in a predicate, but the rule " << parent << " ---> " << rhs 
		  << " involves one with " << max_rvec_dim_sofar << " arguments. varmap="
		  << varmap << std::endl;
	abort();
      }
      incrrule(parent, rhs, weight);
      if (pya != default_pya)
	parent_pya[parent] = pya;
      if (pyb != default_pyb)
	parent_pyb[parent] = pyb;
      rule_priorweight[SSs(parent,rhs)] += weight;
      if (rhs.size() == 2) {
	auto& v = rhs_parent_varmap[rhs];
	assert(v.find(parent) == v.end());
	rhs_parent_varmap[rhs][parent] = varmap; //should (parent,rhs) designate unique varmap
      }
      else if (rhs.size() == 1) {
        auto& v = unarychild_parent_varmap[rhs[0]];
	assert(v.find(parent) == v.end());
	unarychild_parent_varmap[rhs[0]][parent] = varmap;;
      } 
      else {
        std::cerr << "Rules can only be unary or binary!" << std::endl;
        std::cerr << "# " << weight << '\t' << parent << " --> " << rhs << std::endl;
	abort();
      }
     if (debug >= 100000)
        std::cerr << "unarychild_parent_varmap = " << unarychild_parent_varmap << std::endl
	       	<< "rhs_parent_varmap = " << rhs_parent_varmap << std::endl;

      parent_priorweight[parent] += weight;
    }
    assert(max_rvec_dim_sofar <= MAX_RANGEVECTOR_DIM_IMPL);
    assert(max_rvec_dim_sofar > 0);
    max_rangevector_dim = max_rvec_dim_sofar; // lower the limit down to what is actually needed for this grammar
    if (debug >= 1000) std::cerr << "Determined from grammar that max_rvec_dim = " << max_rangevector_dim << "\n";
    return is;
  }  // pycfg_type::read()

  //! write() writes a grammar (implements << )
  //
  std::ostream& write(std::ostream& os) const {
    assert(start.is_defined());
    write_rules(os, start);
    cforeach (S_F, it, parent_weight)
      if (it->first != start) 
	write_rules(os, it->first);
    return os;
  }  // pycfg_type::write()

  std::ostream& write_rules(std::ostream& os, S parent) const {
    rhs_parent_weight.for_each(write_rule(os, parent));
    cforeach (S_S_F, it0, unarychild_parent_weight) {
      S child = it0->first;
      const S_F& parent_weight = it0->second;
      cforeach (S_F, it1, parent_weight)
	if (it1->first == parent)
	  os << it1->second << '\t' << parent 
	     << " --> " << child << std::endl;
    }
    bool old_compact_trees_flag = catcountyldtree_type::compact_trees;  // save old flag
    catcountyldtree_type::compact_trees = false;  // turn off compact_trees
    terms_pytrees.for_each(write_pycache(os, parent));
    catcountyldtree_type::compact_trees = old_compact_trees_flag;
    return os;
  }  // pycfg_type::write_rules()

  //! write_rule{} writes a single rule
  //
  struct write_rule {
    std::ostream& os;
    S parent;

    write_rule(std::ostream& os, symbol parent) : os(os), parent(parent) { }

    template <typename Keys, typename Value>
    void operator() (const Keys& rhs, const Value& parentweights) {
      cforeach (typename Value, pwit, parentweights) 
	if (pwit->first == parent) {
	  os << pwit->second << '\t' << parent << " -->";
	  cforeach (typename Keys, rhsit, rhs)
	    os << ' ' << *rhsit;
	  os << std::endl;
	}
    }  // pycfg_type::write_rule::operator()

  };  // pycfg_type::write_rule{}
  
  //! write_pycache{} writes the cache entries for a category
  //
  struct write_pycache {
    std::ostream& os;
    S parent;
    
    write_pycache(std::ostream& os, S parent) : os(os), parent(parent) { }

    template <typename Words, typename TreePtrs>
    void operator() (const Words& words, const TreePtrs& tps) {
      cforeach (typename TreePtrs, tpit, tps) 
	if ((*tpit)->cat == parent)
	  os << (*tpit) << std::endl;
    }  // pycfg_type::write_pycache::operator()
  };  // pycfg_type::write_pycache{}

  //! logPcorpus() returns the log probability of the corpus trees
  //
  F logPcorpus() const {
    F logP = 0;
    // grammar part
    cforeach (SSs_F, it, rule_priorweight) {
      S parent = it->first.first;
      const Ss& rhs = it->first.second;
      F priorweight = it->second;
      F weight = rule_weight(parent, rhs);
      logP += lgamma(weight) - lgamma(priorweight);
    }
    cforeach (S_F, it, parent_priorweight) {
      S parent = it->first;
      F priorweight = it->second;
      F weight =dfind(parent_weight, parent);
      logP += lgamma(priorweight) - lgamma(weight);
    }
    // PY adaptor part
    cforeach (S_U, it, parent_pyn) {
      S parent = it->first;
      U pyn = it->second;
      U pym = afind(parent_pym, parent);
      F pya = get_pya(parent);
      F pyb = get_pyb(parent);
      logP += lgamma(pyb) - lgamma(pyn+pyb);
      for (U i = 0; i < pym; ++i)
	logP += log(i*pya + pyb);
    }
    terms_pytrees.for_each(logPcache(*this, logP));
    return logP;
  }  // pycfg_type::logPcorpus()

  struct logPcache {
    const pycfg_type& g;
    F& logP;

    logPcache(const pycfg_type& g, F& logP) : g(g), logP(logP) { }

    template <typename Words, typename TreePtrs>
    void operator() (const Words& words, const TreePtrs& tps) {
      cforeach (typename TreePtrs, it, tps) {
	S parent = (*it)->cat;
	U count = (*it)->count;
	F pya = g.get_pya(parent);
	logP += lgamma(count-pya) - lgamma(1-pya);
      }
    }  // pycfg_type::logPcache::operator()
  };  // pycfg_type::logPcache{}

  //! logPrior() returns the prior probability of the PY a and b values
  //
  F logPrior() const {
    F logP = 0;
    if (pyb_gamma_s > 0 && pyb_gamma_c > 0)
      cforeach (S_U, it, parent_pyn) {
	S parent = it->first;
	F pya = get_pya(parent);
	F pyb = get_pyb(parent);
	if (pya_beta_a > 0 && pya_beta_b > 0 && pya > 0)
	  logP += pya_logPrior(pya, pya_beta_a, pya_beta_b);
	logP += pyb_logPrior(pyb, pyb_gamma_c, pyb_gamma_s);
      }
    return logP;
  }  // pycfg_type::logPrior()

  //! pya_logPrior() calculates the Beta prior on pya.
  //
  static F pya_logPrior(F pya, F pya_beta_a, F pya_beta_b) {
    F prior = lbetadist(pya, pya_beta_a, pya_beta_b);     //!< prior for pya
    return prior;
  }  // pycfg_type::pya_logPrior()

  //! pyb_logPrior() calculates the prior probability of pyb 
  //! wrt the Gamma prior on pyb.
  //
  static F pyb_logPrior(F pyb, F pyb_gamma_c, F pyb_gamma_s) {
    F prior = lgammadist(pyb, pyb_gamma_c, pyb_gamma_s);  // prior for pyb
    return prior;
  }  // pcfg_type::pyb_logPrior()

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //                      Resample pyb                                //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////
  
  //! resample_pyb_type{} is a function object that returns the part of log prob that depends on pyb.
  //! This includes the Gamma prior on pyb, but doesn't include e.g. the rule probabilities
  //! (as these are a constant factor)
  //
  struct resample_pyb_type {
    U pyn, pym;
    F pya, pyb_gamma_c, pyb_gamma_s;
    resample_pyb_type(U pyn, U pym, F pya, F pyb_gamma_c, F pyb_gamma_s) 
      : pyn(pyn), pym(pym), pya(pya), pyb_gamma_c(pyb_gamma_c), pyb_gamma_s(pyb_gamma_s)
    { }

    //! operator() returns the part of the log posterior probability that depends on pyb
    //
    F operator() (F pyb) const {
      F logPrior = pyb_logPrior(pyb, pyb_gamma_c, pyb_gamma_s);  //!< prior for pyb
      F logProb = 0;
      logProb += (pya == 0 ? pym*log(pyb) : pym*log(pya) + lgamma(pym + pyb/pya) - lgamma(pyb/pya));
      logProb += lgamma(pyb) - lgamma(pyn+pyb);
      return logProb+logPrior;
    }
  };  // pcfg_type::resample_pyb_type{}

  //! resample_pyb() samples new values for pyb for each adapted nonterminal.
  //! It returns the log prior prob of new b values.
  //
  void resample_pyb() {
    U niterations = 20;   //!< number of resampling iterations
    // std::cerr << "\n## resample_pyb(), initial parent_pya = " << parent_pya << ", parent_pyb = " << parent_pyb << std::endl;
    cforeach (S_U, it, parent_pyn) {
      S parent = it->first;
      U pyn = it->second;
      U pym = afind(parent_pym, parent);
      F pya = get_pya(parent);
      F pyb = get_pyb(parent);
      resample_pyb_type pyb_logP(pyn, pym, pya, pyb_gamma_c, pyb_gamma_s);
      pyb = slice_sampler1d(pyb_logP, pyb, random1, 0.0, std::numeric_limits<F>::infinity(), 0.0, niterations, 100*niterations);
      parent_pyb[parent] = pyb;
      // parent_bap[parent].first += naccepted;
      // parent_bap[parent].second += nproposed;
    }
  }  // pcfg_type::resample_pyb()

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //                   Resample pya and pyb                           //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

  //! resample_pya_type{} calculates the part of the log prob that depends on pya.
  //! This includes the Beta prior on pya, but doesn't include e.g. the rule probabilities
  //! (as these are a constant factor)
  //
  struct resample_pya_type {
    U pyn, pym;
    F pyb, pya_beta_a, pya_beta_b;
    const Ts& trees;
    
    resample_pya_type(U pyn, U pym, F pyb, F pya_beta_a, F pya_beta_b, const Ts& trees) 
      : pyn(pyn), pym(pym), pyb(pyb), pya_beta_a(pya_beta_a), pya_beta_b(pya_beta_b), trees(trees)
    { }

    //! operator() returns the part of the log posterior probability that depends on pya
    //
    F operator() (F pya) const {
      F logPrior = pya_logPrior(pya, pya_beta_a, pya_beta_b);     //!< prior for pya
      F logProb = 0;
      F lgamma1a = lgamma(1-pya);
      cforeach (Ts, it, trees) {
	U count = (*it)->count;
	logProb += lgamma(count-pya) - lgamma1a;
      }
      logProb += (pya == 0 ? pym*log(pyb) : pym*log(pya) + lgamma(pym + pyb/pya) - lgamma(pyb/pya));
      return logPrior + logProb;
    }   // pycfg_type::resample_pya_type::operator()

  };  // pycfg_type::resample_pya_type{}
  
  //! resample_pya() samples new values for pya for each adapted nonterminal
  //
  void resample_pya(const S_Ts& parent_trees) {
    U niterations = 20;   //!< number of resampling iterations
    // std::cerr << "\n## Initial parent_pya = " << parent_pya << ", parent_pyb = " << parent_pyb << std::endl;
    cforeach (S_U, it, parent_pyn) {
      S parent = it->first;
      F pya = get_pya(parent);
      if (pya == 0)   // if this nonterminal has pya == 0, then don't resample
	continue;
      F pyb = get_pyb(parent);
      U pyn = it->second;
      U pym = afind(parent_pym, parent);
      const Ts& trees = afind(parent_trees, parent);
      resample_pya_type pya_logP(pyn, pym, pyb, pya_beta_a, pya_beta_b, trees);
      pya = slice_sampler1d(pya_logP, pya, random1, std::numeric_limits<F>::min(), 1.0, 0.0, niterations, 100*niterations);
      parent_pya[parent] = pya;
    }
  }  // pycfg_type::resample_pya()

  //! resample_pyab_parent_trees_helper{} constructs parent_trees from terms_pytrees.
  //
  struct resample_pyab_parent_trees_helper {
    S_Ts& parent_trees;
    resample_pyab_parent_trees_helper(S_Ts& parent_trees) : parent_trees(parent_trees) { }

    template <typename Words, typename TreePtrs>
    void operator() (const Words& words, const TreePtrs& tps) {
      cforeach (typename TreePtrs, it, tps) {
	S parent = (*it)->cat;
	parent_trees[parent].push_back(*it);
      }
    }  // pycfg_type::resample_pyab_parent_trees_helper::operator()
  };  // pycfg_type::resample_pyab_parent_trees_helper{}

  //! resample_pyab() resamples both pya and pyb for each adapted nonterminal.
  //
  void resample_pyab() {
    const U niterations = 5;  //!< number of alternating samples of pya and pyb
    S_Ts parent_trees;
    terms_pytrees.for_each(resample_pyab_parent_trees_helper(parent_trees));
    for (U i=0; i<niterations; ++i) {
      resample_pyb();
      resample_pya(parent_trees);
    }
    resample_pyb();
  }  // pycfg_type::resample_pyab()

  //! write_adaptor_parameters() writes out adaptor parameters to a file
  //
  std::ostream& write_adaptor_parameters(std::ostream& os) const {
    cforeach (S_F, it, parent_priorweight) {
      S parent = it->first;
      F pya = get_pya(parent);
      if (pya == 1)
	continue;
      U pym = dfind(parent_pym, parent);
      U pyn = dfind(parent_pyn, parent);
      F pyb = get_pyb(parent);
      os << ' ' << parent << ' ' << pym << ' ' << pyn << ' ' << pya << ' ' << pyb;
    }
    return os;
  }  // pycfg_type::write_adaptor_parameters()

};  // pycfg_type{}

//! operator>> (pycfg_type&) reads a pycfg_type g, setting g.start
//! to the parent of the first rule read.
//
std::istream& operator>> (std::istream& is, pycfg_type& g) {
  return g.read(is);
}  // operator>> (pycfg_type&)


std::ostream& operator<< (std::ostream& os, const pycfg_type& g) {
  return g.write(os);
}  // operator<< (pycfg_type&)

namespace std { namespace tr1 {
    template <> struct hash<pycfg_type::Stit> 
      : public std::unary_function<pycfg_type::Stit, std::size_t> {
      size_t operator()(const pycfg_type::Stit t) const
      {
	return size_t(&(*t));
      }  // ext::hash<pycfg_type::Stit>::operator()
    };  // ext::hash<pycfg_type::Stit>{}
  }  } // namespace std::tr1


typedef std::initializer_list<pycfg_type::U> rvec_lt; //range vector 'lite'

std::ostream& ilout (std::ostream& os, const rvec_lt& values) {
  for (auto i = begin(values); i<end(values); ++i)
	os << (i==begin(values) ? "" : ",") << *i;
  return os;
}

static const F unaryclosetolerance = 1e-7;

class pycky {
private:
  //max value to parse with pycky. set upon construction.
  //Having it as a static member makes access simpler.
  static mcfg::ushint max_rangevector_dim; 
public:

  const pycfg_type& g;
  F anneal;         // annealing factor (1 = no annealing)
  
  pycky(const pycfg_type& g, F anneal=1) : g(g), anneal(anneal) {
    max_rangevector_dim = g.max_rangevector_dim;
  }

  typedef pycfg_type::tree tree;
  typedef pycfg_type::U U;
  typedef mcfg::R R;
  typedef pycfg_type::S_S_F S_S_F;
  typedef pycfg_type::St_S_F St_S_F;
  typedef pycfg_type::Stit Stit;

  typedef pycfg_type::S_S_L S_S_L;
  typedef pycfg_type::St_S_L St_S_L;
  typedef pycfg_type::SLtit SLtit;

  typedef std::unordered_map<R, S_F> S_Fs;

  // typedef ext::hash_map<Stit,F> Stit_F;
  typedef tr1::unordered_map<Stit,F> Stit_F;
  typedef std::unordered_map<R, Stit_F> Stit_Fs;

  typedef tr1::unordered_map<Stit,SLtit> Stit2SLtit; //associate pointers into separate tries
  typedef std::unordered_map<R, Stit2SLtit> SLtit_Us;
  typedef pycfg_type::sT sT;

  typedef pycfg_type::St_sT St_sT;
  typedef St_sT::const_iterator StsTit;
  typedef std::unordered_map<R, StsTit> StsTits;

  typedef pycfg_type::Sss Sss;
  //! Creates a hash of the values, which are interpreted as a sequence of ranges
  inline static R index(const rvec_lt& values) {
    return index(values.begin(), values.size());
  }
  inline static R index(const U* values, U size) {
    return mcfg::index_bitvec<U,R,get_MAX_INPUT_LEN>(values, size);
  }

  // Extract the terminals covered by rvec. (May return gapped string, in other words.)
  // Assumes that rvec is sanely ordered
  Sss get_terminals(const rvec_lt& rvec) const{
    Sss ret;
    assert(rvec.size() % 2 == 0);
    for (auto ip=begin(rvec); ip<end(rvec); ip+=2)
    {
	auto i = *ip, j = *(ip+1);
	Ss span; span.reserve(j-1);
	span.insert(span.end(), terminals.begin()+i, terminals.begin()+j);
	ret.push_back(span);
    }
    return ret;
  }

  // --- Parser State - reset upon each new call to inside() ----
  Ss terminals;
  S_Fs inactives;

  Stit_Fs actives;
  // Example
  // actives[4,5] = ((0x2446f48, 0.142857) = {(NP):((VP 1)); (NP NP):((VP 1)); } 
  //                 (0x2446b38, 0.1) = {(N):((NP 1)); })
  // Each line encodes multiple dotted items sharing the same sequence of pre-dot RHS symbols (not stored)
  // The initial pair is (mem-address, prob) and stuff after {} show the trie at that address.
  // prob is the inside-prob so-far for the related items.
  // (NP):((VP 1)) means if the item combines with an NP, we'd get a passive item with cat VP
  //               and the ruleweight of the rule applied is 1
  // (NP NP):((VP 1)) means the item requires two NPs to complete a VP
  //
  
  SLtit_Us active_vms;  //similar structure as above, but each item in a cell is a pair (a,b)
			//where 'a' is a pointer into the rule trie of the grammar (i.e. an "active pointer")
			//and 'b' is a pointer into the varmap trie of the grammar.
			//a trienode in 'b' then stores a map (parent -> varmap-of-implied-rule)

  StsTits pytits; // collection of iterators (one per cell)
	          // each one points into a (symbol-sequence, set-of-trees)
                  // i.e. associates the terminal string of that cell with the trees that yield it

  /////// end of Parser State

  // Structure to aggregate things that are relevant during top-down sampling.
  // Objective is largely to simplify method signatures.
  struct TopDownInfo {
    const S parent;
    const rvec_lt& parent_rv;
    const F rulefactor;
    F& probsofar;
    const F probthreshold;
    tree* const tp;

    const pycky* obj;
    TopDownInfo(S parent, const rvec_lt& parent_rv, F rulefactor,
		  F& probsofar, F probthreshold, tree* const tp, const pycky* obj ) :
	    parent(parent), parent_rv(parent_rv), rulefactor(rulefactor),
	    probsofar(probsofar), probthreshold(probthreshold), tp(tp),
	    obj(obj)
	    {}

    //==Functions to be called by handle_case() when sampling a tree top-down from chart
    //for binary rules
    inline bool operator()(const rvec_lt& active_rv, const rvec_lt& passive_rv) {
      return obj->random_binary_passive_do(*this, active_rv, passive_rv);
    }
    //for concatenating unary rules
    inline bool operator()(const rvec_lt& unarychild_rv) {
      return obj->random_unary_passive(*this, unarychild_rv);
    }
  };

  // Functions to be called by handle_case() when parsing to populate the inside chart
  struct InsideFunctor {
    InsideFunctor(const rvec_lt& parent_rv, pycky* obj) : parent_rv(parent_rv), obj(obj) {}

    const rvec_lt& parent_rv;
    pycky* obj;

    //Both functions wrap void methods, so just return false, because
    //our caller would interpret true as 'break'.

    //for binary rules
    inline bool operator()(const rvec_lt& active_rv, const rvec_lt& passive_rv) {
      obj->inside_visit_factorisation(parent_rv, active_rv, passive_rv);
      return false;
    }
    //for concatenating unary rules
    inline bool operator()(const rvec_lt& unarychild_rv) {
      obj->inside_visit_unary_concat(parent_rv, unarychild_rv); 
      return false;
    }
  };

  //save some typing:
  #define truewrap(expr) if (expr) return true

  // Key logic for handling binary sub-cases for a particular parent range vector.
  // f is a struct whose operator() method is called with arguments
  //  - active_rvec, passive_rvec, for each binary rule case
  //  - unarychild_rvec for each factorisation of a contiguous result set into concatenated ranges & unary rule
  // The first case that f returns true for, causes execution to return
  // to the caller by returning true from this method.
  // If no case triggered that, this method returns false.
  // The relevant sub-cases to consider are deduced from the contents/dim of 
  // f.parent_rv.
  template<typename Functor>
  static bool handle_subcases(Functor& f) {
    // for a rule 'parent -> B C'
    // (rlen := number of ranges in predicate)
    // options to ./permute_ranges.py  -p 'truewrap(f(' -s'));'
    assert(f.parent_rv.size() % 2 == 0);
    int parent_num_ranges = f.parent_rv.size() / 2;
    assert(parent_num_ranges  > 0);
    U a = mcfg::at(f.parent_rv, 0), z = mcfg::at(f.parent_rv, f.parent_rv.size()-1);
    if (parent_num_ranges > MAX_RANGEVECTOR_DIM_IMPL) {
      std::cerr << "pycky::handle_subcases() parent_num_ranges =" << parent_num_ranges  
      	  << " not supported; range vector was ";
      ilout(std::cerr, f.parent_rv) << ", maximum supported is " << MAX_RANGEVECTOR_DIM_IMPL << std::endl;
      abort();
    }
    if (parent_num_ranges  == 1) { //////////////////////////////////////////////////////
      // lines ending with *unary concat handle rules like A(abc) -> X3(a,b,c);
      // only added the ones here for up to 5 rhs ranges. cf get_MAX_RANGEVECTOR_DIM_IMPL.
      // NB unary rules involving a single range are still handled in unary_close
      for (U b = a+1; b < z; ++b) {
	if (max_rangevector_dim >= 2) truewrap(f({a,b,b,z})); //*unary concat
	// Case where both RHSs are single ranges
	// -- rlen(RHS1) = 1, rlen(RHS2) = 1; ./permute... abbz 1
  	truewrap(f({a,b}, {b,z}));

	// Cases where gapped RHSs combine into single parent range
	// == 3 rhs ranges
  	for (U c = b+1; c < z && max_rangevector_dim >= 2; ++c) {
	  if (max_rangevector_dim >= 3) truewrap(f({a,b,b,c,c,z})); //*unary concat
	  // -- rlen(RHS1) = 2, rlen(RHS2) = 1
	  truewrap(f({a,b,c,z}, {b,c}));

	  // == 4 rhs ranges
	  for (U d = c+1; d < z; ++d) {
	    if (max_rangevector_dim >= 4) truewrap(f({a,b,b,c,c,d,d,z})); //*unary concat
	    // -- rlen(RHS1) = 2, rlen(RHS2) = 2
	    truewrap(f({a,b,c,d}, {b,c,d,z}));

	    // == 5 rhs ranges
	    for (U e = d+1; e < z && max_rangevector_dim >= 3; ++e) {
	      if (max_rangevector_dim >= 5) truewrap(f({a,b,b,c,c,d,d,e,e,z})); //*unary concat
	      // -- rlen(RHS1) = 3, rlen(RHS2) = 2
	      truewrap(f({a,b,c,d,e,z}, {b,c,d,e}));
	      //== 6 rhs ranges
	      for (U g = e+1; g < z; ++g) {
	        // -- rlen(RHS1) = 3, rlen(RHS2) = 3
		truewrap(f({a,b,c,d,e,g}, {b,c,d,e,g,z}));

		//== 7 rhs ranges
		for (U h = g+1; h < z && max_rangevector_dim >= 4; ++h) {
	          // -- rlen(RHS1) = 4, rlen(RHS2) = 3
	          truewrap(f({a,b,c,d,e,g,h,z}, {b,c,d,e,g,h}));
		  
		  //== 8 rhs ranges
		  for (U i = h+1; i < z; ++i) {
	            // -- rlen(RHS1) = 4, rlen(RHS2) = 4
		    truewrap(f({a,b,c,d,e,g,h,i}, {b,c,d,e,g,h,i,z}));
		   
		    //== 9 rhs ranges
		    for (U j = i+1; j < z && max_rangevector_dim >= 5; ++j) {
	              // -- rlen(RHS1) = 5, rlen(RHS2) = 4
	              truewrap(f({a,b,c,d,e,g,h,i,j,z}, {b,c,d,e,g,h,i,j}));
		    }
		  }
		}
	      }
	    }
	  }
  	}
      }
    }
    //the oddly nested ifs below is for reusing the values read from parent_rv (std::initializer_list doesn't slice nicely)
    else if (parent_num_ranges >= 2) { ////////////////////////////////////////////////////////
      U b = mcfg::at(f.parent_rv,1);
      U c = mcfg::at(f.parent_rv,2);
      if (parent_num_ranges  == 2) {
        truewrap(f({a,b}, {c,z})); // A(x,y) -> B(x) C(y); rlen(RHS1) = rlen(RHS2) = 1
	// partially concatenating rules
	if (max_rangevector_dim >= 3) {
	  for (U ai=a+1; ai<b; ++ai)  { //subdivide first range
	    truewrap(f({a,ai,ai,b,c,z})); // A(xy,z) -> B(x,y,z)
	  }
	  for (U ci=c+1; ci<z; ++ci)  { //subdivide second range
	    truewrap(f({a,b,c,ci,ci,z})); // A(x,yz) -> B(x,y,z)
	  }
	  if (max_rangevector_dim >= 4) { //subdivide both ranges
	    for (U ai=a+1; ai<b; ++ai)  
	      for (U ci=c+1; ci<z; ++ci)  
	        truewrap(f({a,ai,ai,b,c,ci,ci,z})); // A(xq,yz) -> B(x,q,y,z)
	  }
	}
      }
      else if (parent_num_ranges >= 3) {
        U d = mcfg::at(f.parent_rv,3);
        U e = mcfg::at(f.parent_rv,4);
        if (parent_num_ranges  == 3) {
          truewrap(f({a,b,c,d}, {e,z})); // A(x,y,z) -> B(x,y) C(z); rlen(RHS1) 2; rlen(RHS2) = 1
	  // partially concatenating rules
	  if (max_rangevector_dim >= 4) {
	    for (U ai=a+1; ai<b; ++ai)  { //subdivide first range
	      truewrap(f({a,ai,ai,b,c,d,e,z})); // A(xq,y,z) -> B(x,q,y,z)
	    }
	    for (U ci=c+1; ci<d; ++ci)  { //subdivide second range
	      truewrap(f({a,b,c,ci,ci,d,e,z})); // A(x,qy,z) -> B(x,q,y,z)
	    }
	    for (U ei=e+1; ei<z; ++ei)  { //subdivide last range
	      truewrap(f({a,b,c,d,e,ei,ei,z})); // A(x,q,yz) -> B(x,q,y,z)
	    }
	  }
	}
        else if (parent_num_ranges >= 4) {
          U ff = mcfg::at(f.parent_rv,5);
          U g = mcfg::at(f.parent_rv,6);
          if (parent_num_ranges  == 4) {
            truewrap(f({a,b,c,d,e,ff}, {g,z})); // A(p,x,y,z) -> B(p,x,y) C(z); rlen(RHS1) 3; rlen(RHS2) = 1
	  }
          else if (parent_num_ranges >= 5) {
            U h = mcfg::at(f.parent_rv,7);
            U i = mcfg::at(f.parent_rv,8);
            if (parent_num_ranges  == 5) {
              truewrap(f({a,b,c,d,e,ff,g,h}, {i,z})); // A(q,p,x,y,z) -> B(q,p,x,y) C(z); rlen(RHS1) 4; rlen(RHS2) = 1
	    }
          }
        }
      }
    }
    return false;
  }



  //! inside() constructs the inside table, and returns the probability
  //! of the start symbol rewriting to the terminals.
  //
  template <typename terminals_type>
  F inside(const terminals_type& terminals0, S start) {

    terminals = terminals0;

    if (debug >= 10000)
      std::cerr << "# cky::inside() terminals = " << terminals << std::endl;
    if (debug >= 20000) 
	std::cerr << "Parsing with max_rvec_dim = " << max_rangevector_dim << std::endl;

    U n = terminals.size();

    //one 'inactive' per cell. each inactive maps grammar symbol to inside prob of that sym spanning cell
    inactives.clear();
    //one 'active' per cell. basically dot-based incremental parsing
    
    //one 'active' *associates* each sequence of RHS rule symbols covered so far for the span
    //*with* an { S-> F}, where each S is a next rule symbol and F the prob of the item so far
    actives.clear();

    //holds the var mappings of the rules being worked towards and how far an item is into that mapping
    active_vms.clear();

    //holds cached trees applicable to each cell. indexed by cell then terminal string
    pytits.clear();

    inside_terminals(terminals); //parse single terminals

    //variable 'gap' is an artefact of the original code; it does not relate to gaps
    //in constituent ranges
    for (U gap = 2; gap <= n; ++gap) { // non-terminals. e.g. rule shape A -> B C
      for (U left = 0; left + gap <= n; ++left) {
	U right = left + gap;

        for (U c = left+1; c < right && max_rangevector_dim >= 2 ; ++c) { // c and d are indexing 
          for (U d = c; d < right; ++d) {                                 // for parent_ranges_2 
	    
	    for (U e = d+1; e < right && max_rangevector_dim >= 3; ++e) { // e and f are indexing
	      for (U f = e; f < right; ++f) {	                          // for parent_ranges_3

                for (U g = f+1; g < right && max_rangevector_dim >= 4; ++g) { // g and h are indexing
                  for (U h = g; h < right; ++h) {                             // for parent_ranges_4
  
                    for (U i = h+1; i < right && max_rangevector_dim >= 5; ++i) { // i and j are indexing
                      for (U j = i; j < right; ++j) {                             // for parent_ranges_5
                        //Deal with cases where parent rangevec has 5 ranges
 	                auto parent_ranges_5 = {left,c,d,e,f,g,h,i,j,right};
     	                if (debug >= 20100) ilout(std::cerr<< "++(",parent_ranges_5) << ")++" << std::endl;
 	                InsideFunctor func(parent_ranges_5, this);
 	                handle_subcases(func);
     	                StsTit& pytit = get_pytit<terminals_type>(parent_ranges_5); 
     	                inside_finalise_cell<terminals_type>(parent_ranges_5, pytit);
                      }
                    }

                    //Deal with cases where parent rangevec has 4 ranges
 	            auto parent_ranges_4 = {left,c,d,e,f,g,h,right};
     	            if (debug >= 20100) ilout(std::cerr<< "++(",parent_ranges_4) << ")++" << std::endl;
 	            InsideFunctor func(parent_ranges_4, this);
 	            handle_subcases(func);
     	            StsTit& pytit = get_pytit<terminals_type>(parent_ranges_4); 
     	            inside_finalise_cell<terminals_type>(parent_ranges_4, pytit);
                  }
                }

	       //Deal with cases where parent rangevec has 3 ranges
 	        auto parent_ranges_3 = {left,c,d,e,f,right};
     	        if (debug >= 20100) ilout(std::cerr<< "++(",parent_ranges_3) << ")++" << std::endl;
 	        InsideFunctor func(parent_ranges_3, this);
 	        handle_subcases(func);
     	        StsTit& pytit = get_pytit<terminals_type>(parent_ranges_3); 
     	        inside_finalise_cell<terminals_type>(parent_ranges_3, pytit);
	      }
	    }

	    //Deal with cases where parent rangevec has 2 ranges
	    auto parent_ranges_2 = {left,c,d,right};
    	    if (debug >= 20100) ilout(std::cerr<< "++(",parent_ranges_2) << ")++" << std::endl;
	    InsideFunctor func(parent_ranges_2, this);
	    handle_subcases(func);
    	    StsTit& pytit = get_pytit<terminals_type>(parent_ranges_2); 
    	    inside_finalise_cell<terminals_type>(parent_ranges_2, pytit);
          }
        }
	//Deal with cases where parent rangevec has 1 range
	auto parent_ranges_1 = {left,right};
    	  if (debug >= 20100) ilout(std::cerr<< "++(",parent_ranges_1) << ")++" << std::endl;
	// pytit will keep set of trees yielding the whole result (left,right)
	InsideFunctor func(parent_ranges_1, this);
	handle_subcases(func);
	StsTit& pytit = get_pytit<terminals_type>(parent_ranges_1);
	inside_finalise_cell<terminals_type>(parent_ranges_1, pytit);
      }
    }
    return dfind(inactives[index({0,n})], start);
  }  // pycky::inside()

  template <typename terminals_type>
  F inside(const terminals_type& terminals) {
    return inside(terminals, g.start);
  }

  template <typename terminals_type>
  void inside_terminals(const terminals_type& terminals) {
    U n = terminals.size();
    for (U i = 0; i < n; ++i) {   // terminals (single ones)
      //tree-set for this cell is exactly the cached trees currently yielding terminal token i
      pytits[index({i, i+1})] = g.terms_pytrees.find1({terminals[i]});  // PY cache
      //init inside probability for the terminal ( P(term=term)=1)
      inactives[index({i,i+1})][terminals[i]] = 1;
      StsTit& pytit = pytits[index({i,i+1})]; //pytit points to treeset for this cell
      if (pytit != g.terms_pytrees.end()) 
	//if there are such trees yet, add the discounted counts (nk-a)/(n+b) to inside probs
	add_pycache(pytit->data, inactives[index({i,i+1})]);
      inside_unaryclose(inactives[index({i,i+1})], actives[index({i,i+1})], active_vms[index({i,i+1})],
		        {i,i+1});
      
      debug_print_cell_state({i,i+1});
    }
  }
  
  template<typename terminals_type>
  StsTit& get_pytit(const rvec_lt& parent_rv) {

    StsTit& pytit = pytits[index(parent_rv)]; //this quite likely causes an insert in the chart too - want that.
    
    // There might be a more efficient way than this full lookup.
    // For the CFG case, a cell that will already have been visited by the time this is called,
    // e.g. (left,right-1) when parent_rv={left,right},
    // is used as starting point for a call to find1, with just the terminal[right-1] as arg.
    // Not obvious what the 'one-step' simpler key is with the trie-keys now being Sss.
    pytit = g.terms_pytrees.find(get_terminals(parent_rv));

    return pytit;
  }

  // Performs some final steps for finishing a particular inside cell
  // 1. Update passive probabilities with PY-stats:
  //   a. ...include new table factors for adapted NTs
  //   b. ...include cached counts from trees already yielding the string implied by the cell
  // 2. Perform unary closure
  template<typename terminals_type>
  void inside_finalise_cell(const rvec_lt& parent_rv, StsTit& pytit) {
	S_F& parentinactives = inactives[index(parent_rv)]; //both these
	Stit_F& parentactives = actives[index(parent_rv)];  //get update below during visitation
	Stit2SLtit& parentactives_vms = active_vms[index(parent_rv)];

        if (debug >= 23000) ilout(std::cerr << "pycky::inside_finalise_cell(): ", parent_rv) 
	    << ", inactives=" << parentinactives << std::endl;
  	// PY correction (for a particular passive rvec)
	// multiply in the new-table factor for each passive item
	foreach (S_F, it, parentinactives) {
	  F pya = g.get_pya(it->first);    // PY cache statistics
	  if (pya == 1.0)
	    continue;
	  F pyb = g.get_pyb(it->first);
	  U pym = dfind(g.parent_pym, it->first);
	  U pyn = dfind(g.parent_pyn, it->first);
	  it->second *= power( (pym*pya + pyb)/(pyn + pyb), anneal);
	}

        if (debug >= 23000) ilout(std::cerr << "pycky::inside_finalise_cell(): ", parent_rv) 
	    << ", after new-table factor inactives=" << parentinactives << std::endl;
//	StsTit& pytit = get_pytit<terminals_type>(parent_rv); //TODO get rid of the method arg and do this locally?
	// multiply in the discounted counts if there are trees yielding terminals[{parent_rv}]
	// i.e. updates passive probs
	if (pytit != g.terms_pytrees.end()) {
	  add_pycache(pytit->data, parentinactives);
	  if (debug >= 23000) std::cerr << "\trelevant-cache=" << pytit->data << std::endl;
	}

        if (debug >= 23000) ilout(std::cerr << "pycky::inside_finalise_cell(): ", parent_rv) 
	    << ", after cache factor inactives=" << parentinactives << std::endl;
	
	// Some passive items obtained while visiting the mid splits above
	// could be built upon as new active items. This call creates those new actives.
	inside_unaryclose(parentinactives, parentactives, parentactives_vms, parent_rv);

	debug_print_cell_state(parent_rv);
  }

  // Innermost part of Complete action, subject to variable mapping.
  // Operates on the set of active items covering active_rvec 
  // and a single passive item that covers passive_rvec.
  // Newly obtained items covering result_rvec are put into parentinactives.
  void try_completes( const rvec_lt& active_rvec, const rvec_lt& passive_rvec,
		  const Stit parentactive,
		  S rightinactive,
		  const SLtit& mit,
		  F leftrightprob,
		  S_F& parentinactives,
		  const rvec_lt& result_rvec
		  ) {
    const SLtit parentactive_vms = mit->find1(rightinactive);
    //COMPLETE some cases, resulting in updates/creations of passives
    //visit each dotted rule (they differ only in terms of LHS and rule_weight now)
    cforeach (S_F, itparent, parentactive->data) {
            S parent = itparent->first; //the LHS of the dotted rule
            const L& m = parentactive_vms->data.at(parent);
            if (debug >= 20100) 
          	  ilout(ilout(ilout(std::cerr << "    trying to complete parent " << parent 
          		  << ". resultrange=", result_rvec)
          		  << ". activerange=", active_rvec)
			  <<". passiverange=", passive_rvec)
          		  << ". m=" << m << std::endl;
            
	    if (mcfg::compose_activepassive<U>(active_rvec, passive_rvec, m, result_rvec)) {
          	  if (debug >= 20100) std::cerr << "\t\tyes, allowed according to varmap";
          	  parentinactives[parent] += leftrightprob 
          		  * power(itparent->second/afind(g.parent_weight, parent), anneal);
		  if (debug >= 20100) std::cerr << ", inside(" << parent << ")=" << parentinactives[parent] << std::endl;
            } else 
          	  if (debug >= 20100)  std::cerr << "\t\tno according to varmap\n";
    }	
  }


  // Perform all the inside actions for binary rules 
  // for a particular factorisation of the LHS and RHS.
  // Updates the cell's inactive items. 
  // Successive calls with the same parent_rv are fine as long as the (active_rv, passive_rv) change.
  void  inside_visit_factorisation(const rvec_lt& parent_rv, const rvec_lt& active_rv, const rvec_lt& passive_rv) {
	if (debug >= 20100) {
	  ilout(std::cerr<< "+-(",parent_rv) << ")-+ : ";
       	  ilout(std::cerr, active_rv) << " / ";
	  ilout(std::cerr, passive_rv) << std::endl;
	  if (debug >= 20200) {
	    ilout(std::cerr << "parent inactives[",parent_rv) << "] = " << inactives[index(parent_rv)]  << std::endl;
	    ilout(std::cerr << "parent actives[", parent_rv) << "] = " << actives[index(parent_rv)];
	    ilout(std::cerr << ", parent active_vms[", parent_rv) << "] = " << active_vms[index(parent_rv)] << std::endl;
	  }
	}
	S_F& parentinactives = inactives[index(parent_rv)]; //gets updated with result of COMPLETEs below

	const S_F& rightinactives = inactives.at(index(passive_rv));
	if (debug >= 20100) std::cerr<< "  rightpassives = " << rightinactives << std::endl;
	if (rightinactives.empty())
	  return;
	
	Stit_F& leftactives = actives[index(active_rv)];
	Stit2SLtit& leftactives_vms = active_vms[index(active_rv)];
	if (debug >= 20100) std::cerr<< "  leftactives = " << leftactives << std::endl;
	if (debug >= 20100) std::cerr<< "  leftactives_vms = " << leftactives_vms << std::endl;
	//cross product over active items B and passive items C
	cforeach (Stit_F, itleft, leftactives) { 
	  const Stit leftactive = itleft->first;
	  const F leftprob = itleft->second;
          if (debug >= 20100) std::cerr << "   continuing active " << *leftactive << "," << leftprob  << std::endl;
	  assert(leftactives_vms.find(leftactive) != leftactives_vms.end());
	  const SLtit& mit = leftactives_vms[leftactive]; //mit is the place in the varmap trie
	  cforeach (S_F, itright, rightinactives) {
	      S rightinactive = itright->first;
	      const F rightprob = itright->second;
	      if (debug >= 20100) std::cerr << "   trying with passive " << rightinactive << "," << rightprob  << std::endl;
	      // See if there are actually completions of B with C.
	      // parentactive encodes active items whose next rule symbol is C.
	      const Stit parentactive = leftactive->find1(rightinactive);
	      if (parentactive != leftactive->end()) {//if there are such continuations...
		F leftrightprob = leftprob * rightprob;
		//complete them subject to the variable mapping
		try_completes(active_rv, passive_rv, parentactive, rightinactive, mit, leftrightprob, 
				parentinactives, parent_rv);

		assert (parentactive->key_trie.empty()); //binary rules only!
	    }
	  }
	} 
	return;
  }

  // Update the cell[parent_rv] with unary rules that concatenate child ranges into a single parent range
  void inside_visit_unary_concat(const rvec_lt& parent_rv, const rvec_lt& unarychild_rv) {
	if (debug >= 20100) {
	  ilout(std::cerr<< "+-1(",parent_rv) << ")1-+ : ";
       	  ilout(std::cerr, unarychild_rv) << std::endl;
	  if (debug >= 20200) 
	    ilout(std::cerr << "parent inactives[",parent_rv) << "] = " << inactives[index(parent_rv)]  << std::endl;
	}
	S_F& parentinactives = inactives[index(parent_rv)]; //gets updated below
	const S_F& childinactives = inactives.at(index(unarychild_rv));
	if (debug >= 20100) std::cerr<< "  childpassives = " << childinactives << std::endl;
	cforeach (S_F, itchild, childinactives) {
	  S child = itchild->first;
	  const F childprob = itchild->second;
	  const auto it = g.unarychild_parent_weight.find(child);
	  if (it != g.unarychild_parent_weight.end()) {
	    const auto& parent_weight = it->second;
	    const auto& pa_vm = afind(g.unarychild_parent_varmap, child);
	    cforeach (S_F, itparent, parent_weight) {
	      S parent = itparent->first;
	      const auto& m = afind(pa_vm, parent);
	      if (debug >= 20100) std::cerr << "\t apply " << parent << " -> " << child << " ?" << std::endl;
	      
	      if (mcfg::compose_unary<U>(unarychild_rv, m, parent_rv)) {
		  if (debug >= 20100) std::cerr << "\t\tyes, allowed according to varmap";
		  parentinactives[parent] += childprob
          		  * power(itparent->second/afind(g.parent_weight, parent), anneal);
		  if (debug >= 20100) std::cerr << ", inside(" << parent << ")=" << parentinactives[parent] << std::endl;
	      } else if (debug >= 20100)  std::cerr << "\t\tno according to varmap\n";
	    }
	  }
	}
    
  }

  void debug_print_cell_state(const rvec_lt& rvec) {
    if (debug >= 20000) {
	//these references should ideally be const; but more annoying obtaining that from the maps
	S_F& parentinactives = inactives[index(rvec)];
	Stit_F& parentactives = actives[index(rvec)];
	Stit2SLtit& parentactives_vms = active_vms[index(rvec)];

	if (debug >= 20000)
	  ilout(std::cerr << "# cky::inside() inactives[", rvec) << "] = " 
		    << parentinactives << std::endl;
	if (debug >= 20100)
	  ilout(std::cerr << "# cky::inside() actives[", rvec) << "] = " 
		    << parentactives << std::endl;
	if (debug >= 20150)
	  ilout(std::cerr << "# cky::inside() actives_vms[", rvec) << "] = " 
		    << parentactives_vms << std::endl;
	if (debug >= 20100) {
	  ilout(std::cerr << "# cky::inside() pytits[", rvec) << "] = "; 
	  if (pytits[index( rvec)] == g.terms_pytrees.end())
	    std::cerr << "()" << std::endl;
	  else
	    std::cerr << pytits[index(rvec)]->data << std::endl;
	}
    }
  }


  //! Add cache stats for the trees in tps to inside probs in inactives.
  //  Visits each tree and adds to insideprob[treeroot-category] the cache stat (discounted count) 
  // (All this applies to a particular cell)
  void add_pycache(const sT& tps, S_F& inactives) {
    cforeach (sT, it, tps) {
      symbol cat = (*it)->cat;
      F pya = g.get_pya(cat);    // PY cache statistics
      if (pya == 1.0)
	continue;
      F pyb = g.get_pyb(cat);
      U pyn = dfind(g.parent_pyn, cat); //number of times the NT cat has been expanded
      //(*it)->count is num customers at table labelled by this tree
      inactives[cat] += power( ((*it)->count - pya)/(pyn + pyb), anneal);
    }
  }  // pycky::add_cache()

  //! Constructs all items using the unary rules that apply to the cell.
  //  (inactives and actives refer to a single cell)
  void inside_unaryclose(S_F& inactives, Stit_F& actives, Stit2SLtit& active_vms, const rvec_lt& parent_rv) {
    F delta = 1;
    S_F delta_prob1 = inactives;
    S_F delta_prob0;
    if (debug >= 20200)
	std::cerr << "# inside_unaryclose; inactives = " << inactives 
		<< ", actives = " << actives << std::endl;
    // application of a unary rule may produce an item that is amenable for
    // parsing with another rule.
    // The switching between delta_prob0 and 1 effectively manages an agenda - 
    // During a particular iteration, delta_prob0 holds the 'new' items to deal with,
    // and each rule application to one of those places a new item on delta_prob1, 
    // which becomes the agenda for the next iteration.
    // The while condition seems to simultaneously guard against underflow 
    // and detect when the agenda is empty (in which case delta will equal 0 after the emptying
    // iteration.
    while (delta > unaryclosetolerance) {
      if (debug >= 20200)
	std::cerr << "Another while iteration (delta = " << delta << ", delta_prob0 = "
	       << delta_prob1 << std::endl;
      delta = 0;
      delta_prob0 = delta_prob1;

      delta_prob1.clear();
      cforeach (S_F, it0, delta_prob0) {
	S child = it0->first;
        if (debug >= 20200) std::cerr << "\tchild = " << child << std::endl;
	S_S_F::const_iterator it = g.unarychild_parent_weight.find(child);
	if (it != g.unarychild_parent_weight.end()) {
          if (debug >= 20200)
	    std::cerr << "\t\tdelta_prob1 = " << delta_prob1 
		    << ", inactives = " << inactives << std::endl;
	  const S_F& parent_weight = it->second;
	  const auto& pa_vm = afind(g.unarychild_parent_varmap, child);
	  cforeach (S_F, it1, parent_weight) {
	    S parent = it1->first;
	    const auto& m = afind(pa_vm, parent);
	    if (mcfg::compose_unary<U>(parent_rv, m, parent_rv)) {
	      F prob = it0->second;
	      F pya = g.get_pya(parent);
	      if (pya == 1)
	        prob *= power(it1->second/afind(g.parent_weight, parent), 
	          	    anneal);
	      else {
	        F pyb = g.get_pyb(parent);
	        U pym = dfind(g.parent_pym, parent);
	        U pyn = dfind(g.parent_pyn, parent);
	        prob *= power(it1->second/afind(g.parent_weight, parent)
	          	     * (pym*pya + pyb)/(pyn + pyb), 
	          	    anneal);
	      }
	      delta_prob1[parent] += prob;
	      delta = std::max(delta, prob/(inactives[parent] += prob));
	      if (debug >= 20200)
		    std::cerr << "\t\t\tparent = " << parent << ", delta_prob1 = " << delta_prob1 
			    << ", inactives = " << inactives << std::endl;
	    }
	  }
	}
      }
    }
    // make active items out of rules whose rhs starts with the symbols in the inactives
    cforeach (S_F, it0, inactives) {
      Stit it1 = g.rhs_parent_weight.find1(it0->first);
      if (it1 != g.rhs_parent_weight.end()) {
        if (debug >= 30000) 
	  std::cerr << " Updating active items for passive symbol " << it0->first 
	            << ", ruletrie is " << *it1 << ", pre_increment actives[thisone]=" << actives[it1] << std::endl;
	actives[it1] += it0->second;
        SLtit it2 = g.rhs_parent_varmap.find1(it0->first);
	if (it2 != g.rhs_parent_varmap.end()) {
	  auto it3 = active_vms.find(it1);
	  if (debug >= 30000) 
	    std::cerr << " Updating active's varmap according to grammar's varmap for this rule: " << *it2 << std::endl;
	  if (it3 != active_vms.end()) {
	    std::cerr << "BAD BAD BAD active_vms already associates the ruletriepointer with something! Chart indexing error?"
	  		<< "\tcurrently active_vms[" << *it1 << "] = " << *it3 
	  		<< " = (" << *((*it3).first) << ", " << *((*it3).second) << std::endl;
	  }
	  active_vms[it1] = it2;
        } else {
	  std::cerr << "Fail! Could not find rhs starting with " << it0->first 
		  << " in rhs_parent_varmap even though grammar has such rules" << std::endl;
	  assert(it2 != g.rhs_parent_varmap.end());
	  abort();
	}
      }
    }
    if (debug >= 30000)
	  std::cerr << "cell's actives = " << actives << std::endl <<
		  "cell's active_vms = " << active_vms << std::endl;
  } // pycky::inside_unaryclose()

 
  //! random_tree() returns a random parse tree for terminals
  //
  tree* random_tree(S s) {
    U n = terminals.size();
    return random_inactive(s, afind(inactives[index({0, n})], s), {0, n});
  }  // pycky::random_tree

  tree* random_tree() { return random_tree(g.start); }

  //! random_inactive() returns a random expansion for an inactive edge
  //
  tree* random_inactive(const S parent, F parentprob, 
			const rvec_lt& parent_rv) const { 

    //get extent of rvec
    U left = mcfg::at(parent_rv, 0), right = mcfg::at(parent_rv, parent_rv.size()-1);
    if (debug >= 23000) ilout(std::cerr << "pycky::random_inactive(): ", parent_rv) 
	    << ", inside(" << parent << ")=" << parentprob << std::endl;
    //have we reached the base of the chart? (Relies on the grammar restriction that
    //terminating rules only cover a single span or length 1 (i.e. no gaps))
    if (left+1 == right && parent == terminals[left])
      return new tree(parent);
    
    F probthreshold = parentprob * random1();  
    F probsofar = 0;
    F pya = g.get_pya(parent);
    F rulefactor = 1;

    if (debug >= 23000) ilout(std::cerr << "pycky::random_inactive(): ", parent_rv) 
	    << ", inside(" << parent << ")=" << parentprob << ", using choice threshold=" << probthreshold << std::endl;
    if (pya != 1) {
      
      // get tree from cache
      if (debug >= 23000) ilout(std::cerr << "pycky::random_inactive(): ", parent_rv) << ", considering cached trees" << std::endl;

      F pyb = g.get_pyb(parent);
      U pyn = dfind(g.parent_pyn, parent);
      const StsTit& pytit = pytits.at(index(parent_rv));
      if (pytit != g.terms_pytrees.end())
	cforeach (sT, it, pytit->data) {
	  if ((*it)->cat != parent)
	    continue;
          if (debug >= 23000) std::cerr << "\t cell-tree: " << *it << ", sofar(prev)=" << probsofar;
	  probsofar += power( ((*it)->count - pya)/(pyn + pyb), anneal);
	  if (debug >= 23000) std::cerr << ", sofar(now)=" << probsofar << std::endl;
	  if (probsofar >= probthreshold) { //random choice
	    if (debug >= 23000) std::cerr << "\t\tchosen" << std::endl;
	    return *it; //returns a terminating tree
	  }
	}
      U pym = dfind(g.parent_pym, parent);
      rulefactor = (pym*pya + pyb)/(pyn + pyb);
    }

    rulefactor /=  afind(g.parent_weight, parent);
    if (debug >= 23000) std::cerr << "\t constructing fresh tree: sofar="  << probsofar <<  std::endl;
    // tree won't come from cache, so cons up new node
    tree* tp = new tree(parent, get_terminals(parent_rv));
    TopDownInfo tdi(parent, parent_rv, rulefactor, probsofar, probthreshold, tp, this);
    if (random_unary_passive(tdi, tdi.parent_rv))
	return tp;
    if (debug >= 23000) std::cerr << "\t still constructing fresh tree: "
	    << "after random_unary_passive() sofar="  << probsofar <<  std::endl;

    if (random_binary_passive(tdi))
	return tp;

    ilout(std::cerr << "\n## Error in pycky::random_inactive(), parent = " << parent
	      << ", parent_rv = ", parent_rv)
	      << ", probsofar = " << probsofar 
	      << " still below probthreshold = " << probthreshold 
	      << std::endl;
    return tp;
  }  // pycky::random_inactive()

  //!
  // Randomly generate a tree (according to inside chart) that covers the parent_rv with category parent,
  // using a unary rule next.
  // Note that probsofar is a reference and the value gets mutated by this function.
  // Note that tp is the tree that gets extended.
  // Return value is true iff the tree is complete.
  //   For unary rules involving a single range only, call with second arg = parent_rv.
  //   Calls where child_rv are something else are made via handlecases.
  bool random_unary_passive(TopDownInfo& tdi, const rvec_lt& child_rv) const {
    const S_F& childinactives = inactives.at(index(child_rv));
    // try unary rules // A -> B
    if (debug >= 21000) {
      std::cerr << "random_unary_passive: tdi.parent=" << tdi.parent << ", parent_rv=";
      ilout(std::cerr, tdi.parent_rv);
      ilout(std::cerr << ", child_rv=", child_rv);
      std::cerr << ", childinactives=" << childinactives << std::endl;
    }
    cforeach (S_F, it0, childinactives) { //foreach covering category B
      S child = it0->first; //B
      F childprob = it0->second; //inside[..,B]
      S_S_F::const_iterator it1 = g.unarychild_parent_weight.find(child);
      //are there rules for which child is the RHS?
      if (it1 != g.unarychild_parent_weight.end()) {
	const S_F& parent1_weight = it1->second; //list of different LHSs A & their weights.
	const auto itparent = parent1_weight.find(tdi.parent);
	//does parent feature among those rules?
	if (itparent != parent1_weight.end()) {
	  const auto& pa_vm = afind(g.unarychild_parent_varmap, child);
	  if (debug >= 21000) std::cerr << "\tConsidering " << tdi.parent << " -> " << child << std::endl;
	  const auto& m = afind(pa_vm, tdi.parent);
	  //is it compatible?
	  if (mcfg::compose_unary<U>(child_rv, m, tdi.parent_rv)) {
	    //psofar += inside[..,B] * ruleweight       * rulefactor
	    //                         (prior or theta)   (newtable_factor/production normaliser)
	    tdi.probsofar += childprob 
	      * power(itparent->second*tdi.rulefactor, anneal);
	    if (tdi.probsofar >= tdi.probthreshold) { //random choice
	      if (debug >= 21000) std::cerr << "\t\tChosen." << std::endl;
	      tdi.tp->children.push_back(random_inactive(child, childprob, child_rv)); //recurse down
	      return true;
	    }
	  }
	}
      }
    }
    return false; //indicates no expansion was followed; caller should proceed with sampling

  } //pycky::random_unary_passive()

  //random_binary_passive_do()
  bool random_binary_passive_do(TopDownInfo& tdi, const rvec_lt& leftactive_rv, 
			  const rvec_lt& rightpassive_rv) const {
    const auto active_iter = actives.find(index(leftactive_rv));
    const auto passive_iter = inactives.find(index(rightpassive_rv));
    if (active_iter == actives.end() || passive_iter == inactives.end())
      return false;
    const Stit_F& leftactives = active_iter->second;
    const Stit2SLtit& leftactives_vms = active_vms.at(index(leftactive_rv));
    const S_F& rightinactives = passive_iter->second;
    cforeach (Stit_F, itleft, leftactives) {
      const Stit leftactive = itleft->first;
      const F leftprob = itleft->second;
      assert(leftactives_vms.find(leftactive) != leftactives_vms.end());
      const SLtit& mit = leftactives_vms.find(leftactive)->second; //mit is the place in the varmap trie
      cforeach (S_F, itright, rightinactives) {
        S rightinactive = itright->first;
        const F rightprob = itright->second;
        const Stit parentactive = leftactive->find1(rightinactive);
        if (parentactive != leftactive->end()) {
          S_F::const_iterator it = parentactive->data.find(tdi.parent);
          if (it != parentactive->data.end()) {
            const SLtit parentactive_vms = mit->find1(rightinactive);
            assert(parentactive_vms != mit->end());
            const L& m = parentactive_vms->data.at(tdi.parent);
            if (debug>=21000) 
                ilout(ilout(ilout(std::cerr << "\nConsider expanding (" << tdi.parent << ")", tdi.parent_rv)
      	        << " by ", leftactive_rv) << " / ", rightpassive_rv)
      	        << " using passive " << rightinactive << " with m=" << m <<  std::endl;
            if ( mcfg::compose_activepassive<U>(leftactive_rv, rightpassive_rv, m, tdi.parent_rv)) {
	      if (debug>=23000) std::cerr << "\t\tthreshold=" << tdi.probthreshold << ", sofar(prev)=" << tdi.probsofar;
              tdi.probsofar += leftprob * rightprob 
	     	* power(it->second*tdi.rulefactor, anneal); //it->second is ruleweight
	      if (debug>=23000) std::cerr << ", sofar(now)=" << tdi.probsofar << std::endl;
              if (tdi.probsofar >= tdi.probthreshold) { //random choice
		if (debug>=21000) std::cerr << "\t\t\tchosen\n";
                //expand left child (using active items)
                random_active(leftactive, leftprob, leftactive_rv, tdi.tp->children);
                //expand right child (which is passive)
                tdi.tp->children.push_back(random_inactive(rightinactive, rightprob, rightpassive_rv));
                return true;
              }
            }
          }
        }
      }
    } 
    return false; 
  }

  //random_binary_passive
  // Randomly generate a tree (according to inside chart) that covers the parent_rv with category parent,
  // using a binary rule next.
  // Note that probsofar is a reference and the value gets mutated by this function.
  // Note that tp is the tree that gets extended.
  // (The previous two lines refer to stuff stored in tdi.)
  // Return value is true iff the tree is complete.
  inline bool random_binary_passive(TopDownInfo& tdi) const {
    // try binary rules
    // visit all the relevant incremental completes that could have given parent as consequent
    // and randomly 'follow' one of them further down 

    return handle_subcases(tdi);
  }


  void random_active(const Stit parent, F parentprob, const rvec_lt& parent_rv,
		     tree::ptrs_type& siblings) const {
    F probthreshold = random1() * parentprob;
    F probsofar = 0;

    if (debug>=20200) 
	    ilout(std::cerr <<  "random_active:: over ", parent_rv) 
		<< " with parent stit=" << *parent << std::endl;
    // unary rule - in the sense that starting at the root of the grammar's trie and 
    // moving along a single symbol that is inactive in parent_rv, the parent is obtained
    const S_F& parentinactives = inactives.at(index(parent_rv));
    cforeach (S_F, it, parentinactives)
      if (g.rhs_parent_weight.find1(it->first) == parent) {
	probsofar += it->second;
	if (probsofar >= probthreshold) {
	  if (debug>=20200) 
	     std::cerr << "\tchose unary expand using child " << it->first << std::endl;
	  siblings.push_back(random_inactive(it->first, it->second, parent_rv));
	  return;
	}
	break;  // only one unary child can possibly generate this parent
      }

    ilout(std::cerr << "## Error in pycky::random_active(), parent = " << parent
	      << " parent_rv = ", parent_rv)
	      << ", probsofar = " << probsofar << ", probthreshold = " << probthreshold 
	      << std::endl;
    return;
  }  // pycky::random_active()

}; // pycky{}
//define static member
mcfg::ushint pycky::max_rangevector_dim = MAX_RANGEVECTOR_DIM_IMPL;

////// Helps print 'active' items more explicitly - JAB
struct trie_printer_helper
{ 
  std::ostream& os;
  trie_printer_helper(std::ostream& os) : os(os) {}

  template< class K, class V>
  void operator ()(const K& k, const V& v) {
    os << "" << k << ":" << v << "; ";
  }
};

std::ostream& operator<<(std::ostream& os, const pycky::Stit_F::value_type & o)
{
  os << "(" << o.first << ", " << o.second << ")";
  if (debug >= 20150) {
    os << " = {";
    o.first->for_each(trie_printer_helper(os));
    os << "}";
  }
  return os;
}
/////////////////////////// -JAB

struct resample_pycache_helper {
  typedef catcountyldtree_type tree;

  pycfg_type& g;
  pycky& p;

  resample_pycache_helper(pycfg_type& g, pycky& p) : g(g), p(p) { }

  template <typename Words, typename TreePtrs>
  void operator() (const Words& words, TreePtrs& tps) {

    foreach (typename TreePtrs, tit, tps) {
      tree* tp0 = *tit;
      const Sss& words_ = tp0->terminals();
      assert(words_.size() > 0);
      if (words_.size() > 1) {
	if (debug>=20000)
          std::cerr << "WARN: Resampling of a cached tree with discontinuous yield not supported yet. Skipping tree.\n";
        continue;
      }
      const Ss& words = words_[0]; //TODO need an inside() that works with an Sss (gapped yield).
      S start = tp0->category();
      F old_pya = g.set_pya(start, 1.0);
      F pi0 = g.decrtree(tp0); // P(tree tp0 | g without tp0)
      if (pi0 < 0)
	std::cerr << "## pi0 = " << pi0 << ", tp0 = " << tp0 << std::endl;
      assert(pi0 >= 0);
      F r0 = g.tree_prob(tp0); // q(tree tp0 | ? )
      assert(r0 >= 0);
      F tprob = p.inside(words, start);   // parse string 
      if (tprob <= 0) 
	std::cerr << "## Error in resample_pycache(): words = " << words << ", tprob = " << tprob
		  << ", tp0 = " << tp0 << std::endl
		  << "## g = " << g << std::endl;
      assert(tprob >= 0);
      tree* tp1 = p.random_tree(start); //sample proposal tree from inside chart just acquired
      F r1 = g.tree_prob(tp1);
      assert(r1 >= 0);
      
      if (tp0->general_equal(*tp1)) {  // ignore top count
	g.incrtree(tp0); //just replace old tree if it's same as proposed tree
	tp1->selective_delete();
      }
      else {  // *tp1 != *tp0, do acceptance rejection
	F pi1 = g.incrtree(tp1);
	F accept = power( (pi1 * r0) / (pi0 * r1), p.anneal);
	if (!finite(accept))  // accept if there has been an underflow
	  accept = 2.0;
	if (random1() <= accept) {
	  tp0->general_swap(*tp1);  // don't swap top counts
	  tp1->selective_delete();
	}
	else {  // don't accept
	  g.decrtree(tp1);
	  g.incrtree(tp0);
	  tp1->selective_delete();
	}
      }
      g.set_pya(tp0->category(), old_pya);
    }
  }  // resample_pycache_helper::operator()

};  // resample_pycache_helper{}

//! resample_pycache() resamples the strings associated with each cache
//
inline void resample_pycache(pycfg_type& g, pycky& p) {
  resample_pycache_helper h(g, p);
  p.g.terms_pytrees.for_each(h);
}  // resample_pycache()

namespace mcfg {
//! Create a contiguous sequence of terminals from a gapped one,
//  e.g. a_b__c_d => a_b_a_c_d, positions=0,2,3,5
//  Stores the result in spoofed and the range positions for actual within spoofed in positions.
void spoof_contiguous_yield(const ::Sss & actual, ::Ss& spoofed, std::vector<U>& positions) {
  assert(actual.size() > 0);
  positions.push_back(0);
  for (auto it=actual.begin(); it!=actual.end(); ++it) {
    const auto& yld = *it;
    spoofed.insert(spoofed.end(), yld.begin(), yld.end());
    positions.push_back(positions.back() + yld.size());
    if (it+1 != actual.end()) { //skip this in last iteration
      spoofed.push_back(yld.front()); //dummy symbol; reuse of actual terminal ensures parsability
      positions.push_back(positions.back() + 1); //account for dummy
    }
  }
}
} //namespace

#endif // PY_CKY_H
