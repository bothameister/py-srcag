#ifndef YIELDEXTRACTOR_7KL8TF2J
#define YIELDEXTRACTOR_7KL8TF2J

#include "custom-allocator.h"

#include <map>
#include <set>
#include <vector>

#include "sym.h"
#include "utility.h"
#include "xtree.h"
#include "nested_list.h"
#include "py-cky.h"
#include "consts.h"

namespace analysis {
  typedef ::Sss Sss;
  typedef std::set<symbol> symset_t;
  typedef std::map<Sss, symset_t> yield_symset_t;

  typedef std::vector<std::pair<Sss, double> > yield_prob_t;
  typedef std::map<symbol, yield_prob_t> sym_yield_prob_t;
  typedef std::set<symbol> cats_t;

  typedef pycky::tree tree;
  typedef std::vector<tree*> tps_type;

  void node_inserter(const tree* t, const cats_t& targets, yield_symset_t& found_yields) {
    if (targets.find(t->cat) != targets.end())
      found_yields[t->yld].insert(t->cat);
  }
  //! Extract yields for the categories of interest (targets).
  // Outputs to found_yields
  void get_unique_yields_by_cat(const tps_type& tps, const cats_t& targets,
      yield_symset_t& found_yields) {
    for (tree* t : tps) {
      if (t->children.size() > 0) {
        node_inserter(t, targets, found_yields); //this node
        get_unique_yields_by_cat(t->children, targets, found_yields);//recurse for children
      }
    }
  }

  //! Compute inside probabilities for the yields and their given categories.
  //  Outputs to results.
  void score_cat_yields(const yield_symset_t& found_yields, const pycfg_type& g,
    sym_yield_prob_t& results) {

    pycky p(g);
    std::size_t total = found_yields.size();
    std::size_t every = 10, i=0;
    for (auto& kv : found_yields) {
      const auto& yield = kv.first;
      assert(yield.size() <= g.max_rangevector_dim);
      assert(yield.size() > 0);
      
      const Ss* toparse_p = 0;
      mcfg::R idx = 0; //chart cell index
      if (yield.size() == 1) {
        toparse_p = &(yield[0]);
        idx = pycky::index({0,(mcfg::U)toparse_p->size()});
      } else if (yield.size() > 1) {
        toparse_p = new Ss();
        std::vector<mcfg::U> positions;
        mcfg::spoof_contiguous_yield(yield, const_cast<Ss&>(*toparse_p), positions);
        assert(positions.size()>0);
        idx = pycky::index(&(positions[0]), positions.size());
        if (debug >= 5000)
          std::cerr << "Spoofing discontiguous '" << yield << "' as '" << *toparse_p << "' "
                    << "rvec_lt=" << positions << std::endl;
      }
      assert(toparse_p);
      const auto& toparse = *toparse_p;
      if (debug >= 2000)
        std::cerr << "Parsing " << toparse << std::endl;
      p.inside(toparse);
      for (auto& cat : kv.second) {
        double prob = dfind(p.inactives[idx], cat);
        if (debug >= 5000)
          std::cerr << "\tinside(" << cat << ", " << yield << ") = " << prob << std::endl;
        results[cat].push_back(std::make_pair(yield, prob));
      }
      if (yield.size() > 1)
        delete toparse_p;
      if (++i % every == 0 && debug>=1000)
        std::cerr << i << "/" << total << std::endl;
    }
  }

  void output_triplets(const pycfg_type& g, const tps_type& trees, const symset_t& targets, std::ostream* out = &std::cerr) {
    yield_symset_t found;
    get_unique_yields_by_cat(trees, targets, found);
    
    if (debug >= 1000) { 
      std::cerr << "Yields found:" << std::endl;
      for (auto& f : found)
        std::cerr << f << std::endl;
    }
    sym_yield_prob_t results;
    score_cat_yields(found, g, results);
    
    // write: category yield prob
    for (auto& kv : results) {
      for (auto& pair : kv.second) 
        (*out) << kv.first << " " << pair.first << " " << pair.second << " " << std::endl;
    }
  }

}//namespace analysis

#endif /* end of include guard: YIELDEXTRACTOR_7KL8TF2J */
