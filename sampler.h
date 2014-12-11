#ifndef SAMPLER_PYA4554
#define SAMPLER_PYA4554

#include "custom-allocator.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "pstream.h"
#include "py-cky.h"
#include "sym.h"
#include "xtree.h"
#include "consts.h"

#ifndef NOSERIALIZE
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif

#define HERE   __FILE__ << ":" << __LINE__ << " in " << __func__ 

//typedef unsigned int U;
typedef mcfg::U U;
typedef std::vector<Ss> Sss;

typedef pstream::ostream* Postreamp;
typedef std::vector<Postreamp> Postreamps;

typedef pycky::tree tree;
typedef std::vector<tree*> tps_type;

// Sampler doesn't own the grammar g or its data
struct sampler {

  sampler(
		 U niterations = 100, 
		 F anneal_start = 1, F anneal_stop = 1, U anneal_its = 0,
		 F z_temp = 1.0, U z_its = 0,
		 bool hastings_correction = true, bool random_order = true,
		 bool delayed_initialization = false,
		 U resample_pycache_nits = 0,
		 U nparses_iterations = 0)
     : 
      niterations(niterations), anneal_start(anneal_start), anneal_stop(anneal_stop), anneal_its(anneal_its),
      z_temp(z_temp), z_its(z_its), hastings_correction(hastings_correction), random_order(random_order), 
      resample_pycache_nits(resample_pycache_nits),  
      nparses_iterations(nparses_iterations),
      iteration(0), nwords(0), delayed_initialization(delayed_initialization), pp(0)
  {}

  ~sampler() {
    if (pp) {
      delete pp;
      pp = 0;
    }
  }

  //--- Persistent state variables ---
  //* Constructor arguments
  U niterations;
  F anneal_start;
  F anneal_stop;
  U anneal_its;
  F z_temp;
  U z_its;
  bool hastings_correction;
  bool random_order;
  U resample_pycache_nits;
  U nparses_iterations;

  //* Not constructor arguments
  U iteration;
  U nwords;

  typedef pycky::tree tree;
  typedef std::vector<tree*> tps_type;
  tps_type tps;

  //--- non-persisted fields ----
  bool delayed_initialization;
  pycky* pp;

  void deallocate(pycfg_type& g) {
    bool estimate_theta_flag = g.estimate_theta_flag;
    g.estimate_theta_flag = false;
    for (U i = 0; i < tps.size(); ++i) {
      g.decrtree(tps[i], 1);
      tps[i]->selective_delete();
    }
    g.estimate_theta_flag = estimate_theta_flag;
  }

  F gibbs_estimate(pycfg_type& g, const Sss& trains,
       Postreamps& evalcmds, U eval_every,
		   std::ostream* finalparses_stream_ptr = NULL,
		   std::ostream* grammar_stream_ptr = NULL,
		   std::ostream* trace_stream_ptr = NULL,
       const std::string* archive_name = NULL, U archive_every = 1) {

    pp = new pycky(g, anneal_start);
    pycky& p = *pp;
    if (iteration == 0)
      first_pass(g, trains, trace_stream_ptr);

    U n = trains.size();
    if (trace_stream_ptr)
      *trace_stream_ptr << "# " << nwords << " tokens in " << n << " sentences" << std::endl
        << "#\n"
        << "# It\tTemp\t-logP\ttables\tsame\tchanged\treject\t(parent pym pyn pya pyb)*" << std::endl;

    F sum_log2prob = 0;

    typedef std::vector<U> Us;
    Us index(n);
    U unchanged = 0, rejected = 0;

    for (unsigned i = 0; i < n; ++i) 
      index[i] = i;

    for (; iteration < niterations; ++iteration) {

      if (random_order)
        std::random_shuffle(index.begin(), index.end());

      if (iteration + z_its > niterations) 
        p.anneal = 1.0/z_temp;
      else if (iteration == 0 && anneal_its > 0)
        p.anneal = anneal_start;
      else if (iteration < anneal_its) 
        p.anneal = anneal_start*power(anneal_stop/anneal_start,F(iteration)/F(anneal_its-1));
      else 
        p.anneal = anneal_stop;

      assert(finite(p.anneal));

      if (debug >= 100) {
        std::cerr << "# Iteration " << iteration << ", " 
          << g.sum_pym() << " tables, "
          << -g.logPcorpus()/(log(2)*nwords) << " bits/token, " 
          << unchanged << '/' << n << " unchanged";
        if (hastings_correction)
          std::cerr << ", " << rejected << '/' << n-unchanged 
            << " rejected";
        if (p.anneal != 1)
          std::cerr << ", temp = " << 1.0/p.anneal;
        std::cerr << '.' << std::endl;
      }

      if (trace_stream_ptr && iteration % eval_every == 0) {
        *trace_stream_ptr << iteration << '\t'            // iteration
          << 1.0/p.anneal << '\t'         // temperature
          << -(g.logPcorpus()+g.logPrior()) << '\t' // - log P(corpus)
          << g.sum_pym() << '\t'          // # of tables
          << unchanged << '\t'            // # unchanged parses
          << n-unchanged << '\t'          // # changed parses
          << rejected << std::flush;      // # parses rejected
        if (g.pyb_gamma_s > 0 && g.pyb_gamma_c > 0 && debug >= 10)
          g.write_adaptor_parameters(*trace_stream_ptr);
        *trace_stream_ptr << std::endl;
      }

      if (iteration % eval_every == 0)
        foreach (Postreamps, ecit, evalcmds) {
          pstream::ostream& ec = **ecit;
          for (U i = 0; i < n; ++i) 
            ec << tps[i] << std::endl;
          ec << std::endl;
        }

      if (debug >= 500)
        assert(g.sum_pym() == g.terms_pytrees_size());

      if (debug >= 10000)
        std::cerr << g;

#ifndef NOSERIALIZE
      if (archive_name)
        do_archive(g, archive_name, archive_every);
#endif

      sum_log2prob = 0;
      unchanged = 0;
      rejected = 0;

      for (U i0 = 0; i0 < n; ++i0) {
        U i = index[i0];
        if (debug >= 1000)
          std::cerr << "# trains[" << i << "] = " << trains[i] << std::endl;

        tree* tp0 = tps[i];                // get the old parse for sentence to resample

        F pi0 = g.decrtree(tp0);           // remove the old parse's fragments from the CRPs
        if (pi0 <= 0) 
          std::cerr << "## " << HERE 
            << " Underflow in gibbs_estimate() while computing pi0 = decrtree(tp0):"
            << " pi0 = " << pi0
            << ", iteration = " << iteration 
            << ", trains[" << i << "] = " << trains[i] 
            //      << std::endl << "## tp0 = " << tp0 
            << std::endl;

        F r0 = g.tree_prob(tp0);            // compute old tree's prob under proposal grammar
        if (r0 <= 0) 
          std::cerr << "## " << HERE 
            << " Underflow in gibbs_estimate() while computing r0 = tree_prob(tp0):"
            << " r0 = " << r0
            << ", iteration = " << iteration 
            << ", trains[" << i << "] = " << trains[i] 
            //      << std::endl << "## tp0 = " << tp0 
            << std::endl;

        F tprob = p.inside(trains[i]);       // compute inside CKY table for proposal grammar
        if (tprob <= 0) 
          std::cerr << "## " << HERE
            << " Underflow in gibbs_estimate() while computing tprob = inside(trains[i]):"
            << " tprob = " << tprob
            << ", iteration = " << iteration 
            << ", trains[" << i << "] = " << trains[i] 
            //      << std::endl << "## g = " << g 
            << std::endl;
        assert(tprob > 0);

        if (debug >= 1000)
          std::cerr << ", tprob = " << tprob;
        sum_log2prob += log2(tprob);

        tree* tp1 = p.random_tree();         // sample proposal parse from proposal grammar CKY table
        F r1 = g.tree_prob(tp1);
        // assert(r1 > 0);

        if (*tp0 == *tp1) {                  // don't do anything if proposal parse is same as old parse
          if (debug >= 1000)
            std::cerr << ", tp0 == tp1" << std::flush;
          ++unchanged;
          g.incrtree(tp1, 1);
          tps[i] = tp1;
          tp0->selective_delete();
        }
        else {
          F pi1 = g.incrtree(tp1, 1);        // insert proposal parse into CRPs, compute proposal's true probability
          // assert(pi1 > 0);

          if (debug >= 1000)
            std::cerr << ", r0 = " << r0 << ", pi0 = " << pi0
              << ", r1 = " << r1 << ", pi1 = " << pi1 << std::flush;

          if (hastings_correction) {         // perform accept-reject step
            F accept = (pi1 * r0) / (pi0 * r1); // acceptance probability
            if (p.anneal != 1)
              accept = power(accept, p.anneal);
            if (!finite(accept))  // accept if there has been an underflow
              accept = 2.0;
            if (debug >= 1000)
              std::cerr << ", accept = " << accept << std::flush;
            if (random1() <= accept) {      // do we accept the proposal parse?
              if (debug >= 1000)            //  yes
                std::cerr << ", accepted" << std::flush;
              tps[i] = tp1;                 //  insert proposal parse into set of parses
              tp0->selective_delete();      //  release storage associated with old parse
            }
            else {                          // reject proposal parse
              if (debug >= 1000)
                std::cerr << ", rejected" << std::flush;
              g.decrtree(tp1, 1);           // remove proposal parse from CRPs
              g.incrtree(tp0, 1);           // reinsert old parse into CRPs
              tp1->selective_delete();      // release storage associated with proposal parse
              ++rejected;
            }
          }
          else {                            // no hastings correction
            tps[i] = tp1;                   // save proposal parse
            tp0->selective_delete();        // delete old parse
          }
        }

        if (debug >= 1000)
          std::cerr << ", tps[" << i << "] = " << tps[i] << std::endl;
      }

      if (iteration < resample_pycache_nits)
        resample_pycache(g, p);

      if (iteration > 1 && g.pyb_gamma_s > 0 && g.pyb_gamma_c > 0) {
        if (g.pya_beta_a > 0 && g.pya_beta_b > 0)
          g.resample_pyab();
        else
          g.resample_pyb();
      }

      if (finalparses_stream_ptr && iteration + nparses_iterations >= niterations) {
        for (U i = 0; i < n; ++i) 
          (*finalparses_stream_ptr) << tps[i] << std::endl;
        (*finalparses_stream_ptr) << std::endl;
      }

    } //end for-iteration

    if (trace_stream_ptr) {
      *trace_stream_ptr << niterations << '\t'          // iteration
        << 1.0/p.anneal << '\t'         // temperature
        << -g.logPcorpus() << '\t'      // - log P(corpus)
        << g.sum_pym() << '\t'          // # of tables
        << unchanged << '\t'            // # unchanged parses
        << n-unchanged << '\t'          // # changed parses
        << rejected << std::flush;      // # parses rejected
      if (g.pyb_gamma_s > 0 && g.pyb_gamma_c > 0 && debug >= 10)
        g.write_adaptor_parameters(*trace_stream_ptr);
      *trace_stream_ptr << std::endl;
    }

    foreach (Postreamps, ecit, evalcmds) {
      pstream::ostream& ec = **ecit;
      for (U i = 0; i < n; ++i) 
        ec << tps[i] << std::endl;
      ec << std::endl;
    }

    F logPcorpus = g.logPcorpus();

    if (debug >= 10) {
      std::cerr << "# " << niterations << " iterations, " 
        << g.sum_pym() << " tables, "
        << " log P(trees) = " << logPcorpus << ", "
        << -logPcorpus/(log(2)*nwords) << " bits/token, " 
        << unchanged << '/' << n << " unchanged";
      if (hastings_correction)
        std::cerr << ", " << rejected << '/' << n-unchanged 
          << " rejected";
      std::cerr << '.' << std::endl;
    }

    if (debug >= 10000) 
      std::cerr << "# g.terms_pytrees = " << g.terms_pytrees << std::endl;

    if (grammar_stream_ptr) 
      (*grammar_stream_ptr) << g;

#ifndef NOSERIALIZE
      if (archive_name)
        do_archive(g, archive_name, 1);
#endif

    delete pp;
    pp = 0;

    return logPcorpus;
  }


#ifndef NOSERIALIZE
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & niterations;
    ar & anneal_start;
    ar & anneal_stop;
    ar & anneal_its;
    ar & z_temp;
    ar & z_its;
    ar & hastings_correction;
    ar & random_order;
    ar & resample_pycache_nits;
    ar & nparses_iterations;
    ar & iteration;
    ar & nwords;
    ar & tps;

//    std::cerr << "---sampler::serialize()---\n";
//    std::cerr << "niterations" << niterations << "\n\n";
//    std::cerr << "anneal_start" << anneal_start << "\n\n";
//    std::cerr << "anneal_stop" << anneal_stop << "\n\n";
//    std::cerr << "anneal_its" << anneal_its << "\n\n";
//    std::cerr << "z_temp" << z_temp << "\n\n";
//    std::cerr << "z_its" << z_its << "\n\n";
//    std::cerr << "hastings_correction" << hastings_correction << "\n\n";
//    std::cerr << "random_order" << random_order << "\n\n";
//    std::cerr << "resample_pycache_nits" << resample_pycache_nits << "\n\n";
//    std::cerr << "nparses_iterations" << nparses_iterations << "\n\n";
//    std::cerr << "iteration" << iteration << "\n\n";
//    std::cerr << "nwords" << nwords << "\n\n";
//    std::cerr << "tps" << tps << "\n\n";
  }

  static void save(pycfg_type& g, sampler& model, std::ostream& os)
  {
    boost::archive::text_oarchive oa(os);
    oa << g;
    oa << model;
  }

  static void load(pycfg_type& g, sampler& model, std::istream& is)
  {
    boost::archive::text_iarchive ia(is);
    ia >> g;
    ia >> model;
    if (debug >=1000) {
      std::cerr << "Loaded." << std::endl 
        << "g.logPcorpus = " << g.logPcorpus() << std::endl
        << "g.terms_pytrees.size = " << g.terms_pytrees_size() << std::endl;
      float f=0;
      for (std::size_t i = 0; i < model.tps.size(); ++i) {
        f += log(g.tree_prob(model.tps[i])); 
      }
      std::cerr << "sum_log_treeprobs = " << f << std::endl;
      std::cerr << "sum_pym=" << g.sum_pym() << std::endl;
    }
  }
private:
  void do_archive(pycfg_type& g, const std::string* archive_name, U archive_every = 1) {
    if (archive_every>0 && iteration % archive_every == 0) {
      std::stringstream filename;
      filename << *archive_name << iteration;
      if (debug>=1000)
        std::cerr << "# Archiving to '" << filename.str() << std::endl;
      std::ofstream ofs(filename.str().c_str());
      save(g, *this, ofs);
    }
  }
#endif
private:
  void first_pass(pycfg_type& g, const Sss& trains, std::ostream* trace_stream_ptr = NULL ) {
    assert(iteration == 0); //should only be called the first time
    if (tps.size() > 0) //don't allow this; complicates thinking re delayed_init etc.
    {
      std::cerr << "Reloading the state of iteration 0 is not allowed. Exiting." << std::endl;
      abort();
    }

    assert(pp); //pp must be allocated
    pycky& p = *pp;

    for (std::size_t i=0; i<trains.size(); ++i)
      tps.push_back(tps_type::value_type());

    U n = trains.size();
    F sum_log2prob = 0;

    // initialize tps with trees, but don't learn during first pass

    for (unsigned i = 0; i < n; ++i) {
      if (debug >= 1000)
        std::cerr << "# trains[" << i << "] = " << trains[i] << std::endl;
      nwords += trains[i].size();

      F tprob = p.inside(trains[i]);

      if (debug >= 1000)
        std::cerr << ", tprob = " << tprob;

      if (tprob <= 0) 
        std::cerr << "\n## " << HERE << " Error in py-cfg::gibbs_estimate(), tprob = " << tprob
          << ", trains[" << i << "] = " << trains[i] << " failed to parse." << std::endl;
      assert(tprob > 0);

      sum_log2prob += log2(tprob);
      tps[i] = p.random_tree();

      if (debug >= 1000)
        std::cerr << ", tps[" << i << "] = " << tps[i] << std::endl;

      if (!delayed_initialization)
        g.incrtree(tps[i]);
    }

    if (delayed_initialization)    // collect statistics from the random trees
      for (unsigned i = 0; i < n; ++i) 
        g.incrtree(tps[i]);

  }

};

void read_sentences(Sss& sents, std::istream& in)
{
  Ss terminals;
  while (readline_symbols(in, terminals)) 
    if (terminals.empty())
      std::cerr << "## Error in input data: sentence " 
        << sents.size()+1 << " is empty"
        << std::endl;
    else
      sents.push_back(terminals);
}

void write_tps_to_stream(std::ostream& finalparses_stream, const tps_type& tps) {
  for (U i = 0; i < tps.size(); ++i) 
    finalparses_stream << tps[i] << std::endl;
  finalparses_stream << std::endl;
}

//! Only parse the data with the basic grammar, no learning
double static_parse(pycfg_type& g, const Sss& trains, 
     tps_type& returned_trees,
		 std::ostream* finalparses_stream_ptr = NULL) {

  U n = trains.size();
  U nwords = 0;
  tps_type& tps = returned_trees;
  pycky p(g, 1);

  for (unsigned i = 0; i < n; ++i) {
    if (debug >= 1000)
      std::cerr << "# trains[" << i << "] = " << trains[i] << std::endl;
    nwords += trains[i].size();
    
    F tprob = p.inside(trains[i]);

    if (debug >= 1000)
      std::cerr << ", tprob = " << tprob;

    if (tprob <= 0) 
      std::cerr << "\n## " << HERE << " Error in py-cfg::static_parse(), tprob = " << tprob
		<< ", trains[" << i << "] = " << trains[i] << " failed to parse." << std::endl;
    assert(tprob > 0);

    tps.push_back(p.random_tree());
    
    if (debug >= 1000)
      std::cerr << ", tps[" << i << "] = " << tps.back() << std::endl;

  }
  if (finalparses_stream_ptr) {
    write_tps_to_stream(*finalparses_stream_ptr, tps);
  }
	
  double f=0;
  for (std::size_t i = 0; i < tps.size(); ++i) {
    f += log(g.tree_prob(tps[i])); 
  }
  return f;
}  //static_parse


#endif /* end of include guard: SAMPLER_PYA4554 */
