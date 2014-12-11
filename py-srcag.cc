// vim: tabstop=8:softtabstop=8:shiftwidth=8:noexpandtab
// py-srcag.cc
//
// py-srcag runs a Pitman-Yor process for each nonterminal
// to estimate an Adaptor Grammar

const char usage[] =
"py-srcag version of 01/2013 by Jan Botha,\n"
"  based on py-cfg by Mark Johnson (18th October, 2009\n"
"\n"
"(This help output does not document the SRCG-specific extensions\n"
"and accompanying restrictions. See README.md.)\n"
"\n"
"py-srcag [-d debug]\n"
"       [-A parses-file] [-C] [-D] [-E] [-F trace-file] [-G grammar-file]\n"
"       [-H] [-I] [-R nr]\n"
"       [-r rand-init] [-n niterations]\n"
"       [-a a] [-b b] [-w weight]\n"
"       [-e pya-beta-a] [-f pya-beta-b] [-g pyb-gamma-s] [-h pyb-gamma-c]\n"
"       [-T anneal-temp-start] [-t anneal-temp-stop] [-m anneal-its]\n"
"       [-Z ztemp] [-z zits]\n"
"       [-S model-file] [-s store-every] [-L]\n"
"       [-X eval-cmd] [-x eval-every]\n"
"       grammar.lt < train.yld\n"
"\n"
" -d debug        -- debug level\n"
" -A parses-file  -- print analyses of training data at termination\n"
" -N nanal-its    -- print analyses during last nanal-its iterations\n"
" -C              -- print compact trees omitting uncached categories\n"
" -D              -- delay grammar initialization until all sentences are parsed\n"
" -E              -- estimate rule prob (theta) using Dirichlet prior\n"
" -F trace-file   -- file to write trace output to\n"
" -G grammar-file -- print grammar at termination\n"
" -H              -- skip Hastings correction of tree probabilities\n"
" -I              -- parse sentences in order (default is random order)\n"
" -R nr           -- resample PY cache strings during first nr iterations (-1 = forever)\n"
" -n niterations  -- number of iterations\n"
" -r rand-init    -- initializer for random number generator (integer)\n"
" -a a            -- default PY a parameter\n"
" -b b            -- default PY b parameter\n"
" -e pya-beta-a   -- if non-zero, parameter of Beta prior on pya\n"
" -f pya-beta-b   -- parameter of Beta prior on pya\n"
" -g pyb-gamma-s  -- if non-zero, parameter of Gamma prior on pyb\n"
" -h pyb-gamma-c  -- parameter of Gamma prior on pyb\n"
" -w weight       -- default value of theta (or Dirichlet prior) in generator\n"
" -T tstart       -- start with this annealing temperature\n"
" -t tstop        -- stop annealing at this annealing temperature\n"
" -m anneal-its   -- anneal for this many iterations\n"
" -Z ztemp        -- temperature used just before stopping\n"
" -z zits         -- perform zits iterations at temperature ztemp at end of run\n"
" -S model-file   -- save sampler state to model-file (suffixed with iteration number)\n"
" -s store-every  -- if -S given, stores model every store-every iterations (default=1)\n"
" -L              -- continue sampling from loaded state. Specify model-file instead of grammar.lt\n"
" -c		  -- in conjunction with -L, but start doing parses-file immediately\n"
" -X eval-cmd     -- pipe each run's parses into this command (empty line separates runs)\n"
" -x eval-every   -- pipe trees into the eval-cmd every eval-every iterations\n"
"\n"
"The grammar consists of a sequence of rules, one per line, in the\n"
"following format:\n"
"\n"
"   [theta [a [b]]] Parent --> Child1 Child2 ...\n"
"\n"
"where theta is the rule's probability (or, with the -E flag, the Dirichlet prior\n"
"            parameter associated with this rule) in the generator, and\n"
"      a, b (0<=a<=1, 0<b) are the parameters of the Pitman-Yor adaptor process.\n"
"\n"
"If a==1 then the Parent is not adapted.\n"
"\n"
"If a==0 then the Parent is sampled with a Chinese Restaurant process\n"
"           (rather than the more general Pitman-Yor process).\n"
"\n"
"If theta==0 then we use the default value for the rule prior (given by the -w flag).\n"
"\n"
"The start category for the grammar is the Parent category of the\n"
"first rule.\n"
"\n"
"If you specify the -C flag, these trees are printed in \"compact\" format,\n"
"i.e., only cached categories are printed.\n"
"\n"
"If you don't specify the -C flag, cached nodes are suffixed by a \n"
"'#' followed by a number, which is the number of customers at this\n"
"table.\n"
"\n"
"The -A parses-file causes it to print out analyses of the training data\n"
"for the last few iterations (the number of iterations is specified by the\n"
"-N flag).\n"
"\n"
"The -X eval-cmd causes the program to run eval-cmd as a subprocess\n"
"and pipe the current sample trees into it (this is useful for monitoring\n"
"convergence).  Note that the eval-cmd is only run _once_ per training run.\n"
"Trees belonging to different iterations are separated by blank lines.\n"
"\n"
"The program can now estimate the Pitman-Yor hyperparameters a and b for each\n"
"adapted nonterminal.  To specify a uniform Beta prior on the a parameter, set\n"
"\n"
"   -e 1 -f 1\n"
"\n"
"and to specify a vague Gamma prior on the b parameter, set\n"
"\n"
"   -g 10 -h 0.1\n"
"or\n"
"   -g 100 -h 0.01\n"
"\n"
"If you want to estimate the values for a and b hyperparameters, their\n"
"initial values must be greater than zero.  The -a flag may be useful here.\n"
"\n"
"If a nonterminal has an a value of 1, this means that the nonterminal\n"
"is not adapted.\n"
"\n"
"The sampler state saved through -S includes the grammar and most other parameters.\n"
"It excludes parameters that do not affect the model itself, including:\n"
"* output files/streams (-F, -G, -X) and -A, which implies -N too.\n"
"* debug level -d\n"
"* persistence (-S and -s can be specified afresh when doing a reload)\n"
"* random seed (neither -r nor the current seed state)\n"
"* randomisation of order in which sentences are visited (-I)\n\n"
"To reload, pass -L: ./pycfg -L [-d 100 ...] modelfile < train.yld\n"
"* train.yld should be exactly what was used for the original invocation.\n"
"* Sampler state will be loaded from 'modelfile' (i.e. don't specify grammar.lt again)\n"
"* Sampling then continues from the stored iteration index till whatever original -n was\n"
"When -c is passed, the stored values of niterations, nanal-its and internal itercounter\n"
"are ignored, and overriden by -N. (For generating posterior samples; gibbs sampling continues)\n"
"\n";

#include "custom-allocator.h"
         
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <unistd.h>
#include <vector>

int debug = 0;

#include "sampler.h"

struct S_F_incrementer {
  const F increment;
  S_F_incrementer(F increment) : increment(increment) { }

  template <typename arg_type>
  void operator() (const arg_type& arg, S_F& parent_weights) const
  {
    foreach (S_F, it, parent_weights)
      it->second += increment;
  }
};


int main(int argc, char** argv) {

  typedef std::string Str;
  typedef std::vector<Str> Strs;

  pycfg_type g;
  sampler model; 
  bool hastings_correction = true;
  bool random_order = true;
  bool delayed_initialization = false;
  U niterations = 100;
  F anneal_start = 1;
  F anneal_stop = 1;
  U anneal_its = 100;
  int resample_pycache_nits = 0;
  F z_temp = 1;
  U z_its = 0;
  unsigned long rand_init = 0;
  Str parses_filename = "", grammar_filename = "", trace_filename = "", archive_filename = "";
  Strs evalcmdstrs;
  Postreamps evalcmds;
  U eval_every = 1;
  U archive_every = 1;
  U nparses_iterations = 1;
  bool only_parse = false;
  bool load_model = false;
  bool make_output_samples_now = false;
  int chr;
  while ((chr = getopt(argc, argv, "A:CDEF:G:H:I:LN:R:S:T:X:Z:a:b:cd:e:f:g:h:m:n:r:s:t:w:x:z:q")) 
	 != -1)
    switch (chr) {
    case 'c':
      make_output_samples_now = true;
      break;
    case 'A':
      parses_filename = optarg;
      break;
    case 'C':
      catcounttree_type::compact_trees = true;
      catcountyldtree_type::compact_trees = true;
      break;
    case 'D':
      delayed_initialization = true;
      break;
    case 'E':
      g.estimate_theta_flag = true;
      break;
    case 'F':
      trace_filename = optarg;
      break;
    case 'G':
      grammar_filename = optarg;
      break;
    case 'H':
      hastings_correction = false;
      break;
    case 'I':
      random_order = false;
      break;
    case 'L':
      load_model = true;
      break;
    case 'S':
      archive_filename = optarg;
      break;
    case 'N':
      nparses_iterations = atoi(optarg);
      break;
    case 'R':
      resample_pycache_nits = atoi(optarg);
      break;
    case 'T':
      anneal_start = 1.0/atof(optarg);
      break;
    case 'X':
      evalcmdstrs.push_back(std::string(optarg));
      evalcmds.push_back(new pstream::ostream(optarg));
      break;
    case 'Z':
      z_temp = atof(optarg);
      break;
    case 'a':
      g.default_pya = atof(optarg);
      break;
    case 'b':
      g.default_pyb = atof(optarg);
      break;
    case 'd':
      debug = atoi(optarg);
      break;
    case 'e':
      g.pya_beta_a = atof(optarg);
      break;
    case 'f':
      g.pya_beta_b = atof(optarg);
      break;
    case 'g':
      g.pyb_gamma_s = atof(optarg);
      break;
    case 'h':
      g.pyb_gamma_c = atof(optarg);
      break;
    case 'm':
      anneal_its = atoi(optarg);
      break;
    case 'n':
      niterations = atoi(optarg);
      break;
    case 'r':
      rand_init = strtoul(optarg, NULL, 10);
      break;
    case 's':
      archive_every = atoi(optarg);
      break;
    case 't':
      anneal_stop = 1.0/atof(optarg);
      break;
    case 'w':
      g.default_weight = atof(optarg);
      break;
    case 'x':
      eval_every = atoi(optarg);
      break;
    case 'q':
      only_parse = true;
      break;
    case 'z':
      z_its = atoi(optarg);
      break;    default:
      std::cerr << "# Error in " << argv[0] 
		<< ": can't interpret argument -" << char(chr) << std::endl;
      std::cerr << usage << std::endl;
      exit(EXIT_FAILURE);
    }

  if (argc - optind != 1) 
    std::cerr << "# Error in " << argv[0] 
	      << ", argc = " << argc << ", optind = " << optind << '\n' 
	      << usage << std::endl;

  if (debug >= 1000) 
    std::cerr << "# eval_cmds = " << evalcmdstrs << std::endl;

#ifdef NOSERIALIZE
  if (archive_filename != "" || load_model) {
	std::cerr << "Error: Cannot specify loading/saving for binary compiled"
		<< " with -DNOSERIALIZE." << std::endl;
	return 0;
  }
#endif

  Sss trains;
  read_sentences(trains, std::cin);

  if (debug >= 1000) 
    std::cerr << "# trains.size() = " << trains.size() << std::endl;

  { 
    std::ifstream is(argv[optind]);
    if (!is)
      std::cerr << "# Error in " << argv[0] 
		<< ", can't open " << (load_model ? "model" : "grammar") << " file " << argv[optind] << std::endl;
    assert(is);

    if (load_model) {
#ifndef NOSERIALIZE
      sampler::load(g, model, is);
      if (make_output_samples_now) {
	model.nparses_iterations = nparses_iterations;
	model.niterations = model.iteration + nparses_iterations;
      }
#endif
      assert(only_parse || model.tps.size() == trains.size()); //must be same if going to sample further
    }
    else 
      is >> g;
  }

  if (rand_init == 0)
    rand_init = time(NULL);

  mt_init_genrand(rand_init);
   
  std::string* archive_name_ptr = NULL;
  if (!archive_filename.empty())
    archive_name_ptr = &archive_filename;

  std::ostream* trace_stream_ptr = NULL;
  if (!trace_filename.empty()) 
    trace_stream_ptr = new std::ofstream(trace_filename.c_str());

  if (trace_stream_ptr) {
   if (load_model)
    *trace_stream_ptr << "# Loaded model, E = " << g.estimate_theta_flag
		      << ", I = " << model.random_order 
		      << ", R = " << model.resample_pycache_nits
		      << ", n = " << model.niterations
		      << ", w = " << g.default_weight
		      << ", m = " << model.anneal_its
		      << ", Z = " << model.z_temp
		      << ", z = " << model.z_its
		      << ", T = " << 1.0/model.anneal_start
		      << ", t = " << model.anneal_stop
		      << ", a = " << g.default_pya
		      << ", b = " << g.default_pyb
		      << ", e = " << g.pya_beta_a
		      << ", f = " << g.pya_beta_b
		      << ", g = " << g.pyb_gamma_s
		      << ", h = " << g.pyb_gamma_c
		      << ", r = " << rand_init
		      << std::endl;
   else
    *trace_stream_ptr << "# D = " << delayed_initialization 
		      << ", E = " << g.estimate_theta_flag
		      << ", I = " << random_order 
		      << ", R = " << resample_pycache_nits
		      << ", n = " << niterations
		      << ", w = " << g.default_weight
		      << ", m = " << anneal_its
		      << ", Z = " << z_temp
		      << ", z = " << z_its
		      << ", T = " << 1.0/anneal_start
		      << ", t = " << anneal_stop
		      << ", a = " << g.default_pya
		      << ", b = " << g.default_pyb
		      << ", e = " << g.pya_beta_a
		      << ", f = " << g.pya_beta_b
		      << ", g = " << g.pyb_gamma_s
		      << ", h = " << g.pyb_gamma_c
		      << ", r = " << rand_init
		      << std::endl;
  }
  std::ostream* finalparses_stream_ptr = NULL;
  if (!parses_filename.empty())
    finalparses_stream_ptr = new std::ofstream(parses_filename.c_str());

  std::ostream* grammar_stream_ptr = NULL;
  if (!grammar_filename.empty())
    grammar_stream_ptr = new std::ofstream(grammar_filename.c_str());
  
  if (debug >= 1000)
    std::cerr << "# py-srcag Initial grammar = " << std::endl << g << std::endl;

  sampler* active_model;
  if (load_model)
    active_model = &model;
  else
    active_model = new sampler(niterations, anneal_start, anneal_stop, anneal_its, z_temp, z_its,
       hastings_correction, random_order, delayed_initialization,
       static_cast<U>(resample_pycache_nits), nparses_iterations);

  if (only_parse) {
    tps_type trees;
    static_parse(g, trains, trees, finalparses_stream_ptr);
    //deallocate trees
    for (auto& t : trees)
      t->selective_delete();

  } else {
    active_model->gibbs_estimate(g, trains, evalcmds, eval_every,
         finalparses_stream_ptr, grammar_stream_ptr, trace_stream_ptr, archive_name_ptr, archive_every);

    active_model->deallocate(g);
  }
  if (active_model != &model) // delete if one was allocated in the preceding block
    delete active_model;

  
  if (finalparses_stream_ptr)
    delete finalparses_stream_ptr;

  if (grammar_stream_ptr)
    delete grammar_stream_ptr;

  if (trace_stream_ptr)
    delete trace_stream_ptr;

  foreach (Postreamps, it, evalcmds)
    delete *it;

}
