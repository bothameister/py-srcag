const char usage[] =
"yieldextractor -- does various things with a loaded model\n"
"\n"
"yieldextractor [-d debug]\n"
"       [-A parses-file] [-G grammar-file]\n"
"       [-i cats]\n"
"       modelfile\n"
"\n"
" -d debug        -- debug level\n"
" -A parsed-file  -- print the parse trees used\n"
" -i cats         -- grammar categories of interest\n"
" -G grammar-file -- print grammar at termination\n"
"\n"
"Loads 'modelfile'. Extracts the set of unique yields dominated by each cat in cats\n"
"across the analyses contained in the model file.\n"
"Then outputs triples <cat> <yield> <prob> for each one, where <prob>\n"
"is the inside probability of that yield under that category\n"
"according to a static snapshot of the PY-AG.\n"
"\n"
"Alternative invocation:\n"
"  yieldextractor [options] -p modelfile < data_to_parse\n"
"In this mode, it ignores the analyses in the modelfile and instead parses the data\n"
"passed via stdin.\n\n"
"When -i is not given, it simply outputs via -A and -G if they're specified\n"
"\n";

#include <iostream>
#include <unistd.h>
#include <time.h>

#include "yieldextractor.h"
#include "sampler.h"

int debug = 0;

int main(int argc, char** argv)
{
  using namespace analysis;

  pycfg_type g;
  sampler model; 
  std::string parses_filename = "", grammar_filename = "", targets_str = "";
  bool parse_input = false;
  unsigned long rand_init = 0;

  int chr;
  while ((chr = getopt(argc, argv, "A:pr:i:d:G:")) 
	 != -1)
    switch (chr) {
    case 'G':
      grammar_filename = optarg;
      break;
    case 'd':
      debug = atoi(optarg);
      break;
    case 'A':
      parses_filename = optarg;
      break;
    case 'p':
      parse_input = true;
      break;
    case 'i':
      targets_str = optarg;
      break;
    case 'r':
      rand_init = strtoul(optarg, NULL, 10);
      break;
    default:
      std::cerr << "# Error in " << argv[0] 
		<< ": can't interpret argument -" << char(chr) << std::endl;
      std::cerr << usage << std::endl;
      exit(EXIT_FAILURE);
    }

  if (argc - optind != 1) 
    std::cerr << "# Error in " << argv[0] 
	      << ", argc = " << argc << ", optind = " << optind << '\n' 
	      << usage << std::endl;

  std::ostream* finalparses_stream_ptr = NULL;
  if (!parses_filename.empty())
    finalparses_stream_ptr = new std::ofstream(parses_filename.c_str());

  std::ostream* grammar_stream_ptr = NULL;
  if (!grammar_filename.empty())
    grammar_stream_ptr = new std::ofstream(grammar_filename.c_str());

  if (rand_init == 0)
    rand_init = time(NULL);
  mt_init_genrand(rand_init);

  //load model
  std::ifstream is(argv[optind]);
  if (!is)
    std::cerr << "# Error in " << argv[0] 
      << ", can't open model file " << argv[optind] << std::endl;
  assert(is);
  sampler::load(g, model, is);

  //find trees
  tps_type* tps = 0;
  if (parse_input) {
    Sss trains;
    read_sentences(trains, std::cin);
    tps = new tps_type();
    static_parse(g, trains, *tps);
  } else {
    tps = &model.tps;
  }
  assert(tps);
  if (finalparses_stream_ptr)
    write_tps_to_stream(*finalparses_stream_ptr, *tps);

  if (grammar_stream_ptr) 
    (*grammar_stream_ptr) << g;

  // get categories of interest
  std::istringstream istream(targets_str);
  Ss targets_v;
  readline_symbols(istream, targets_v);
  symset_t targets;
  targets.insert(targets_v.begin(), targets_v.end());
  if (targets.size() > 0) {
    output_triplets(g, *tps, targets, &std::cout);
  }

  //teardown
  if (tps != &model.tps) {
    for (auto& t : (*tps))
      t->selective_delete();
    delete tps;
  }

  model.deallocate(g);

  if (finalparses_stream_ptr)
    delete finalparses_stream_ptr;
  if (grammar_stream_ptr)
    delete grammar_stream_ptr;

  return 0;
}

