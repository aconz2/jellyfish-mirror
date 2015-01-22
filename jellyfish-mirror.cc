#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>
#include <vector>

#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/large_hash_iterator.hpp>

using jellyfish::mer_dna;
typedef jellyfish::stream_manager<char**> stream_manager;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**>> sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, mer_dna> mer_iterator;
typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**>> read_parser;
 
/* ========== mer_counter ========== */
class mer_counter : public jellyfish::thread_exec {
  bool            canonical_;
  int             nb_threads_;
  mer_hash&       hash_;
  stream_manager  streams_;
  sequence_parser parser_;

public:
  mer_counter(int nb_threads, mer_hash& hash, char** file_begin, char** file_end, bool canonical) :
    canonical_(canonical),
    nb_threads_(nb_threads),
    hash_(hash),
    streams_(file_begin, file_end, 1), // 1: parse one file at a time
    parser_(mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
  { }

  virtual void start(int thid) {
    mer_dna tmp;
    for(mer_iterator mers(parser_, canonical_) ; mers; ++mers) {
      hash_.update_add(*mers, 1, tmp);
    }
    hash_.done();
  }

};

/* ========== main ========== */
namespace po = boost::program_options;

int main(int argc, char *argv[]) {

  /* provided by Guillaume */
  const bool canonical = false;

  /* default values */
  int num_threads = 1;
  int out_counter_length = 4;

  /* manadatory arguments */
  std::string in_file;
  std::string out_file;
  std::string jf_file;

  po::options_description desc("Count k-mers which appear in the jellyfish file, 'jf' This uses the same hash function (and other parameters) so that the dumped order is identical. Does not support auto merging of files right now");
  desc.add_options()
    ("help,h", "Show help message")
    ("jf,j", po::value< std::string>(&jf_file)->required(), ".jf file")
    ("input,i", po::value< std::string>(&in_file)->required(), "fast[a|q] file to count")
    ("output,o", po::value< std::string>(&out_file)->required(), "File to write to")
    ("out-counter", po::value<int>(&out_counter_length)->default_value(out_counter_length), "Number of bytes to output counts as")
    ("threads,t", po::value<int>(&num_threads)->default_value(num_threads), "Number of threads to use");

  /* ---------- parse arguments ---------------- */
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help")) {
      std::cout << desc;
      exit(0);
    } 
    po::notify(vm);
  } catch(po::error &e) {
    std::cerr << "ERROR:" << e.what() << std::endl << desc;
    exit(1);
  }

  /* ---------- Get header from jf file ----------- */
  std::ifstream jf_file_stream(jf_file);
  if(!jf_file_stream.good()) {
    std::cerr << "Error opening " << jf_file << std::endl;
    exit(1);
  }
  jellyfish::file_header header;
  header.read(jf_file_stream);

  mer_dna::k(header.key_len() / 2);
  // create the hash table with same parameters and hash function as given
  mer_hash hash(header.size(), header.key_len(), header.val_len(), num_threads, header.max_reprobe());
  hash.ary()->matrix(header.matrix());

  /* ---------- prime the hash table with the k-mers from given file ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << " === Priming hash ===" << std::endl;
    binary_reader reader(jf_file_stream, &header);
    while(reader.next()) {
      hash.set(reader.key());
    }
  }

  
  std::auto_ptr<jellyfish::dumper_t<mer_array> > dumper;
  dumper.reset(new binary_dumper(out_counter_length, header.key_len(), num_threads, out_file.c_str(), &header));
  hash.dumper(dumper.get());

  /* ---------- Count k-mers from input ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Counting k-mers ===" << std::endl;
    // jellyfish likes c style strings 'n stuff
    char **reads_c = new char*[1];
    reads_c[0] = (char *) in_file.c_str();

    mer_counter counter(num_threads, hash, reads_c, reads_c + 1, canonical);
    counter.exec_join(num_threads);

    delete[] reads_c;
  }

  /* ---------- Dump to output ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Dumping jellyfish hash to " << out_file << " ===" << std::endl;
    dumper->one_file(true);
    dumper->dump(hash.ary());

  }

  return 0;
}
