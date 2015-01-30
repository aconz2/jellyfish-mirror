#pragma GCC diagnostic ignored "-Wunused-local-typedefs"

#include <string>
#include <fstream>
#include <vector>
#include <cmath>

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
  unsigned int    limit_;
  int             nb_threads_;
  mer_hash&       hash_;
  stream_manager  streams_;
  sequence_parser parser_;

public:
  mer_counter(int nb_threads, mer_hash& hash, char** file_begin, char** file_end, bool canonical, int limit) :
    canonical_(canonical),
    limit_(limit),
    nb_threads_(nb_threads),
    hash_(hash),
    streams_(file_begin, file_end, 1), // 1: parse one file at a time
    parser_(mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
  { }

  virtual void start(int thid) {
    if(limit_ > 0) {
      mer_array* ary = hash_.ary();
      uint64_t val;
      for(mer_iterator mers(parser_, canonical_) ; mers; ++mers) {
        if(ary->get_val_for_key(*mers, &val)) {
           if(val < limit_) { 
             hash_.add(*mers, 1);
           } // else the val is already at limit_
        } // else the key doesn't exist, no op 
      }
      
    } else {
      mer_dna tmp;
      for(mer_iterator mers(parser_, canonical_) ; mers; ++mers) {
        // this will only add 1 if the mer is already in the hash
        // tmp is some optimization thing
        hash_.update_add(*mers, 1, tmp);
      }
    }
    hash_.done();
  }

};

/* ========== main ========== */
namespace po = boost::program_options;

int main(int argc, char *argv[]) {


  /* default values */
  bool canonical = false;
  int num_threads = 1;
  int out_counter_length = 4;
  int count_limit = 0;

  /* manadatory arguments */
  std::vector<std::string> in_files;
  std::string out_file;
  std::string jf_file;

  po::options_description desc("Count k-mers which appear in the jellyfish file, 'jf' This uses the same hash function (and other parameters) so that the dumped order is identical. Does not support auto merging of files right now");
  desc.add_options()
    ("help,h", "Show help message")
    ("jf,j", po::value< std::string>(&jf_file)->required(), ".jf file")
    ("input,i", po::value<std::vector<std::string> >(&in_files)->multitoken()->required(), "fast[a|q] file(s) to count")
    ("output,o", po::value< std::string>(&out_file)->required(), "File to write to")
    ("limit,l", po::value<int>(&count_limit)->default_value(count_limit), "Limit counts to no greater than limit (0 means no limit)")
    ("canonical,c", po::bool_switch(&canonical)->default_value(canonical), "Count k-mers canonically")
    //("out-counter", po::value<int>(&out_counter_length)->default_value(out_counter_length), "Number of bytes to output counts as")
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

  //std::cerr << "The input header file says size: " << header.size() << std::endl;

  mer_dna::k(header.key_len() / 2);
  int val_len = count_limit == 0 ? header.val_len() : (((int) std::log2(count_limit)) + 1);
  mer_hash hash(header.size(), header.key_len(), val_len, num_threads, header.max_reprobe());
  hash.ary()->matrix(header.matrix());

  /* ---------- prime the hash table with the k-mers from given file ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Priming hash ===" << std::endl;
    binary_reader reader(jf_file_stream, &header);
    while(reader.next()) {
      hash.set(reader.key());
    }
  }

  /* ---------- create the dumper ----------------- */
  std::auto_ptr<jellyfish::dumper_t<mer_array> > dumper;
  out_counter_length = count_limit ? (int) std::ceil(std::log2(count_limit) / 8) : out_counter_length;
  dumper.reset(new binary_dumper(out_counter_length, header.key_len(), num_threads, out_file.c_str(), &header));
  hash.dumper(dumper.get());

  /* ---------- Count k-mers from input ---------- */
  {
    boost::timer::auto_cpu_timer t(std::cerr, 2);
    std::cerr << "=== Counting k-mers ===" << std::endl;
    // jellyfish likes c style strings 'n stuff
    char **in_files_c = new char*[in_files.size()];
    for(size_t i = 0; i < in_files.size(); i++){
      in_files_c[i] = new char[in_files[i].size() + 1];
      strcpy(in_files_c[i], in_files[i].c_str());
    } 

    mer_counter counter(num_threads, hash, in_files_c, &in_files_c[in_files.size() - 1], canonical, count_limit);
    counter.exec_join(num_threads);

    delete[] in_files_c;
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
