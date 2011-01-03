// References
// "A Revealing Introduction to Hidden Markov Models"
// http://www.cs.sjsu.edu/faculty/stamp/RUA/SimpleHMM.pdf

// #define _GLIBCXX_DEBUG
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include "SimpleHMM.hpp"
DEFINE_int32(N, 10, "Number of hidden state");
DEFINE_string(DATA, "data/wiki-sample.word", "Data file");

int main(int argc, char *argv[]) {
  google::InitGoogleLogging(*argv);
  google::ParseCommandLineFlags(&argc, &argv, true);


  SimpleHMM hmm;
  SimpleHMM::Observation o;
  
  hmm.read(FLAGS_DATA, o, true);
  hmm.train(FLAGS_N);

  hmm.forward_backward(o);
  hmm.viterbi(o);
  
  //  std::cout << "backward: " << hmm.backward(o) << std::endl;
  //  hmm.viterbi();
	
  return 0;
};
