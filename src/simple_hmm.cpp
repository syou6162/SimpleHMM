// References
// "A Revealing Introduction to Hidden Markov Models"
// http://www.cs.sjsu.edu/faculty/stamp/RUA/SimpleHMM.pdf

// #define _GLIBCXX_DEBUG
#include <set>
#include <algorithm>
#include <boost/foreach.hpp>
#include "SimpleHMM.hpp"

int main(int argc, char *argv[]) {
  google::InitGoogleLogging(*argv);
  google::ParseCommandLineFlags(&argc, &argv, true);


  SimpleHMM hmm;
  SimpleHMM::Observation o;
  
  hmm.read("nonparametric-exercise-ja-1/wiki-sample.word", o, true);

  hmm.train(10);

  hmm.forward_backward(o);
  hmm.viterbi(o);
  
  //  std::cout << "backward: " << hmm.backward(o) << std::endl;
  //  hmm.viterbi();
	
  return 0;
};
