// References
// "A Revealing Introduction to Hidden Markov Models"
// http://www.cs.sjsu.edu/faculty/stamp/RUA/SimpleHMM.pdf

#ifndef SIMPLE_HMM_HPP
#define SIMPLE_HMM_HPP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <numeric>
#include <limits>
#include <fstream>
#include <sstream>
#include <tr1/unordered_map>
#include <boost/random.hpp>
#include <gflags/gflags.h>
#include <glog/logging.h>
// #include "util.hpp"

class SimpleHMM {
public:
  typedef std::vector<std::vector<double> > TransitionMatrix;
  typedef std::vector<std::vector<double> > EmissionMatrix;
  typedef std::vector<int> Observation;
  typedef std::vector<double> InitialProbVector;

  SimpleHMM();
  ~SimpleHMM();

  void read(const std::string& filename, Observation& o, const bool train);
  void init_alpha_and_c(const Observation& o); 
  void init_beta(const Observation& o);
  void set_a(TransitionMatrix& a);
  void set_b(EmissionMatrix& b);
  void set_pi(InitialProbVector& pi);
  void train(const int N);
	
  void forward(const Observation& o, bool init = true);
  void backward(const Observation& o, bool init = true);
  void forward_backward(const Observation& o);
  void viterbi(const Observation& o);
  int getID(const std::string& str, const bool train);
protected:
  std::map<std::string, int> _word2id;
  std::vector<std::string>   _id2word;
  ////////////////////////////////////////
  // parameters to be estimated
  ////////////////////////////////////////
  TransitionMatrix _a;
  EmissionMatrix _b;
  InitialProbVector _pi;

  ////////////////////////////////////////
  // tables for DP
  ////////////////////////////////////////
  std::vector<std::vector<double> > _alpha;
  std::vector<std::vector<double> > _beta;
	
  std::vector<double> _c; // scaling coefficient

  std::vector<std::vector<std::tr1::unordered_map<int, double> > > _tmp;

  int _N; // Number of Hidden State
  int _M; // Number of Emission Symbol

  void randomize();
  void EStep(const Observation& o, 
			 std::vector<std::vector<std::vector<double> > >& gamma,
			 std::vector<std::vector<double> >& sum_of_gamma);

  void MStep(const Observation& o, 
			 const std::vector<std::vector<std::vector<double> > >& gamma,
			 const std::vector<std::vector<double> >& sum_of_gamma);
  double likelihood(const Observation& o);
  double log_likelihood(const Observation& o);
  void check_parameters_consistency();
};

#endif
