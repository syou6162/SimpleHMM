#include "SimpleHMMExample.hpp"

SimpleHMMExample::SimpleHMMExample() {
  {
	std::vector<double> tmp;
	tmp.push_back(0.3);
	tmp.push_back(0.0);
	tmp.push_back(0.4);
	tmp.push_back(0.1);
	tmp.push_back(0.2);
	_a.push_back(tmp);
  }
	
  {
	std::vector<double> tmp;
	tmp.push_back(0.7);
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	tmp.push_back(0.3);
	tmp.push_back(0.0);
	_a.push_back(tmp);
  }
	
  {
	std::vector<double> tmp;
	tmp.push_back(0.3);
	tmp.push_back(0.2);
	tmp.push_back(0.1);
	tmp.push_back(0.2);
	tmp.push_back(0.2);
	_a.push_back(tmp);
  }
	
  {
	std::vector<double> tmp;
	tmp.push_back(0.5);
	tmp.push_back(0.1);
	tmp.push_back(0.0);
	tmp.push_back(0.4);
	tmp.push_back(0.0);
	_a.push_back(tmp);
  }
	
  {
	std::vector<double> tmp;
	tmp.push_back(0.6);
	tmp.push_back(0.3);
	tmp.push_back(0.0);
	tmp.push_back(0.1);
	tmp.push_back(0.0);
	_a.push_back(tmp);
  }

  {
	std::vector<double> tmp;
	tmp.push_back(0.6);
	tmp.push_back(0.1);
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	tmp.push_back(0.3);
	_b.push_back(tmp);
  }
	
  {
	std::vector<double> tmp;
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	tmp.push_back(1.0);
	tmp.push_back(0.0);
	_b.push_back(tmp);
  }
	
  {
	std::vector<double> tmp;
	tmp.push_back(0.1);
	tmp.push_back(0.2);
	tmp.push_back(0.7);
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	_b.push_back(tmp);
  }
	
  {
	std::vector<double> tmp;
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	tmp.push_back(1.0);
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	_b.push_back(tmp);
  }
	
  {
	std::vector<double> tmp;
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	tmp.push_back(1.0);
	tmp.push_back(0.0);
	tmp.push_back(0.0);
	_b.push_back(tmp);
  }

  _pi.push_back(0.6); 
  _pi.push_back(0.4);
  _pi.push_back(0.0);
  _pi.push_back(0.0);
  _pi.push_back(0.0);

  _M = _b[0].size();
  _N = _pi.size();
};

void SimpleHMMExample::run() {
  SimpleHMMExample::Observation o;
  o.push_back(0);
  o.push_back(1);
  o.push_back(2);
  o.push_back(3);
  o.push_back(4);
  
  viterbi(o);
};
