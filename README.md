# SimpleHMM

# Objective
To practice writing codes about Forward-Backward Algorithm.

# Current Features
Implemented three algorithm including viterbi, forward(backward), and forward-backward algorithm.

# Usage
	./wscript configure
	./wscript build
  
# Example
	# simple viterbi example
	./build/default/src/simple_hmm_example 
	wget http://www.phontron.com/data/nonparametric-exercise-ja-1.tar.gz
	tar zxvf nonparametric-exercise-ja-1.tar.gz
	mv nonparametric-exercise-ja-1 data
	cd data
	ruby -e 'output_file = File.open("gold.txt", "w"); File.open("wiki-sample.wordpart") {|file| file.each{|line| output_file.print line.chomp! + " "}}; output_file.puts'
	# see command line option
	 ./build/default/src/simple_hmm -helppackage
	 ./build/default/src/simple_hmm -ITER=10 > data/result.txt
	 ./grade-bayes-hmm.pl gold.txt result.txt
# TODO
* Treatment with unkown words
* 