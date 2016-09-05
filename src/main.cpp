/* 
 * File:                        main.cpp
 * Author:                      Aline Castendiek
 * Student ID:                  768297
 * Date:                        30/10/15
 * 1st operating system:        Linux [Ubuntu 3.13.0-37-generic]
 * 2nd operating system:        Mac OS X [El Capitan 10.11]       
 * 1st Compiler:                clang [3.4]
 * 2nd Compiler:                g++ [4.8.4]
 * Doxygen version:             1.8.6          
 */


#ifndef __MAIN_CPP__
#define __MAIN_CPP__

#include <iostream>
#include "../include/HMM.hpp"


int main (int argc, char* argv[])
{
  if (argc == 3) {

    /// HMM instance constructed from file in command line. Debugging is enabled
    HMM hmm(argv[1],true);

    // Computes probability of observation sequence read in from file in command line
    hmm.get_forward_prob(argv[2]);

    // Computes most likely sequence of hidden states for observation sequence 
    // read in from file in command line
    //hmm.get_most_likely_seq(argv[2]);
  }

  else {
    std::cerr << 
      "USAGE: main.cpp hmm_file test_file\n"
      "hmm_file: has to match the following syntax and order:\n"
      "1. List the transitions including transition probabilities:\n" 
      "state  state  probability. e.g.: 0  0  0.3.\n"
      "Use a line of tildes (~~~) to end the section.\n"
      "2. List the states and their string representations:\n"
      "state  representation. e.g.: 1  NN \n"
      "Use a line of tildes (~~~) to end the section.\n"
      "3. List the observations including observation probabilities:\n"
      "state  observation  probability. e.g.: 0  otter  0.4 \n"
      "Use a line of tildes (~~~) to end the section.\n"
      "4. List the initial probabilities for a state:\n"
      "state  probability. e.g.: 1  0.1 \n"
      "Give only ONE transition, observation, initial probability or "
      "state-string representation per line.\n" 
      "Separate fields from each other using tabs.\n" 
      "In addition, make sure that all states are natural numbers.\n"
      "An example HMM file is included in /examples.\n\n"

      "test_file: has to contain one word per line. "
      "Each beginning and end of a sentence in the data must be marked "
      "with an additional <BOS> and respectively <EOS> tag. "
      "Some example files are included in /data.\n";
    exit(1);
  }
  return 0;
}

#endif 