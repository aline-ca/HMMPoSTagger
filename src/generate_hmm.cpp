/* 
 * File:                        generate_hmm.hpp
 * Author:                      Aline Castendiek
 * Student ID:                  768297
 * Date:                        16/10/16
 * 1st operating system:        Mac OS X [El Capitan 10.11.5]  
 * 2nd operating system:        Linux [Ubuntu 3.13]     
 * 1st Compiler:                clang [3.4]
 * 2nd Compiler:                g++ [4.8.4]
 * Doxygen version:             1.8.11          
 */


#ifndef __GENERATE_HMM_CPP__
#define __GENERATE_HMM_CPP__

#include <iostream>
#include "../include/HMMGenerator.hpp"

int main (int argc, char* argv[])
{
   if (argc == 3) {

      // HMM generator instance that will compute probabilities without smoothing
      // HMMGenerator hmm_generator(argv[1],argv[2],false);

      // HMM generator instance that will perform standard add-one smoothing
      //HMMGenerator hmm_generator(argv[1], argv[2]);

      // HMM generator instance that will perform smoothing with own specified summand
      HMMGenerator hmm_generator(argv[1],argv[2],true,0.75);
   }

   else {
      std::cerr <<
              "USAGE: generate_hmm training_file output_file \n"
              "training_file: Training file that contains one word and its pos tag per line, seperated by a tabulator.\n"
              "output_file: File in which the generated model shall be saved.\n\n";
      exit(1);
   }
   return 0;
}

#endif