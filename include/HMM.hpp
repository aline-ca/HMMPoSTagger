/* 
 * File:                        HMM.hpp
 * Author:                      Aline Castendiek
 * Student ID:                  768297
 * Date:                        30/10/15
 * 1st operating system:        Linux [Ubuntu 3.13.0-37-generic]
 * 2nd operating system:        Mac OS X [El Capitan 10.11]       
 * 1st Compiler:                clang [3.4]
 * 2nd Compiler:                g++ [4.8.4]
 * Doxygen version:             1.8.6          
 */

/* 
  TODO: 
    - Implement logarithmic probabilities for higher precision and to prevent numeric overflow
    - Escape line separator sign (~)
*/

#ifndef __HMM_HPP__
#define __HMM_HPP__

#include <boost/numeric/ublas/matrix.hpp>   // boost::matrix
#include <boost/numeric/ublas/io.hpp>       // for printing matrices
#include <boost/tokenizer.hpp>              // boost::tokenizer
#include <string>                           // std:string
#include <map>                              // std:map
#include <tuple>                            // std::tuple (requires C++11)
#include <iostream>                         // std::cout
#include <fstream>                          // std::fstream
#include <cstdlib>                          // abort(), exit()

#define UNSEENWORDPROB                      0.01

using namespace boost;

/** 
  @brief The HMM class is used to read in a trained HMM as tsv file and a file containing test data. 
         Then either the viterbi algorithm for finding the most probable hidden state sequence 
         or the forward algorithm for getting the forward probability of a sentence can be executed.
  @note The class can either process an HMM generated from the HMMGenerator class or a hand-written
        HMM.
  @note Although add-one is used to smooth unseen transitions for all words in the training data, the
        HMM can NOT yet handle completely unseen words in the test data. Unseen words will result in
        a zero probability for the complete test data.
  @note The tilde sign is used as sezparator sign between the different phases of reading in the HMM.
        It was specifically chosen because it did not occur in any training data, including the whole
        section 02-21 of the Wall Street Journal. If a tilde occurs in the training data, the 
        HMMGenerator class will be able to integrate it in the HMM correctly, but the file cannot be
        read in by the HMM class.
*/
class HMM
{
 public: // Types 
  bool                                      debug_mode;           ///< Prints out intermediate steps. Set to false by default.

 private: // Typedefs
  typedef std::string                       Word;                 ///< Observed word as string
  typedef std::string                       HiddenSymbol;         ///< Hidden symbol/state in which word was emitted
  typedef int                               Index;                ///< Index used for observation and state mapping 
  typedef unsigned                          State;                ///< Representation of hidden state as number   
  typedef std::map<Word,Index>              WordIndexMap;         ///< For mapping observed word to an index
  typedef std::map<State,HiddenSymbol>      StateSymbolMap;       ///< Maps a state number to its string representation 
  typedef std::map<HiddenSymbol,State>      SymbolStateMap;       ///< Maps a state as string representation to its number
  typedef std::vector<std::string>          StringVector;         ///< String vector (mostly for saving observations)
  typedef std::vector<StringVector>         StringVectorVector;
  typedef std::vector<State>                StateVector;          ///< Vector of state numbers
  typedef std::vector<double>               DoubleVector;         ///< Vector of double values
  typedef std::tuple<unsigned,unsigned>     UnsignedTuple;        ///< Represents a transition from one state to another
  typedef std::vector<UnsignedTuple>        TupleVector;          ///< Vector of state-tuples
  typedef std::map<State,double>            StateDoubleMap;       ///< Maps a state to a probability (init probs)
  
  typedef numeric::ublas::matrix<double>    Matrix;               ///< Matrix of double values (probabilities)
  typedef numeric::ublas::matrix<State>     StateMatrix;          ///< Matrix of states (for backpointer trellis)
  
  typedef tokenizer< char_separator<char> > Tokenizer;            ///< Boost tokenizer for reading in data


 private: // Types
  Matrix          transition_matrix;    ///< Maps two states to their transition probability
  Matrix          observation_matrix;   ///< Maps a state and a word index to their observation probability
  WordIndexMap    word_index_map;       ///< Maps observed word to an index (used to represent word in the observation matrix)
  Index           current_index;        ///< Generated index for mapping words to word_index_map
  StateSymbolMap  state_symbol_map;     ///< Maps the unsigned representation of a state to its string representation   
  SymbolStateMap  symbol_state_map;     ///< Maps the string representation of a state to its unsigned representation
  StateDoubleMap  init_probs;           ///< Maps a state to its initial probability
  unsigned        num_of_states;        ///< Number of states
  unsigned        num_of_observations;  ///< Number of observations

  
 public: // Functions
  /**
    @brief Constructor that reads in a HMM from a tsv file.
    @param file Reference to the file that contains the HMM.
    @param d Debug mode that can be enabled if desired. Is set to false by default.
  */
  HMM(const std::string& file, bool d=false) {
    
    debug_mode = d; 
    // Ifstream object constructed from file 
    std::ifstream in(file.c_str()); 
    
    if (in) {
      init();                              // Initialize counters
      read_in(in);                         // Read in HMM from file and build internal data structures
      
      std::cout << "HMM BUILD FROM FILE '" << file << "': SUCCESSFUL.\n";         
      if (debug_mode) { print(); }
    }
    else std::cerr << "Error: Unable to open '" << file << "'\n";  
  }   


  /***************************************
  *          Public Functions            *
  ***************************************/
  
  
  /**
    @brief Gets most likely sequence of hidden states.
    @param file Reference to the file that contains the test data.
    @return String that contains the most likely state sequence.
  */
  const std::string get_most_likely_seq(const std::string& file) 
    { return viterbi(file); }

  /**
    @brief Gets a vector of forward probabilities for all the sentences in the test data.
    @param file Reference to the file that contains the test data.
    @return Vector of forward probabilities. 
  */
  const DoubleVector get_forward_probs(const std::string& file) 
    {   return compute_forward_probs(file); }


  private: // Functions

  /***************************************
  *        Initializer Functions         *
  ***************************************/
  
  /**
    @brief Initializes an empty matrix by writing zero into all its cells.
    @param m Reference to the matrix that shall be initialized.
  */
  void initialize_matrix(Matrix& m) {
    for (unsigned i = 0; i < m.size1(); ++i) {
      for (unsigned j = 0; j < m.size2(); ++j) {
        m(i,j) = 0;
      }
    }
  }


  /**
    @brief Initializes counter variables.
  */
  void init () {
    current_index = 0;
    num_of_states = 0;
    num_of_observations = 0;
  }


 /*********************************************************************
  *          Functions for Forward and Viterbi Calculations           *
  ********************************************************************/

  /***************************************
  *          Forward Algorithm           *
  ***************************************/
  /**
    @brief Computes the forward probability for a single observed sentence.
    @param observations Reference to the sentence saved in a string vector.
    @return Probability of the observed sequence.
  */
  const double forward(const StringVector& observations) {
      
      const unsigned& N = transition_matrix.size1();   // Number of states 
      const unsigned& T = observations.size();         // Number of time steps in the observation
      Matrix forward_trellis(N,T);                     // Forward trellis
      
      std::cout << "----------------------------------------------------------------------" << "\n";
      std::cout << "COMPUTING FORWARD PROBABILITY FOR FOLLOWING OBSERVATION SEQUENCE:\n";
      for (auto it = observations.begin(); it != observations.end(); ++it) {
  	std::cout << *it << ' ';
      }
      std::cout << "\n\n";
      
      // Initalization: Probability to start in state i * probability to emit 
      // first observation (observations[0]) in state i.  
      for (unsigned i = 0; i < N; ++i) {
        // If we get the invalid word index -1, we know that the word did not occur in the training data.
        // So the probability to start in it is zero. Right now, this will result in a zero probability for
        // all the test data which is not very good...
        if (get_word_index(observations[0]) == -1) { forward_trellis(i,0) = init_probs[i] * UNSEENWORDPROB; }     
        else { forward_trellis(i,0) = init_probs[i] * observation_matrix(i,get_word_index(observations[0])); }
      }
      // Recursion Step (Iterative solution): 
      for (unsigned t = 1; t < T; ++t) {              // t represents current time step
        for (unsigned s = 0; s < N; ++s) {            // s represents current state
          double sum = 0;                             // Value that will be written in current trellis cell
          for (unsigned j = 0; j < N; ++j) {          // j represents previous state
            // Again, test whether word has an entry in the HMM, and if not, use a zero probability.
            if (get_word_index(observations[t]) == -1) {
              sum += forward_trellis(j,t-1) * transition_matrix(j,s) * UNSEENWORDPROB;
            }
            else {                                    // Word has an entry in the HMM
              // Compute previous sum + forward probability at previous trellis position * transition probability from previous
              // to current state * probability to emit observations[t] in current state.   
              sum += forward_trellis(j,t-1) * transition_matrix(j,s) * observation_matrix(s,get_word_index(observations[t]));
            }
          }
          // Write calculated value into current trellis cell. Afterwards, set sum variable back to zero.
          forward_trellis(s,t) = sum;
        }
      }
      if (debug_mode) { std::cout << "Complete trellis:\n" << forward_trellis << "\n";}
      
      // Termination Step: Sum up all final forward probabilties (last cells in trellis).
      double result = 0;                              // Probability that will be returned
      for (unsigned i = 0; i < N; ++i) {
        result += forward_trellis(i,T-1);
      }
      std::cout << "THE FORWARD PROBABILITY IS:\n" << result << "\n";
      return result;
  }

  /***************************************
  *          Viterbi Algorithm           *
  ***************************************/

  /**
    @brief Computes the best path for a given observation sequence and the trained HMM.
           The test data gets read in from a file.
    @param file Reference to the file that contains the test data.
    @return Most probable hidden state sequence as string.
  */
  /// Compute the best path for a given observation sequence and HMM.
  const std::string viterbi(const std::string& file) {
    std::ifstream in(file.c_str());   ///< Ifstream object created from file
    if (in) {
      // Vector containing the sequence of observations
      const StringVector& observations = parse_test_data(in);
      const unsigned& N = transition_matrix.size1();  // Number of states 
      const unsigned& T = observations.size();        // Number of time steps in the observation
      Matrix viterbi(N,T);                            // Viterbi trellis
      StateMatrix backpointer(N,T);                   // Backpointer trellis
      
      std::cout << "----------------------------------------------------------------------" << "\n";
      std::cout << "COMPUTING THE BEST PATH FOR FOLLOWING OBSERVATION SEQUENCE:\n";
      for (auto it = observations.begin(); it != observations.end(); ++it) {
  	std::cout << *it << ' ';
      }
      std::cout << "\n\n";

      // Initalization: Probability to start in state i * probability to emit
      // first observation (observations[0]) in state i.
      for (unsigned i = 0; i < N; ++i) {
        // If we get the invalid word index -1, we know that the word did not occur in the training data.
        // The actual probability to start in it would be 0, but then we will always multiply by this 0
        // and only get zero probabilities. Because of this, use constant UNSEENWORDPROB for multiplication instead.
        if (get_word_index(observations[0]) == -1) { viterbi(i,0) = init_probs[i] * UNSEENWORDPROB; } 
        else {
          viterbi(i,0) = init_probs[i] * observation_matrix(i,get_word_index(observations[0]));
        }
        // (This step is actually not necessary because we will never use the first row
        // of the backpointer trellis because it does not point to any previous state.)
        backpointer(i,0) = 0;
      }

      // Recursion Step (Iterative solution): 
      for (unsigned t = 1; t < T; ++t) {              // t represents current time step
        for (unsigned s = 0; s < N; ++s) {            // s represents current state
          double maximum = 0;                         // Value that will be written in current viterbi trellis cell
          double best_state_prob = 0;                 // Probability of current best state (for backpointer trellis)
          State best_state;                           // Value that will be written in current backpointer trellis cell
          for (unsigned j = 0; j < N; ++j) {          // j represents previous state
            // Again, test whether word has an entry in the HMM.
            if (get_word_index(observations[t]) == -1) { 
              maximum = std::max(maximum, viterbi(j,t-1) * transition_matrix(j,s) * UNSEENWORDPROB);
            }
            else {
              // Compute previous sum + forward probability at previous trellis position * transition probability 
              // from previous to current state * probability to emit observations[t] in current state.
              maximum = std::max(maximum, viterbi(j,t-1) * transition_matrix(j,s) * observation_matrix(s,get_word_index(observations[t])));
            }

            // Write calculated value into current viterbi trellis cell. Afterwards, set sum variable back to zero.
            viterbi(s,t) = maximum;
            // Saves the current probability (necessary for filling in the backpointer trellis)
            double current_prob = viterbi(j,t-1) * transition_matrix(j,s);
            // If the current probability is higher than the best probability, make it the new best probability
            // and set the current state as best state:
            if (current_prob > best_state_prob) {
              best_state_prob = current_prob;
              best_state = j;
            }
            // Then write best state into current backpointer trellis cell:
            backpointer(s,t) = best_state;
          }
        }
      }

      if (debug_mode) { 
        std::cout << "Viterbi trellis:\n" << viterbi << "\n";
        std::cout << "Backpointer trellis:\n" << backpointer << "\n";
      }

      // Get best last state by iterating over last column of viterbi matrix. 
      // Necessary because viterbi and backpointer trellis have the same size, but position [s,1]
      // in backpointer trellis refers to state position [s,0] in viterbi trellis, so we have to save the 
      // last state seperately.  
      State best_last_state;                      // Best state in last row of viterbi matrix
      float best_last_prob = 0;                   // Best last probability
      for (unsigned i = 0; i < N; ++i) {
        if (viterbi(i,T-1) > best_last_prob) {
          best_last_prob = viterbi(i,T-1);
          best_last_state = i;
        }
      }

      // Backtrace: Follow backpointers in matrix bottum-up to reconstruct best path.
      StateVector path;                           // Vector of states that stores the best path 
      path.push_back(best_last_state);            // First append the last best state we just computed
      State current_state = best_last_state;     // Current position in backpointer matrix
      // For all time steps, starting in last column and ignoring the very first time step
      // (because it does not point to any cell):

      for (int t = T-1; t > 0; --t) {
        // Set next current state and save it in path vector: 
        current_state = backpointer(current_state,t);
        path.push_back(current_state);
      }

      // Backtrace path is returned backwards, therefore reverse it:
      std::reverse(path.begin(),path.end());

      // Get the result sequence by concatenating all string representations of the state path:   
      std::string result = get_string_repr(path);
      
      std::cout << "THE MOST PROBABLE HIDDEN STATE SEQUENCE IS:\n"
                <<  result << "\n"
                << "----------------------------------------------------------------------" << "\n";
      return result;
    }
    else {
      std::cerr << "Error: Unable to open '" << file << "'\n";
      exit(1);
    } 
  }


 /*********************************************************
  *             Functions for reading in HMM              *
  ********************************************************/

  /**
    @brief Reads in HMM from an input stream.
    @param in Reference to the instream object that will be read from.
  */
  void read_in(std::istream& in) {    
    StringVector vec;                       // Vector that reads in lines
    std::string line;                       // Current line
    std::cout << "Call read_in Function\n";
    unsigned separator_counter = 0;         // Counts occurred lines of number signs (separators)
    while (getline(in,line)) {              // Read in each line seperately
          std::cout << line << "\n";

        char_separator<char> sep("\t");       // Define tab field separator
      /// Tokenizer constructed from current line and defined separator:     
      Tokenizer tok(line, sep); 
      
      // Get metadata from first line (The structure of the line is: # States Observations #)
      if(line[0] == '#' && line[line.length()-1] == '#') {
        vec.assign(tok.begin(),tok.end());  // Assign current line to vector    
        build_hmm_matrices(vec);            // Build and initialize transition and observation matrix 
        continue;
      }    
      
      if (line[0] != '~') {                 // If current line has no separator
        vec.assign(tok.begin(),tok.end());  // Assign current line to vector    
        // Ignore all lines that have less than two fields:
        if (vec.size() < 2) continue;
        switch(separator_counter) {
          case 0:                           // No separator occured. Read in transitions.
            get_transitions(vec);
            break;
          case 1:                           // One separator occured. Read in symbol state representations.
            get_symbol_state_repr(vec); 
            break;
          case 2:                           // Two separators occured. Read in observations.
            get_observations(vec);
            break;
          case 3:                           // Three separators occured. Read in initial probabilities.
            get_init_probs(vec);
            break;
          default:
            print_read_in_error();
            exit(1);
        }
      }
      // Else: We found a separator, so we switch to the next case and increment counter.
      else {
        ++separator_counter;
      }
    }
  }

  /**
    @brief Reads in metadata section at the beginning of the file to get the number of states and observations.
           Resizes the transition and observation matrix to this size and initializes them.
    @param v Reference to the vector.
  */
  void build_hmm_matrices(const StringVector& v) {

    // Get number of states and number of observations: 
    num_of_states = atoi(v[1].c_str());         // Convert strings to ints   
    num_of_observations = atoi(v[2].c_str());

    // Set the size of transition and observation matrix:
    transition_matrix.resize(num_of_states,num_of_states);
    observation_matrix.resize(num_of_states,num_of_observations);
    
    // Initialize matrices:
    initialize_matrix(transition_matrix);
    initialize_matrix(observation_matrix);   
  }
  
  
  /**
    @brief Reads in transitions and transition probabilities from a vector.
    @param v Reference to the vector.
  */
  void get_transitions(const StringVector& v) {
    const int& state1 = atoi(v[0].c_str());           // Convert strings to ints
    const int& state2 = atoi(v[1].c_str());
    const double& prob = ::atof(v[2].c_str());        // Convert string to double
    transition_matrix(state1,state2) = prob;
  }


  /**
    @brief Reads in observations and observation probabilities from vector.
    @param v Reference to the vector.
  */
  void get_observations(const StringVector& v) {
    auto it = word_index_map.find(v[1]);              // Pointer to map index of observed word.
    // If observed word is not in map yet, add it and increase index by one;    
    if (it == word_index_map.end()) {
      word_index_map.insert(std::make_pair(v[1],current_index));
      ++current_index;
    }
    const double& prob = ::atof(v[2].c_str());          // Convert string to double
    const State& hidden_state = get_hidden_state(v[0]);
    const Index word_index = get_word_index(v[1]);
    // Write observation into observation matrix:
    observation_matrix(hidden_state,get_word_index(v[1])) = prob;
  }


  /**
    @brief Reads in initial probabilities from vector.
    @param v Reference to the vector.
  */
  void get_init_probs(const StringVector& v) {
    const int& state = atoi(v[0].c_str());            // Convert string to int
    const double& prob = ::atof(v[1].c_str());        // Convert string to double
    init_probs[state] = prob;
  }


  /**
    @brief Reads in state symbol map from vector.
    @param v Reference to the vector.
  */
  void get_symbol_state_repr(const StringVector& v) {
    const int& state = atoi(v[0].c_str());            // Convert string to int
    state_symbol_map[state] = v[1];
    symbol_state_map[v[1]] = state;
  }


  /**
    @brief Prints error message and instructions on how to use the program.
           Gets called if something goes wrong while reading in HMM.
  */
  void print_read_in_error() const 
  {
    std::cerr << 
      "An error occurred while trying to read in the HMM from a tsv file.\n"
      "Please either use the HMMGenerator class to generate an HMM from raw data or "
      "make sure that your hand-written HMM matches the following syntax and order:\n\n"

      "1. List the transitions including transition probabilities:\n" 
      "state  state  probability. e.g.: 0  0  0.3.\n"
      "Use a line of tildes (~~~) to end the section.\n\n"

      "2. List the states and their string representations:\n"
      "state  representation. e.g.: 1  NN \n"
      "Use a line of tildes (~~~) to end the section.\n\n"

      "3. List the observations including observation probabilities:\n"
      "state  observation  probability. e.g.: 0  otter  0.4 \n"
      "Use a line of tildes (~~~) to end the section.\n\n"

      "4. List the initial probabilities for a state:\n\n"
      "state  probability. e.g.: 1  0.1 \n"

      "Give only ONE transition, observation, initial probability or "
      "state-string representation per line.\n" 
      "Separate fields from each other using tabs.\n" 
      "In addition, make sure that all states are natural numbers.\n"
      "An example of a hand-written HMM is included in /examples.\n\n";
  }



  /**************************************************
   *           Private Getter Functions              *
   *************************************************/

  /**
    @brief Gets index for an observed word.
    @param word Reference to observed word.
    @return Index for word. Returns the invalid index -1 if a word has no index.
  */
  const Index get_word_index(const Word& word) const
  {
    auto f = word_index_map.find(word);               // Pointer to index of current word
    if (f != word_index_map.end()) return f->second;  
    else {
      //std::cout << "Note: No HMM entry for '" <<  word << "' found. \n";
      return -1;
    }
  }


  /**
    @brief Gets string representation for a hidden state.
    @param s Reference to the hidden state as number.
    @return Reference to string representation of the hidden state.
  */
  const HiddenSymbol& get_hidden_symbol(const State& s) const
  {
    auto f = state_symbol_map.find(s);
    if (f != state_symbol_map.end()) return f->second;  
    else {
      std::cerr << "ERROR: No string representation for hidden state '" <<  s << "' found. \n";
      exit(1);
    }
  }


  /**
    @brief Get hidden state from its string representation.
    @param str Reference to hidden state as string.
    @return Reference to the hidden state as number.
  */
  const State& get_hidden_state(const HiddenSymbol& str) const
  {
    auto f = symbol_state_map.find(str);
    if (f != symbol_state_map.end()) return f->second;  
    else {
      std::cerr << "ERROR: No hidden state for string representation '" << str << "'' found.\n";
      exit(1);
    }
  }


  /**
    @brief Takes a vector of hidden states as argument and returns the sequence of its
           string representations.
    @param vec Reference to vector of hidden states as numbers.
    @return The same sequence as string representation.
  */
  const std::string get_string_repr(const StateVector& vec) const
  {
    std::string str;              // String that will be returned
    for (auto it = vec.begin(); it != vec.end(); ++it) {
      // If string is empty, write first hidden state into it:
      if (str.empty()) {
        str = get_hidden_symbol(*it);
      }
      // Else append all other hidden states to it:
      else {
        str += (" " + get_hidden_symbol(*it)); 
      }
    }
    return str;
  }


  /**************************************************
   *               Other Functions                 *
   *************************************************/

  /**
    @brief Prints out all properties of the HMM.
  */
  void print() const
  {
    std::cout << "----------------------------------------------------------------------" << "\n";
    std::cout << "HMM PROPERTIES:\n"
          << "Transition Matrix: " << transition_matrix << "\n\n"  
          << "Observation Matrix: " << observation_matrix << "\n\n"
          << "Word Index Map:" << "\n";
    for (auto e = word_index_map.begin(); e != word_index_map.end(); ++e) {
      std::cout << e->first << " = " << e->second << "\n"; 
    }
    std::cout << "\nInitial Probabilities:" << "\n";
    for (auto e = init_probs.begin(); e != init_probs.end(); ++e) {
      std::cout << e->first << " = " << e->second << "\n"; 
    }
    std::cout << "\nState Symbol Map:" << "\n";
    for (auto e = state_symbol_map.begin(); e != state_symbol_map.end(); ++e) {
      std::cout << e->first << " = " << e->second << "\n";
    }
    std::cout << "----------------------------------------------------------------------" << "\n";
  }


  /**
    @brief Parses the test data by reading it in from a file line by line, saving the observed
           sequence in a vector and returning it.
    @param in Reference to the instream object that will be read from.
    @return Vector containing the observation sequence.

  */
  const StringVector parse_test_data(std::istream& in) const
  {
    StringVector vec;                       // Vector that will contain all sentences
    std::string line;                       // Current line
    while (getline(in,line)) {              
      vec.push_back(line);
    }
    return vec;
  }

  
    /**
    @brief Takes a test file as input and splits it up into sentences 
           by looking for the <EOS> tag. Saves the single sentences as a vector
           of string vectors.
    @param file Reference to the file that contains the data that shall be 
           split into sentences.
    @return Vector of string vectors that contains the sentences.
  */
  const StringVectorVector get_sentences_from_file(const std::string& file)
  {
    std::ifstream in(file.c_str());
    if (in) {
        StringVectorVector all_sentences;         // Vector of all sentences
        StringVector sentence;                    // Current sentence
        std::string line;                         // Current line
      
        while (getline(in,line)) {      
            if(line != "<EOS>") {
                sentence.push_back(line);
            }
            else {
                sentence.push_back(line);
                all_sentences.push_back(sentence);
                sentence.clear();
            }        
        }
        return all_sentences;
    }
    else {
        std::cerr << "Error: Unable to open '" << file << "'\n";
        exit(1);
    } 
  }
 
  /**
    @brief Gets vector of forward probabilites by calling the forward() function 
           for every single sentence in the vector.
    @param sentences Vector of string vectors where every string vector represents a single sentence.
    @return Vector of forward probabilities with one entry for every sentence.
  */
    const DoubleVector get_forward_probs_for_vector(StringVectorVector sentences)  
    {       
        DoubleVector all_forward_probs;
        
        for (auto it = sentences.begin(); it != sentences.end(); it++) {
            double current_forward_prob = forward(*it);
            all_forward_probs.push_back(current_forward_prob);
        }        
        return all_forward_probs;
      
    }  
    
    
    /**
    @brief Executes all necessary steps for calculating the forward probability of sentences 
           in a file. Extracts the sentences from the file first and then calculates 
           the forward probability of every sentence seperately.
    @param file Reference to the file that contains the sentences.
    @return Vector of forward probabilities with one entry for every sentence.
  */
    const DoubleVector compute_forward_probs(const std::string& file) 
    {
        StringVectorVector sentences = get_sentences_from_file(file);
        DoubleVector forward_probs = get_forward_probs_for_vector(sentences);
        return forward_probs; 
    }
  
}; // HMM

#endif