#include <Rcpp.h>
#include <string>
#include <iostream>
#include <sstream>
//#include <regex>
#include <vector>
#include <iterator>
#include <map>
#include <exception>
#include <stdexcept>

using namespace std;

//internal C++ functions. not exported.
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}
 
std::vector<std::string> split(const std::string &s, char delim) {
  //#Function to split a string into a vector of strings.
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

bool cpp_eval_bool(vector<string> arule, vector<string> irule, vector<string> mvar, vector<bool> val) {
  //#Function to evaluate the Boolean rule of one gene with one rule at a time. Return a Boolean value for that gene.
  //Note that this function assumes the elements in mvar and val are all ordered.
  
  //Construct a map object between mvar and val.
  std::map<std::string, bool> mmap;
  
  mmap.insert(std::make_pair("0", false)); //handles 0 in the rule.
  for(unsigned int i=0; i<mvar.size(); i++) {
    mmap.insert(std::make_pair(mvar[i], val[i]));
  }
  
  //Replace rules with values.
  std::vector<bool> arule_bool; 
  for(unsigned int  i=0; i<arule.size(); i++) {
    if(arule[i].find("&") == std::string::npos) { //string::npos=-1, which represents a non-position, and is equivalent to not found.
      //if not containing & in the string. 
      arule_bool.push_back(mmap[arule[i]]);
      
    } else {
      //if containing & in the string. evaluate both term and return only a single value.
      std::vector<std::string> tmp = split(arule[i], '&'); //split up the A&B term into A, B. must be single quote, as the function only take single char as delimiter.
      
      bool tmp_ans = true; //initialised var.
      for(unsigned int j=0;j<tmp.size(); j++) {
        tmp_ans = tmp_ans && mmap[tmp[j]];
      }
      
      arule_bool.push_back(tmp_ans);
    }
  }
  
  std::vector<bool> irule_bool; 
  for(unsigned int i=0; i<irule.size(); i++) {
    if(irule[i].find("&") == std::string::npos) { //string::npos=-1, which represents a non-position, and is equivalent to not found.
      //if not containing & in the string. 
      irule_bool.push_back(mmap[irule[i]]);
      
    } else {
      //if containing & in the string. evaluate both term and return only a single value.
      std::vector<std::string> tmp = split(irule[i], '&'); //split up the A&B term into A, B. must be single quote, as the function only take single char as delimiter.
      
      bool tmp_ans = true; //initialised var.
      for(unsigned int j=0;j<tmp.size(); j++) {
        tmp_ans = tmp_ans && mmap[tmp[j]];
      }
      
      irule_bool.push_back(tmp_ans);
    }
  }
  
  //evaluate both arule and irule to give a single final boolean value.
  bool final_ans = false;
  for(unsigned int i=0; i<arule_bool.size(); i++) {
    if(arule_bool[i] == true) { //if any activator is true, then final ans will be true.
      final_ans = true;
      break;
    }
  }
  
  for(unsigned int i=0; i<irule_bool.size(); i++) {
    if(irule_bool[i] == true) { //if any inhibitor is true, then final ans will be false.
      final_ans = false;
      break;
    }
  }
  
  return final_ans;
}

//external C++ functions. exported.
//' @title Simulate a Boolean model.
//' 
//' @description
//' (&&&Not for public use&&&)This function simulates the Boolean model using an initial state. For use within simulate_model(). Returns a matrix of full asynchronous state space.
//' 
//' @param bmodel list. A list of 4 lists created by decreate_model().
//' @param fstate data frame. It must have been initialised by initialise_data(), and has gene names as column names. Must contain only 1 row.
//' @param verbose logical. Indicates whether to output progress.
// [[Rcpp::export]]
Rcpp::List rcpp_simulate(Rcpp::List bmodel, Rcpp::LogicalVector fstate, bool verbose=false) {
  //(1) First step must be to convert all R data structures into C++ data structures for easier manipulation.
  std::vector<bool> first_state = Rcpp::as<std::vector<bool> >(fstate);
  std::vector<std::string> model_gene = Rcpp::as<std::vector<std::string> >(bmodel["gene"]);
  std::vector<std::string> model_var = Rcpp::as<std::vector<std::string> >(bmodel["var"]);

  std::vector< std::vector<std::string> > model_actrule; 
  std::vector< std::vector<std::string> > model_inhrule;
  
  std::vector<std::string> tmp_string;
  //fill in model_act_rule vector of vector.
  for(unsigned int  i=0; i<model_gene.size(); i++) {
    Rcpp::List tmp_list = bmodel["act_rule"];
    tmp_string = Rcpp::as<std::vector<std::string> >(tmp_list(i));
    model_actrule.push_back(tmp_string);
  }
  
  //fill in model_inh_rule vector of vector.
  for(unsigned int  i=0; i<model_gene.size(); i++) {
    Rcpp::List tmp_list = bmodel["inh_rule"];
    tmp_string = Rcpp::as<std::vector<std::string> >(tmp_list(i));
    model_inhrule.push_back(tmp_string);
  }
  
  //(2) Start simulation cycle.
	int num_gene = first_state.size();
  std::vector< std::vector<bool> > next_loop;
  next_loop.push_back(first_state);
  std::vector< std::vector<bool> > all_state;
	all_state.push_back(first_state);
	std::vector< std::vector<bool> > steady_state;
	int lvl = 0;
	
  //std::copy(all_state[0].begin(), all_state[0].end(), std::ostream_iterator<bool>(Rcpp::Rcout, ","));
  
	while(true) //for debug only. change back to true.
	{
		std::vector< std::vector<bool> > part_state;
    lvl += 1;

    //Rcpp::Rcout << "no\n";
    //Rcpp::Rcout << "Next_loop size :" << next_loop.size();
    //Rcpp::Rcout << "\n";

		if(next_loop.size() != 0) {
			for(unsigned int i=0; i<next_loop.size(); i++) { //take each state from the previous level.
				std::vector<bool> old_state = next_loop[i];
        
        if(verbose == true) {
          Rcpp::Rcout << "\rSearching at level " << lvl << ", cell state: ";
          std::copy(old_state.begin(), old_state.end(), std::ostream_iterator<bool>(Rcpp::Rcout, ","));
        }
        
        int steady_count = 0; //for checking steady state.
        for(signed int j=0; j<num_gene; j++) { //check for each gene if it changes in this state.
          std::vector<bool> new_state = old_state;
          
          bool ans = false;
          //Checking for the case where both arule and irule are '0'.
          //In this case, the val should remain the same, and should not change. Only change the value if both arule and irule are not '0'.
          //Equal to 0 indicates a match. A mismatch is -1.
          if(model_actrule[j][0].compare("0") != 0 || model_inhrule[j][0].compare("0") != 0) {
            ans = cpp_eval_bool(model_actrule[j], model_inhrule[j], model_var, old_state);
            
            if(old_state[j] != ans) {
              new_state[j] = ans;
              part_state.push_back(new_state);
            } else {
              steady_count += 1;
            }
            
            if(steady_count == num_gene) //if all genes stay the same, then it is a point steady state.
            {
              steady_state.push_back(old_state);
            } 
          } 
        }
			}
		}
    
    //To replace the R unique() called on part_state. all 3 std functions called below are smart to work on vector of vectors.
    std::sort(part_state.begin(), part_state.end()); //NOTE that std::unique only works on consecutive duplicates, so sorting is a must!
		std::vector< std::vector<bool> >::iterator unique_point;
		unique_point = std::unique(part_state.begin(), part_state.end()); //unique does not change the data. it creates a new copy of sorted data with all the duplicates placed at the end.
    part_state.erase(unique_point, part_state.end());
    
    //To concatenate two vectors of vectors.
    std::vector< std::vector<bool> > tmp_df = all_state;
    tmp_df.insert(tmp_df.end(), part_state.begin(), part_state.end());
    
    //Check for duplicated vector in tmp_df.
    std::sort(tmp_df.begin(), tmp_df.end());
    std::vector< std::vector<bool> >::iterator duplicated_point;
    duplicated_point = std::adjacent_find(tmp_df.begin(), tmp_df.end());
    
    //Rcpp::Rcout << "Size of before part_state : " << part_state.size() << "\n";
    
    if(duplicated_point != tmp_df.end()) { //if there is duplication
      for(unsigned int j=0; j<all_state.size(); j++) {
        std::vector< std::vector<bool> >::iterator x;
        x = std::find(part_state.begin(), part_state.end(), all_state[j]);
        if(x != part_state.end()) { //if found. std::end(part_state)
          part_state.erase(x);
        }
      }
    }

    //Rcpp::Rcout << "Size of after part_state : " << part_state.size() << "\n";
/*
    //debugging line.
    for(int k=0; k<part_state.size(); k++) {
      std::copy(part_state[k].begin(), part_state[k].end(), std::ostream_iterator<bool>(Rcpp::Rcout, ","));
      Rcpp::Rcout << "\n";
    }
*/

    //Pass for next iteration, and save the unique states into all states.
    next_loop = part_state;
    all_state.insert(all_state.end(), next_loop.begin(), next_loop.end());
    
    //Rcpp::Rcout << "Size of all_state : " << all_state.size() << "\n";
    
    //Final breakout condition. Break out of loop when there is no more unique state to go.
    if(next_loop.size() == 0) {
      if(verbose) {
        Rcpp::Rcout << "End of simulation.\n";
      }
      break;
    }
	}
  if(verbose == true) {
    Rcpp::Rcout << "\n"; //to make the output cleaner.
  }
  
  //Convert C++ objects back to R object.
  int state_list_nrow = all_state.size();
  int steady_list_nrow = steady_state.size();
  Rcpp::List state_list(state_list_nrow);
  Rcpp::List steady_list(steady_list_nrow);
  
  //Rcpp::Rcout << "Size of all_state : " << all_state.size() << "\n";
  //Rcpp::Rcout << "Size of steady_state : " << steady_state.size() << "\n";
  
  for(int i=0; i<state_list_nrow; i++) {
    state_list[i] = Rcpp::wrap(all_state[i].begin(), all_state[i].end());
  }
  
  for(int i=0; i<steady_list_nrow; i++) {
    steady_list[i] = Rcpp::wrap(steady_state[i].begin(), steady_state[i].end());
  }
  
  //better to use List rather than matrix for returning.
  return Rcpp::List::create(Rcpp::Named("statespace") = state_list,
                            Rcpp::Named("steady") = steady_list);
}
