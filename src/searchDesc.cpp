#include <Rcpp.h>
using namespace Rcpp;

// Date:      9/13/2015
// Creator:   Sven Nelson
// Function:  searchDesc
// C++ function for R called searchDesc takes searchTerms and newFrame
// searchTerms is
// newFrame is
// source with sourceCpp("~/Dropbox/00. Code files/searchDesc.cpp")
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// Define contains function: Make sure this ignores case

// [[Rcpp::export]]
bool contains(const std::string & str, const std::string substr) {
  if(str.size()<substr.size()) {
    return false;
  }

  for(int i=0; i<str.size(); i++) {
    if(str.size()-i < substr.size()) {
      return false;
    }

    bool match = true;
    for(int j=0; j<substr.size(); j++) {
      if(str.at(i+j) != substr.at(j)) {
        match = false;
        break;
      }
    }
    if(match) {
      return true;
    }
  }
  return false;
}

// [[Rcpp::export]]
int searchDesc(std::vector< std::string > desc, std::vector< std::string > searchTerms) { // types fixed
  int searchTermHits = 0;
  //Rcpp::DataFrame newFrame // doesn't work for some reason
  //std::vector< std::string > desc // works

  //std::vector< std::string > desc = newFrame["DESC"];

  int n = desc.size();
  int o = searchTerms.size();
  bool hit;
  for (int i = 0; i < n; ++i) {
    hit = false;
    for (int j = 0; j < o; ++j) {
    //while (hit == false) {
      if (hit == true) {
        break; // break exits loop (continue goes to next number in loop)
      } else if (contains(desc[i], searchTerms[j])) { // if a searh term in contained in the DESC column of i^th row
        ++searchTermHits;
        hit = true;
      }
    }
  }
  return searchTermHits;
}


// # newFrame is geneFrame of those not already categorized based on AGI match
// # rowname is probe_id, then ACCNUM, SYMBOL, DESC
// # searchTermHits <- searchDesc(newFrame, searchTerms) # returns searchTermHits
// for (j in 1:nrow(newFrame)) { # data.frame[row, column]
//   if (i!=27 & !is.null(searchTerms) & any(sapply(searchTerms, grepl, newFrame$DESC, ignore.case=FALSE)[j,])) {
//     searchTermHits <- searchTermHits + 1
// #print("Entered if statement, Counting hits")
// #print(searchTermHits)
//   } else if (i==27) {
//     searchTermHits <- searchTermHits + sum(newFrame$DESC == "expressed protein") + sum(newFrame$DESC == "hypothetical protein")
//   }
// }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
searchDesc(myAnnot[,3], c("mitochon", "integrase","ribosom","gibberellin"))
*/
