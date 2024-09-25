// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;
using namespace std;

void printVector(IntegerVector inputVector) {
  Rcpp::IntegerVector::iterator it;
  for (it = inputVector.begin(); it != inputVector.end(); ++it) {
    Rcpp::Rcout << *it << " ";
  }
  Rcpp::Rcout << std::endl;
}

// void printVector(NumericVector inputVector) {
//   Rcpp::NumericVector::iterator it;
//   for (it = inputVector.begin(); it != inputVector.end(); ++it) {
//     Rcpp::Rcout << *it << " ";
//   }
//   Rcpp::Rcout << std::endl;
// }

List createListFromIndex(IntegerVector index) {
  int n = index.size()-1; // Number of elements in the index vector
  List myList(n); // Create a list with n elements
  
  for (int i = 0; i < n; i++) {
    int len = index[i+1]-index[i]; // Length of the current integer vector
    IntegerVector intVector = seq_len(len+1)-1; // Create an integer vector of values from 1 to len
    myList[i] = intVector; // Assign the integer vector to the i-th element of the list
  }
  
  return myList; // Return the list
}

IntegerVector in_range2(NumericVector x, NumericVector y, double lower, double upper, IntegerVector indices, double thr = 0) {
//  Rcout<<indices[0]<<": "<<lower<<" "<<upper<<" ";
  return indices[x >= lower & x <= upper  & y > thr];
}

IntegerVector in_range(NumericVector x, NumericVector y, double lower, double upper, double thr = 0) {
  // Use lower_bound to find the starting point for values within the range
  auto start_itr = std::lower_bound(x.begin(), x.end(), lower);
  // Use the iterator from lower_bound as the starting point for upper_bound
  auto end_itr = std::upper_bound(start_itr, x.end(), upper);
  // Initialize vector to store the indices
  IntegerVector indices;
  for (auto itr = start_itr; itr != end_itr; ++itr) {
    // Convert to 0-based index for C++ and check if y[index] > 0
    int index = std::distance(x.begin(), itr);
    if (y[index] > thr) {
      indices.push_back(index);
    }
  }
//  Rcout<<lower<<" "<<upper<<" "<<indices<<endl;
  return indices;
}

// [[Rcpp::export]]
std::vector<List> find_peak_targets(NumericVector v,
                                    NumericVector mz,
                                    NumericVector mzd,
                                    IntegerVector scan_index,
                                    IntegerVector ov,
                                    int num,
                                    double thr = 0,
                                    int n_skip = 0) {
  int nv = v.size();
  int n_gap = 0;
  int ns = scan_index.size() - 1;
  
  List indices = createListFromIndex(scan_index);
  NumericVector mz_min = mz - mzd;
  NumericVector mz_max = mz + mzd;
  
  
  NumericVector temp_mz(ns);
  NumericVector temp_intensity(ns);
  temp_mz.fill(0);
  temp_intensity.fill(0);
  
  std::vector<List> res;
  int num_targets = 0;
//  int nc = 0;

  Progress p(nv, true); // create a progress bar

  // Rcpp::Timer timer;
  for (IntegerVector::iterator iov = ov.begin(); iov != ov.end(); ++iov) {
    if (Progress::check_abort() ) return res; // check for user interrupts
    p.increment(); // update progress

    if (v[*iov] == 0) {
      continue;
    }
//    nc += 1;
    NumericVector target_mzs = clone(temp_mz);
    NumericVector target_ints = clone(temp_intensity);
    int i_start=0;
    int i_end=ns-1;
    
    int max_idx = *iov;
    int cnt = 0;
    
    // Rcout<<max_idx<< " \n";
    
    double mz_min_i = mz_min[max_idx];
    double mz_max_i = mz_max[max_idx];
    double mz_t = mz[max_idx];
    
    int tarteg_scan = sum(scan_index < max_idx) - 1;
    
    target_mzs[tarteg_scan] = mz_t;
    target_ints[tarteg_scan] = v[max_idx];
    
    // Rcout<<*iov << " " << v[*iov] << endl;
    
    // Rcout<<v[max_idx] << "----" << max(v) << "===" << max(v[Rcpp::Range(scan_index[tarteg_scan], scan_index[tarteg_scan+1])]) << "\n";
    // Rcout<<mz_t<<": "<<mz_min_i<<"----"<<mz_max_i<<endl;
    
    for (int i=tarteg_scan+1; i<ns; i++) {
      // Rcout<< i << "---" << scan_index[i] << "  " << scan_index[i+1] << "--" << scan_index.size() << "\n";
      
      Range irg = Rcpp::Range(scan_index[i], scan_index[i+1]);
      NumericVector mzv = mz[irg];
      NumericVector vv = v[irg];
//     printVector(indices[i]);
      IntegerVector idx_in_range = in_range(mzv, vv, mz_min_i, mz_max_i, thr);
//      IntegerVector idx_in_range = in_range2(mzv, vv, mz_min_i, mz_max_i, indices[i], thr);
//      Rcout << idx_in_range2 << endl;

      NumericVector mzv2 = mzv[idx_in_range];
      NumericVector vv2 = vv[idx_in_range];
      
      // printVector(mzv2);
      
      if (idx_in_range.size() > 0) {
        int i_target = Rcpp::which_min(Rcpp::abs(mzv2 - mz_t));
        target_mzs[i] = mzv2[i_target];
        target_ints[i] = vv2[i_target];
        // Rcout << "mz:" <<mzv2[i_target] <<" - "<< mz[scan_index[i] + idx_in_range[i_target]] << "ints:" <<vv2[i_target] <<" - "<< v[scan_index[i] + idx_in_range[i_target]]<< "\n";
        v[scan_index[i] + idx_in_range[i_target]] = 0;
        n_gap = 0;
        // cnt+=1;
      } else{
        n_gap += 1;
        if (n_gap > n_skip) {
          i_end = i-n_gap;
          break;
        }
      }
    }
    
    // Rcout<<"++++++++++++++++++++++++++++++++++++\n";    
    for (int i=tarteg_scan; i>=0; i--) {
      // Rcout<< i << "---" << scan_index[i] << "  " << scan_index[i+1] << "--" << scan_index.size() << "\n"; 
      
      Range irg = Rcpp::Range(scan_index[i], scan_index[i+1]);
      NumericVector mzv = mz[irg];
      NumericVector vv = v[irg];
      
      IntegerVector idx_in_range = in_range(mzv, vv, mz_min_i, mz_max_i, thr);
//      IntegerVector idx_in_range = in_range2(mzv, vv, mz_min_i, mz_max_i, indices[i], thr);
//      Rcout << idx_in_range2 << "--" << endl;

      NumericVector mzv2 = mzv[idx_in_range];
      NumericVector vv2 = vv[idx_in_range];
      
      if (idx_in_range.size() > 0) {
        int i_target = Rcpp::which_min(Rcpp::abs(mzv2 - mz_t));
        target_mzs[i] = mzv2[i_target];
        target_ints[i] = vv2[i_target];
        // Rcout << "mz:" <<mzv2[i_target] <<" - "<< mz[scan_index[i] + idx_in_range[i_target]] << "ints:" <<vv2[i_target] <<" - "<< v[scan_index[i] + idx_in_range[i_target]]<< "\n";
        v[scan_index[i] + idx_in_range[i_target]] = 0;
        n_gap = 0;
        // cnt+=1;
      } else{
        n_gap += 1;
        if (n_gap > n_skip) {
          i_start = i+n_gap;
          break;
        }
      }
    }
    v[max_idx] = 0;
    
    // printVector(target_mzs[Rcpp::Range(i_start, i_end)]);
    // printVector(target_ints[Rcpp::Range(i_start, i_end)]);
    
    // Rcout<<i_start<<"  "<<i_end << " " << cnt <<" "<< target_mzs[Rcpp::Range(i_start, i_end)].size() <<"\n";
    
//    if (nc % 1000 == 0) {
//      // String str_nc = std::to_string(nc);
//      // timer.step(str_nc);
//
//      NumericVector x = v[v>0];
//      Rcout<<nc<<"\t"<<num_targets << "\t" << x.size() << "\n";
//      // NumericVector tx(timer);
//      // Rcout<<tx;
////      break;
//    }
    
    if(i_end - i_start + 1 >= num){
      List peak_target;
      IntegerVector scan_range = {i_start, i_end};
      Range trg = Rcpp::Range(i_start, i_end);
      peak_target["mzs"] = target_mzs[trg];
      peak_target["ints"] = target_ints[trg];
      peak_target["rg"] = scan_range+1;
      res.push_back(peak_target);
      num_targets += 1;
    }

  }
  
  return res;
}
