#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector find_greater_equal_than(NumericVector data, NumericVector values) {
    int value_size = values.length();
    int in_size = data.length();

    int idx = 0;

    NumericVector res(value_size);
    cout << "value_size: " << value_size << endl;
    cout << "in_size: " << in_size << endl;
    for (int i = 0; i < value_size; i++) {
        while (idx < in_size && data[idx] < values[i])
            idx++;
        res[i] = idx + 1;
    }
    return(res);
}

// [[Rcpp::export]]
NumericVector find_dense_min(NumericVector values, int i_max) {
    int i, num_values = values.length();
    NumericVector res(2);

    for (i = i_max; i > 0; i--)
        if (values[i-1] > values[i])
            break;
    res[0] = i + 1;

    for (i = i_max; i < num_values-1; i++)
        if (values[i+1] > values[i])
            break;
    res[1] = i + 1;

    return(res);
}

// [[Rcpp::export]]
NumericVector rect_unique(NumericMatrix data, NumericVector order, NumericVector diff) {
    int i, j, io, jo;
    int num_row = data.nrow();
    int num_diff = diff.length();
    bool is_unique;

    if (2*num_diff != data.ncol()) {
        Rcout << "Error: number of columns in data must be twice the number of differences" << endl;
        throw std::exception();
        return NumericVector();
    }

    NumericVector keep(num_row);

    for (i = 0; i < num_row; i++) {
        io = order[i];
        keep[io] = 1;
        for (j = 0; j < i; j++) {
            jo = order[j];
            is_unique = 0;
            for (int k = 0; k < num_diff; k++) {
                is_unique = is_unique || data(io, 2*k) - data(jo, 2*k+1) > diff[k] || data(jo, 2*k) - data(io, 2*k+1) > diff[k];
            }
            if (keep[jo] && !is_unique) {
                keep[io] = 0;
                break;
            }
        }
    }

    return(keep);
}

// [[Rcpp::export]]
IntegerVector findRangeIndices(NumericVector x, NumericVector bounds) {
  // Check that bounds has exactly two elements
  if (bounds.size() != 2) {
    stop("Bounds vector must have exactly two elements.");
  }

  // Assign the first element as lower and the second as upper
  double lower = bounds[0];
  double upper = bounds[1];

  // Initialize vector to store the indices
  IntegerVector indices;

  // Use lower_bound to find the starting point for values within the range
  auto start_itr = std::lower_bound(x.begin(), x.end(), lower);

  // Use the iterator from lower_bound as the starting point for upper_bound
  auto end_itr = std::upper_bound(start_itr, x.end(), upper);

  // Calculate indices based on iterator positions
  for (auto itr = start_itr; itr != end_itr; ++itr) {
    // Convert to 1-based index for R
    indices.push_back(std::distance(x.begin(), itr) + 1);
  }

  return indices;
}

// [[Rcpp::export]]
IntegerVector findRangeIndicesDecreasing(NumericVector x, NumericVector bounds) {
  if (bounds.size() != 2) {
    stop("Bounds vector must have exactly two elements.");
  }

  double lower = bounds[1]; // Note the order is swapped
  double upper = bounds[0]; // because the vector is decreasing

  IntegerVector indices;

  // Use std::lower_bound with a greater-than comparison to find the starting point
  auto start_itr = std::lower_bound(x.begin(), x.end(), upper, std::greater<double>());

  // Use the iterator from lower_bound as the starting point for upper_bound
  auto end_itr = std::upper_bound(start_itr, x.end(), lower, std::greater<double>());

  // Calculate indices based on iterator positions
  for (auto itr = start_itr; itr != end_itr; ++itr) {
    indices.push_back(std::distance(x.begin(), itr) + 1); // Convert to 1-based index for R
  }

  return indices;
}

/*** R
# Test code in R

# Create a sorted numeric vector to test
test_vector <- sort(runif(100, min = 0, max = 100), decreasing = TRUE)

# Define the bounds for the range you want to find indices for
bounds <- c(60, 30)  # Lower and upper bounds

# Use the function to find indices of values within the specified range
indices <- findRangeIndicesDecreasing(test_vector, bounds)

# Display the test vector, the bounds, and the resulting indices
print(test_vector)
print(paste("Lower bound:", bounds[1], "Upper bound:", bounds[2]))
print(indices)

# Verify by showing the actual values from the test vector that are in the range
print(test_vector[indices])
*/
