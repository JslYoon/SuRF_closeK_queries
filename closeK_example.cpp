#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <chrono>
#include <set>
#include <ctime>
#include <random>

#include "include/surf.hpp"

using namespace surf;
using namespace std;


// Function to find if a query performs a false positive or not; it checks if a value exists within a certain range in source_vec
bool find_key_in(const std::vector<uint64_t>& source_vec, uint64_t left_end, uint64_t right_end) {
    auto lower = std::lower_bound(source_vec.begin(), source_vec.end(), left_end);

    if (lower != source_vec.end() && *lower <= right_end) {
        return true; // Found a value within the range
    }

    return false; // No values found within the range
}


std::vector<std::string> intsToStrings(const std::vector<uint64_t>& intVector) {
    std::vector<std::string> stringVector;
    stringVector.reserve(intVector.size());  // Reserve space to improve efficiency
    std::transform(intVector.begin(), intVector.end(), std::back_inserter(stringVector),
                   [](uint64_t num) { return std::to_string(num); });
    
    return stringVector;
}


std::string intToString(uint64_t intNumber) {
    return std::to_string(intNumber);
}

// To get normal distribution
vector<uint64_t> get_normal_distribution(uint64_t N, double mean, double stddev, uint64_t range_min, uint64_t range_max) {

    std::vector<uint64_t> results;
    results.reserve(N); 

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> dist(mean, stddev);

    for (size_t i = 0; i < N; ++i) {
        double number;
        do {
            number = dist(gen);
            number = (number - mean) / (4 * stddev) * (range_max - range_min) + (range_min + range_max) / 2.0;
        } while (number < range_min || number > range_max); // Repeat if the number is outside the range

        results.push_back(static_cast<uint64_t>(number));
    }
    return results;
}



// To get uniform distribution
std::vector<uint64_t> get_uniform_distribution(uint64_t N, uint64_t range_min, uint64_t range_max) {
    std::vector<uint64_t> v_keys(N, 0);
    std::random_device rd;
    std::mt19937_64 gen(rd()); 
    std::uniform_int_distribution<uint64_t> dist(range_min, range_max);

    for (uint64_t i = 0; i < N; ++i) {
        v_keys[i] = dist(gen);
    }

    return v_keys;
}

// To get exponential distribution
std::vector<uint64_t> get_exponential_distribution(uint64_t N, double lambda, uint64_t range_min, uint64_t range_max) {
    std::vector<uint64_t> v_keys(N, 0);
    std::random_device rd;
    std::mt19937 gen(rd()); 
    std::exponential_distribution<> dist(lambda);

    double max_exp_value = std::log(range_max) / lambda;
    
    for (uint64_t i = 0; i < N; ++i) {
        double exp_value = dist(gen);
        
        
        exp_value = std::min(max_exp_value, exp_value); 
        double scale = exp_value / max_exp_value;
        v_keys[i] = range_min + static_cast<uint64_t>((range_max - range_min) * scale);
    }
    return v_keys;
}




// Helper function for interface display
template<typename T>
uint64_t display_select_vec(vector<T>& vec) {
  uint64_t num = 0;
  bool is_valid = false;
  for(int i = 1; i <= vec.size(); i++) { // display vec
    cout << i << ". " << vec[i-1] << endl;
  }

  while(!is_valid) { // input loop
    cout << "Choose an option (1 to " << vec.size() << "): ";  
    if (cin >> num && num > 0 && num <= vec.size()) { 
      is_valid = true; 
    } else {
      cout << "Invalid input." << endl;
      cin.clear(); 
      cin.ignore(numeric_limits<streamsize>::max(), '\n'); 

    }
  }
  return num;
}

// Until valid number
uint64_t until_number_input(uint64_t min, uint64_t max) {
  
  uint64_t number;
    std::cout << "Enter a number between " << min << " and " << max << ": ";
    while (true) {
        if (std::cin >> number) {
            if (number >= min && number <= max) {
                break;  
            } else {
                std::cout << "Please enter a number within the range " << min << " to " << max << ": ";
            }
        } else {
            cout << "Invalid input." << endl;
            cin.clear();  
            cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
        }
    }
    return number;
}
        



// Function to test surf
void test_surf(uint64_t batch_size, string key_distribution, string query_distribution, uint64_t test_num, uint64_t N) {

  //----------------------------------------
  //GENERATING DATA
  //----------------------------------------

  vector<uint64_t> v_keys_temp;

  if (key_distribution == "uniform") { // uniform distribution
    v_keys_temp = get_uniform_distribution(N, 0, static_cast<uint64_t>(pow(2, 50))-1);

  } else if (key_distribution == "normal") { // normal distribution
    v_keys_temp = get_normal_distribution(N, 100.0, 20.0, 0, static_cast<uint64_t>(pow(2, 50))-1);
  }

  // Sort v_keys for faster verification of false positiviry
  vector<uint64_t> sorted_v_keys = v_keys_temp;
  sort(sorted_v_keys.begin(), sorted_v_keys.end());
  std::vector<std::string> keys = intsToStrings(sorted_v_keys);
  
  //----------------------------------------
  //SuRF CONSTRUCTION
  //---------------------------------------- 
  for(int i = 0; i < keys.size(); i++) {
    cout << "asfd " << keys[i] << endl;
  }
  SuRF* surf = new SuRF(keys);

  //----------------------------------------
  //BUILDING WORKLOAD
  //----------------------------------------

  vector<uint64_t> test_queries;
  if (query_distribution == "uniform") { // uniform distribution
    test_queries = get_uniform_distribution(N, 0, static_cast<uint64_t>(pow(2, 50))-1);

  } else if (query_distribution == "normal") { // normal distribution
    test_queries = get_normal_distribution(N, 100.0, 20.0, 0, static_cast<uint64_t>(pow(2, 50))-1);

  } else if (query_distribution == "exponential") { // exponential distribution
    test_queries = get_exponential_distribution(N, 10.0, 10, static_cast<uint64_t>(pow(2, 50))-1);
  }

  //----------------------------------------
  //QUERYING SuRF
  //----------------------------------------

  vector<uint64_t> rq_ranges({0, 16, 64, 256});
  uint64_t fp;
  uint64_t tn;
  uint64_t tp;
  double all_rate = 0;
  for(int i = 0; i < rq_ranges.size(); i++) {
    fp = 0;
    tn = 0;
    tp = 0;
    for (int j = 0; j < test_queries.size(); j++) {
      uint64_t lower_bound = test_queries[j];
      uint64_t upper_bound = lower_bound + rq_ranges[i];
    //   cout << "lower bound: " << intToString(lower_bound) << " upper bound: " << intToString(upper_bound) << endl;

      if(surf->lookupRange(intToString(lower_bound), true, intToString(upper_bound), true)) {
        if(find_key_in(sorted_v_keys, lower_bound, upper_bound)) {
          tp++;
        } else {
          fp++;
        }
      } else {
        tn++;
      }
    }    
    double rate = static_cast<double>(fp) / (fp + tn);
    all_rate += rate;
  }
  all_rate = all_rate / rq_ranges.size();
  cout << "    The false positive rate overall for mixed range query " << key_distribution << " keys and " << query_distribution << " is " << all_rate << endl;


  //----------------------------------------
  //SuRF WITH kEY K WE QUERY FROM K+(TEST_NUM)
  //----------------------------------------
  uint64_t TEST_NUM = test_num;
  all_rate = 0;
  for(int i = 0; i < rq_ranges.size(); i++) {
    fp = 0;
    tn = 0;
    tp = 0;
    for(int j = 0; j < sorted_v_keys.size(); j++) {
      uint64_t lower_bound = sorted_v_keys[j] + TEST_NUM;
      uint64_t upper_bound = lower_bound + rq_ranges[i];
      if(surf->lookupRange(intToString(lower_bound), true, intToString(upper_bound), true)) {

        if(find_key_in(sorted_v_keys, lower_bound, upper_bound)) {
          tp++;
        } else {
          fp++;
        }

      } else {
        tn++;
      }
    }
    double rate = static_cast<double>(fp) / (fp + tn);
    all_rate += rate;
  }
  all_rate = all_rate / rq_ranges.size();
  cout << "    The false positive rate for K, K+" << test_num << " queries is " << all_rate << endl;
}


int main() {  

  // Testable options
  vector<string> key_dists({"normal", "uniform"});
  vector<string> query_dists({"normal", "uniform", "exponential"});
  vector<string> interface_options({"Start test", "Choose key distribution", "Choose query distribution", "Choose K, K+n", "Choose number of tests", "Exit test"});

  string key_dist = "normal";
  string query_dist = "normal";
  uint64_t test_num = 1;  
  uint64_t N=10000000; 

 
  cout << "Welcome to SuRF test!" << endl;
  while(true) {

    cout << endl << "----------------------------------------------" << endl
      << "Current options" << endl
      << "  Number of tests: " << N << endl
      << "  Key distribution: " << key_dist << endl
      << "  Query distribution: " << query_dist << endl
      << "  Testing for K, K+" << test_num << endl
      << "----------------------------------------------" << endl << endl;

    switch(display_select_vec(interface_options)) {
      case 1: // Start test
        cout << endl;
        test_surf(100.0,  key_dist,query_dist, test_num, N);
        break;

      case 2: // Choose key distribution
        key_dist = key_dists[display_select_vec(key_dists)-1];
        break;

      case 3: // Choose query distribution
        query_dist = query_dists[display_select_vec(query_dists)-1];
        break;

      case 4: // Choose K, K+n
        test_num = until_number_input(1, N);
        break;

      case 5: // Choose N
        N = until_number_input(0, 100000000);
        break;

      case 6: // Exit
        cout << "Goodbye!" << endl;
        return 0;

    }
  }
  
  return 0;
}


