#include <limits>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <csignal>

#include <omp.h>
#include <stdio.h>
#define DEFAULT_STRIPES_IN_RESULTS 5
#define DEFAULT_PRINT_STRIPES 20

using namespace std;
using namespace boost::random;

const int G = -1;
const int E = -2;

bool stop_flag = false;

void sigint_handler(int) { stop_flag = true; }

unsigned long long __rdtsc()
{
	unsigned long lo, hi;
		/**
		 * __asm__ __volatile__ (      // serialize
		 * "xorl %%eax,%%eax \n        cpuid"
		 * ::: "%rax", "%rbx", "%rcx", "%rdx");
		**/
		/* We cannot use "=A", since this would use %rax on x86_64 */
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return (unsigned long long)hi << 32 | lo;
}

class Stripe {
 public:
  Stripe(int groups_count, int group_len, int disks_count, int seed)
    : generator_(seed) {
    v_.push_back(E);
    v_.push_back(G);
    for (int group = 1; group <= groups_count; group++) {
      for (int i = 0; i < group_len; i++) {
        v_.push_back(group);
      }
    }

    realSize_ = v_.size();
    distribution_ = uniform_int_distribution<>(0, realSize_ - 1);
    v_.resize(disks_count);
  }

  void next() {
    for (size_t i = 0; i < realSize_; i++) {
      swap(v_[i], v_[distribution_(generator_)]);
    }
  }

  int operator[](size_t index) const { return v_[index]; }

  vector<int> v() const { return v_; }

 private:
  size_t realSize_;
  vector<int> v_;
  mt19937_64 generator_;
  uniform_int_distribution<> distribution_;
};

struct Result {
  int maxMinDiff;
  vector<vector<int> > sums;
  vector<vector<int> > stripes;
};

unsigned long long calc_possible_stripes(unsigned long long stripe_len,
                                         unsigned long long groups_count,
                                         unsigned long long group_len) {
  unsigned long long result = 1;
  unsigned long long group_len_fac = 1;
  for (unsigned long long i = 2; i <= group_len; i++) {
    group_len_fac *= i;
  }
  for (unsigned long long i = 2; i <= stripe_len; i++) {
    result *= i;
    if (groups_count > 0 && result % group_len_fac == 0) {
      groups_count--;
      result /= group_len_fac;
    }
  }
  return result;
}

template <int disks_count, int groups_count, int group_len>
Result bruteforce(size_t stripes_in_result) {
  Result result;
  const size_t stripe_len = groups_count * group_len + 2;
  result.maxMinDiff = disks_count;


  const size_t total_items = disks_count * stripe_len;
  vector<int> sum(disks_count);

  unsigned long long ts;
  ts = __rdtsc();
  Stripe stripe(groups_count, group_len, disks_count, ts);
  unsigned long long stripe_counter = 0;

  do {
    if (stripe_counter % (50 * 1000 * 1000 / stripe_len / disks_count) == 0) {
      cerr << '\r' << "Diff: " << setw(3) << result.maxMinDiff
           << "  Good stripes: " << setw(15) << result.stripes.size()
           << "  Stripes: " << setw(15) << stripe_counter;
    }

    size_t pos = 0;
    fill(sum.begin(), sum.end(), 0);
    while (pos < total_items) {
      size_t failed_index =
          disks_count - ((pos - 1 + disks_count) % disks_count + 1);
      int failed = stripe[failed_index];
      for (size_t i = 0; i < stripe_len; i++) {
        int source = stripe[i];
        sum[(pos + i) % disks_count] +=
            (i != failed_index) &&
            ((failed > 0 && source == failed) || (failed == G && source != E));
      }
      pos += stripe_len;
    }

    int new_min = sum[1];
    int new_max = sum[1];
    for (size_t i = 1; i < sum.size(); i++) {
      new_min = min(new_min, sum[i]);
      new_max = max(new_max, sum[i]);
    }
    int new_max_min_diff = new_max - new_min;
    if (new_max_min_diff < result.maxMinDiff) {
      result.maxMinDiff = new_max_min_diff;
      result.sums.clear();
      result.stripes.clear();
      result.sums.push_back(sum);
      result.stripes.push_back(stripe.v());
    } else if (result.maxMinDiff == new_max_min_diff) {
      result.sums.push_back(sum);
      result.stripes.push_back(stripe.v());
    }

    stripe_counter++;
    stripe.next();
  } while (!(result.maxMinDiff <= 2 && result.stripes.size() == stripes_in_result) &&
           !stop_flag);

  cerr << '\r' << "Diff: " << setw(3) << result.maxMinDiff
       << "  Good stripes: " << setw(15) << result.stripes.size()
       << "  Stripes: " << setw(15) << stripe_counter << endl;

  return result;
}

int CheckStripeLen(size_t disks_count, size_t groups_count, size_t group_len)
{
	const size_t stripe_len = groups_count * group_len + 2;
	if (stripe_len > disks_count)
	{
		cerr << "Error: stripe is shorter than disks count" << endl;
	    exit(1);
	}
	unsigned long long possible_stripes =
	      calc_possible_stripes(stripe_len, groups_count, group_len);

	cerr << "Possible stripes: " << possible_stripes << endl;
	return 0;
}
string stripe2string(vector<int> stripe) {
  ostringstream oss;
  for (size_t i = 0; i < stripe.size() && stripe[i] != 0; i++) {
    oss << setw(3);
    switch (stripe[i]) {
      case G:
        oss << 'G';
        break;
      case E:
        oss << 'E';
        break;
      default:
        oss << stripe[i];
    }
  }
  return oss.str();
}

int printResults(vector<Result> results, int list_len)
{
	int i;
	Result result;
	for(i = 0; i < list_len; i++)
	{
    result = results[i];
	  cout << result.sums.size() << " with diff=" << result.maxMinDiff << endl;
	  int sums_counter = 0;
	  for (size_t i = 0; i < result.sums.size(); i++) {
	    for (size_t j = 0; j < result.sums[i].size(); j++) {
	      cout << result.sums[i][j] << ' ';
	    }
	    cout << '\t' << stripe2string(result.stripes[i]);
	    cout << endl;

	    sums_counter++;
	    if (sums_counter > DEFAULT_PRINT_STRIPES) {
	      break;
	    }
	  }
	}
	return 0;
}

#define disks_count 30
#define groups_count 3
#define group_len 9

int main(int argc, char *argv[])
{
  signal(SIGINT, sigint_handler);
  int thread_num, i;
  size_t stripes_in_result;

  CheckStripeLen(disks_count, groups_count, group_len);

  if(argc == 3)
  {
	  thread_num = atoi(argv[1]);
	  stripes_in_result = atoi(argv[2]);
  }
  else
  {
	  thread_num = omp_get_max_threads();
	  stripes_in_result = DEFAULT_STRIPES_IN_RESULTS;
  }
  printf("thread_num = %d, stripes_in_result = %lu\n", thread_num, stripes_in_result);
  omp_set_num_threads(thread_num);

  vector<Result> results(thread_num);

#pragma omp parallel
  {
	  Result result;
#pragma omp for
	  for(i = 0; i < thread_num; i++)
	  {
      results[i] = bruteforce<disks_count, groups_count, group_len>(stripes_in_result);
	  }

  }

  printResults(results, thread_num);

  return 0;
}
