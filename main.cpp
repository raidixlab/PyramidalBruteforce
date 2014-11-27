#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <sstream>

using namespace std;

const int G = -1;
const int E = -2;

class Stripe {
 public:
  struct Node {
    Node(int v, Node* n) : value(v), next(n) {}

    void fillVector(vector<int>& result) {
      Node* node = this;
      int i = 0;
      while (node != NULL) {
        result[i] = node->value;
        node = node->next;
        i++;
      }
    }

    int value;
    Node* next;
  };

  Stripe(int groups_count, int group_len) : justConstructed_(true) {
    afteri_ = new Node(E, NULL);
    head_ = i_ = new Node(G, afteri_);
    for (int group = 1; group <= groups_count; group++) {
      for (int i = 0; i < group_len; i++) {
        head_ = new Node(group, head_);
      }
    }
  }

  bool next(vector<int>& result) {
    if (justConstructed_) {
      justConstructed_ = false;
      head_->fillVector(result);
      return true;
    } else if (afteri_->next != NULL || afteri_->value < head_->value) {
      Node* beforek;
      Node* k;
      if (afteri_->next != NULL && i_->value >= afteri_->next->value) {
        beforek = afteri_;
      } else {
        beforek = i_;
      }
      k = beforek->next;
      beforek->next = k->next;
      k->next = head_;
      if (k->value < head_->value) {
        i_ = k;
      }
      afteri_ = i_->next;
      head_ = k;
      head_->fillVector(result);
      return true;
    } else {
      return false;
    }
  }

 private:
  bool justConstructed_;
  Node* head_;
  Node* i_;
  Node* afteri_;
};

struct Result {
  int max;
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
Result bruteforce() {
  Result result;
  const size_t stripe_len = groups_count * group_len + 2;
  result.max = stripe_len;

  const size_t total_items = disks_count * stripe_len;
  vector<int> table(total_items);
  vector<int> sum(disks_count);

  Stripe stripe(groups_count, group_len);
  vector<int> stripe_v(disks_count);
  unsigned long long stripe_counter = 0;

  unsigned long long possible_stripes =
      calc_possible_stripes(stripe_len, groups_count, group_len);

  while (stripe.next(stripe_v)) {
    if (stripe_counter % (1 * 1000 * 1000) == 0) {
      cerr << '\r' << 100.0 * stripe_counter / possible_stripes << "%\t\t";
    }
    size_t pos = 0;
    fill(sum.begin(), sum.end(), 0);
    while (pos < total_items) {
      size_t failed_index =
          disks_count - ((pos - 1 + disks_count) % disks_count + 1);
      int failed = stripe_v[failed_index];
      for (size_t i = 0; i < stripe_len; i++) {
        int source = stripe_v[i];
        sum[(pos + i) % disks_count] +=
            (i != failed_index) &&
            ((failed > 0 && source == failed) || (failed == G && source != E));
      }
      pos += stripe_len;
    }

    int new_max = *max_element(sum.begin(), sum.end());
    if (result.max > new_max) {
      result.max = new_max;
      result.sums.clear();
      result.stripes.clear();
      result.sums.push_back(sum);
      result.stripes.push_back(stripe_v);
    } else if (result.max == new_max) {
      result.sums.push_back(sum);
      result.stripes.push_back(stripe_v);
    }

    stripe_counter++;
  }

  cerr << "\r100%\t\t" << endl;

  return result;
}

string stripe2string(vector<int> stripe) {
  ostringstream oss;
  for (size_t i = 0; i < stripe.size() && stripe[i] != 0; i++) {
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

int main() {
  Result result = bruteforce<15, 3, 4>();
  cout << result.sums.size() << " with max=" << result.max << endl;
  int sums_counter = 0;
  for (size_t i = 0; i < result.sums.size(); i++) {
    for (size_t j = 0; j < result.sums[i].size(); j++) {
      cout << result.sums[i][j] << ' ';
    }
    cout << '\t' << stripe2string(result.stripes[i]);
    cout << endl;

    sums_counter++;
    if (sums_counter == 100) {
      break;
    }
  }
}
