#include <rplus_utils.hpp>

#define COST double

//DATA, DIMENSIONS, MAX NUMBER OF ENTRIES PER NODE
template<typename DATA_TYPE, size_t N, size_t M, size_t ff = M / 2>
class RPlus {
private:
  struct Node;

  struct Entry {
    shared_ptr<DATA_TYPE> data;
    shared_ptr<Node> child;

    Entry();
    Entry(shared_ptr<Node> &child);
    Entry(shared_ptr<DATA_TYPE> &data, shared_ptr<Node> &child);
    HyperRectangle<N> get_mbr();
    bool is_in_leaf();
  };

  template<size_t axis>
  struct comparator_ENTRYSINGLEDIM {
    bool operator() (const Entry A, const Entry B) {
      return (A.get_mbr().get_boundaries().first[axis] < B.get_mbr().get_boundaries().first[axis]);
    }
  };

  struct ENTRYDIST {
    double distance;//Priority criteria
    Entry entry;//Object
    ENTRYDIST(HyperPoint &p, Entry &md_obj){
      entry = md_obj;
      distance = MINDIST(p, md_obj.get_mbr());//Roussopoulos said that a R tree with minimun overlap don't need MINMAXDIST, and R+ ... don't have overlap so...
    }
  };

  struct comparator_ENTRYDIST {
    bool operator()(const ENTRYDIST &A, const ENTRYDIST &B) {
      return A.distance > B.distance;
    }
  };

  struct Node {
    HyperRectangle<N> mbr;
    shared_ptr<Node> parent;
    ENTRY_GROUP entries;
    bool is_leaf();

    Node();
    void add(Entry &IR);
    void add(ENTRY_GROUP &S);
    size_t get_size();

    size_t size;
  };

  shared_ptr<Node> root;

  void hidden_insert(shared_ptr<Node> &R, Entry &IR);
  vector<DATA_TYPE> hidden_search(shared_ptr<Node> &R, const HyperRectangle<N> &W);
  void split_node(shared_ptr<Node> &R);
  tuple<ENTRY_GROUP, ENTRY_GROUP, ENTRY_GROUP> partition(ENTRY_GROUP *S);
  pair<COST, double> sweep(size_t axis, ENTRY_GROUP *S);
  COST calculate_cost(ENTRY_GROUP *test_set);//min coverage
  double MINDIST(DATA_TYPE p, const HyperRectangle<N> r);//necessary
  double MINMAXDIST(DATA_TYPE p, const HyperRectangle<N> r);//not necessary but could be if overlap would exist
  void push_node_in_queue(HyperPoint refdata, shared_ptr<Node> &current, priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> &q_NN);
  //pack algorithm discarded because it is not really necessary, only for better partition

public:
  RPlus();
  ~RPlus();
  void insert(const DATA_TYPE data);
  vector<DATA_TYPE> search(const HyperRectangle<N> &W);
  vector<DATA_TYPE> kNN_query(HyperPoint refdata, size_t k);
};

//===============================R-PLUS-TREE-IMPLEMENTATION============================================

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
RPlus<DATA_TYPE, N, M, ff>::RPlus() {
  try {
    ENTRY_GROUP temp;
    if (N < 2 || M < 2 || M > temp.max_size()) {
      throw runtime_error(ERROR_M_N_VALUES);
    }
    else {
      root = make_shared<Node>();
    }
  }
  catch (const exception &error) {
    ALERT(error.what());
  }
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
RPlus<DATA_TYPE, N, M, ff>::~RPlus() {
  root.reset();
}

//INSERTION PUBLIC
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::insert(const DATA_TYPE data) {
  shared_ptr<Node> R = root;
  Entry new_entry;
  new_entry.data = data;
  hidden_insert(R, new_entry);
}

//INSERTION PRIVATE
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::hidden_insert(shared_ptr<Node> &R, Entry &IR) {
  if (!R->is_leaf()) {
    for (Entry &entry_ : R->entries) {
      if (entry_.get_mbr().contains(IR.data)) {
        hidden_insert(entry_.child, IR);
      }
    }
  }
  else {
    R->add(IR);
    if (R->get_size() > M) {
      split_node(R);
    }
  }
}

//RANGE QUERY SEARCH PUBLIC
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
vector<DATA_TYPE> RPlus<DATA_TYPE, N, M, ff>::search(const HyperRectangle<N> &W) {
  shared_ptr<Node> R = root;
  return hidden_search(R, W);
}

//RANGE QUERY SEARCH PRIVATE
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
vector<DATA_TYPE> RPlus<DATA_TYPE, N, M, ff>::hidden_search(shared_ptr<Node> &R, const HyperRectangle<N> &W) {
  vector<DATA_TYPE> data_found;
  if (!R->is_leaf()) {
    for (Entry &entry : R->entries) {
      if (entry.get_mbr().isOverlaping(W)) {
        hidden_search(entry.child, W);
      }
    }
  }
  else {
    for (Entry &entry : R->entries) {
      if (entry.entry_rect.isOverlaping(W)) {
        data_found.push_back(entry.data);
      }
    }
  }
  return data_found;
}

//KNN QUERY SEARCH PUBLIC AND ONLY THAT, BRANCH AND BOUND ALGORITHM
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
vector<DATA_TYPE> RPlus<DATA_TYPE, N, M, ff>::kNN_query(HyperPoint refdata, size_t k) {
  priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> best_branchs_queue;
  vector<DATA_TYPE> kNN(k);
  push_node_in_queue(refdata, root, best_branchs_queue);
  size_t i = size_t(0);
  while (i < k && !best_branchs_queue.empty()) {
    ENTRYDIST closest_entry = best_branchs_queue.top();
    if (!closest_entry.entry.is_in_leaf())
      push_node_in_queue(refdata, closest_entry.entry.child, best_branchs_queue);
    else {
      kNN[i++] = *(closest_entry.entry.data);
      best_branchs_queue.pop();
    }
  }
  return kNN;
}

//PUSH EACH ENTRY OF A NODE IN THE PRIORITY QUEUE
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::push_node_in_queue(HyperPoint refdata, shared_ptr<Node> &current, priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> &q_NN) {
  for (size_t i = size_t(0); i < current->size; ++i) {
    ENTRYDIST packed_entry(refdata, current->entries[i]);
    q_NN.push(packed_entry);
  }
}

//SPLIT NODE METHOD
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::split_node(shared_ptr<Node> &R) {
  tuple<ENTRY_GROUP, ENTRY_GROUP, ENTRY_GROUP> partitioned_group = partition(&R->entries);
  shared_ptr<Node> N2 = make_shared<Node>();
                   N2->add(get<1>(partitioned_group));
                   R = make_shared<Node>();
                   R->add(get<0>(partitioned_group));
  for (Entry &entry : get<2>(partitioned_group)) {
    if (R->is_leaf()) {
      R->add(entry);
      N2->add(entry);
    }
    else {
      split_node(entry.child);
      //add the two newnodes to entry.child
    }
  }
  if (R == root) {
    shared_ptr<Node> new_root = make_shared<Node>();
    Entry en1(R), en2(N2);
    new_root->add(en1); new_root->add(en2);
    root = new_root;
  }
  else {
    shared_ptr<Node> PR = R->parent;
    //replace R with R' and N2
    if (PR->get_size() > M) {
      split_node(PR);
    }
  }
}

//PARTITION METHOD | 3 GROUPS(ONLY FIRST NODE, ONLY SECOND NODE, ENTRIES THAT ARE IN BOTH NODES)
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
tuple<vector<typename RPlus<DATA_TYPE, N, M, ff>::Entry>, vector<typename RPlus<DATA_TYPE, N, M, ff>::Entry>,
    vector<typename RPlus<DATA_TYPE, N, M, ff>::Entry>> RPlus<DATA_TYPE, N, M, ff>::partition(ENTRY_GROUP *S) {
  if (S->size() <= ff) {
    ENTRY_GROUP empty_group;
    return make_tuple(*S, empty_group, empty_group);
  }
  COST cheapest_cost = DBL_MAX;
  size_t optimal_dim = size_t(0);
  double optimal_cutline = 0.0;

  for (size_t current_dim = size_t(0); current_dim < N; ++current_dim) {
    pair<COST, double> cost_and_cutline = sweep(current_dim, S);//sweep in any dimension
    COST temp_min_cost = cheapest_cost;
    cheapest_cost = min(cost_and_cutline.first, cheapest_cost);

    if (cheapest_cost != temp_min_cost) {
      optimal_cutline = cost_and_cutline.second;
      optimal_dim = current_dim;
    }
  }
  ENTRY_GROUP s1, s2, s3;
  for (Entry &entry : *S) {
    if (entry.get_mbr().get_boundaries().second[optimal_dim] <= optimal_cutline)//only in the first group
      s1.push_back(entry);
    else if (entry.get_mbr().get_boundaries().first[optimal_dim] > optimal_cutline)//only in the second group
      s2.push_back(entry);
    else //entry match in the two groups
      s3.push_back(entry);
  }
  return make_tuple(s1, s2, s3);
}

//SWEEP DIMENSIONS METHOD IN ONE DIMENSION
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
pair<COST, double> RPlus<DATA_TYPE, N, M, ff>::sweep(size_t axis, ENTRY_GROUP *S) {
  comparator_ENTRYSINGLEDIM<axis> comparator;
  partial_sort(S->begin(), S->begin() + ff, S->end(), comparator);//sort the first ff entries to "sweep" -> O(n log m)
  ENTRY_GROUP test_set(ff);
  for (size_t id_entry = 0; id_entry < ff; test_set[id_entry] = (*S)[id_entry++]) {}//picking the ff firsts entries of the sorted set
  COST total_cost = calculate_cost(&test_set);
  return make_pair(total_cost, test_set[ff - 1].get_mbr().get_boundaries().first[axis]);
}

//MINIMAL COVERAGE COST METHOD
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
COST RPlus<DATA_TYPE, N, M, ff>::calculate_cost(ENTRY_GROUP *test_set) {
  COST group_cost = 0.0;
  HyperRectangle<N> coverage_mbr;
  for (Entry &entry : *test_set) {
    coverage_mbr.adjust_with_hrect(entry.get_mbr());
    group_cost += entry.get_mbr().get_hypervolume();
  }
  return abs(coverage_mbr.get_hypervolume() - group_cost);
}

//DIST TO NEAREST SIDE OF HYPERETANGLE
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
double RPlus<DATA_TYPE, N, M, ff>::MINDIST(DATA_TYPE p, const HyperRectangle<N> r) {
  double sum = 0.0;
  for (size_t i = size_t(0); i < N; ++i) {
    if (p[i] < r.get_boundaries().first[i])
      sum += pow(p[i] - r.get_boundaries().first[i], 2);
    else if (p[i] > r.get_boundaries().second[i])
      sum += pow(p[i] - r.get_boundaries().second[i], 2);
  }
  return sqrt(sum);
}

//DIST TO FURTHEST CORNER - NEAREST SIDE OF HYPERETANGLE
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
double RPlus<DATA_TYPE, N, M, ff>::MINMAXDIST(DATA_TYPE p, const HyperRectangle<N> r) {
  double minmax = DBL_MAX;
  for (size_t i = size_t(0); i < N; ++i) {
    double RMsum = 0.0;
    for (size_t j = size_t(0); j < N; ++j) {
      if (j != i) {
        if (p[j] >= (r.get_boundaries().first[j] + r.get_boundaries().second[j]) / 2)
          RMsum += pow(p[j] - r.get_boundaries().first[j], 2);
        else
          RMsum += pow(p[j] - r.get_boundaries().second[j], 2);
      }
    }
    if (p[i] <= (r.get_boundaries().first[i] + r.get_boundaries().second[i]) / 2)
      minmax = min(pow(p[i] - r.get_boundaries().first[i], 2) + RMsum, minmax);
    else
      minmax = min(pow(p[i] - r.get_boundaries().second[i], 2) + RMsum, minmax);
  }
  return sqrt(minmax);
}

//========================================NODE-IMPLEMENTATION==========================================

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
RPlus<DATA_TYPE, N, M, ff>::Node::Node() {
  entries.resize(M + 2);
  size = size_t(0);
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
bool RPlus<DATA_TYPE, N, M, ff>::Node::is_leaf() {
  return entries[0].is_in_leaf();
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::Node::add(Entry &IR) {
  entries[size++] = IR;
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::Node::add(ENTRY_GROUP &S) {
  for (Entry &entry : S) {
    add(entry);
  }
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
size_t RPlus<DATA_TYPE, N, M, ff>::Node::get_size() {
  return size;
}

//=======================================ENTRY-IMPLEMENTATION==========================================

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
RPlus<DATA_TYPE, N, M, ff>::Entry::Entry() {
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
RPlus<DATA_TYPE, N, M, ff>::Entry::Entry(shared_ptr<Node> &child) {
  this->child = child;
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
RPlus<DATA_TYPE, N, M, ff>::Entry::Entry(shared_ptr<DATA_TYPE> &data, shared_ptr<Node> &child) {
  this->data = data;
  this->child = child;
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
HyperRectangle<N> RPlus<DATA_TYPE, N, M, ff>::Entry::get_mbr() {
  if (child) {
    return child->mbr;
  }
  return HyperRectangle(data, data);
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
bool RPlus<DATA_TYPE, N, M, ff>::Entry::is_in_leaf() {
  return !child;
}