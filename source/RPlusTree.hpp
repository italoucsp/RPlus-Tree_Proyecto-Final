#include <rplus_utils.hpp>

#define COST double

/*TPLT: data type | number of dimensions | max entries per node | fill factor(by default = M / 2)
  Contains: Node, Entry, comparators(ENTRYSINGLEDIM, ENTRYDIST)
  Operations that you are able to do: insertion, range query, k-nearest neighbors query.*/
template<typename T, size_t N, size_t M, size_t ff = M / 2>
class RPlus {
private:
  struct Node;

  struct Entry {
    shared_ptr<HyperRectangle<T, N>> data;
    shared_ptr<Node> child;

    Entry();
    Entry(shared_ptr<Node> &child);
    Entry(shared_ptr<HyperRectangle<T, N>> &data, shared_ptr<Node> &child);
    HyperRectangle<T, N> get_mbr();
    bool is_in_leaf();
    void show_entry(size_t index);
  };

  struct comparator_ENTRYSINGLEDIM {
    size_t axis;
    comparator_ENTRYSINGLEDIM(size_t axis) {
      this->axis = axis;
    }
    bool operator() (Entry &A, Entry &B) {
      return (A.get_mbr().get_boundaries().first[axis] < B.get_mbr().get_boundaries().first[axis]);
    }
  };

  struct ENTRYDIST {
    double distance;//Priority criteria
    Entry entry;//Object
    ENTRYDIST(HyperPoint<T, N> &p, Entry &md_obj){
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
    HyperRectangle<T, N> mbr;
    shared_ptr<Node> parent;
    ENTRY_GROUP entries;
    bool is_leaf();

    Node();
    Entry& operator[](size_t index);
    void add(Entry &IR);
    void add(ENTRY_GROUP &S);
    size_t get_size();
    void print_node();
  private:
    size_t size;
  };

  shared_ptr<Node> root;

  void hidden_insert(shared_ptr<Node> &R, Entry &IR);
  vector<HyperPoint<T, N>> hidden_search(shared_ptr<Node> &R, const HyperRectangle<T, N> &W, vector<HyperPoint<T, N>> &data_found);
  void split_node(shared_ptr<Node> &R);
  tuple<ENTRY_GROUP, ENTRY_GROUP, ENTRY_GROUP> partition(ENTRY_GROUP *S);
  pair<COST, double> sweep(size_t axis, ENTRY_GROUP *S);
  COST calculate_cost(ENTRY_GROUP *test_set);//min coverage
  double MINDIST(HyperPoint<T, N> p, const HyperRectangle<T, N> r);
  void push_node_in_queue(HyperPoint<T, N> refdata, shared_ptr<Node> &current, priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> &q_NN);
  //pack algorithm discarded because it is not really necessary, only for better partition

public:
  RPlus();
  ~RPlus();
  void insert(HyperRectangle<T, N> &packed_data);
  vector<HyperPoint<T, N>> search(const HyperRectangle<T, N> &W);
  vector<HyperPoint<T, N>> kNN_query(HyperPoint<T, N> refdata, size_t k);
  void read_tree();
};

//===============================R-PLUS-TREE-IMPLEMENTATION============================================

template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::RPlus() {
  try {
    ENTRY_GROUP temp;
    if (N < 2 || M < 2 || M > temp.max_size()) {
      throw runtime_error(ERROR_M_N_VALUES);
    }
    else if (ff > M) {
      throw runtime_error(ERROR_FF_VALUE);
    }
    else {
      root = make_shared<Node>();
    }
  }
  catch (const exception &error) {
    ALERT(error.what());
  }
}

template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::~RPlus() {
  root.reset();
}

//INSERTION PUBLIC
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::insert(HyperRectangle<T, N> &packed_data) {
  try {
    ENTRY_GROUP temp;
    if (!packed_data.data_container()) {
      throw runtime_error(ERROR_NO_RECTDATA);
    }
    else {
      shared_ptr<Node> R = root;
      Entry new_entry;
      new_entry.data = make_shared<HyperRectangle<T, N>> (packed_data);
      hidden_insert(R, new_entry);
    }
  }
  catch (const exception &error) {
    ALERT(error.what());
  }
}

//INSERTION PRIVATE
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::hidden_insert(shared_ptr<Node> &R, Entry &IR) {
  if (!R->is_leaf()) {
    for (size_t i(0); i < R->get_size(); ++i) {
      if ((*R)[i].get_mbr().overlaps(*IR.data))
        hidden_insert((*R)[i].child, IR);
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
template<typename T, size_t N, size_t M, size_t ff>
vector<HyperPoint<T, N>> RPlus<T, N, M, ff>::search(const HyperRectangle<T, N> &W) {
  shared_ptr<Node> R = root;
  vector<HyperPoint<T, N>> data_found;
  return hidden_search(R, W, data_found);
}

//RANGE QUERY SEARCH PRIVATE
template<typename T, size_t N, size_t M, size_t ff>
vector<HyperPoint<T, N>> RPlus<T, N, M, ff>::hidden_search(shared_ptr<Node> &R, const HyperRectangle<T, N> &W, vector<HyperPoint<T, N>> &data_found) {
  for (size_t i(0); i < R->get_size(); ++i) {
    if ((*R)[i].get_mbr().overlaps(W)) {
      if (!R->is_leaf()) {
        hidden_search((*R)[i].child, W);
      }
      else {
        data_found.push_back(entry.data);
      }
    }
  }
  return data_found;
}

//KNN QUERY SEARCH PUBLIC AND ONLY THAT, BRANCH AND BOUND ALGORITHM
template<typename T, size_t N, size_t M, size_t ff>
vector<HyperPoint<T, N>> RPlus<T, N, M, ff>::kNN_query(HyperPoint<T, N> refdata, size_t k) {
  priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> best_branchs_queue;
  vector<HyperPoint<T>> kNN(k);
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
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::push_node_in_queue(HyperPoint<T, N> refdata, shared_ptr<Node> &current, priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> &q_NN) {
  for (size_t i = size_t(0); i < current->size; ++i) {
    ENTRYDIST packed_entry(refdata, current->entries[i]);
    q_NN.push(packed_entry);
  }
}

//SPLIT NODE METHOD
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::split_node(shared_ptr<Node> &R) {
  ALERT("__split_init__");
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
    Entry right_node_after_split(N2);
    PR->add(right_node_after_split);
    if (PR->get_size() > M) {
      split_node(PR);
    }
  }
  ALERT("__split_done__");
}

//PARTITION METHOD | 3 GROUPS(ONLY FIRST NODE, ONLY SECOND NODE, ENTRIES THAT ARE IN BOTH NODES)
template<typename T, size_t N, size_t M, size_t ff>
tuple<vector<typename RPlus<T, N, M, ff>::Entry>, vector<typename RPlus<T, N, M, ff>::Entry>,
    vector<typename RPlus<T, N, M, ff>::Entry>> RPlus<T, N, M, ff>::partition(ENTRY_GROUP *S) {
  ALERT("__partition_init__");
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
  ALERT("__partition_done__");
  return make_tuple(s1, s2, s3);
}

//SWEEP DIMENSIONS METHOD IN ONE DIMENSION
template<typename T, size_t N, size_t M, size_t ff>
pair<COST, double> RPlus<T, N, M, ff>::sweep(size_t axis, ENTRY_GROUP *S) {
  ALERT("__sweep_init__");
  comparator_ENTRYSINGLEDIM comparator(axis);
  partial_sort(S->begin(), S->begin() + ff, S->end(), comparator);//sort the first ff entries to "sweep" -> O(n log m)
  ENTRY_GROUP test_set(ff);
  for (size_t id_entry = 0; id_entry < ff; test_set[id_entry] = (*S)[id_entry], ++id_entry) {}//picking the ff firsts entries of the sorted set
  COST total_cost = calculate_cost(&test_set);
  ALERT("__sweep_done__");
  return make_pair(total_cost, test_set[ff - 1].get_mbr().get_boundaries().first[axis]);
}

//MINIMAL COVERAGE COST METHOD
template<typename T, size_t N, size_t M, size_t ff>
COST RPlus<T, N, M, ff>::calculate_cost(ENTRY_GROUP *test_set) {
  ALERT("__ccost_init__");
  COST group_cost = 0.0;
  HyperRectangle<T, N> coverage_mbr;
  for (Entry &entry : *test_set) {
    coverage_mbr.adjust(entry.get_mbr());
    group_cost += entry.get_mbr().get_hypervolume();
  }
  ALERT("__ccost_done__");
  return abs(coverage_mbr.get_hypervolume() - group_cost);
}

//DIST TO NEAREST SIDE OF HYPERECTANGLE
template<typename T, size_t N, size_t M, size_t ff>
double RPlus<T, N, M, ff>::MINDIST(HyperPoint<T, N> p, const HyperRectangle<T, N> r) {
  double sum = 0.0;
  for (size_t i = size_t(0); i < N; ++i) {
    if (p[i] < r.get_boundaries().first[i])
      sum += pow(p[i] - r.get_boundaries().first[i], 2);
    else if (p[i] > r.get_boundaries().second[i])
      sum += pow(p[i] - r.get_boundaries().second[i], 2);
  }
  return sqrt(sum);
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::read_tree() {
  root->print_node();
}

//========================================NODE-IMPLEMENTATION==========================================

template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::Node::Node() {
  entries.resize(M + 1);
  size = size_t(0);
}

template<typename T, size_t N, size_t M, size_t ff>
bool RPlus<T, N, M, ff>::Node::is_leaf() {
  return entries[0].is_in_leaf();
}

template<typename T, size_t N, size_t M, size_t ff>
typename RPlus<T, N, M, ff>::Entry& RPlus<T, N, M, ff>::Node::operator[](size_t index) {
  try {
    if (index >= size) {
      throw runtime_error(ERROR_NODE_OFR);
    }
    else {
      return entries[index];
    }
  }
  catch (const exception &error) {
    ALERT(error.what());
  }
  return entries[0];
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::add(Entry &IR) {
  entries[size++] = IR;
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::add(ENTRY_GROUP &S) {
  for (Entry &entry : S)
    add(entry);
}

template<typename T, size_t N, size_t M, size_t ff>
size_t RPlus<T, N, M, ff>::Node::get_size() {
  return size;
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::print_node() {
  cout << "\tNODE : size(" << size << ") = [" << endl;
  cout << "\t\tID : " << this << endl;
  cout << "\t\tType : " << ((is_leaf())?"LEAF":"INTERNAL") << endl;
  for (size_t i(0); i < size; ++i) {
    (*this)[i].show_entry(i + 1);
  }
  cout << "\t]\n";
}

//=======================================ENTRY-IMPLEMENTATION==========================================

template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::Entry::Entry() {
  //smart pointers did the job
}

template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::Entry::Entry(shared_ptr<Node> &child) {
  this->child = child;
}

template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::Entry::Entry(shared_ptr<HyperRectangle<T, N>> &data, shared_ptr<Node> &child) {
  this->data = data;
  this->child = child;
}

template<typename T, size_t N, size_t M, size_t ff>
HyperRectangle<T, N> RPlus<T, N, M, ff>::Entry::get_mbr() {
  if (child) {
    return child->mbr;
  }
  return *data;
}

template<typename T, size_t N, size_t M, size_t ff>
bool RPlus<T, N, M, ff>::Entry::is_in_leaf() {
  return !child;
}

template<typename T, size_t N, size_t M, size_t ff>

void RPlus<T, N, M, ff>::Entry::show_entry(size_t index) {
  cout << "\t\tEntry<" << index << ">{\n";
  data->show_data();
  cout << "\t\tChild : " << child << endl;
  cout << "\t\t}\n";
}