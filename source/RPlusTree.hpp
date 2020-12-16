#include <rplus_utils.hpp>

#define GET_BOUNDARIES(entry) entry.get_mbr().get_boundaries()

/*TEMPLATE PARAMETERS: (1)data type | (2)number of dimensions | (3)max entries per node | (4)fill factor(by default = 2)
  Contains: Node, Entry, comparators(ENTRYSINGLEDIM, ENTRYDIST)
  Approach: Packed Static - Point R+ Tree
  Operations that you are able to do: insertion(insert,"massive"), range query(search), k-nearest neighbors query(kNN_query).*/
template<typename T, size_t N, size_t M, size_t ff = 2>
class RPlus {
private:
  struct Node;

  struct Entry {
    HyperPoint<T, N> data;
    shared_ptr<Node> child;

    Entry();
    Entry(shared_ptr<Node> &child);
    Entry(HyperPoint<T, N> &data);
    HyperRectangle<T, N> get_mbr();
    Entry& operator=(const Entry &other);
    bool is_in_leaf();
    void show_entry(size_t index);
  };

  struct comparator_ENTRYSINGLEDIM {
    size_t axis;
    comparator_ENTRYSINGLEDIM(size_t axis) {
      this->axis = axis;
    }
    bool operator() (Entry &A, Entry &B) {
      return (GET_BOUNDARIES(A).first[axis] < GET_BOUNDARIES(B).first[axis]);
    }
  };

  struct ENTRYDIST {
    double distance;//Priority criteria
    Entry entry;//Object
    ENTRYDIST(HyperPoint<T, N> &p, Entry &md_obj){
      entry = md_obj;
      if(md_obj.get_mbr().get_hypervolume() > 0)
        distance = RPlus::MINDIST(p, md_obj.get_mbr());//Roussopoulos said that a R tree with minimun overlap don't need MINMAXDIST, and R+ ... don't have overlap so...
      else
        distance = RPlus::EUCDIST(p, md_obj.data);
    }
  };

  struct comparator_ENTRYDIST {
    bool operator()(const ENTRYDIST &A, const ENTRYDIST &B) {
      return A.distance > B.distance;
    }
  };

  struct Node {
    HyperRectangle<T, N> mbr;
    vector<Entry> entries;
    bool is_leaf();

    Node();
    Entry& operator[](size_t index);
    void add(Entry &new_entry);
    void add(vector<Entry> &S);
    size_t get_size();
    void resize(size_t new_size);
    void print_node(bool rp_root = false);
  private:
    size_t size;
  };

  shared_ptr<Node> root;

  shared_ptr<Node> split_node(shared_ptr<Node> &A, size_t axis, T optimal_cutline);
  vector<Entry> partition(vector<Entry> &S, vector<Entry> &s_to_split, size_t &optimal_dim, T &optimal_cutline);
  void pack(vector<Entry> &S);
  inline pair<double, T> sweep(size_t axis, vector<Entry> &S, vector<Entry> &grouped_test, vector<Entry> &rest);
  inline int min_number_splits(vector<Entry> &test_set, size_t axis, T optimal_cutline);//min splits
  static double MINDIST(HyperPoint<T, N> p, HyperRectangle<T, N> r);
  static double EUCDIST(HyperPoint<T, N> p1, HyperPoint<T, N> p2);
  inline void push_node_in_queue(HyperPoint<T, N> refdata, shared_ptr<Node> &current, priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> &q_NN);

public:
  RPlus();
  virtual ~RPlus();
  void insert(vector<HyperPoint<T, N>> &unpacked_data);
  vector<HyperPoint<T, N>> search(const HyperRectangle<T, N> &W);
  vector<HyperPoint<T, N>> kNN_query(HyperPoint<T, N> refdata, size_t k);
  void read_tree();
};

//===============================R-PLUS-TREE-IMPLEMENTATION============================================

//BUILDER RPLUS
template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::RPlus() {
  try {
    vector<Entry> temp;
    if (N < 2 || M < 2 || M > temp.max_size()) {
      throw runtime_error(ERROR_M_N_VALUES);
    }
    else if (ff > M || ff < 2) {
      throw runtime_error(ERROR_FF_VALUE);
    }
    else {
      root = make_shared<Node>();
    }
  }
  catch (const exception &error) {
    ALERT(error.what());
    exit(0);
  }
}

//ERASER RPLUS
template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::~RPlus() {
  root.reset();
}

//INSERTION METHOD
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::insert(vector<HyperPoint<T, N>> &unpacked_data) {
  vector<Entry> unpacked_entries(unpacked_data.size());
  size_t i(0);
  for (HyperPoint<T, N> &udata : unpacked_data) {
    unpacked_entries[i++] = Entry(udata);
  }
  pack(unpacked_entries);
}

//RANGE QUERY SEARCH PUBLIC
template<typename T, size_t N, size_t M, size_t ff>
vector<HyperPoint<T, N>> RPlus<T, N, M, ff>::search(const HyperRectangle<T, N> &W) {
  try {
    if (!root) {
      throw runtime_error(ERROR_EMPTY_TREE);
    }
    else {     
      vector<HyperPoint<T, N>> range_query;
      unordered_set<string> songs_names;
      stack<shared_ptr<Node>> dfs_s;
      dfs_s.push(root);
      while (!dfs_s.empty()) {
        shared_ptr<Node> current = dfs_s.top();
        dfs_s.pop();
        for (size_t i(0); i < current->get_size(); ++i) {
          if ((*current)[i].get_mbr().overlaps(W)) {
            if (!current->is_leaf())
              dfs_s.push((*current)[i].child);
            else {
              size_t temp_songs_names_size = songs_names.size();
              //try pick data
              songs_names.insert((*current)[i].data.get_songs_name());
              if (temp_songs_names_size != songs_names.size())
                range_query.push_back((*current)[i].data);
              //end try
            }
          }
        }
      }
      return range_query;
    }
  }
  catch (const exception &error) {
    ALERT(error.what())
    exit(0);
  }
}

/*KNN METHOD: k-Nearest Neighbors query using branch and bound algorithm with MINDIST function
  PAPER: A. Papadopoulos, Y. Manolopoulos, "Performance of Nearest Neighbor Queries in R-trees *",
         Department of Informatics Aristotle University - 54006 Thessaloniki , Greece */
template<typename T, size_t N, size_t M, size_t ff>
vector<HyperPoint<T, N>> RPlus<T, N, M, ff>::kNN_query(HyperPoint<T, N> refdata, size_t k) {
  try {
    if (!root) {
      throw runtime_error(ERROR_EMPTY_TREE);
    }
    else {
      priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> best_branchs_queue;
      vector<HyperPoint<T, N>> kNN(k);
      unordered_set<string> songs_names;
      push_node_in_queue(refdata, root, best_branchs_queue);
      size_t i = size_t(0);
      while (i < k && !best_branchs_queue.empty()) {
        ENTRYDIST closest_entry = best_branchs_queue.top();
        best_branchs_queue.pop();
        if (!closest_entry.entry.is_in_leaf())
          push_node_in_queue(refdata, closest_entry.entry.child, best_branchs_queue);
        else {
          size_t temp_songs_names_size = songs_names.size();
          //try pick k(i) nn
          songs_names.insert(closest_entry.entry.data.get_songs_name());
          if(temp_songs_names_size != songs_names.size())
            kNN[i++] = closest_entry.entry.data;
          //end try
          best_branchs_queue.pop();
        }
      }
      return kNN;
    }
  }
  catch (const exception &error) {
    ALERT(error.what())
    exit(0);
  }
}

//PUSH EACH ENTRY OF A NODE IN THE PRIORITY QUEUE
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::push_node_in_queue(HyperPoint<T, N> refdata, shared_ptr<Node> &current, priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> &q_NN) {
  for (size_t i = size_t(0); i < current->get_size(); ++i) {
    ENTRYDIST packed_entry(refdata, (*current)[i]);
    q_NN.push(packed_entry);
  }
}

/*SPLIT NODE METHOD: Division of a node in given axis and optimal cutline,
                     then do downward propagation of the split.*/
template<typename T, size_t N, size_t M, size_t ff>
shared_ptr<typename RPlus<T, N, M, ff>::Node> RPlus<T, N, M, ff>::split_node(shared_ptr<Node> &A, size_t axis, T optimal_cutline) {
  SAY("SPLIT PROCESS")
  shared_ptr<Node> B = make_shared<Node>();
  vector<Entry> set_A, set_B;
  for (size_t i(0); i < A->get_size(); ++i) {
    if (GET_BOUNDARIES((*A)[i]).second[axis] < optimal_cutline)
      set_A.push_back((*A)[i]);
    else if (GET_BOUNDARIES((*A)[i]).first[axis] > optimal_cutline)
      set_B.push_back((*A)[i]);
    else {
      if (!A->is_leaf()) {
        set_B.emplace_back(split_node((*A)[i].child, axis, optimal_cutline));
        set_A.emplace_back((*A)[i].child);
      }//if A were a leaf, would be impossible to split a point
    }
  }
  A->resize(0); A->add(set_A);
  B->resize(0); B->add(set_B);
  return B;
}

//PARTITION METHOD
template<typename T, size_t N, size_t M, size_t ff>
vector<typename RPlus<T, N, M, ff>::Entry> RPlus<T, N, M, ff>::partition(vector<Entry> &S, vector<Entry> &s_to_split, size_t &optimal_dim, T &optimal_cutline) {
  if (S.size() <= ff) {
    vector<Entry> nn_set;
    nn_set.swap(S);
    return nn_set;
  }
  double cheapest_cost = DBL_MAX;
  optimal_dim = size_t(0);
  optimal_cutline = 0;
  //sweep and find the best partition cutline/axis
  vector<Entry> test_group, test_remainder, possible_remainder_entries, only_new_group, remainder_entries;
  for (size_t current_dim(0); current_dim < N; ++current_dim) {
    pair<double, T> cost_and_cutline = sweep(current_dim, S, test_group, test_remainder);
    double temp_min_cost = cheapest_cost;
    cheapest_cost = min(cost_and_cutline.first, cheapest_cost);
    if (cheapest_cost != temp_min_cost) {
      optimal_cutline = cost_and_cutline.second;
      optimal_dim = current_dim;
      if (!only_new_group.empty()) only_new_group.clear();
      only_new_group.assign(test_group.begin(), test_group.end());
      test_group.clear();
      if (!possible_remainder_entries.empty()) possible_remainder_entries.clear();
      possible_remainder_entries.assign(test_remainder.begin(), test_remainder.end());
      test_remainder.clear();
    }
  }
  for (Entry entry : possible_remainder_entries) {
    if (GET_BOUNDARIES(entry).first[optimal_dim] > optimal_cutline)
      remainder_entries.push_back(entry);
    else
      s_to_split.push_back(entry);
  }
  S.clear(); S.swap(remainder_entries);
  return only_new_group;
}

//PACK METHOD
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::pack(vector<Entry> &S) {
  if (S.size() <= ff) {
    root->add(S);
    return;
  }
  vector<Entry> next_level_upward, s_to_split;
  size_t axis;
  T cutline;
  while (!S.empty()) {
    vector<Entry> current_tff_set(partition(S, s_to_split, axis, cutline));//new set S, given cutline and axis
    if (!s_to_split.empty()) {
      for (Entry entry : s_to_split) {
        if (!entry.is_in_leaf()) {
          S.emplace_back(split_node(entry.child, axis, cutline));
          current_tff_set.emplace_back(entry.child);
        }
        else 
          S.emplace_back(entry.data);//to have at least ff entries in the new node
      }
      s_to_split.clear();
    }
    shared_ptr<Node> new_node = make_shared<Node>();
    new_node->add(current_tff_set);
    next_level_upward.emplace_back(new_node);
  }
  pack(next_level_upward);
}

//SWEEP DIMENSIONS METHOD IN ONE DIMENSION
template<typename T, size_t N, size_t M, size_t ff>
pair<double, T> RPlus<T, N, M, ff>::sweep(size_t axis, vector<Entry> &S, vector<Entry> &grouped_test, vector<Entry> &rest) {
  comparator_ENTRYSINGLEDIM comparator(axis);
  partial_sort(S.begin(), S.begin() + ff, S.end(), comparator);//sort the first ff entries to "sweep" -> O(n log ff)
  grouped_test.assign(S.begin(), S.begin() + ff);//picking the ff first entries of the sorted set
  rest.assign(S.begin() + ff, S.end());//picking the rest entries for the S set if the group test is selected as the best in partition
  T optimal_cutline = GET_BOUNDARIES(grouped_test.back()).second[axis];
  return make_pair(min_number_splits(grouped_test, axis, optimal_cutline), optimal_cutline);
}

//MINIMAL NUMBER OF SPLITS METHOD
template<typename T, size_t N, size_t M, size_t ff>
int RPlus<T, N, M, ff>::min_number_splits(vector<Entry> &test_set, size_t axis, T optimal_cutline) {
  int cost = 0;
  for (Entry &entry : test_set) {
    if (GET_BOUNDARIES(entry).first[axis] < optimal_cutline && GET_BOUNDARIES(entry).second[axis] > optimal_cutline)
      ++cost;
  }
  return cost;
}

//DIST TO NEAREST SIDE OF AN HYPERECTANGLE
template<typename T, size_t N, size_t M, size_t ff>
double RPlus<T, N, M, ff>::MINDIST(HyperPoint<T, N> p, HyperRectangle<T, N> r) {
  double sum = 0.0;
  for (size_t i = size_t(0); i < N; ++i) {
    if (p[i] < r.get_boundaries().first[i])
      sum += pow(p[i] - r.get_boundaries().first[i], 2);
    else if (p[i] > r.get_boundaries().second[i])
      sum += pow(p[i] - r.get_boundaries().second[i], 2);
  }
  return sqrt(sum);
}

//DIST BETWEEN TWO POINTS
template<typename T, size_t N, size_t M, size_t ff>
double RPlus<T, N, M, ff>::EUCDIST(HyperPoint<T, N> p1, HyperPoint<T, N> p2) {
  double d = 0.0;
  for (size_t i(0); i < N; ++i) {
    d += pow(p1[i] - p2[i], 2);
  }
  return sqrt(d);
}

//READ R+ TREE METHOD: Using bfs to read the levels of the tree since the root.
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::read_tree() {
  if (root) {
    root->print_node(true);
    queue<shared_ptr<Node>> bfs_q;
    for (size_t i(0); i < root->get_size(); ++i) {
      if ((*root)[i].child)
        bfs_q.push((*root)[i].child);
    }
    while (!bfs_q.empty()) {
      for (size_t i(0); i < bfs_q.front()->get_size(); ++i) {
        if ((*bfs_q.front())[i].child)
          bfs_q.push((*bfs_q.front())[i].child);
      }
      bfs_q.front()->print_node();
      bfs_q.pop();
    }
  }
  else
    cout << "The r+ tree is empty.\n";
}

//========================================NODE-IMPLEMENTATION==========================================

template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::Node::Node() {
  entries.resize(M);
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
    ALERT(error.what())
    exit(0);
  }
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::add(Entry &new_entry) {
  if (size == 0)
    mbr = new_entry.get_mbr();
  else
    mbr.adjust(new_entry.get_mbr());
  entries[size++] = new_entry;
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::add(vector<Entry> &S) {
  if (max(S.size(), M) == S.size())
    entries.resize(S.size());
  else
    if (entries.size() == 0)
      entries.resize(M);
  for (Entry &entry : S) {
    add(entry);
  }
}

template<typename T, size_t N, size_t M, size_t ff>
size_t RPlus<T, N, M, ff>::Node::get_size() {
  return size;
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::resize(size_t new_size) {
  size = new_size;
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::print_node(bool rp_root) {
  cout << "\tNODE : size(" << size << ") = [" << endl;
  cout << "\t\tA. ID : " << this << endl;
  cout << "\t\tB. Type : " << ((rp_root)? "ROOT" :((is_leaf())?"LEAF":"INTERNAL")) << endl;
  cout << "\t\tC. Boundaries : \n";
  mbr.show_rect();
  cout << "\t\tD. [Entries] : \n";
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

//Entry for root and internal nodes
template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::Entry::Entry(shared_ptr<Node> &child) {
  this->child = child;
}

//Entry for leaves
template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::Entry::Entry(HyperPoint<T, N> &data) {
  this->data = data;
}

template<typename T, size_t N, size_t M, size_t ff>
typename RPlus<T, N, M, ff>::Entry& RPlus<T, N, M, ff>::Entry::operator=(const Entry &other) {
  child = other.child;
  data = other.data;
  return *this;
}

template<typename T, size_t N, size_t M, size_t ff>
HyperRectangle<T, N> RPlus<T, N, M, ff>::Entry::get_mbr() {
  if (!is_in_leaf()) {
    return child->mbr;
  }
  return make_hyper_rect(data);
}

template<typename T, size_t N, size_t M, size_t ff>
bool RPlus<T, N, M, ff>::Entry::is_in_leaf() {
  return !child;
}

template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Entry::show_entry(size_t index) {
  cout << "\t\t" << char(192) << "->Entry<" << index << ">{\n";
  if (is_in_leaf())
    cout << "\t\t" << char(175) << " Data : ", data.show_data(), cout << "\t\t" << char(175) << " Song\'s name : " << data.get_songs_name() << endl;
  else
    cout << "\t\t" << char(175) << " Boundaries(child) : \n", get_mbr().show_rect();
  cout << "\t\t" << char(175) << " Child : " << child << endl;
  cout << "\t\t}\n";
}