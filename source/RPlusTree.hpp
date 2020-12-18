#include <rplus_utils.hpp>

#define GET_BOUNDARIES(entry) entry.get_mbr().get_boundaries()

/*Comment VISUALIZE_INSERT_COUNT if you don't want to see the insert counter ex.: [STEP] : 151623
                                                                                  [STEP] : 151624
                                                                                  [STEP] : 151625
                                                                                  [STEP] : 151626
                                                                                  [STEP] : 151627
                                                                                  ...*/
#define VISUALIZE_INSERT_COUNT

//Comment NON_REPEATED_SONGS if you want repeated songs by the id(this case is "name"), by default commented because this is a R+Tree for points, not for shapes with volume

//#define NON_REPEATED_SONGS

/*TEMPLATE PARAMETERS: (1)data type | (2)number of dimensions | (3)max entries per node | (4)fill factor(by default = 2)
  Contains: Node, Entry, comparators(ENTRYSINGLEDIM, ENTRYDIST).
  Approach: P R+ Tree (Point R+ Tree non packed) - insertion 1x1 - knn query and range query using queues and stacks.
  Features: No overlap (geometric and by saturation propagated split), structure to store hyperpoints, non repeatable data (because this structure store points).
  Link: https://github.com/italoucsp/RPlus-Tree_Proyecto-Final.
  Why not pack algorithm?: too (a lot) slow at first for entries more than 10k, Time Complexity: O(n^2/k log ff) aprox. 
                           But samely I have the code with pack algorithm ("garbage.txt").
  Operations that you are able to do: insertion(insert,"1x1"), range query(search), k-nearest neighbors query(kNN_query).*/
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
    Entry entry;//Object for the queue
    ENTRYDIST(HyperPoint<T, N> &p, Entry &md_obj){
      entry = md_obj;
      if(md_obj.get_mbr().get_hypervolume() > 0)
        distance = RPlus::MINDIST(p, md_obj.get_mbr());//page 5, rule 3, line 5 to 7 - Roussopoulos et al. suggest that when the overlap is small...
      else
        distance = RPlus::EUCDIST(p, md_obj.data);//If the entry is in leaf, then use euclidean distance because the distance is now between hyperpoints
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

  void insert(Entry &entry);
  shared_ptr<Node> choose_leaf(Entry &entry, stack<shared_ptr<Node>> &parents);
  shared_ptr<Node> split_by_parent_cut(shared_ptr<Node> &A, size_t axis, T optimal_cutline);
  shared_ptr<Node> split_by_saturation(shared_ptr<Node> &A);
  inline void partition(shared_ptr<Node> &danger_node, size_t &optimal_dim, T &optimal_cutline);
  inline pair<double, T> sweep(size_t axis, vector<Entry> &S);
  inline int min_number_splits(vector<Entry> &test_set, size_t axis, T optimal_cutline);
  inline void push_node_in_queue(HyperPoint<T, N> refdata, shared_ptr<Node> &current, priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> &q_NN);
  static double MINDIST(HyperPoint<T, N> p, HyperRectangle<T, N> r);
  static double EUCDIST(HyperPoint<T, N> p1, HyperPoint<T, N> p2);

public:
  RPlus();
  virtual ~RPlus();
  void assign(vector<HyperPoint<T, N>> &unpacked_data);
  vector<HyperPoint<T, N>> search(const HyperRectangle<T, N> &W);
  vector<HyperPoint<T, N>> kNN_query(HyperPoint<T, N> refdata, size_t k);
  void read_tree();
};

//===============================R-PLUS-TREE-IMPLEMENTATION============================================

//BUILDER RPLUS: Create empty root
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
    exit(1);
  }
}

//DESTROYER RPLUS: Simple class destroyer
template<typename T, size_t N, size_t M, size_t ff>
RPlus<T, N, M, ff>::~RPlus() {
  root.reset();
}

//RANGE QUERY METHOD: Give an hyperrectangle W and get the entries that overlaps with it.
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
#ifdef NON_REPEATED_SONGS
              size_t temp_songs_names_size = songs_names.size();
              //try pick data
              songs_names.insert((*current)[i].data.get_songs_name());
              if (temp_songs_names_size != songs_names.size())
                range_query.push_back((*current)[i].data);
              //end try
#else
              range_query.push_back((*current)[i].data);
#endif // NON_REPEATED_SONGS
            }
          }
        }
      }
      return range_query;
    }
  }
  catch (const exception &error) {
    ALERT(error.what())
      exit(1);
  }
}

/*KNN METHOD: k-Nearest Neighbors query using branch and bound algorithm with MINDIST function.
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
#ifdef NON_REPEATED_SONGS
          size_t temp_songs_names_size = songs_names.size();
          //try pick k(i) nn
          songs_names.insert(closest_entry.entry.data.get_songs_name());
          if (temp_songs_names_size != songs_names.size())
            kNN[i++] = closest_entry.entry.data;
          //end try
#else
          kNN[i++] = closest_entry.entry.data;
#endif // NON_REPEATED_SONGS
        }
      }
      return kNN;
    }
  }
  catch (const exception &error) {
    ALERT(error.what())
      exit(1);
  }
}

//--PUSH EACH ENTRY OF A NODE IN THE PRIORITY QUEUE--
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::push_node_in_queue(HyperPoint<T, N> refdata, shared_ptr<Node> &current, priority_queue<ENTRYDIST, vector<ENTRYDIST>, comparator_ENTRYDIST> &q_NN) {
  for (size_t i = size_t(0); i < current->get_size(); ++i) {
    ENTRYDIST packed_entry(refdata, (*current)[i]);
    q_NN.push(packed_entry);
  }
}

//ASSIGN METHOD: "Massive" insertion(1x1x(size of unpacked_data vector)). Give a list of hyperpoints (data) to insert in the R+
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::assign(vector<HyperPoint<T, N>> &unpacked_data) {
  size_t step_insert(1);
  for (HyperPoint<T, N> &hp : unpacked_data) {
    Entry data_entry(hp);
    insert(data_entry);
#ifdef VISUALIZE_INSERT_COUNT
    SAY(step_insert)
      ++step_insert;
#endif // VISUALIZE_INSERT_COUNT

  }
}

//INSERTION METHOD: Single insertion (1x1), need assign method to be called because it is private.
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::insert(Entry &entry) {
  stack<shared_ptr<Node>> parents;
  shared_ptr<Node> &candidate_node = choose_leaf(entry, parents);
  //if saturated node ->split, else -> simple insert
  candidate_node->add(entry);
  if (candidate_node->get_size() > M) {
    parents.push(candidate_node);
    while (parents.top()->get_size() > M) {
      shared_ptr<Node> &current_to_split = parents.top();
      parents.pop();
      Entry new_entry(split_by_saturation(current_to_split));

      if (!parents.empty()) {//parent is an internal node
        shared_ptr<Node> &currents_parent = parents.top();
        currents_parent->add(new_entry);
      }
      else {//no more parents?? -> new root = grow up the tree
        shared_ptr<Node> new_root = make_shared<Node>();
        Entry root_entry(root);
        new_root->add(root_entry); new_root->add(new_entry);
        root = new_root;
        return;
      }
    }
  }
}

//CHOOSE LEAF METHOD: Search the node to place the new entry and build a parent's path for split upward propagation
template<typename T, size_t N, size_t M, size_t ff>
shared_ptr<typename RPlus<T, N, M, ff>::Node> RPlus<T, N, M, ff>::choose_leaf(Entry &entry, stack<shared_ptr<Node>> &parents) {
  shared_ptr<Node> candidate_node = root;
  while (!candidate_node->is_leaf()) {
    parents.push(candidate_node);
    shared_ptr<Node> temp = candidate_node;
    for (size_t i(0); i < temp->get_size(); ++i) {
      candidate_node = (*temp)[i].child;
      if ((*temp)[i].get_mbr().contains(entry.data))
        break;
    }
  }
  return candidate_node;
}

/*SPLIT BY PARENT'S CUT METHOD: Division of a node A in given axis and optimal cutline,
                                then do downward propagation of the split by parent's cut.*/
template<typename T, size_t N, size_t M, size_t ff>
shared_ptr<typename RPlus<T, N, M, ff>::Node> RPlus<T, N, M, ff>::split_by_parent_cut(shared_ptr<Node> &A, size_t axis, T cutline) {
  shared_ptr<Node> B = make_shared<Node>();
  vector<Entry> set_A, set_B;
  for (size_t i(0); i < A->get_size(); ++i) {
    Entry &entry = (*A)[i];
    if (A->is_leaf()) {
      if(entry.data[axis] < cutline)
        set_A.push_back(entry);
      else if (entry.data[axis] > cutline)
        set_B.push_back(entry);
      else {
        if(set_A.size() > set_B.size())
          set_B.push_back(entry);
        else
          set_A.push_back(entry);
      }
    }
    else {
      if (GET_BOUNDARIES(entry).second[axis] <= cutline)
        set_A.push_back(entry);
      else if (GET_BOUNDARIES(entry).first[axis] >= cutline)
        set_B.push_back(entry);
      else {
        set_B.emplace_back(split_by_parent_cut(entry.child, axis, cutline));
        set_A.emplace_back(entry.child);
      }
    }
  }
  A->resize(0); A->add(set_A);
  B->resize(0); B->add(set_B);
  return B;
}

/*SPLIT BY SATURATION METHOD: When a parent node was affected by split, this could be
                              saturated (size of the node > M), so is neccessary a split
                              with a new partition line.*/
template<typename T, size_t N, size_t M, size_t ff>
shared_ptr<typename RPlus<T, N, M, ff>::Node> RPlus<T, N, M, ff>::split_by_saturation(shared_ptr<Node> &A) {
  size_t axis;
  T cutline;
  partition(A, axis, cutline);
  return split_by_parent_cut(A, axis, cutline);
}

//PARTITION METHOD: Returns the best(min. cost) cutline and axis to split a saturated node using sweep.
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::partition(shared_ptr<Node> &danger_node, size_t &optimal_dim, T &optimal_cutline) {
  double cheapest_cost = numeric_limits<double>::max();
  optimal_dim = size_t(0);
  optimal_cutline = 0;
  //sweep and find the best partition cutline/axis
  vector<Entry> S(danger_node->entries);
  for (size_t current_dim(0); current_dim < N; ++current_dim) {
    pair<double, T> cost_and_cutline = sweep(current_dim, S);
    double temp_min_cost = cheapest_cost;
    cheapest_cost = min(cost_and_cutline.first, cheapest_cost);
    if (cheapest_cost != temp_min_cost) {
      optimal_cutline = cost_and_cutline.second;
      optimal_dim = current_dim;
    }
  }
}

/*SWEEP METHOD: Using sweep line method, this algorithm returns the cost and cutline for a given axis and an entry set.
                ALG: partial sort(sweep line) + pick first ff entries + take the ff entry's max bound in the given axis as cutline.*/
template<typename T, size_t N, size_t M, size_t ff>
pair<double, T> RPlus<T, N, M, ff>::sweep(size_t axis, vector<Entry> &S) {
  comparator_ENTRYSINGLEDIM comparator(axis);
  partial_sort(S.begin(), S.begin() + ff, S.end(), comparator);//sort the first ff entries to "sweep" -> O((M + 1) log ff)
  vector<Entry> grouped_test(S.begin(), S.begin() + ff);//picking the ff first entries of the sorted set
  T optimal_cutline = GET_BOUNDARIES(grouped_test.back()).second[axis];
  return make_pair(min_number_splits(grouped_test, axis, optimal_cutline), optimal_cutline);
}

//MIN NUMBER OF SPLITS METHOD: Counts how many entries for a given set intersecs (need split) with a given cutline in given axis.
template<typename T, size_t N, size_t M, size_t ff>
int RPlus<T, N, M, ff>::min_number_splits(vector<Entry> &test_set, size_t axis, T optimal_cutline) {
  int cost = 0;
  if (!test_set[0].is_in_leaf()) {
    for (Entry &entry : test_set) {
      if (GET_BOUNDARIES(entry).first[axis] < optimal_cutline && GET_BOUNDARIES(entry).second[axis] > optimal_cutline)
        ++cost;
    }
  }
  else {
    for (Entry &entry : test_set) {
      if (entry.data[axis] == optimal_cutline)
        ++cost;
    }
  }
  return cost;
}

//MINDIST METHOD: Calculates the distance between an hyperpoint and the nearest side of an hyperrectangle.
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

//EUCDIST METHOD: Simple function to calculate euclidean distance between two hyperpoints
template<typename T, size_t N, size_t M, size_t ff>
double RPlus<T, N, M, ff>::EUCDIST(HyperPoint<T, N> p1, HyperPoint<T, N> p2) {
  double d = 0.0;
  for (size_t i(0); i < N; ++i) {
    d += pow(p1[i] - p2[i], 2);
  }
  return sqrt(d);
}

//READ TREE METHOD: Using bfs, read the levels of the tree since the root.
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

//Access to the entries of a node by index
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
    exit(1);
  }
}

//add single entry
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::add(Entry &new_entry) {
  if (size == 0) {
    entries.resize(M);
    mbr = new_entry.get_mbr();
    entries[size++] = new_entry;
  }
  else {
    if (size >= M) {
      entries.resize(size + 1);
      entries[size++] = new_entry;//saturated - temporaly break the rule : M entries per node as max
    }
    else {
      mbr.adjust(new_entry.get_mbr());
      entries[size++] = new_entry;
    }
  }
}

//add many entries
template<typename T, size_t N, size_t M, size_t ff>
void RPlus<T, N, M, ff>::Node::add(vector<Entry> &S) {
  for (Entry &entry : S) {
    add(entry);
  }
}

//Returns how many active entries has the node
template<typename T, size_t N, size_t M, size_t ff>
size_t RPlus<T, N, M, ff>::Node::get_size() {
  return size;
}

//Change how many active entries has the node
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

//Copy for entry
template<typename T, size_t N, size_t M, size_t ff>
typename RPlus<T, N, M, ff>::Entry& RPlus<T, N, M, ff>::Entry::operator=(const Entry &other) {
  child = other.child;
  data = other.data;
  return *this;
}

//If the entry is in a leaf node -> returns data made hyperrectangle (0 volume), else -> returns MBR of its child
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