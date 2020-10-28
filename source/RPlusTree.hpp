#include <rplus_utils.hpp>

#define COST double

//DATA, DIMENSIONS, MAX NUMBER OF ENTRIES PER NODE
template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
class RPlus {
  struct Node;

  struct Entry {
    shared_ptr<DATA_TYPE> data;
    shared_ptr<Node> child;

    HyperRectangle<N> get_mbr();
  };

  struct Node {
    HyperRectangle<N> mbr;
    ENTRY_GROUP entries;
    bool isLeaf();

    Node();
    void add(Entry &IR);
    size_t get_size();

    size_t size;
  };

  shared_ptr<Node> root;
  
  void hidden_insert(shared_ptr<Node> &R, Entry &IR);
  void hidden_remove(shared_ptr<Node> &R, Entry &IR);
  vector<DATA_TYPE> hidden_search(shared_ptr<Node> &R, const HyperRectangle<N> &W);
  void split_node(shared_ptr<Node> &R);
  pair<shared_ptr<Node>, ENTRY_GROUP> partition(ENTRY_GROUP &S);
  pair<COST, double> sweep(size_t axis, DATA_TYPE Okd);
  shared_ptr<Node> pack(ENTRY_GROUP &S);

public:
  RPlus();
  ~RPlus();
  void insert(const DATA_TYPE data);
  void remove(const DATA_TYPE data);
  vector<DATA_TYPE> search(const HyperRectangle<N> &W);
  vector<DATA_TYPE> kNN_query(DATA_TYPE refdata, size_t k);
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

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::insert(const DATA_TYPE data) {
  shared_ptr<Node> R = root;
  Entry new_entry;
  new_entry.data = data;
  hidden_insert(R, new_entry);
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::remove(const DATA_TYPE data) {
  shared_ptr<Node> R = root;
  //hidden_remove(R, IR);
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
vector<DATA_TYPE> RPlus<DATA_TYPE, N, M, ff>::search(const HyperRectangle<N> &W) {
  shared_ptr<Node> R = root;
  return hidden_search(R, W);
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
vector<DATA_TYPE> RPlus<DATA_TYPE, N, M, ff>::kNN_query(DATA_TYPE refdata, size_t k) {

}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
vector<DATA_TYPE> RPlus<DATA_TYPE, N, M, ff>::hidden_search(shared_ptr<Node> &R, const HyperRectangle<N> &W) {
  vector<DATA_TYPE> data_found;
  if (!R->isLeaf()) {
    for (Entry &entry : R->entries) {
      if (entry.entry_rect.isOverlaping(W)) {
        hidden_search(entry.child,W);
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

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::hidden_insert(shared_ptr<Node> &R, Entry &IR) {
  if (!R->isLeaf()) {
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

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::hidden_remove(shared_ptr<Node> &R, Entry &IR) {

}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::split_node(shared_ptr<Node> &R) {
  /*SN1.[Find a Partition]
      Partition R using the Partition routine of the
      Pack algorithm(see next section).Let RECT and
      p be the rectangle and pointer respectively associated with node R.Also, let S1 and S2 denote the
      two sub - regions resulting after the partition.
      Create n1 = (p1, RECT1) and n2 = (p2, RECT2), the
      two nodes resulting from the split of R, where
      RECTi = RECT.overlap(Si), for i = 1, 2.*/
  shared_ptr<Node> N1 = make_shared<Node>();
  shared_ptr<Node> N2 = make_shared<Node>();
  /*SN2.[Populate New Nodes]
      Put in ni all nodes(pk, RECTk) of R such that RECTk
      lies completely in RECTi, for i = 1, 2. For
      those nodes that RECTk n RECTi no= RECTk
      (i.e.they just overlap with the sub - region)
      a) if R is a leaf node, then put RECTk in both new nodes
      b) Otherwise, use SplitNode to recursively split the children nodes along the partition.Let
      (pk 1, RECTk 1) and (pk 2, RECTk 2) be the two
      nodes after splitting(pk, RECTk), where RECTki
      lies completely in RECTi, i = 1, 2. Add those two
      nodes to the corresponding node ni.*/
  
  /*SN3.[Propagate Node Split Upward]
      If R is the root, create a new root with only two
      children, n1 and n2.
      Otherwise, let PR be R’s parent node.Replace R
      in PR with n1 and n2.If PR has now more than
      M entries, invoke SplitNode(PR).*/
  if (R == root) {
    shared_ptr<Node> new_root = make_shared<Node>();
  }
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
pair<shared_ptr<RPlus<DATA_TYPE, N, M, ff>::Node>, vector<RPlus<DATA_TYPE, N, M, ff>::Entry>> RPlus<DATA_TYPE, N, M, ff>::partition(ENTRY_GROUP &S) {
  if (S.size() <= ff) {
    shared_ptr<Node> newnode = make_shared<Node>();
    ENTRY_GROUP empty_group;
    return make_pair(newnode, empty_group);
  }
  vector<COST> cost_for_all_axises;
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
pair<COST, double> RPlus<DATA_TYPE, N, M, ff>::sweep(size_t axis, DATA_TYPE Okd) {

}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
shared_ptr<RPlus<DATA_TYPE, N, M, ff>::Node> RPlus<DATA_TYPE, N, M, ff>::pack(ENTRY_GROUP &S) {
  if (S.size() <= ff) {
    shared_ptr<Node> newnode = make_shared<Node>();
    for (Entry &entry : S) {
      newnode.add(entry);
    }
    return newnode;
  }
  ENTRY_GROUP AN;
  while (!S.empty()) {
    partition(S).second;
  }
  return pack(AN);
}

//========================================NODE-IMPLEMENTATION==========================================

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
RPlus<DATA_TYPE, N, M, ff>::Node::Node() {
  entries.resize(M + 1);
  size = size_t(0);
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
bool RPlus<DATA_TYPE, N, M, ff>::Node::isLeaf() {
  return child == nullptr;
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
void RPlus<DATA_TYPE, N, M, ff>::Node::add(Entry &IR) {
  entries[size++] = IR;
}

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
size_t RPlus<DATA_TYPE, N, M, ff>::Node::get_size() {
  return size;
}

//=======================================ENTRY-IMPLEMENTATION==========================================

template<typename DATA_TYPE, size_t N, size_t M, size_t ff>
HyperRectangle<N> RPlus<DATA_TYPE, N, M, ff>::Entry::get_mbr() {
  return child->mbr;
}