#include <rplus_utils.hpp>

#define COST double

//DATA, DIMENSIONS, MAX NUMBER OF ENTRIES PER NODE
template<typename DATA_TYPE, size_t N, size_t M>
class RPlus {
  struct Node;

  struct Entry {
    DATA_TYPE data;
    shared_ptr<Node> child;

    HyperRectangle<N> get_mbr();
  };

  struct Node {
    HyperRectangle<N> mbr;
    ENTRY_GROUP entries;
    bool isLeaf();

    Node();
  };

  shared_ptr<Node> root;

  
  void hidden_insert(shared_ptr<Node> &R, Entry &IR);
  void hidden_remove(shared_ptr<Node> &R, Entry &IR);
  vector<DATA_TYPE> hidden_search(shared_ptr<Node> &R, const HyperRectangle<N> &W);
  void split_node(shared_ptr<Node> &R);
  pair<shared_ptr<Node>, ENTRY_GROUP> partition(ENTRY_GROUP &S, unsigned int ff);
  pair<COST, double> sweep(size_t axis, DATA_TYPE Okd, unsigned int ff);
  void pack(ENTRY_GROUP &S, unsigned int ff);

public:
  RPlus();
  ~RPlus();
  void insert(const DATA_TYPE data);
  void remove(const DATA_TYPE data);
  vector<DATA_TYPE> search(const HyperRectangle<N> &W);
  vector<DATA_TYPE> kNN_query(DATA_TYPE refdata, size_t k);
};

//===============================R-PLUS-TREE-IMPLEMENTATION============================================

template<typename DATA_TYPE, size_t N, size_t M>
RPlus<DATA_TYPE, N, M>::RPlus() {
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

template<typename DATA_TYPE, size_t N, size_t M>
RPlus<DATA_TYPE, N,  M>::~RPlus() {
  root.reset();
}

template<typename DATA_TYPE, size_t N, size_t M>
void RPlus<DATA_TYPE, N,  M>::insert(const DATA_TYPE data) {
  shared_ptr<Node> R = root;
  Entry new_entry;
  new_entry.data = data;
  hidden_insert(R, new_entry);
}

template<typename DATA_TYPE, size_t N, size_t M>
void RPlus<DATA_TYPE, N,  M>::remove(const DATA_TYPE data) {
  shared_ptr<Node> R = root;
  hidden_remove(R, IR);
}

template<typename DATA_TYPE, size_t N, size_t M>
vector<DATA_TYPE> RPlus<DATA_TYPE, N,  M>::search(const HyperRectangle<N> &W) {
  shared_ptr<Node> R = root;
  return hidden_search(R, W);
}

template<typename DATA_TYPE, size_t N, size_t M>
vector<DATA_TYPE> RPlus<DATA_TYPE, N,  M>::hidden_search(shared_ptr<Node> &R, const HyperRectangle<N> &W) {
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

//========================================NODE-IMPLEMENTATION==========================================

template<typename DATA_TYPE, size_t N, size_t M>
RPlus<DATA_TYPE, N,  M>::Node::Node() {
  entries.resize(M);
}

template<typename DATA_TYPE, size_t N, size_t M>
bool RPlus<DATA_TYPE, N,  M>::Node::isLeaf() {
  return true;
}