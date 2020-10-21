#include <rplus_utils.hpp>

//DIMENSIONS, MAX NUMBER OF ENTRIES PER NODE
template<typename DATA_TYPE, int N_,int M_>
class RPlus {
  size_t N, M, m;

  struct Node;

  struct Entry {
    HyperRectangle<size_t(N_)> entry_rect;
    DATA_TYPE data;
    shared_ptr<Node> child;
  };

  struct Node {
    HyperRectangle<size_t(N_)> node_rect;
    ENTRY_GROUP entries;
    bool isLeaf();

    Node(size_t M);
  };

  shared_ptr<Node> root;

  void hidden_insert(shared_ptr<Node> &R, const HyperRectangle<size_t(N_)> &IR);
  void hidden_remove(shared_ptr<Node> &R, const HyperRectangle<size_t(N_)> &IR);
  vector<DATA_TYPE> hidden_search(shared_ptr<Node> &R, const HyperRectangle<size_t(N_)> &W);
  void split_node(shared_ptr<Node> &R);
  pair<shared_ptr<Node>, RECT_GROUP> partition(RECT_GROUP &S, unsigned int ff);
  pair<double, size_t> sweep(size_t axis, size_t Okd, unsigned int ff);
  void pack(RECT_GROUP &S, unsigned int ff);
public:
  RPlus();
  ~RPlus();
  void insert(const HyperRectangle<size_t(N_)> &IR);
  void remove(const HyperRectangle<size_t(N_)> &IR);
  vector<DATA_TYPE> search(const HyperRectangle<size_t(N_)> &W);
  void kNN_query(size_t k);
};

//===============================R-PLUS-TREE-IMPLEMENTATION============================================

template<typename DATA_TYPE, int N_, int M_>
RPlus<DATA_TYPE, N_, M_>::RPlus() {
  try {
    if (N_ < 2 || M_ < 2) {
      throw runtime_error(ERROR_M_N_VALUES);
    }
    else {
      M = size_t(M_);
      N = size_t(N_);
      m = M / 2;
      root = make_shared<Node>(M);
    }
  }
  catch(const exception &error){
    ALERT(error.what());
  }
}

template<typename DATA_TYPE, int N_, int M_>
RPlus<DATA_TYPE, N_, M_>::~RPlus() {
  root.reset();
}

template<typename DATA_TYPE, int N_, int M_>
void RPlus<DATA_TYPE, N_, M_>::insert(const HyperRectangle<size_t(N_)> &IR) {
  shared_ptr<Node> R = root;
  hidden_insert(R, IR);
}

template<typename DATA_TYPE, int N_, int M_>
void RPlus<DATA_TYPE, N_, M_>::remove(const HyperRectangle<size_t(N_)> &IR) {
  shared_ptr<Node> R = root;
  hidden_remove(R, IR);
}

template<typename DATA_TYPE, int N_, int M_>
vector<DATA_TYPE> RPlus<DATA_TYPE, N_, M_>::search(const HyperRectangle<size_t(N_)> &W) {
  shared_ptr<Node> R = root;
  return hidden_search(R, W);
}

template<typename DATA_TYPE, int N_, int M_>
vector<DATA_TYPE> RPlus<DATA_TYPE, N_, M_>::hidden_search(shared_ptr<Node> &R, const HyperRectangle<size_t(N_)> &W) {
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

template<typename DATA_TYPE, int N_, int M_>
RPlus<DATA_TYPE, N_, M_>::Node::Node(size_t M) {
  entries.resize(M);
}

template<typename DATA_TYPE, int N_, int M_>
bool RPlus<DATA_TYPE, N_, M_>::Node::isLeaf() {
  return true;
}