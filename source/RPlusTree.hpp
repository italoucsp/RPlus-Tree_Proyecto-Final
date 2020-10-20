#include <rplus_utils.hpp>

//DIMENSIONS, MAX NUMBER OF ENTRIES PER NODE
template<typename DATA_TYPE, int N_,int M_>
class RPlus {
  size_t N, M, m;

  struct Node;

  struct Entry {
    HyperRectangle<N_> entry_rect;
    DATA_TYPE data;
    shared_ptr<Node> child;
  };

  struct Node {
    HyperRectangle<N_> node_rect;
    ENTRY_GROUP entries;

    Node(size_t M);
  };

  shared_ptr<Node> root;

  void hidden_search(shared_ptr<Node> &R, const HyperRectangle<N_> &W);
  void split_node(shared_ptr<Node> &R);
  pair<shared_ptr<Node>, RECT_GROUP> partition(RECT_GROUP &S, unsigned int ff);
  pair<double,DATA_TYPE> sweep(size_t axis, DATA_TYPE Okd, unsigned int ff);
  void pack(RECT_GROUP &S, unsigned int ff);
public:
  RPlus();
  ~RPlus();
  void search(const HyperRectangle<N_> &W);
  void kNN_query(size_t k);
  void insert(shared_ptr<Node> &R, const HyperRectangle<N_> &IR);
  void remove(shared_ptr<Node> &R, const HyperRectangle<N_> &IR);
};

//===============================R-PLUS-TREE-IMPLEMENTATION============================================

template<typename T, int N_, int M_>
RPlus<T, N_, M_>::RPlus() {
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

template<typename T, int N_, int M_>
RPlus<T, N_, M_>::~RPlus() {
  root.reset();
}

//========================================NODE-IMPLEMENTATION==========================================

template<typename T, int N_, int M_>
RPlus<T, N_, M_>::Node::Node(size_t M) {
  entries.resize(M);
}