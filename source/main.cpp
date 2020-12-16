#include <RPlusTree.hpp>

//TIMER CLASS - Copyright 2020 Roger Peralta Aranibar Advanced Data Estructures
template <typename>
class Timer;

template <typename R, typename... T>
class Timer<R(T...)> {
public:
  typedef R(*function_type)(T...);
  function_type function;

  explicit Timer(function_type function, std::string process_name = "")
    : function_(function), process_name_(process_name) {}

  R operator()(T... args) {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    start = std::chrono::high_resolution_clock::now();

    R result = function_(std::forward<T>(args)...);

    end = std::chrono::high_resolution_clock::now();
    int64_t duration =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start)
      .count();

    std::cout << std::setw(10) << process_name_ << std::setw(30)
      << "Duration: " + std::to_string(duration) + " ns\n";
    return result;
  }

private:
  function_type function_;
  std::string process_name_;
};

//#######################################################################

RPlus<float, 14, 128, 128> demo;

int build_data_structure() {
  return 0;
}

int query_knn() {
  return 0;
}

int main() {
  RPlus<double, 2, 3, 3> demo_1;
  HyperPoint<double, 2> p1(array<double, 2>{1.2, 5.3}, "dwe"), p2(array<double, 2>{4.9, 2.7}, "ade"), p3(array<double, 2>{4.8, 25.2}, "dde"), p4(array<double, 2>{8.8, 12.4}, "zde");
  HyperPoint<double, 2> p5(array<double, 2>{15.2, 6.3}, "dfe"), p6(array<double, 2>{102.5, 77.7}, "bde"), p7(array<double, 2>{-6.8, -2.45}, "ede"), p8(array<double, 2>{-18.9, -17.12}, "xde");
  HyperPoint<double, 2> p9(array<double, 2>{5.2, 5.4}, "det"), p10(array<double, 2>{4.8, 33.7}, "cde"), p11(array<double, 2>{4.5, 25.01}, "gde"), p12(array<double, 2>{7.8, 12.5}, "yde");
  HyperRectangle<double, 2> rect1(p1, 1);
  HyperRectangle<double, 2> rect2(p2, 1);
  HyperRectangle<double, 2> rect3(p3, 1);
  HyperRectangle<double, 2> rect4(p4, 1);
  HyperRectangle<double, 2> rect5(p5, 1);
  HyperRectangle<double, 2> rect6(p6, 1);
  HyperRectangle<double, 2> rect7(p7, 1);
  HyperRectangle<double, 2> rect8(p8, 1);
  HyperRectangle<double, 2> rect9(p9, 1);
  HyperRectangle<double, 2> rect10(p10, 1);
  HyperRectangle<double, 2> rect11(p11, 1);
  HyperRectangle<double, 2> rect12(p12, 1);
  vector<HyperRectangle<double, 2>> pd = { rect1, rect2, rect3, rect4,rect5, rect6, rect7, rect8,rect9, rect10, rect11, rect12 };
  demo_1.insert(pd);
  SAY("insertion_done");
  demo_1.read_tree();
  SAY("done");
  HyperPoint<double, 2> pr(array<double, 2>{1.3, 5.5});
  HyperRectangle<double, 2> W(p11, 5);
  for (auto res : demo_1.search(W)) {
    res.show_data();
  }
  for (auto res : demo_1.kNN_query(pr, 5)) {
    res.show_data();
  }
  /*
  vector<HyperRectangle<double, KUSED_DIMENSIONS>> DB_CONTAINER;
  vector<string> considered_features = { "name", "acousticness", "danceability", "duration_ms", "energy",
                                         "explicit", "instrumentalness", "key", "liveness", "loudness", "mode",
                                         "popularity", "speechiness", "tempo", "valence" };
  read_data_from_file("data_set_example.txt", considered_features, "name", DB_CONTAINER);
  for (auto i : DB_CONTAINER) {
    i.show_rect();
  }*/
  return 0;
}