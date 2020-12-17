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
  vector<HyperPoint<double, KUSED_DIMENSIONS>> DB_CONTAINER;
  vector<string> considered_features = { "name", "acousticness", "danceability", "duration_ms", "energy",
                                         "explicit", "instrumentalness", "key", "liveness", "loudness", "mode",
                                         "popularity", "speechiness", "tempo", "valence" };
  read_data_from_file("data_set.csv", considered_features, "name", DB_CONTAINER);
  SAY("extraction_done")
  RPlus<double, KUSED_DIMENSIONS, 16, 8> demo_2;
  demo_2.assign(DB_CONTAINER);
  SAY("insertion_done")
  //demo_2.read_tree();
  /*for (auto i : DB_CONTAINER) {
    i.show_data();
    cout << i.get_songs_name() << endl;
  }*/
  return 0;
}