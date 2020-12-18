#include <RPlusTree.hpp>

/*
-----------------------[ADVANCED DATA STRUCTURES]------------------------
-----------------SPOTIFY DATASET 1921-2020, 160k+TRACKS------------------
-------------------------DATA-STRUCT.: R+ Tree---------------------------
------------------------Italo Mamani Huaricallo--------------------------
---------------------------Arequipa-Perú-2020----------------------------
*/

//TIMER CLASS - Copyright 2020 Roger Peralta Aranibar Advanced Data Structures
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

//##########################GLOBAL-RPLUS-TREE#############################
//Go to the file "RPlusTree.hpp to see the code for R+ and documentation
//Go to the file "rplus_utils.hpp to see the code for R+ tools
//Best try tested: M = 16, ff = 8
RPlus<double, KUSED_DIMENSIONS, 16, 8> demo;

vector<HyperPoint<double, KUSED_DIMENSIONS>> DB_CONTAINER;//To store the dataset for insertion

//#########################FUNCTIONS-FOR-TIMER############################

int build_data_structure() {
  demo.assign(DB_CONTAINER);
  return 0;
}

vector<HyperPoint<double, KUSED_DIMENSIONS>> query_knn(HyperPoint<double, KUSED_DIMENSIONS> query, size_t k) {
  return demo.kNN_query(query, k);
}

int main() {
  //Vector that contains the 14 USED DIMENSIONS + 1(ID = "name")
  vector<string> considered_features = { "name", "acousticness", "danceability", "duration_ms", "energy",
                                         "explicit", "instrumentalness", "key", "liveness", "loudness",
                                         "mode", "popularity", "speechiness", "tempo", "valence" };

  //WARNING[!] : The delimiter in the csv file must be ';'.

    read_data_from_file("data_set.csv", considered_features, "name", DB_CONTAINER);
    SAY("Data extraction done")

  //=============================TESTING=CAMP==============================

    Timer<int()> timed_built(build_data_structure, "R+ Insertion");
    timed_built();

    Timer<vector<HyperPoint<double, KUSED_DIMENSIONS>>(HyperPoint<double, KUSED_DIMENSIONS>, size_t)> timed_query(query_knn, "R+ Query kNN");
    //EXAMPLE
    //HyperPoint<double, KUSED_DIMENSIONS> query(array<double, KUSED_DIMENSIONS>{1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 5, 6, 8, 10});
    //OR HyperPoint<double, KUSED_DIMENSIONS> query; cin >> query; input: 1 2 3 4 5 6 7 8 9 0 5 6 8 10
    int option;
    while (true) {
      SAY("Write 1 if you want to do a knn query or 0 if you want to exit:")
        cin >> option;
        if (option) {
          SAY("Write the content of the query separated by space:")
          HyperPoint<double, KUSED_DIMENSIONS> query; cin >> query;
          SAY("Write how many neigbors do you want(k):")
          size_t k; cin >> k;
          vector<HyperPoint<double, KUSED_DIMENSIONS>> result = timed_query(query, k);
          for (auto &r : result) {
            r.show_data();
            cout << "song\'s name : " << r.get_songs_name() << endl;
          }
        }
        else {
          exit(0);
        }
    }

  return 0;
}