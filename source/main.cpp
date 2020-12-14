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

int main() {
  RPlus<double, 2, 3, 1> demo_1;
  HyperPoint<double, 2> p1(array<double, 2>{1.2,5.3}), p2(array<double, 2>{4.9, 2.7}), p3(array<double, 2>{4.8, 25.2}), p4(array<double, 2>{8.8, 12.4});
  HyperRectangle<double, 2> rect1(p1, 0.1);
  HyperRectangle<double, 2> rect2(p2, 0.1);
  HyperRectangle<double, 2> rect3(p3, 0.1);
  HyperRectangle<double, 2> rect4(p4, 0.1);
  vector<HyperRectangle<double, 2>> pd = {rect1, rect2, rect3, rect4};
  demo_1.insert(pd);
  ALERT("insertion_done");
  demo_1.read_tree();
  ALERT("done");
  return 0;
}