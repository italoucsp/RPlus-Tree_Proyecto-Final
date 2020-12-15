#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <queue>
#include <sstream>
#include <stdarg.h>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace std;

#define ALERT(message) cout << "[ERROR] : " << message << endl;
#define SAY(message) cout << "[STEP] : " << message << endl;
#define ERROR_M_N_VALUES "The number of dimensions and/or the number of entries per node should be greater than 1."
#define ERROR_FF_VALUE "The value for fill factor should be between 2 and M."
#define ERROR_NODE_OFR "The index is out of range in the node."
#define ERROR_NO_RECTDATA "This rectangle is not a data container."

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, size_t N>
struct HyperPoint {
  HyperPoint();
  HyperPoint(array<T, N> data);
  HyperPoint(array<T, N> data, string sg_name);
  bool operator<(const HyperPoint<T, N> &other) const;
  bool operator>(const HyperPoint<T, N> &other) const;
  HyperPoint<T, N>& operator=(const HyperPoint<T, N> &other);
  T& operator[](size_t index);
  T operator[](size_t index) const;
  void show_data();

private:
  array<T, N> multidata;
  string songs_name;
};

template<typename T, size_t N>
HyperPoint<T, N>::HyperPoint() {
  multidata.fill(T(0));
}

template<typename T, size_t N>
HyperPoint<T, N>::HyperPoint(array<T, N> data) {
  for (size_t i(0); i < N; ++i) {
    multidata[i] = data[i];
  }
}

template<typename T, size_t N>
HyperPoint<T, N>::HyperPoint(array<T, N> data, string sg_name) {
  for (size_t i(0); i < N; ++i) {
    multidata[i] = data[i];
  }
  songs_name = sg_name;
}

template<typename T, size_t N>
bool HyperPoint<T, N>::operator<(const HyperPoint<T, N> &other) const{
  for (size_t i(0); i < N; ++i) {
    if (multidata[i] >= other.multidata[i]) {
      return false;
    }
  }
  return true;
}

template<typename T, size_t N>
bool HyperPoint<T, N>::operator>(const HyperPoint<T, N> &other) const {
  return !(*this < other);
}

template<typename T, size_t N>
HyperPoint<T, N>& HyperPoint<T, N>::operator=(const HyperPoint<T, N>& other) {
  songs_name = other.songs_name;
  for (size_t i(0); i < N; ++i)
    multidata[i] = other.multidata[i];
  return *this;
}

template<typename T, size_t N>
T& HyperPoint<T, N>::operator[](size_t index) {
  return multidata[index];
}

template<typename T, size_t N>
T HyperPoint<T, N>::operator[](size_t index) const{
  return multidata[index];
}

template<typename T, size_t N>
void HyperPoint<T, N>::show_data() {
  cout << "\tHyperPoint<" << N << "> : {"; for (size_t i(0); i < N; ++i) cout << setprecision(10) << multidata[i] << ((i != N - 1) ? "," : ""); cout << "}" << endl;
}

template<typename T, size_t N>
HyperPoint<T, N> get_min(const HyperPoint<T, N> A, const HyperPoint<T, N> B) {
  return (A < B) ? A : B;
}

template<typename T, size_t N>
HyperPoint<T, N> get_max(const HyperPoint<T, N> A, const HyperPoint<T, N> B) {
  return (A > B) ? A : B;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, size_t N>
struct HyperRectangle {
  HyperRectangle();
  HyperRectangle(HyperPoint<T, N> &DATA, T fatness);
  HyperRectangle(HyperPoint<T, N> &A, HyperPoint<T, N> &B);
  HyperRectangle<T, N>& operator=(const HyperRectangle<T, N>& other);
  const HyperPoint<T, N>& get_data();
  bool overlaps(const HyperRectangle<T, N> &other);
  bool contains(const HyperPoint<T, N> &point);
  void adjust(const HyperRectangle<T, N> &other);
  pair<HyperRectangle<T, N>, HyperRectangle<T, N>> cut(size_t axis, T cutline);
  pair<HyperPoint<T, N>, HyperPoint<T, N>> get_boundaries();
  double get_hypervolume();
  bool data_container();
  void show_rect();

private:
  HyperPoint<T, N> bottom_left, top_right;
  HyperPoint<T, N> data;
  bool is_data;
};

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle() {
  is_data = false;
}

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle(HyperPoint<T, N> &DATA, T fatness) {
  is_data = true;
  data = DATA;
  for (size_t i(0); i < N; ++i) {
    bottom_left[i] = DATA[i] - fatness;
    top_right[i] = DATA[i] + fatness;
  }
}

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle(HyperPoint<T, N> &A, HyperPoint<T, N> &B) {
  is_data = false;
  for (size_t i(0); i < N; ++i) {
    bottom_left[i] = min(A[i], B[i]);
    top_right[i] = max(A[i], B[i]);
  }
}

template<typename T, size_t N>
HyperRectangle<T, N>& HyperRectangle<T, N>::operator=(const HyperRectangle<T, N>& other) {
  bottom_left = other.bottom_left;
  top_right = other.top_right;
  if (is_data == other.is_data && is_data) {
    data = other.data;
  }
  return *this;
}

template<typename T, size_t N>
const HyperPoint<T, N>& HyperRectangle<T, N>::get_data() {
  try {
    if (!data_container())
      throw runtime_error(ERROR_NO_RECTDATA);
    else
      return data;
  }
  catch (const exception &error) {
    ALERT(error.what());
    exit(0);
  }
}

template<typename T, size_t N>
bool HyperRectangle<T, N>::overlaps(const HyperRectangle<T, N> &other) {
  for (size_t i(0); i < N; ++i) {
    if (!(bottom_left[i] <= other.top_right[i] &&
          top_right[i] >= other.bottom_left[i])) {
      return false;
    }
  }
  return true;
}

template<typename T, size_t N>
bool HyperRectangle<T, N>::contains(const HyperPoint<T, N> &point) {
  for (size_t i(0); i < N; ++i) {
    if (!(bottom_left[i] <= point[i] && point[i] <= top_right[i]))
      return false;
  }
  return true;
}

template<typename T, size_t N>
void HyperRectangle<T, N>::adjust(const HyperRectangle<T, N> &other) {
  for (size_t i(0); i < N; ++i) {
    bottom_left[i] = min(other.bottom_left[i], bottom_left[i]);
    top_right[i] = max(other.top_right[i], top_right[i]);
  }
}

template<typename T, size_t N>
pair<HyperRectangle<T, N>, HyperRectangle<T, N>> HyperRectangle<T, N>::cut(size_t axis, T cutline) {
  HyperPoint<T, N> f_top_right = top_right, s_bottom_left = bottom_left;
  f_top_right[axis] = s_bottom_left[axis] = cutline;
  HyperRectangle<T, N> half_part_f(bottom_left, f_top_right), half_part_s(s_bottom_left, top_right);
  if (is_data) {
    half_part_f.is_data = half_part_s.is_data = true;
    half_part_f.data = half_part_s.data = data;
  }
  return make_pair(half_part_f, half_part_s);
}

template<typename T, size_t N>
pair<HyperPoint<T, N>, HyperPoint<T, N>> HyperRectangle<T, N>::get_boundaries() {
  return make_pair(bottom_left, top_right);
}

template<typename T, size_t N>
double HyperRectangle<T, N>::get_hypervolume() {
  double hypervolume = 1.0;
  for (size_t i(0); i < N; ++i) {
    hypervolume *= double(top_right[i] - bottom_left[i]);
  }
  return hypervolume;
}

template<typename T, size_t N>
bool HyperRectangle<T, N>::data_container() {
  return is_data;
}

template<typename T, size_t N>
void HyperRectangle<T, N>::show_rect() {
  cout << "\t\tHyperRectangle<" << N << "> :   bottom_left: "; bottom_left.show_data(); cout << "\t\t\t\t\ttop_right: "; top_right.show_data();
  if (is_data) {
    cout << "\t\t\tDATA:\n" << "\t\t\t";
    data.show_data();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, size_t N>
void read_data_from_file(string file_path, vector<HyperRectangle<T, N>> &db_container) {
  ifstream data_set_file(file_path);
  string single_line;
  while (!data_set_file.eof()) {

  }
}