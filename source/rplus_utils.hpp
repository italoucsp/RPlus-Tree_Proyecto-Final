#include <algorithm>
#include <array>

#include <fstream>

#include <iomanip>
#include <iostream>

#include <math.h>
#include <memory>

#include <queue>

#include <sstream>
#include <stack>
#include <stdarg.h>
#include <stdexcept>
#include <string>

#include <tuple>

#include <unordered_set>
#include <utility>

#include <vector>

using namespace std;

#define ALERT(message) cerr << "[ERROR] : " << message << endl;
#define SAY(message) cout << "[STEP] : " << message << endl;
#define ERROR_M_N_VALUES "The number of dimensions and/or the number of entries per node should be greater than 1."
#define ERROR_FF_VALUE "The value for fill factor should be between 2 and M."
#define ERROR_NODE_OFR "The index is out of range in the node."
#define ERROR_EMPTY_TREE "This R+ Tree is empty."

const size_t T_DIMENSIONS_NUM = 19;
const size_t KUSED_DIMENSIONS = 14;

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
  string get_songs_name();
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
string HyperPoint<T, N>::get_songs_name() {
  return songs_name;
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
  HyperRectangle(HyperPoint<T, N> &A, HyperPoint<T, N> &B);
  HyperRectangle<T, N>& operator=(const HyperRectangle<T, N>& other);
  bool overlaps(const HyperRectangle<T, N> &other);
  bool contains(const HyperPoint<T, N> &point);
  void adjust(const HyperRectangle<T, N> &other);
  pair<HyperRectangle<T, N>, HyperRectangle<T, N>> cut(size_t axis, T cutline);
  pair<HyperPoint<T, N>, HyperPoint<T, N>> get_boundaries();
  double get_hypervolume();
  void show_rect();

  template<typename T, size_t N>
  friend HyperRectangle<T, N> make_hyper_rect(HyperPoint<T, N> &h_point);

private:
  HyperPoint<T, N> bottom_left, top_right;
};

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle() {
  //...
}

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle(HyperPoint<T, N> &A, HyperPoint<T, N> &B) {
  for (size_t i(0); i < N; ++i) {
    bottom_left[i] = min(A[i], B[i]);
    top_right[i] = max(A[i], B[i]);
  }
}

template<typename T, size_t N>
HyperRectangle<T, N>& HyperRectangle<T, N>::operator=(const HyperRectangle<T, N>& other) {
  bottom_left = other.bottom_left;
  top_right = other.top_right;
  return *this;
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
void HyperRectangle<T, N>::show_rect() {
  cout << "\t\tHyperRectangle<" << N << "> :   bottom_left: "; bottom_left.show_data();
  cout << "\t\t\t\t\ttop_right: "; top_right.show_data();
}

template<typename T, size_t N>
HyperRectangle<T, N> make_hyper_rect(HyperPoint<T, N> &h_point) {
  HyperRectangle<T, N> conv;
  conv.bottom_left = conv.top_right = h_point;
  return conv;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//CSV file reader : path of the file | features that were considered | container for the data in hyperrectangles
template<typename T, size_t N>
void read_data_from_file(string file_path, vector<string> &features, string id, vector<HyperPoint<T, N>> &db_container) {
  string cols_labels, row_data_line;
  ifstream data_set_file(file_path);
  if (!data_set_file.is_open()) {
    ALERT("Couldn\'t open the file " + file_path)
    exit(1);
  }
  if (N + 1 < features.size()) {
    ALERT("The value of N is not enough for the given features.")
    exit(1);
  }

  getline(data_set_file, cols_labels);
  string label;
  istringstream iss_cols_labels(cols_labels);
  array<bool, T_DIMENSIONS_NUM> table_check_features;
  table_check_features.fill(false);
  size_t tcfi(0), id_id(0);
  while (!iss_cols_labels.eof()) {
    getline(iss_cols_labels, label, ',');
    for (size_t fi(0); fi < features.size(); ++fi) {
      if (features[fi] == label) {
        table_check_features[tcfi] = true;
        if (label == id)
          id_id = tcfi;
      }
    }
    ++tcfi;
  }

  while (getline(data_set_file, row_data_line)) {
    istringstream iss_cols_data(row_data_line);
    size_t i(0), ii(0);
    string data_col_row, songs_name;
    array<T, KUSED_DIMENSIONS> raw_multidimensional_data;
    while (!iss_cols_data.eof()) {
      getline(iss_cols_data, data_col_row, ',');
      if (i < T_DIMENSIONS_NUM && table_check_features[i]) {
        if (id_id == i)
          songs_name = data_col_row;
        else {
          stringstream ss_raw_datacol(data_col_row);
          ss_raw_datacol >> raw_multidimensional_data[ii];
          ++ii;
        }
      }
      ++i;
    }
    HyperPoint<T, KUSED_DIMENSIONS> multidimensional_data(raw_multidimensional_data, songs_name);
    db_container.push_back(multidimensional_data);
  }

  data_set_file.close();
}