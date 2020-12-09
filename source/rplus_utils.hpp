#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <array>
#include <utility>
#include <stdarg.h>
#include <iomanip> 
#include <queue>
#include <tuple>
#include <math.h>

using namespace std;

#define ENTRY_GROUP vector<Entry>
#define ALERT(message) cout << message << endl;
#define ERROR_M_N_VALUES "ERROR: The number of dimensions and/or the number of entries per node should be greater than 1."
#define ERROR_FF_VALUE "ERROR: The value for fill factor should be between 1 and M."
#define ERROR_NODE_OFR "ERROR: The index is out of range in the node."
#define ERROR_NO_RECTDATA "ERROR: This rectangle is not a data container."

template<typename T, size_t N>
struct HyperPoint {
  HyperPoint();
  HyperPoint(array<T, N> data);
  bool operator<(const HyperPoint<T, N> &other) const;
  bool operator>(const HyperPoint<T, N> &other) const;
  HyperPoint<T, N>& operator=(const HyperPoint<T, N> &other);
  T& operator[](size_t index);
  T operator[](size_t index) const;
  void show_data();

private:
  array<T, N> multidata;
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
  for (size_t i(0); i < N; multidata[i] = other.multidata[i], ++i) {}
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
  HyperRectangle(HyperPoint<T, N> &A, HyperPoint<T, N> &B, bool is_data = false);
  HyperPoint<T, N>& operator[](size_t index);
  HyperRectangle<T, N>& operator=(const HyperRectangle<T, N>& other);
  void insert_data(HyperPoint<T, N> *point);
  bool overlaps(const HyperRectangle<T, N> &other);
  bool contains(const HyperPoint<T, N> &point);
  void adjust(const HyperRectangle<T, N> &other);
  pair<HyperPoint<T, N>, HyperPoint<T, N>> get_boundaries();
  double get_hypervolume();
  bool data_container();
  void show_data();

private:
  HyperPoint<T, N> bottom_left, top_right;
  vector<HyperPoint<T, N>*> data;
  bool is_data;
};

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle() {
  is_data = false;
}

template<typename T, size_t N>
HyperRectangle<T, N>::HyperRectangle(HyperPoint<T, N> &A, HyperPoint<T, N> &B, bool is_data) {
  for (size_t i(0); i < N; ++i) {
    bottom_left[i] = min(A[i], B[i]);
    top_right[i] = max(A[i], B[i]);
  }
  this->is_data = is_data;
  if (is_data) {
    data.resize(2);
    data[0] = &A;
    data[1] = &B;
  }
}

template<typename T, size_t N>
HyperPoint<T, N>& HyperRectangle<T, N>::operator[](size_t index) {
  return data[index];
}

template<typename T, size_t N>
HyperRectangle<T, N>& HyperRectangle<T, N>::operator=(const HyperRectangle<T, N>& other) {
  bottom_left = other.bottom_left;
  top_right = other.top_right;
  is_data = other.is_data;
  if (is_data) {
    for (size_t i(0); i < other.data.size(); ++i)
      data[i] = other.data[i];
  }
  return *this;
}

template<typename T, size_t N>
void HyperRectangle<T, N>::insert_data(HyperPoint<T, N> *point) {
  try {
    if (!is_data) {
      throw runtime_error(ERROR_NO_RECTDATA);
    }
    else {
      for (size_t i(0); i < N; ++i) {
        bottom_left[i] = min((*point)[i], bottom_left[i]);
        top_right[i] = max((*point)[i], top_right[i]);
      }
      data.push_back(point);
    }
  }
  catch (const exception &error) {
    ALERT(error.what());
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
pair<HyperPoint<T, N>, HyperPoint<T, N>> HyperRectangle<T, N>::get_boundaries() {
  return make_pair(bottom_left, top_right);
}

template<typename T, size_t N>
double HyperRectangle<T, N>::get_hypervolume() {
  double hypervolume = 1.0;
  for (size_t i(0); i < N; ++i) {
    hypervolume *= abs(top_right[i] - bottom_left[i]);
  }
  return hypervolume;
}

template<typename T, size_t N>
bool HyperRectangle<T, N>::data_container() {
  return is_data;
}

template<typename T, size_t N>
void HyperRectangle<T, N>::show_data() {
  cout << "\t\tHyperRectangle<" << N << "> :   bottom_left: "; bottom_left.show_data(); cout << "\t\t\t\t\ttop_right: "; top_right.show_data();
  if (is_data) {
    cout << "\t\t\tDATA:\n";
    for (HyperPoint<T, N>* &d : data) {
      cout << "\t\t\t";  d->show_data();
    }
  }
}