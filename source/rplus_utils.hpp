#include <iostream>
#include <memory>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <array>
#include <utility>
#include <stdarg.h>
#include <iomanip> 
#include <queue>
#include <math.h>

using namespace std;

#define ENTRY_GROUP vector<Entry>
#define RECT_GROUP vector<HyperRectangle<size_t(N_)>>
#define ALERT(message) cout << message << endl;
#define ERROR_M_N_VALUES "ERROR: El n\243mero de dimensiones y el n\243mero de espacios por nodo deben ser mayor a 1."
#define ERROR_BND_VALUES "ERROR: El valor de maxBound dado debe ser mayor al de minBound."
#define ERROR_NEQ_DIMS "ERROR: Las dimensiones de ambos puntos deben ser iguales."

struct HyperPoint {
  HyperPoint();
  explicit HyperPoint(size_t N, ...);

  void reset(size_t N);
  bool operator<(const HyperPoint &other) const;
  HyperPoint& operator=(const HyperPoint& other);
  double& operator[](size_t index);
  void show_data();

private:
  vector<double> dimensions;
};

HyperPoint::HyperPoint() {
  dimensions.resize(size_t(2));
  for (auto &d : dimensions) d = 0.0;
}

HyperPoint::HyperPoint(size_t N, ...) {
  va_list list;
  va_start(list, N);
  dimensions.resize(N);
  for (double &d : dimensions) { d = va_arg(list, double); }
  va_end(list);
}

void HyperPoint::reset(size_t N) {
  dimensions.resize(N);
  for (size_t i = 0; i < N; dimensions[i] = 0.0, ++i) {}
}

bool HyperPoint::operator<(const HyperPoint &other) const{
  try {
    if (this->dimensions.size() != other.dimensions.size()) {
      throw runtime_error(ERROR_NEQ_DIMS);
    }
    else {
      for (size_t i = 0; i < dimensions.size(); ++i) {
        if (this->dimensions[i] >= other.dimensions[i]) {
          return false;
        }
      }
      return true;
    }
  }
  catch (const exception &error) {
    ALERT(error.what());
  }
  return false;
}

HyperPoint& HyperPoint::operator=(const HyperPoint& other) {
  this->dimensions.resize(other.dimensions.size());
  for (size_t i = 0; i < dimensions.size(); this->dimensions[i] = other.dimensions[i], ++i) {}
  return *this;
}

double& HyperPoint::operator[](size_t index) {
  return dimensions[index];
}

void HyperPoint::show_data() {
  cout << '['; for (auto d : dimensions) {
    cout << setprecision(10) << d << ',';
  }cout << "//]" << endl;
}

template<size_t N>
struct HyperRectangle {
  HyperRectangle();
  HyperRectangle(HyperPoint minBound, HyperPoint maxBound);
  bool isOverlaping(const HyperRectangle<N> &other);
  bool contains(const HyperPoint &point);
  void adjust_with_hrect(const HyperRectangle<N> &other);
  void adjust_with_hpoint(const HyperPoint &point);
  pair<HyperPoint, HyperPoint> get_boundaries();
  double get_hypervolume();

private:
  HyperPoint bottom_left, top_right;
};

template<size_t N>
HyperRectangle<N>::HyperRectangle() {
  bottom_left.reset(N);
  top_right.reset(N);
}

template<size_t N>
HyperRectangle<N>::HyperRectangle(HyperPoint minBound, HyperPoint maxBound) {
  try {
    if (maxBound < minBound) {
      throw runtime_error(ERROR_BND_VALUES);
    }
    else {
      bottom_left = minBound;
      top_right = maxBound;
    }
  }
  catch (const exception &error) {
    ALERT(error.what());
  }
}

template<size_t N>
bool HyperRectangle<N>::isOverlaping(const HyperRectangle<N> &other) {
  for (size_t i = 0; i < N; ++i) {
    if (!(bottom_left[i] <= other.top_right[i] &&
          top_right[i] >= other.bottom_left[i])) {
      return false;
    }
  }
  return true;
}

template<size_t N>
bool HyperRectangle<N>::contains(const HyperPoint &point) {

}

template<size_t N>
void HyperRectangle<N>::adjust_with_hrect(const HyperRectangle<N> &other) {
  for (size_t i = 0; i < N; ++i) {
    bottom_left[i] = min(other.bottom_left[i], bottom_left[i]);
    top_right[i] = max(other.top_right[i], top_right[i]);
  }
}

template<size_t N>
void HyperRectangle<N>::adjust_with_hpoint(const HyperPoint &point) {

}

template<size_t N>
pair<HyperPoint, HyperPoint> HyperRectangle<N>::get_boundaries() {
  return make_pair(bottom_left, top_right);
}

template<size_t N>
double HyperRectangle<N>::get_hypervolume() {
  double hypervolume = 1.0;
  for (size_t i = 0; i < N; ++i) {
    hypervolume *= abs(top_right[i] - bottom_left[i]);
  }
  return hypervolume;
}