#include <iostream>
#include <memory>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <array>
#include <utility>

using namespace std;

#define ENTRY_GROUP vector<Entry>
#define RECT_GROUP vector<HyperRectangle<size_t(N_)>>
#define ALERT(message) cout << message << endl;
#define ERROR_M_N_VALUES "ERROR: El n\243mero de dimensiones y el n\243mero de espacios por nodo deben ser mayor a 1."
#define ERROR_BND_VALUES "ERROR: El valor de maxBound dado debe ser mayor al de minBound."

template<size_t N>
struct HyperPoint {
  HyperPoint();
  explicit HyperPoint(const array<double, N> &dimensions_);
  bool operator<=(const HyperPoint<N> &hp2) const;

  array<double, N> dimensions;
};

template<size_t N>
HyperPoint<N>::HyperPoint() {
  for (auto &d : dimensions) d = 0.0;
}

template<size_t N>
HyperPoint<N>::HyperPoint(const array<double, N> &dimensions_) {
  dimensions = dimensions_;
}

template<size_t N>
bool HyperPoint<N>::operator<=(const HyperPoint<N> &hp2) const {
  for (size_t i = 0; i < N; ++i) {
    if (this->dimensions[i] > hp2->dimensions[i]) {
      return false;
    }
  }
  return true;
}

template<size_t N>
struct HyperRectangle {
  HyperRectangle();
  HyperRectangle(HyperPoint<N> minBound, HyperPoint<N> maxBound);
  bool isOverlaping(const HyperRectangle &other);
  bool contains(const HyperPoint<N> &point);

private:
  HyperPoint<N> bottom_left, top_right;
};

template<size_t N>
HyperRectangle<N>::HyperRectangle() {}

template<size_t N>
HyperRectangle<N>::HyperRectangle(HyperPoint<N> minBound, HyperPoint<N> maxBound) {
  try {
    if (maxBound <= minBound) {
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
bool HyperRectangle<N>::isOverlaping(const HyperRectangle &other) {
  for (size_t i = 0; i < N; ++i) {
    if (!(bottom_left.dimensions[i] <= other.top_right.dimensions[i] &&
          top_right.dimensions[i] >= other.bottom_left.dimensions[i])) {
      return false;
    }
  }
  return true;
}
