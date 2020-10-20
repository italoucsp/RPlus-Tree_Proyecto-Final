#include <iostream>
#include <memory>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <utility>

using namespace std;

#define ENTRY_GROUP vector<Entry>
#define RECT_GROUP vector<HyperRectangle<N_>>
#define ALERT(message) cout << message << endl;
#define ERROR_M_N_VALUES "ERROR: El n\243mero de dimensiones y el n\243mero de espacios por nodo deben ser mayor a 1."


template<int N>
struct HyperRectangle {
  HyperRectangle();
  bool isOverlaping(const HyperRectangle &other);
};

template<int N>
HyperRectangle<N>::HyperRectangle() {

}

template<int N>
bool HyperRectangle<N>::isOverlaping(const HyperRectangle &other) {
  return true;
}
