#include <RPlusTree.hpp>

int main() {
  RPlus<string, 5, 8> demo;
  HyperPoint point(5,12.0,1987.0,15.0,26.0,12548.0256);
  point.show_data();
  point.reset(14);
  point.show_data();
  return 0;
}