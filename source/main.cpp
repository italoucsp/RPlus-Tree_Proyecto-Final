#include <RPlusTree.hpp>

int main() {
  RPlus<double, 2, 3> demo;
  HyperPoint<double, 2> p1(array<double, 2>{1.2,5.3}), p2(array<double, 2>{4.9, 2.7});
  HyperRectangle<double, 2> rect1(p1,p2, true);
  HyperPoint<double, 2> p3(array<double, 2>{4.8, 25.2}), p4(array<double, 2>{8.8, 12.4});
  HyperRectangle<double, 2> rect2(p3, p4, true);
  HyperPoint<double, 2> p5(array<double, 2>{130.1, 42.51}), p6(array<double, 2>{8.7, 1.04});
  HyperRectangle<double, 2> rect3(p5, p6, true);
  HyperPoint<double, 2> p7(array<double, 2>{0.25, 30.2}), p8(array<double, 2>{75.5, 42.42});
  HyperRectangle<double, 2> rect4(p7, p8, true);
  vector<HyperRectangle<double, 2>> pd = {rect1, rect2, rect3, rect4};
  demo.insert(pd);
  demo.read_tree();
  return 0;
}