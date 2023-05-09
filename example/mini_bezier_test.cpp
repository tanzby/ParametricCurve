#include <iostream>
#include "bezier.hpp"

using namespace std;

int main() {
  // a line
  Bezier<1> line;

  using PointType = Bezier<1>::PointType;

  // line.set({{1,1},{2,2}});
  line.set({PointType{1}, PointType{1}, PointType{2}, PointType{2}});
  for (double i = 0; i < 1.0; i += 0.1) {
    cout << line.at(i).transpose() << ", " << line.at(i, 1).transpose() << endl;
  }

  // cubic bezier
  Bezier cubic_bezier{{1, 1}, {3, 1}, {4, 2}, {6, 3}};
  cout << "first order derivative in t=0.5 should be [4.5, 2.25]^T, the output "
          "is:\n"
       << cubic_bezier.at(0.5, 1) << endl;

  return 0;
}