#include "bezier.hpp"
#include <iomanip>
#include <iostream>

using namespace std;

int main() {
  using namespace std::placeholders;

  // define an Bezier
  Bezier bezier{{1, 1}, {3, 1}, {4, 2}, {6, 3}};

  // as $L = \int_{0}^{1}|f'(t)|dt$, the norm of the first derivative of the
  // curve is defined.
  auto df_norm = [&](double t) -> double { return bezier.at(t, 1).norm(); };

  cout << setprecision(16);

  // naive method 1
  double t = 0, step = 0.5, length = 0;
  auto prev = bezier.at(0);
  while (t <= 1) {
    auto p = bezier.at(t);
    length += (p - prev).norm();
    prev = p;
    t += step;
  }
  cout << "naive method with step=0.5:     \t" << length << endl;

  // naive method 2
  t = 0, step = 0.001, length = 0;
  prev = bezier.at(0);
  while (t <= 1) {
    auto p = bezier.at(t);
    length += (p - prev).norm();
    prev = p;
    t += step;
  }
  cout << "naive method with step=0.001:   \t" << length << endl;

  cout << "result by using simpson's rule: \t" << bezier.length() << endl;

  return 0;
}
