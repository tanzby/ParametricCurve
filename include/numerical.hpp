#pragma once

#include <complex>

struct NumericalQuadrature {
  template <typename F>
  static double simpson_1_3(F &&derivative_func, const double &L,
                            const double &R) {
    double mid = (L + R) / 2.0;
    return (derivative_func(L) + 4.0 * derivative_func(mid) +
            derivative_func(R)) *
           (R - L) / 6.0;
  }

  template <typename F>
  static double adaptive_simpson_1_3(F &&derivative_func, const double &L,
                                     const double &R,
                                     const double &eps = 0.0001) {
    const double mid = (L + R) / 2.0;
    double ST = simpson_1_3(derivative_func, L, R),
           SL = simpson_1_3(derivative_func, L, mid),
           SR = simpson_1_3(derivative_func, mid, R);
    double ans = SL + SR - ST;
    if (std::abs(ans) <= 15.0 * eps)
      return SL + SR + ans / 15.0;
    return adaptive_simpson_1_3(derivative_func, L, mid, eps / 2.0) +
           adaptive_simpson_1_3(derivative_func, mid, R, eps / 2.0);
  }

  template <typename F>
  static double simpson_3_8(F &&derivative_func, const double &L,
                            const double &R) {
    double mid_L = (2 * L + R) / 3.0, mid_R = (L + 2 * R) / 3.0;
    return (derivative_func(L) + 3.0 * derivative_func(mid_L) +
            3.0 * derivative_func(mid_R) + derivative_func(R)) *
           (R - L) / 8.0;
  }

  template <typename F>
  static double adaptive_simpson_3_8(F &&derivative_func, const double &L,
                                     const double &R,
                                     const double &eps = 0.0001) {
    const double mid = (L + R) / 2.0;
    double ST = simpson_3_8(derivative_func, L, R),
           SL = simpson_3_8(derivative_func, L, mid),
           SR = simpson_3_8(derivative_func, mid, R);
    double ans = SL + SR - ST;
    if (fabs(ans) <= 15.0 * eps)
      return SL + SR + ans / 15.0;
    return adaptive_simpson_3_8(derivative_func, L, mid, eps / 2.0) +
           adaptive_simpson_3_8(derivative_func, mid, R, eps / 2.0);
  }
};
