#include "cubic_spline.hpp"

auto CubicInterpolation::at(
    const double& t, const int& derivative_order ) const -> CubicInterpolation::PointType   {
  switch (derivative_order) {
    case 0: {
      PointType point(0.0, 0.0);
      double tpow = 1.0;
      for (const auto & coeff : coeffs_) {
        point += coeff * tpow;
        tpow *= t;
      }
      return point;
    }
    case 1:
      return coeffs_[1] + 2.0 * coeffs_[2] * t + 3.0 * coeffs_[3] * t * t;
    case 2:
      return 2.0 * coeffs_[2] + 6.0 * coeffs_[3] * t;
    case 3:
      return 6.0 * coeffs_[3];
    default:
      return PointType::Zero();
  }
}

void CubicInterpolation::interpolation(const PointType& p1,
                                       double theta1,
                                       const PointType& p2,
                                       double theta2) {
  auto dp = p2 - p1;

  double k = dp.norm();

  double kc1 = k * cos(theta1);
  double ks1 = k * sin(theta1);
  double kc2 = k * cos(theta2);
  double ks2 = k * sin(theta2);

  coeffs_[0] = p1;
  coeffs_[1] << kc1, ks1;
  coeffs_[2] << 3.0 * dp(0) - 2.0 * kc1 - kc2, 3.0 * dp(1) - 2.0 * ks1 - ks2;

  coeffs_[3] << -2.0 * dp(0) + kc1 + kc2, -2.0 * dp(1) + ks1 + ks2;

  computeLength();
}

void CubicInterpolation::print(std::ostream& out, const std::string& s ) const   {
  out << s;
  for (const auto& c : coeffs_) {
    out << c;
  }
  out << '\n';
}
