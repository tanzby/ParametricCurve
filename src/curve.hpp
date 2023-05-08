#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "Eigen/Dense"

#include "numerical.hpp"

using namespace std;

/**
 * @brief 曲线基类，定义了必要的接口
 * @tparam N 曲线阶数
 * @tparam PointDim 数据点所在的维度
 */
template <int N = 3, int PointDim = 2>
class Curve {
 protected:
  using PointsType = Eigen::Matrix<double, PointDim, N + 1>;

  double length_{-1};

  /**
   * @brief 该方法计算 length_
   */
  virtual void computeLength() {
    auto df = [&](double t) -> double { return this->at(t, 1).norm(); };
    this->length_ = NumericalQuadrature::adaptive_simpson_3_8(df, 0, 1);
  }

  constexpr static const double EPS = 0.000001;

 public:
  typedef Eigen::Matrix<double, PointDim, 1> PointType;

  /**
   * @brief 返回该曲线的长度，而长度的信息应在 computeLength() 中计算完成
   * @return 派生类中曲线的长度
   */
  [[nodiscard]] virtual auto length() const -> const double { return length_; }

  /**
   * @brief 返回参数 $t$ 位置曲线上对应的点的位置，或$N$阶导数
   * @param t 参数位置，定义域为 $[0,1]$
   * @param derivative_order 求解的导数阶次，默认为$0$阶，即原函数值
   * @return $t$参数位置，$N$阶导数的值
   */
  virtual auto at(const double& t, const int& derivative_order = 0) const -> PointType = 0;

  /**
   * @brief 返回参数 $t$ 位置曲线上对应的点的位置，或$N$阶导数
   * @param t 参数位置，定义域为 $[0,1]$
   * @param derivative_order 求解的导数阶次，默认为$0$阶，即原函数值
   * @return $t$参数位置，$N$阶导数的值
   */
  virtual auto findClosestParameter(const PointType& point,
                                      double init_param,
                                      const int& max_iter = 20) const -> double {
    /**
     * 找到参数u使得离点p是最近的，即求解
     * f = (C(u)-p)*C'(u) = 0
     * f'= C'(u)*C'(u) + (C(u)-p)*C''(u)
     * 当给定u_{n+1} = u_{n} - f/f'
     */
    assert(max_iter > 0);
    double numerator;
    double denominator;
    for (int iter = 0; iter < max_iter; ++iter) {
      auto d = this->at(init_param) - point;
      auto first_order = this->at(init_param, 1);
      auto second_order = this->at(init_param, 2);

      numerator = d.dot(first_order);
      denominator =
          first_order.dot(first_order) + d.dot(second_order) + std::numeric_limits<double>::min();
      init_param = init_param - numerator / denominator;
    }

    return init_param;
  }

  /**
   * @brief 将曲线的参数内容传送到输出流对象中，默认必须末尾携带 `\n`
   * @param out 流对象，例如std::cout, std::stringsteam
   * @param s 设定的前缀
   */
  virtual void print(std::ostream& out, const std::string& s = "") const = 0;

  /**
   * @brief 近似地以$delta$间距均匀采集点
   * @param delta 参数间距或者距离间距，由$for_arc_length$指定
   * @param for_arc_length 如果为`true`，则$delta$以长度为单位，否则以$t\in[0,1]$参数为单位
   * @param max_iter_time 最大迭代次数
   * @return 参数化采样点数据
   */
  virtual auto sampleWithArcLengthParameterized(
      const double& delta,
      bool arc_length_base = true,
      const int& max_iter_time = 4,
      vector<double>* arc_length_t = nullptr) -> vector<PointType> {
    const double avg_distance = arc_length_base ? delta : this->length_ * delta;

    const double avg_t = arc_length_base ? avg_distance / this->length_ : delta;

    assert(avg_distance >= 0 && avg_distance <= this->length_ / 2.0);

    const int n(this->length_ / avg_distance);

    vector<double> t_array(n, 0);
    vector<double> dists(n, 0);
    vector<PointType> ret(n + 1);  // 实际上总共有 n+1 个点
    ret.front() = this->at(0.0);
    ret.back() = this->at(1.0);

    // 第一个点与最后一个点保持在开始与结尾处
    for (int i = 1; i < n; ++i) {
      t_array[i] = i * avg_t;
      ret[i] = this->at(t_array[i]);
    }
    t_array.emplace_back(1.0);

    double prev_offset = -1;

    for (int iter = 0; iter < max_iter_time; ++iter) {
      // 1. 计算上一次迭代确定的 t 参数下，每一个分段的近似长度
      for (int j = 1; j < n; j++) { dists[j] = (ret[j] - ret[j - 1]).norm();
}

      double offset = 0;
      for (int j = 1; j < n; j++) {
        // 2. 累计近似弧长并计算误差
        offset += dists[j] - avg_distance;

        // 3. Newton's method
        double first_order = this->at(t_array[j], 1).norm();
        double second_order = this->at(t_array[j], 2).norm();
        double numerator = offset * first_order;
        double denominator = offset * second_order + first_order * first_order;

        t_array[j] = t_array[j] - numerator / denominator;

        ret[j] = this->at(t_array[j]);
      }

      if (offset < EPS || abs(offset - prev_offset) < EPS) {
        break;
      }
      prev_offset = offset;
    }

    if (arc_length_t != nullptr) {
      arc_length_t->swap(t_array);
    }

    return ret;
  }

  /**
   * @brief 在参数$t\in[0, 1]$取与弧长线性相关点
   * @param t 弧长参数化下的参数值，$t\in [0, 1]$
   * @param arc_length_t 返回最终计算得到的参数位置
   * @param derivative_order 给定的阶次, 默认为0阶
   * @param max_iter_time 最大迭代次数
   * @return 参数化点位置
   */
  virtual auto atWithArcLengthParameterized(const double& t,
                                                 const int& derivative_order = 0,
                                                 const int& max_iter_time = 4,
                                                 double* arc_length_t = nullptr) const -> PointType {
    assert(t >= 0.0 && t <= 1.0);

    double approx_t = t;
    double target_length = t * this->length_;
    double prev_approx_t = approx_t;

    const auto df = [&](double t) -> double { return this->at(t, 1).norm(); };

    for (int iter = 0; iter < max_iter_time; ++iter) {
      double approx_length = NumericalQuadrature::adaptive_simpson_3_8(df, 0, approx_t);
      double d = approx_length - target_length;
      if (abs(d) < EPS) { break;
}

      // Newton's method
      double first_order = this->at(approx_t, 1).norm();
      double second_order = this->at(approx_t, 2).norm();
      double numerator = d * first_order;
      double denominator = d * second_order + first_order * first_order;

      approx_t = approx_t - numerator / denominator;

      if (abs(approx_t - prev_approx_t) < EPS) { {
        break;
      } } else { {
        prev_approx_t = approx_t;
}
}
    }
    if (arc_length_t != nullptr) {
      *arc_length_t = approx_t;
    }

    return this->at(approx_t, derivative_order);
  }
};
