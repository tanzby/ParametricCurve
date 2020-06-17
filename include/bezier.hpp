#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "curve.hpp"

using namespace std;

template <int N = 3, int PointDim = 2>
class Bezier: public Curve<N, PointDim>
{
public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using PointType     = typename Curve<N, PointDim>::PointType;

private:

    using PointsType    = typename Curve<N, PointDim>::PointsType;

    PointsType derivative_precal_points_[N+1];
    int        derivative_precal_prefix_[N+1];

    /**
     * @brief 计算导数缓存
     */
    void derivativePrecompute()
    {
        derivative_precal_prefix_[0] = 1;
        derivative_precal_prefix_[1] = N;
        for(int k = 1; k <= N; ++k)
        {
            int I_range = N - k;
            for(int i = 0; i <= I_range; ++i)
            {
                derivative_precal_points_[k].col(i) = derivative_precal_points_[k-1].col(i+1) - derivative_precal_points_[k-1].col(i);
            }

            derivative_precal_prefix_[k] = derivative_precal_prefix_[k-1]*(N-k+1);
        }
    }

    /**
     * @brief 基于数值积分，计算曲线的长度
     */
    void computeLength() override
    {
        derivativePrecompute();
        Curve<N, PointDim>::computeLength();
    }

public:

    Bezier()
    {
        derivative_precal_points_[0] = PointsType::Zero();
    }

    Bezier(const std::initializer_list<PointType>& list)
    {
        set(std::forward<decltype(list)>(list));
    }

    template <typename T = array<double, PointDim>>
    Bezier(const std::initializer_list<T>& list)
    {
        set(std::forward<decltype(list)>(list));
    }

    void set(const std::initializer_list<PointType>& list)
    {
        int i = 0;
        for(const auto & it : list)
        {
            derivative_precal_points_[0].col(i) = it;
            if (++i > N) break;
        }
        computeLength();
    }

    template <typename T = array<double, PointDim>>
    void set(const std::initializer_list<PointType>& list)
    {
        int i = 0;
        for(const auto & it : list)
        {
            derivative_precal_points_[0].col(i) = it;
            if (++i > N) break;
        }
        computeLength();
    }

    PointType at(const double& t, const int& derivative_order = 0)  const override
    {
        double t_ = std::clamp(t,0.0,1.0);

        if (derivative_order > N) return PointType::Zero();

        /// De Casteljau’s Algorithm ///
        PointsType temp = derivative_precal_points_[derivative_order];

        int I_range =  N - derivative_order;
        for (int i = 0; i < I_range; ++i)
        {
            int J_range = I_range - i;
            for (int j = 0; j < J_range; ++j)
                temp.col(j) = (1.0 - t_) * temp.col(j) + t_ * temp.col(j + 1);
        }
        return derivative_precal_prefix_[derivative_order] * temp.col(0);
    }

    PointType param(const int& i) const
    {
        if (i > N)
        {
            throw out_of_range("Bezier::coeff: Current curve is in N=" + std::to_string(N) +
            " order, but the given i is more than N+1="+std::to_string(N));
        }
        return derivative_precal_points_[0].col(i);
    }

    void print(std::ostream &out, const string &s = "") const override
    {
        out << s << derivative_precal_points_[0] << '\n';
    }
};


template <int PointDim = 2>
class PiecewiseBezierCurve: public Curve<3, PointDim>
{
public:

    using BezierType = Bezier<3, PointDim>;
    using PointType  = typename Curve<3, PointDim>::PointType;

private:

    vector<BezierType> beziers_;
    vector<double>     length_table_;

    mutable const vector<PointType>* fitting_points_;

    std::vector<double> reParameterize(
            const int& begin, const int& end, const std::vector<double>& parameter, BezierType& bezier)
    {
        std::vector<double> uPrime(end - begin + 1); //  New parameter values

        for (int i = begin; i <= end; i++) {
            uPrime[i - begin] = bezier.findClosestParameter(fitting_points_->at(i), parameter[i - begin], 1);
        }
        return uPrime;
    }
    
    BezierType generateBezier(const int& begin, const int& end, std::vector<double>& parameters, 
            const PointType& left_tangent, const PointType& right_tangent)
    {
        std::vector<std::array<PointType, 2>> A(parameters.size());

        for(int i = 0; i < parameters.size(); ++i)
        {
            auto& u = parameters[i];
            A[i][0] = left_tangent  * 3.0 * (1.0 - u) * (1.0 - u) * u;
            A[i][1] = right_tangent * 3.0 * (1.0 - u) * u * u;
        }

        double C[2][2] = {0, 0, 0, 0};
        double X[2] = {0, 0};

        PointType tmp;

        for(int i = 0; i < parameters.size(); ++i)
        {
            C[0][0] += A[i][0].dot(A[i][0]);
            C[0][1] += A[i][0].dot(A[i][1]);
            C[1][0] =  C[0][1];
            C[1][1] += A[i][1].dot(A[i][1]);

            BezierType b{fitting_points_->at(begin),fitting_points_->at(begin),fitting_points_->at(end),fitting_points_->at(end)};

            tmp = fitting_points_->at(i + begin) - b.at(parameters[i]);

            X[0] += A[i][0].dot(tmp);
            X[1] += A[i][1].dot(tmp);
        }

        double det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
        double det_C0_X  = C[0][0] * X[1] - C[1][0] * X[0];
        double det_X_C1  = X[0] * C[1][1] - X[1] * C[0][1];

        double alpha_l =  det_C0_C1 == 0?0: det_X_C1 / det_C0_C1;
        double alpha_r =  det_C0_C1 == 0?0: det_C0_X / det_C0_C1;

        // Checks for "dangerous" points, meaning that the alpha_l or alpha_r are abnormally large
        // from here http://newsgroups.derkeiler.com/Archive/Comp/comp.graphics.algorithms/2005-08/msg00419.html
        // This is a common problem with this algorithm.

        double segLength = (fitting_points_->at(begin) - fitting_points_->at(end)).norm();

        bool danger = false;
        if ((alpha_l > segLength * 2) || (alpha_r > segLength * 2)) {
            // begin  += 0;
            danger = true;
        }

        //  If alpha negative, use the Wu/Barsky heuristic (see text)
        //  (if alpha is 0, you get coincident control points that lead to
        //  divide by zero in any subsequent NewtonRaphsonRootFind() call.
        PointType c1, c2;
        double epsilon = Curve<>::EPS * segLength;
        if (alpha_l < epsilon || alpha_r < epsilon || danger)
        {
            // fall back on standard (probably inaccurate) formula, and subdivide further if needed.
            c1 = fitting_points_->at(begin) + left_tangent  * (segLength / 3.0);
            c2 = fitting_points_->at(end)   + right_tangent * (segLength / 3.0);
        } else {
            c1 = fitting_points_->at(begin) + left_tangent  * alpha_l;
            c2 = fitting_points_->at(end)   + right_tangent * alpha_r;
        }

        return {fitting_points_->at(begin), c1, c2, fitting_points_->at(end)};
    }
    
    std::pair<double, int> computeMaxError(const int& begin, const int& end,
                                           BezierType & bezier,  const  std::vector<double>& parameters)
    {
        double max_dist = 0.0;
        int    split_point = (end - begin + 1) / 2;

        for(int i = begin+1; i < end; ++i)
        {
            // increase error weight.
            double dist = (bezier.at(parameters[i-begin]) - fitting_points_->at(i)).squaredNorm();
            if(dist>max_dist)
            {
                max_dist = dist;
                split_point = i;
            }
        }
        return {max_dist, split_point};
    }
    
    void fitCubic(const int& begin, const int& end,
            const PointType& left_tangent,const PointType& right_tangent, double max_error)
    {
        // Use heuristic if region only has two points in it
        if (end - begin + 1 == 2)
        {
            double dist = (fitting_points_->at(begin) - fitting_points_->at(end)).norm()/3.0;
            BezierType temp_bezier
            {
                fitting_points_->at(begin),
                fitting_points_->at(begin) + left_tangent  * dist,
                fitting_points_->at(end)   + right_tangent * dist,
                fitting_points_->at(end)
            };
            beziers_.emplace_back(temp_bezier);
            return;
        }

        auto u = chordLengthParameterize(begin, end);
        auto bez = generateBezier(begin, end, u, left_tangent, right_tangent);
        auto res = computeMaxError(begin,end,bez,u);
        double error = res.first;
        int    split_idx = res.second;

        if (error < max_error)
        {
            beziers_.emplace_back(bez);
            return;
        }

        //  If error not too large, try some reparameterization
        //  and iteration
        int max_iteration = 20;
        if (error < max_error * max_error)
        {
            for (int  i = 0; i < max_iteration; ++i)
            {
                auto uPrime = reParameterize(begin, end, u, bez);
                auto _bez = generateBezier(begin, end, uPrime, left_tangent, right_tangent);
                auto _res = computeMaxError(begin,end,bez,u);
                double _error = res.first;
                int    _split_idx = res.second;
                if (_error < error)
                {
                    beziers_.emplace_back(_bez);
                    return;
                }
                u = uPrime;
            }
        }

        PointType center_tangent = (fitting_points_->at(split_idx-1)-fitting_points_->at(split_idx+1)).normalized();
        fitCubic(begin, split_idx, left_tangent, center_tangent, max_error);
        fitCubic(split_idx, end, -center_tangent, right_tangent, max_error);
    }

    std::vector<double> chordLengthParameterize(int begin, int end)
    {
        std::vector<double> u(end - begin + 1);
        u[0] = 0.0;

        for (int i = begin + 1; i <= end; ++i)
        {
            u[i - begin] = u[i - begin - 1] + (fitting_points_->at(i)-fitting_points_->at(i-1)).norm();
        }

        for (int i = begin + 1; i <= end; ++i)
        {
            u[i - begin] = u[i - begin] / u[end - begin];
        }
        return u;
    }

protected:

    void computeLength() override
    {
        if (beziers_.empty()) return;

        for(BezierType& b: beziers_)
        {
            this->length_ +=  b.length();
        }

        length_table_.resize(beziers_.size());
        length_table_[0] = beziers_.front().length();

        for(int i = 1; i < beziers_.size(); ++i)
        {
            length_table_[i] = length_table_[i-1] + beziers_[i].length();
        }

        for (int i = 0; i < beziers_.size(); ++i)
        {
            length_table_[i] /= length_table_.back();
        }
    }

public:

    PointType at(const double &t, const int &derivative_order = 0) const override
    {
        if (beziers_.empty()) return {};
        if (t > 1.0) return {};

        auto bidx = std::upper_bound(length_table_.begin(), length_table_.end(), t);

        assert(bidx != length_table_.end()); // 理应任何时候不会等于end

        double t_s = bidx==length_table_.begin()?0:*(bidx-1);
        double t_b = (t-t_s) / (*bidx - t_s);

        return beziers_.at(bidx-length_table_.begin()).at(t_b, derivative_order);
    }

    void print(std::ostream &out, const string &s = "") const override
    {
        out << s;
        for(const auto& b: beziers_)
        {
            b.print(out);
        }
    }

    void setExternal(vector<BezierType> bezier_list)
    {
        fitting_points_ = nullptr;
        beziers_ = bezier_list;
        computeLength();
    }

    /**
     * @brief 给定有序的数据点，使用数值优化方法在给定误差范围内进行自动拟合与分段，拟合完成后，自动调用 computeLength()
     *        完成长度的计算。所有拟合结果被保存在PiecewiseBezierCurve结构中
     * @param points 有序的数据点
     * @param max_error 限定最大误差
     * @return
     */
    const std::vector<BezierType>& fit(const std::vector<PointType>& points, const double& max_error = 10.0)
    {
        beziers_.clear();

        if (points.size() > 1)
        {
            fitting_points_ = &points;
            auto left_tangent  = (fitting_points_->at(1) - fitting_points_->at(0)).normalized();
            auto right_tangent = (fitting_points_->at(points.size()-2) - fitting_points_->at(points.size()-1)).normalized();

            fitCubic(0, points.size()-1, left_tangent, right_tangent, max_error);
        }
        computeLength();
        return beziers_;
    }

    const std::vector<BezierType>& getPiecewiseBeziers() const
    {
        return beziers_;
    }

};