#pragma once

# include <Eigen/Sparse>
# include <Eigen/SparseLU>

#include "curve.hpp"

/**
 * @brief 一个简单的三次样条，只有两个数据点
 */
class CubicInterpolation: public Curve<3, 2>
{

public:

    using PointType = typename Curve<3, 2>::PointType;

private:

    PointType coeffs_[4];

public:

    CubicInterpolation()
    {
        for(int i = 0; i < 4; ++i) coeffs_[i].setZero();
    }

    PointType at(const double &t, const int &derivative_order = 0) const override
    {
        switch(derivative_order)
        {
            case 0:
            {
                PointType point(0.0, 0.0);
                double tpow = 1.0;
                for (int i = 0; i < 4; ++i)
                {
                    point += coeffs_[i] * tpow;
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

    /**
     * @brief 使用两个点，及起点与终点的方向，进行三次样条插值
     * @param p1 起点坐标
     * @param theta1 起点方向，值域为 $[0, 2 \pi]$
     * @param p2 终点坐标
     * @param theta2 终点方向，值域为 $[0, 2 \pi]$
     */
    void interpolation(const PointType& p1, double theta1,
             const PointType& p2, double theta2)
    {
        auto dp = p2 - p1;

        double k = dp.norm();

        double kc1 = k * cos(theta1);
        double ks1 = k * sin(theta1);
        double kc2 = k * cos(theta2);
        double ks2 = k * sin(theta2);

        coeffs_[0] = p1;
        coeffs_[1] << kc1, ks1;
        coeffs_[2] <<  3.0 * dp(0) - 2.0 * kc1 - kc2,
                       3.0 * dp(1) - 2.0 * ks1 - ks2;

        coeffs_[3] << -2.0 * dp(0) + kc1 + kc2,
                      -2.0 * dp(1) + ks1 + ks2;

        computeLength();
    }

    void print(std::ostream& out, const std::string& s = "") const override
    {
        out << s;
        for(const auto& c: coeffs_){
            out << c;
        }
        out <<'\n';
    }

};

/**
 * @brief 分段三次样条插值
 */
template <int PointDim = 2>
class PieceWiseCubicSpline: public Curve<3, PointDim>
{
public:

    using PointType = typename Curve<3, PointDim>::PointType;

    enum BoundaryCondition
    {
        FirstDeriv = 0,
        SecondDeriv
    };

protected:

    vector<double> t_;

    BoundaryCondition left_condition_, right_condition_;
    double  left_cond_value_, right_cond_value_;
    bool    force_linear_extrapolation_;

    struct _Spline
    {
        mutable const vector<double>*  t_ptr;
        Eigen::VectorXd y;
        Eigen::VectorXd a, b, c;
        double b0, c0;

        _Spline(){};

        void interp(PieceWiseCubicSpline<>& host, const vector<double>& t, const Eigen::VectorXd& y)
        {
            assert(t.size() == y.size());
            const int N = y.size();
            this->t_ptr = &t;
            this->y = y;

            // build Am=d
            Eigen::VectorXd d = Eigen::VectorXd::Zero(N);
            vector<Eigen::Triplet<double>> triplets;
            triplets.reserve(3*N);

            // build h
            for (int i = 1; i < N-1; ++i)
            {
                d[i] = (y[i+1]-y[i])/(t[i+1]-t[i]) - (y[i]-y[i-1])/(t[i]-t[i-1]);
                triplets.emplace_back(i, i-1,    (t[i]   - t[i-1]) /3.0);
                triplets.emplace_back(i, i,  2 * (t[i+1] - t[i-1]) /3.0);
                triplets.emplace_back(i, i+1,    (t[i+1] - t[i])   /3.0);
            }

            // boundary conditions
            if(host.left_condition_ == SecondDeriv)
            {
                // 2*b[0] = f''
                triplets.emplace_back(0, 0, 2.0);
                triplets.emplace_back(0, 1, 0.0);

                d[0] = host.left_cond_value_;
            }
            else if(host.left_condition_ == FirstDeriv)
            {
                // c[0] = f', needs to be re-expressed in terms of b:
                // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
                triplets.emplace_back(0, 0, 2.0*(t[1]-t[0]));
                triplets.emplace_back(0, 1, 1.0*(t[1]-t[0]));

                d[0]=3.0*((y[1]-y[0])/(t[1]-t[0])-host.left_cond_value_);
            } else{
                assert(false);
            }
            if(host.right_condition_ == SecondDeriv)
            {
                // 2*b[n-1] = f''
                triplets.emplace_back(N-1, N-1, 2.0);
                triplets.emplace_back(N-1, N-2, 0.0);

                d[N-1]=host.right_cond_value_;

            }
            else if(host.right_condition_ == FirstDeriv)
            {
                // c[n-1] = f', needs to be re-expressed in terms of b:
                // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
                // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
                triplets.emplace_back(N-1, N-1, 2.0*(t[N-1]-t[N-2]));
                triplets.emplace_back(N-1, N-2, 1.0*(t[N-1]-t[N-2]));

                d[N-1]=3.0*(host.right_cond_value_-(y[N-1]-y[N-2])/(t[N-1]-t[N-2]));
            } else{
                assert(false);
            }

            // solve the equation system to obtain the parameters b[]
            Eigen::SparseMatrix<double> A(N, N);
            A.setFromTriplets(triplets.begin(), triplets.end());

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            b = solver.solve(d);

            a.resize(N);
            c.resize(N);
            for(int i=0; i < N-1; i++)
            {
                a[i] = (b[i+1]-b[i])/(t[i+1]-t[i]) / 3.0;
                c[i] = (y[i+1]-y[i])/(t[i+1]-t[i]) - (2.0*b[i]+b[i+1])*(t[i+1]-t[i])/3.0;
            }

            // for left extrapolation coefficients
            b0 = (host.force_linear_extrapolation_==false) ? b[0] : 0.0;
            c0 = c[0];

            // for the right extrapolation coefficients
            // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y{n-1}
            double h = t[N-1]-t[N-2];
            // b[N-1] is determined by the boundary condition
            a[N-1]=0.0;
            c[N-1]=3.0*a[N-2]*h*h+2.0*b[N-2]*h+c[N-2];   // = f'_{n-2}(x_{n-1})
            if(host.force_linear_extrapolation_)
                b[N-1]=0.0;
        }

        double at(const int& idx, const double& t, const int derivative_order = 0) const
        {
            const vector<double>&  t_ = *t_ptr;
            size_t n   = t_.size();

            double h = t - t_[idx];
            double interpol = y [0];

            if(t<t_[0])
            {
                // extrapolation to the left
                switch(derivative_order)
                {
                    case 0:
                        interpol=(b0*h + c0)*h + y[0];
                        break;
                    case 1:
                        interpol=2.0*b0*h+c0;
                        break;
                    case 2:
                        interpol=2.0*b0*h;
                        break;
                    default:
                        interpol=0.0;
                        break;
                }
            }
            else if(t>t_[n-1])
            {
                // extrapolation to the right
                switch(derivative_order)
                {
                    case 0:
                        interpol=(b[n-1]*h + c[n-1])*h + y [n-1];
                        break;
                    case 1:
                        interpol=2.0*b[n-1]*h + c[n-1];
                        break;
                    case 2:
                        interpol=2.0*b[n-1];
                        break;
                    default:
                        interpol=0.0;
                        break;
                }
            }
            else
            {
                // interpolation
                switch(derivative_order)
                {
                    case 0:
                        interpol=((a[idx]*h + b[idx])*h + c[idx])*h + y[idx];
                        break;
                    case 1:
                        interpol=(3.0*a[idx]*h + 2.0*b[idx])*h + c[idx];
                        break;
                    case 2:
                        interpol=6.0*a[idx]*h + 2.0*b[idx];
                        break;
                    case 3:
                        interpol=6.0*a[idx];
                        break;
                    default:
                        interpol=0.0;
                        break;
                }
            }
            return interpol;
        }
    };

    _Spline spline_objects_[PointDim];

    void computeLength() override
    {

        int idx = 0;
        const auto df = [&](const double& t)->double
        {
            double acc = 0;
            for (int d = 0; d < PointDim; ++d)
            {
                const double temp = spline_objects_[d].at(idx, t, 1);
                acc += temp * temp;
            }
            return sqrt(acc);
        };

        for (idx = 0; idx < t_.size() - 1; ++idx)
        {
            this->length_ += NumericalQuadrature::adaptive_simpson_3_8(df, t_[idx], t_[idx+1]);
        }
//        auto df = [&](const double& t) -> double
//        {
//            return this->at(t, 1).norm();
//        };
//        this->length_ = NumericalQuadrature::adaptive_simpson_3_8(df, t_.front(), t_.back());

    }

public:

    PieceWiseCubicSpline()
    {
        left_cond_value_ = right_cond_value_ = 0;
        left_condition_ = right_condition_ = SecondDeriv;
        force_linear_extrapolation_ = false;
    }

    void setBoundary(BoundaryCondition left,  double left_value,
                     BoundaryCondition right, double right_value,
                     bool force_linear_extrapolation)
    {
        assert(t_.size()==0);          // set_points() must not have happened yet
        left_condition_=left;
        right_condition_=right;
        left_cond_value_=left_value;
        right_cond_value_=right_value;
        force_linear_extrapolation_=force_linear_extrapolation;
    }


    void interpolation(const vector<double>& t, const vector<PointType>& points)
    {
        assert(points.size()>2);

        const unsigned N = points.size();

        t_ = t;
        // TODO: maybe sort x and y, rather than returning an error
        for(int i=0; i<N-1; i++) {
            assert(t_[i]<t_[i+1]);
        }

        Eigen::Map<const Eigen::Matrix<double, PointDim, -1>> points_mat(points.data()->data(), PointDim, points.size());

        for(int d = 0; d < PointDim; ++d)
        {
            spline_objects_[d].interp(*this, t_, points_mat.row(d));
        }

        computeLength();
    }

    /**
     * @brief 给定$t$, 返回三次样条插值点
     * @param t 参数$t$，对于三次样条来说，$t$取值范围为插值数据所在范围
     * @param derivative_order 给定的阶次
     * @return 返回三次样条插值点，类型为PieceWiseCubicSpline::PointType
     */
    PointType at(const double &t, const int &derivative_order = 0) const override
    {
        auto begin = t_.data();
        auto it    = std::lower_bound(begin, begin + t_.size(), t);
        int idx = std::max(int(it - begin)-1, 0);

        PointType ret {};
        for(int d=0; d < PointDim; ++d)
        {
            ret[d] = this->spline_objects_[d].at(idx, t, derivative_order);
        }

        return ret;
    }

    void print(std::ostream& out, const std::string& s = "") const override
    {
        /*
        out << s;
        out << "t:\n";
        out << t_ << '\n';
        out << "y:\n";
        out << y_ << '\n';
        out << "a:\n";
        out << a_ << '\n';
        out << "b:\n";
        out << b_ << '\n';
        out << "c:\n";
        out << c_ << '\n';
        out << "b0:\n";
        out << b0_ << '\n';
        out << "c0:\n";
        out << c0_ << '\n';
        out << "boundary condition:\n";
        out << left_condition_ << " "<< left_cond_value_ << " " << right_condition_
        << " " << right_cond_value_ << " is_linear_extrap:" << force_linear_extrapolation_ <<'\n';
         */
    }

};