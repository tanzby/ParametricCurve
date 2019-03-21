#include <vector>
#include <array>
#include <cmath>

namespace bf {
    
    struct Vec2
    {
        double x, y;
        
        Vec2()
        {
            x = 0;
            y = 0;
        }
        Vec2(const double&x , const double&  y)
        {
            this->x = x;
            this->y = y;
        }
        
        template <typename Tp>
        Vec2(Tp& t)
        {
            this->x = t.x;
            this->y = t.y;
        }
        
        Vec2 operator - (const Vec2& _vec)
        {
            return {x - _vec.x, y - _vec.y};
        }
        
        Vec2 operator - (const double& m)
        {
            return {x - m, y - m};
        }
        
        Vec2 operator - (void)
        {
            return {-x, -y};
        }
        
        template <typename Tp>
        Vec2 operator - (const Tp& m)
        {
            return {x - m.x, y - m.y};
        }
        
        friend Vec2 operator-(const double& m, const Vec2& _vec)
        {
            return {m-_vec.x,m-_vec.y};
        }
        
        Vec2 operator + (const Vec2& _vec)
        {
            return {x + _vec.x, y + _vec.y};
        }
        
        Vec2 operator + (const double& m)
        {
            return {x + m, y + m};
        }
        
        friend Vec2 operator+(const double& m, const Vec2& _vec)
        {
            return {m+_vec.x,m+_vec.y};
        }
        
        Vec2 operator / (const double& m)
        {
            return {x/m, y/m};
        }
        
        Vec2 operator * (const double& m)
        {
            return {x*m, y*m};
        }
        
        Vec2 operator * (const double& m) const
        {
            return {x*m, y*m};
        }
        
        Vec2 operator * (const Vec2& _vec)
        {
            return {x*_vec.x, y* _vec.y};
        }
        
        friend Vec2 operator*(const double& m,  const Vec2& _vec)
        {
            return {m*_vec.x,m*_vec.y};
        }
        
        double dot(const Vec2& _vec)
        {
            return x*_vec.x+ y*_vec.y;
        }
        
        double norm(const double& n = 2)
        {
            return std::pow(std::pow(x,n) + std::pow(y,n), 1.0/n);
        }
        
        double sum()
        {
            return x+y;
        }
        
        Vec2 pow(const double& m)
        {
            return {std::pow(x,m), std::pow(y,m)};
        }
        
        Vec2 normalize()
        {
            double _m = this->norm();
            return {x/_m,y/_m};
        }
    };
    
    
    struct Bezier
    {
        Vec2 p1;
        Vec2 c1; // ctrl point 1
        Vec2 c2; // ctrl point 2
        Vec2 p2;
        
        Vec2 q(double t)
        {
            return std::pow(1.0-t,3.0) * p1 + 3.0*std::pow(1.0-t,2.0) * t * c1 + 3.0*(1.0-t)* std::pow(t, 2.0) * c2 + std::pow(t, 3.0) * p2;
        }
        
        Vec2 qprime(double t)
        {
            return 3.0*std::pow(1.0-t, 2.0) * (c1-p1) + 6.0*(1.0-t) * t * (c2-c1) + 3.0*std::pow(t, 2.0) * (p2-c2);
        }
        
        Vec2 qprimeprime(double t)
        {
            return 6.0*(1.0-t) * (c2-2.0*c1+p1) + 6.0*t* (p2-2.0*c2+c1);
        }
        
    };
    
    class Curve
    {
        double newtonRaphsonRootFind(Bezier& bezier, const Vec2& point, const double& u)
        {
            Vec2 d = bezier.q(u) - point;
            
            double numerator = (d * bezier.qprime(u)).sum();
            double denominator = (bezier.qprime(u).pow(2.0) + d * bezier.qprimeprime(u)).sum();
            
            return denominator == 0.0? u : (u - numerator/denominator);
        }
        
        template <typename T, typename E>
        std::vector<double> reparameterize(const std::vector<T, E>& points, const int& begin, const int& end, const std::vector<double>& parameter, Bezier& bezier)
        {
            std::vector<double> uPrime(end - begin + 1); //  New parameter values
            
            for (int i = begin; i <= end; i++) {
                uPrime[i - begin] = newtonRaphsonRootFind(bezier, points[i], parameter[i - begin]);
            }
            return uPrime;
        }
        
        template <typename T, typename E>
        Bezier generateBezier(const std::vector<T, E>& points, const int& begin, const int& end,
                              std::vector<double>& parameters, const Vec2& left_tangent, const Vec2& right_tangent)
        {
            std::vector<std::array<Vec2,2>> A(parameters.size());
            
            for(int i = 0; i < parameters.size(); ++i)
            {
                auto& u = parameters[i];
                A[i][0] = left_tangent  * 3.0 * std::pow(1.0 - u, 2.0) * u;
                A[i][1] = right_tangent * 3.0 * (1.0 - u) * u * u;
            }
            
            double C[2][2] = {0, 0, 0, 0};
            double X[2] = {0, 0};
            
            Vec2 tmp;
            
            for(int i = 0; i < parameters.size(); ++i)
            {
                C[0][0] += A[i][0].dot(A[i][0]);
                C[0][1] += A[i][0].dot(A[i][1]);
                C[1][0] = C[0][1];
                C[1][1] += A[i][1].dot(A[i][1]);
                
                Bezier b{points[begin],points[begin],points[end],points[end]};
                
                tmp = Vec2(points[i + begin]) - b.q(parameters[i]);
                
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
            // This is a common problem with this algoithm.
            
            double segLength = Vec2(points[begin].x-points[end].x,points[begin].y-points[end].y).norm();
            
            bool danger = false;
            if ((alpha_l > segLength * 2) || (alpha_r > segLength * 2)) {
                // begin  += 0;
                danger = true;
            }
            
            //  If alpha negative, use the Wu/Barsky heuristic (see text)
            //  (if alpha is 0, you get coincident control points that lead to
            //  divide by zero in any subsequent NewtonRaphsonRootFind() call.
            Bezier bez {points[begin], Vec2(), Vec2(), points[end]};
            double epsilon = 1.0e-6 * segLength;
            if (alpha_l < epsilon || alpha_r < epsilon || danger)
            {
                // fall back on standard (probably inaccurate) formula, and subdivide further if needed.
                bez.c1 = bez.p1 + left_tangent *  (segLength / 3.0);
                bez.c2 = bez.p2 + right_tangent * (segLength / 3.0);
            } else{
                bez.c1 = bez.p1 + left_tangent * alpha_l;
                bez.c2 = bez.p2 + right_tangent * alpha_r;
            }
            
            return bez;
        }
 
        template <typename T, typename E>
        std::pair<double, int> computeMaxError(const std::vector<T,E> &points, const int& begin, const int& end,
                                               Bezier& bezier,  const  std::vector<double>& parameters)
        {
            double max_dist = 0.0;
            int    split_point = (end - begin + 1) / 2;
            
            for(int i = begin+1; i < end; ++i)
            {
                // increase error weight.
                double dist = std::pow((bezier.q(parameters[i-begin])-Vec2(points[i])).norm(),2);
                if(dist>max_dist)
                {
                    max_dist = dist;
                    split_point = i;
                }
            }
            return {max_dist, split_point};
        }
        
        template <typename T, typename E>
        void fitCubic(const std::vector<T, E>& points,const int& begin,
                      const int& end,const Vec2& left_tangent,const Vec2& right_tangent, double max_error)
        {
            // Use heuristic if region only has two points in it
            if (end - begin + 1 == 2)
            {
                double dist = Vec2(points[begin].x-points[end].x,points[begin].y-points[end].y).norm()/3.0;
                Bezier _bezier {
                    points[begin],
                    Vec2(points[begin])+ left_tangent  * dist,
                    Vec2(points[end])  + right_tangent * dist,
                    points[end]
                };
                beziers.emplace_back(_bezier);
                return;
            }
            
            auto u = chordLengthParameterize(points, begin, end);
            auto bez = generateBezier(points, begin, end, u, left_tangent, right_tangent);
            auto res = computeMaxError(points,begin,end,bez,u);
            double error = res.first;
            int    split_idx = res.second;
            
            if (error < max_error)
            {
                beziers.emplace_back(bez);
                return;
            }
            
            //  If error not too large, try some reparameterization
            //  and iteration
            int max_iteration = 20;
            if (error < max_error * max_error)
            {
                for (int  i = 0; i < max_iteration; ++i) {
                    auto uPrime = reparameterize(points, begin, end, u, bez);
                    Bezier _bez = generateBezier(points, begin, end, uPrime, left_tangent, right_tangent);
                    auto _res = computeMaxError(points,begin,end,bez,u);
                    double _error = res.first;
                    int    _split_idx = res.second;
                    if (_error < error) {
                        beziers.emplace_back(_bez);
                        return;
                    }
                    u = uPrime;
                }
            }
            
            Vec2 center_tangent = (Vec2(points[split_idx-1])-Vec2(points[split_idx+1])).normalize();
            fitCubic(points, begin, split_idx, left_tangent, center_tangent, max_error);
            fitCubic(points, split_idx, end, -center_tangent, right_tangent, max_error);
            
        }
        
        template <typename T, typename E>
        std::vector<double> chordLengthParameterize(const std::vector<T, E>& points, int begin, int end)
        {
            std::vector<double> u(end - begin + 1);
            
            u[0] = 0.0;
            for (int i = begin + 1; i <= end; ++i)
            {
                u[i - begin] = u[i - begin - 1] + Vec2(points[i].x-points[i-1].x,points[i].y-points[i-1].y).norm();
            }
            
            for (int i = begin + 1; i <= end; ++i)
            {
                u[i - begin] = u[i - begin] / u[end - begin];
            }
            return u;
        }
        
    public:
        
        std::vector<Bezier> beziers;
        
        
        /**
         fit a set of 2D points to one or more Bezier Curve.

         @param points set of data, element from which must include class member x and y.
         @param max_error tolerance for the error that curves fit the point set.
         @return one or more Bezier model.
         */
        template <typename T, typename E> const std::vector<Bezier>& fitCurve(const std::vector<T,E>& points, const double& max_error)
        {
            beziers.clear();
            
            if (points.size() > 1)
            {
                Vec2 left_tangent  = (Vec2(points[1])-Vec2(points[0])).normalize();
                Vec2 right_tangent = (Vec2(points[points.size()-2])-Vec2(points[points.size()-1])).normalize();
                
                fitCubic(points, 0, points.size()-1, left_tangent, right_tangent, max_error);
            }
            return beziers;
        }
        
    };
    
}

