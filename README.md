# ParametricCurve

ParametricCurve implements two kinds of curve, including

* **N order Bezier Curve**
* **Cubic PieceWise Spline**

Both of them inherit from `Curve`, which provides

* Compute the n-order derivative at a point
* Compute the parameter position (t) closest to a point on the curve
* Sample points with equal arc distance by numerical method

For more details, read the codes of examples, a simple description about N-orders Bezier Curve can be found in [https://zhuanlan.zhihu.com/p/130247362](https://zhuanlan.zhihu.com/p/130247362) 

## usage

ParametricCurve only relies on `Eigen`, just copy the `include` files to your project.

## reference

* https://github.com/volkerp/fitCurves
* https://kluge.in-chemnitz.de/opensource/spline/


