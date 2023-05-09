#include <iostream>
#include <opencv2/opencv.hpp>
#include "cubic_spline.hpp"

using namespace std;
using namespace cv;

int main() {
  //    vector<double> param {60, 90, 170, 230, 250};
  //    vector<Eigen::Vector2d> points{
  //        {60,160},{90,220},{170,210},{230,260},{250,240}
  //    };

  vector<double> param{60, 90, 170, 230, 250};
  vector<Eigen::Vector2d> points{{60, 80}, {90, 100}, {170, 180}, {230, 240}, {250, 240}};

  PieceWiseCubicSpline s;
  s.interpolation(param, points);

  int width = 300;
  Mat img(width, width, CV_8UC3, Scalar(0));

  // draw curve
  vector<Point> vec;
  double t = 0, step = 1, end = width;
  while (t <= end) {
    auto p = s.at(t);
    vec.emplace_back(int(p[0]), int(p[1]));
    t += step;
  }
  polylines(img, vec, false, {0, 0, 255}, 1, LINE_AA);

  // draw fitted point
  for (auto p : points) circle(img, {int(p[0]), int(p[1])}, 4, {255, 255, 255}, -1);

  // draw derivative
  t = 0, step = 50, end = width;
  while (t <= end) {
    Eigen::Vector2d p = s.at(t, 0);
    Eigen::Vector2d dir = s.at(t, 1) * 50;
    Eigen::Vector2d v = p + dir;
    cout << s.at(t, 1) << endl;
    cv::arrowedLine(img, {int(p.x()), int(p.y())}, {int(v.x()), int(v.y())}, {255, 255, 255});
    t += step;
  }

  // compute length
  double acc_length = 0;
  for (int i = 1; i < vec.size(); ++i) {
    auto p = (vec[i] - vec[i - 1]);
    acc_length += sqrt(p.x * p.x + p.y * p.y);
  }

  s.print(cout);

  cout << "simpson's rule:\t" << s.length() << endl;
  cout << "naive method:   " << acc_length << endl;

  /*
  // find closest point
  double guss = 170;
  cv::Point2i target {200,200};
  guss = s.findClosestParameter({target.x,target.y}, guss);
  cv::line(img, target, {int(guss),int(s.at(guss).y())}, {255,255,255},2);
  */

  imshow("piece wise cubic spline", img);
  waitKey(0);

  return 0;
}
