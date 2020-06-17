#include <iostream>
#include <opencv2/opencv.hpp>
#include "cubic_spline.hpp"

using namespace cv;
using namespace std;

int main()
{
    // fit a cubic spline
    CubicInterpolation cubic;
    cubic.interpolation({100,100},0,{300,300},0);

    int width = 400;
    Mat img(width, width, CV_8UC3, Scalar(0));

    vector<Point> vec;
    double t = 0.0, step = 0.01;
    while(t <= 1.0)
    {
        auto p = cubic.at(t);
        vec.emplace_back(int(p[0]),int(p[1]));
        t += step;
    }
    polylines(img, vec, false, {0,0,255},1, LINE_AA);
    double acc_length = 0;
    for(int i = 1;i<vec.size();++i)
    {
        auto p = (vec[i]-vec[i-1]);
        acc_length += sqrt(p.x*p.x+p.y*p.y);
    }
    cout  << "naive length: " << acc_length << endl;
    cout  << "length: " << cubic.length() << endl;

    // Uniform sampling
    for(auto p: cubic.sampleWithArcLengthParameterized(20))
    {
        cv::circle(img, {int(p.x()),int(p.y())},4,{255,255,255},1);
    }

    imshow("Cubic Interpolation", img);
    waitKey(0);

    return 0;
}
