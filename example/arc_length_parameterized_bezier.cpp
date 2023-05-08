#include <iostream>
#include "opencv2/opencv.hpp"
#include "bezier.hpp"

using namespace cv;
using namespace std;

int main()
{
    // define a bezier
    Bezier bezier {{100,200},{400,100},{320,300},{600,400}};

    Mat img(500,700, CV_8UC3, Scalar(0));


    Eigen::Vector2d prev_point{0,0};
    vector<Point> vec;

    // draw curve
    double t = 0.0, step = 0.01;
    while(t<=1.0)
    {
        auto ps4 = bezier.atWithArcLengthParameterized(t,0,4);
        vec.emplace_back(int(ps4[0]),int(ps4[1]));
        t+=step;
    }
    polylines(img, vec, false, {0,0,255},1,LINE_AA);

    // draw equal arc distance point
    t=0, step = 0.25;
    while(t<=1.0)
    {
        auto ps4 = bezier.atWithArcLengthParameterized(t,0,4);
        circle(img,{int(ps4[0]),int(ps4[1])},4,{0,255,0},-1);
        if (prev_point != Eigen::Vector2d::Zero())
        {
            line(img, {int(prev_point[0]),int(prev_point[1])}, {int(ps4[0]),int(ps4[1])},{255,255,255},1, LINE_AA);
        }
        prev_point = ps4;
        t+=step;
    }

    for(auto& p: bezier.sampleWithArcLengthParameterized(step,false,4))
    {
        circle(img,{int(p[0]),int(p[1])},2,{255,0,0},-1);
    }

    imshow("arc-length parameterized of bezier curve", img);
    waitKey(0);

    return 0;
}

