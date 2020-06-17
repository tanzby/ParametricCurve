#include <iostream>
#include <sstream>
#include <opencv2/opencv.hpp>
#include "bezier.hpp"

using namespace std;
using namespace cv;

class BezierFittingDemo
{
    string window_name = "test curve fitting";

    mutex points_mutex_;
    vector<Eigen::Vector2d> points;
    int max_error {100};
    PiecewiseBezierCurve<2> curve;
    Mat img {800, 800, CV_8UC3, Scalar(0)};

    void static onMouseCallback(int event, int x, int y, int flags, void* userdata)
    {
        static int prev_event = 0;
        if (prev_event == EVENT_LBUTTONDBLCLK)
        {
            if (event == EVENT_LBUTTONDOWN)      return;
            else if (event == EVENT_LBUTTONUP)   prev_event = 0;
            return;
        }

        BezierFittingDemo* demo = (BezierFittingDemo*)(userdata);

        switch (event)
        {
            case EVENT_LBUTTONUP:
            {
                lock_guard<mutex> lock(demo->points_mutex_);
                demo->points.emplace_back(x, y);
            }break;
            case EVENT_LBUTTONDBLCLK:
            {
                lock_guard<mutex> lock(demo->points_mutex_);
                if (demo->points.size()>1)
                {
                    demo->points.pop_back();
                    demo->points.pop_back();
                }
            }break;
        }
        prev_event = event;
    }

    void static onTrackerCallback(int pos, void* userdata)
    {
        BezierFittingDemo* demo = (BezierFittingDemo*)(userdata);
        demo->curve.fit(demo->points, max(double(demo->max_error),0.001));
        demo->curve.print(cout, "curve\n");
        cout << endl;
    }

public:

    void run()
    {
        namedWindow(window_name);
        setMouseCallback(window_name, &BezierFittingDemo::onMouseCallback, this);
        createTrackbar("max_error", window_name, &max_error,5000, &BezierFittingDemo::onTrackerCallback, this);

        int key = 0;
        while((key = waitKey(10)) != 'q')
        {
            lock_guard<mutex> lock(points_mutex_);
            // draw
            img = 0;
            for (auto& p: points)
            {
                circle(img,{int(p.x()),int(p.y())},4, {0,255,255},2);
            }

            if (key == 'f')
            {
                curve.fit(points, max_error);

                curve.print(cout, "curve\n");
                cout << endl;
            }

            for(const auto& b: curve.getPiecewiseBeziers())
            {
                vector<Point> draw_points;
                double t = 0.0, step = 0.01;
                while(t<=1.0)
                {
                    auto p = b.at(t);
                    draw_points.emplace_back(p[0], p[1]);
                    t+=step;
                }
                polylines(img, draw_points, false,
                              Scalar(int(b.length())*134214%255, int(b.length())*523622%255,int(b.length())*357326%255), 2, CV_AA);

            }
            printHelp();
            imshow(window_name, img);

        }
        destroyAllWindows();
    }
    
    void printHelp()
    {
        stringstream help;
        help << "\nThis is a demo for showing Bezier Curve Fitting\n"
             << "1. Click left button of mouse to add a new point\n"
             << "2. Double click the left mouse button quickly to cancel the previous point\n"
             << "3. Press 'f' to fit a piecewise bezier curve within the specified error threshold\n"
             << "4. Press 'q' to exit the application\n"
             << "5. you can drag the button bar to adjust the error rate\n";
        static once_flag flag;
        call_once(flag,[&](){cout << help.str();});
    }
};

int main()
{
    BezierFittingDemo demo;
    demo.printHelp();
    demo.run();
    return 0;
}