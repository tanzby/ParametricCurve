#include <iostream>

#include "Curve.hpp"

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <opencv2/opencv.hpp>


void test_diy_class()
{
    struct _point
    {
        double x,y;
    };

    std::vector<_point> points{
        {45,105},
        {63,163},
        {95,205},
        {121,220},
        {146,256},
        {178,341}
    };

    bf::Curve curve;

    curve.fitCurve(points, 4.0);

    puts("\nfrom diy class");
    for(auto b: curve.beziers)
    {
        std::cout <<b.p1.x <<" "<< b.p1.y<< std::endl;
        std::cout <<b.c1.x <<" "<< b.c1.y<< std::endl;
        std::cout <<b.c2.x <<" "<< b.c2.y<< std::endl;
        std::cout <<b.p2.x <<" "<< b.p2.y<< std::endl;
        putchar('\n');
    }

}

void test_pcl_class()
{
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_ptr(new pcl::PointCloud<pcl::PointXYZI>);

    pcl::PointXYZI pointXYZI;
    pointXYZI.x = 45;
    pointXYZI.y = 105;
    cloud_ptr->push_back(pointXYZI);
    pointXYZI.x = 63;
    pointXYZI.y = 163;
    cloud_ptr->push_back(pointXYZI);
    pointXYZI.x = 95;
    pointXYZI.y = 205;
    cloud_ptr->push_back(pointXYZI);
    pointXYZI.x = 121;
    pointXYZI.y = 220;
    cloud_ptr->push_back(pointXYZI);
    pointXYZI.x = 146;
    pointXYZI.y = 256;
    cloud_ptr->push_back(pointXYZI);
    pointXYZI.x = 178;
    pointXYZI.y = 341;
    cloud_ptr->push_back(pointXYZI);

    bf::Curve curve;

    curve.fitCurve(cloud_ptr->points, 4.0);

    puts("\nfrom pcl class");
    for(auto b: curve.beziers)
    {
        std::cout <<b.p1.x <<" "<< b.p1.y<< std::endl;
        std::cout <<b.c1.x <<" "<< b.c1.y<< std::endl;
        std::cout <<b.c2.x <<" "<< b.c2.y<< std::endl;
        std::cout <<b.p2.x <<" "<< b.p2.y<< std::endl;
        putchar('\n');
    }
}

void test_opencv_class()
{
    std::vector<cv::Point2i> points{
        {45,105},
        {63,163},
        {95,205},
        {121,220},
        {146,256},
        {178,341}
    };

    bf::Curve curve;

    curve.fitCurve(points, 4.0);

    puts("\nfrom opencv class");
    for(auto b: curve.beziers)
    {
        std::cout <<b.p1.x <<" "<< b.p1.y<< std::endl;
        std::cout <<b.c1.x <<" "<< b.c1.y<< std::endl;
        std::cout <<b.c2.x <<" "<< b.c2.y<< std::endl;
        std::cout <<b.p2.x <<" "<< b.p2.y<< std::endl;
        putchar('\n');
    }
}

int main() {

    test_diy_class();
    test_pcl_class();
    test_opencv_class();

}