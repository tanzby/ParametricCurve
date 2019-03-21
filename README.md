# Bezier Fitting

Implement of Philip J. Schneider's "An Algorithm for Automatically Fitting Digitized Curves". Adaptively fit a set of 2D points with one or more cubic Bezier curves.



And the C++ and Python implementations refer to [fitCurves](https://github.com/volkerp/fitCurves) and the original C code.

- [fitCurves](https://github.com/volkerp/fitCurves) has a example gui application (demo.py) using Tkinter
- The original C code is available on <http://graphicsgems.org/> as well as in <https://github.com/erich666/GraphicsGems>



## usage



```c++
#include <iostream>
#include <vector>

#include "Curve.h"

using namespace std;

int main()
{
    std::vector<T> points; // class member of T must include x and y. 
    
    // do something such as setting data.
  
    bf::Curve curve;
    curve.fitCurve(points, 4.0);

    for(auto b: curve.beziers)
    {
        std::cout <<b.p1.x <<" "<< b.p1.y<< std::endl;
        std::cout <<b.c1.x <<" "<< b.c1.y<< std::endl;
        std::cout <<b.c2.x <<" "<< b.c2.y<< std::endl;
        std::cout <<b.p2.x <<" "<< b.p2.y<< std::endl;
        putchar('\n');
    }
}
```

