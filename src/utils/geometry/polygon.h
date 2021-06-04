#ifndef POLYGON_H
#define POLYGON_H

#define SIGN(x) ((0<x)-(x<0))
#define XOR(l1, l2)  (l1 && (!l2)) || ((!l1) && l2)

#include <vector>
#include <iostream>
using namespace std;

// https://docs.microsoft.com/en-us/cpp/cpp/
// namespaces-cpp?view=msvc-160#:~:text=A%20namespace%20is%20a%20declarative,code%20base%20includes%20multiple%20libraries.

namespace polygon
{
//! a Point in the 2D space
struct Point2D
{
    double x,y;
    Point2D(double x0 = 0, double y0 = 0):x(x0),y(y0){;}
    Point2D(const Point2D &rhs){ x=rhs.x; y = rhs.y; }
    bool operator == (const Point2D &rhs) { if(x==rhs.x && y==rhs.y) return true; return false; }
    Point2D operator = (const Point2D &rhs) {x=rhs.x; y = rhs.y; return *this; }
};

//! a Point in the 3D space
struct Point
{
    double x, y, z;
    Point(double x0=0, double y0=0, double z0=0):x(x0),y(y0),z(z0){;}
    Point(const Point &other) {x = other.x; y = other.y; z = other.z; }
    Point operator = (const Point &other) { x = other.x; y = other.y; z = other.z; return *this; }
    bool operator == (const Point &other) { if(x == other.x && y == other.y && z == other.z) return true; return false; }
};

//! triangleArea
double triangleArea(const polygon::Point &P0, const polygon::Point &P1, const polygon::Point &P2);

//! areCollinear
bool areCollinear(const polygon::Point &P0, const polygon::Point &P1, const polygon::Point &P2);

//! area3D_Polygon
double area3D_Polygon(const std::vector<polygon::Point> &V);

//! getNormal
std::vector<double> getNormal(const std::vector<polygon::Point> &points);

//! isPointInside
bool isPointInside(const polygon::Point &P,
                   const std::vector<polygon::Point> &V,
                   bool normalizeCrossProducts,
                   double tolerance = std::numeric_limits<double>::epsilon());

//! planeCoefficients
bool planeCoefficients(const std::vector<polygon::Point> &points, double &a, double &b, double &c, double &d);

//! pointPlaneDistance
double pointPlaneDistance(const polygon::Point &P, double a, double b,double c, double d);

//! getRotationMatrix
bool getRotationMatrix(const std::vector<polygon::Point> &aPolygon,
                       std::vector<double> &rotationMatrix);

//! getPolygonCenter
bool getPolygonCenter(const std::vector<polygon::Point> &aPolygon,
                      polygon::Point &faceCenter, double &faceArea);

//! testPolygonPlaneIntersection
bool testPolygonPlaneIntersection(const std::vector<polygon::Point> &aPolygon, double a, double b, double c, double d/*, int &NbPoints*/);

//! polygon_triangulate
bool polygon_triangulate(const std::vector<polygon::Point> &aPolygon, std::vector<std::vector<polygon::Point> > &triangles);

//! polygon_triangulate_planar
bool polygon_triangulate_planar(const std::vector<double> &x, const std::vector<double> &y, std::vector<int> &triangles);

//! triangle area
double triangle_area (double xa, double ya, double xb, double yb, double xc, double yc);

//! diagonal
bool diagonal (int im1, int ip1, const std::vector<int> &prev_node, const std::vector<int> &next_node,
               const std::vector<double> &x, const std::vector<double> &y);

//! in_cone
bool in_cone (int im1, int ip1, const std::vector<int> &prev_node, const std::vector<int> &next_node,
              const std::vector<double> &x, const std::vector<double> &y);

//! diagonalie
bool diagonalie (int im1, int ip1, const std::vector<int> &next_node,
                 const std::vector<double> &x, const std::vector<double> &y);

//! intersect
bool intersect(double xa, double ya, double xb, double yb, double xc,
               double yc, double xd, double yd);

//! collinear
bool collinear(double xa, double ya, double xb, double yb, double xc, double yc);

//! between
bool between (double xa, double ya, double xb, double yb, double xc, double yc);

//! intersect_prop
bool intersect_prop(double xa, double ya, double xb, double yb, double xc, double yc, double xd, double yd);

// end namespace
}

#endif // POLYGON_H
