#ifndef TETQUALITYCLASS_H
#define TETQUALITYCLASS_H

#include <Eigen/Dense>

#include <mesh.h>
#include <polygon.h>

#include <iostream>
using namespace std;

class TetQualityClass
{
public:

    TetQualityClass();
    TetQualityClass(const std::vector<polygon::Point> &aTet);
    TetQualityClass(const std::vector<mesh::meshPoint> &aTet);

    TetQualityClass(double ax0, double ay0, double az0,
                    double ax1, double ay1, double az1,
                    double ax2, double ay2, double az2,
                    double ax3, double ay3, double az3);

    //! set points
    void setPoints(const std::vector<mesh::meshPoint> &aTet);

private:

    //! ------------------------------
    //! definition of the tetrahedron
    //! ------------------------------
    double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;

    double h(double sigma);

public:

    double Volume();
    void invertTetForCalculation();
    void getQualityMeasure(double &q0, double &q1, double &q2, double &V);
};

#endif // TETQUALITYCLASS_H
