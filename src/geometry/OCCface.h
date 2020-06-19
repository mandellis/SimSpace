#ifndef OCCFACE_H
#define OCCFACE_H

#include "geometryface.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>

//! ----
//! C++
//! ----
#include <vector>
#include <map>

class OCCFace: public geometryFace
{
public:

    OCCFace(const TopoDS_Face &aFace = TopoDS_Face());
    OCCFace(const std::vector<TopoDS_Face> &listOfFaces);

    double area() override;
    bool pointProjection(double *aPoint, double *theProjectedPoint, double eps=0.0) override;

    void setGeometry(const TopoDS_Face &aFace);
    void setGeometry(const std::vector<TopoDS_Face> &listOfFaces);

    double deflection();

private:

    std::vector<TopoDS_Face> myFaces;
    double myMaxDeflection;
    std::vector<double> myDeflections;

private:

    void computeDeflection();
};

#endif // OCCFACE_H
