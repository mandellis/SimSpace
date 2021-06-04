#ifndef GEOMETRYFACE_H
#define GEOMETRYFACE_H
#include <iostream>
using namespace std;

class geometryFace
{
public:

    geometryFace();

    //! return the area of the face
    virtual double area() = 0;

    //! project the point aPoint onto the face, returning the coordinates
    //! of the projected point; return bool is the operation has been
    //! succesfull, false if not
    virtual bool pointProjection(double *aPoint, double *theProjectedPoint, double eps=0.0) = 0;

protected:

};

#endif // GEOMETRYFACE_H
