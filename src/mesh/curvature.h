#ifndef CURVATURE_H
#define CURVATURE_H

//! redefinition of opencascade Handle() macro
#include "occhandle.h"

//! custom includes
#include <ng_meshvs_datasourceface.h>

//! OCC
#include <TopoDS_Face.hxx>
#include <GeomLib_Tool.hxx>
#include <GeomLib.hxx>

//! Qt
#include <QMap>

class curvature
{
public:

    //! constructor
    curvature(const TopoDS_Face& aFace, double *P);

    //! constructor
    curvature();

    //! set the face
    void setFace(const TopoDS_Face &aFace){ myFace = aFace; }

    //! set the point
    void setPoint(double *P) {
        myPoint[0] = P[0];
        myPoint[1] = P[1];
        myPoint[2] = P[2];
    }

    //! type of curvature
    enum typeOfCurvature
    {
        typeOfCurvature_gaussian,
        typeOfCurvature_min,
        typeOfCurvature_max,
        typeOfCurvature_average
    };

    bool getCurvature(double &curvature, typeOfCurvature type);

    static bool getCurvature_static(const TopoDS_Face &aFace,
                                    double *aPoint,
                                    double &curvature,
                                    typeOfCurvature type,
                                    double *normal);

private:

    TopoDS_Face myFace;
    double myPoint[3];

    bool parameters(double *uv);
};

#endif // CURVATURE_H
