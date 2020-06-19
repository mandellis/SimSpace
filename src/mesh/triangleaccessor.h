#ifndef TRIANGLEACCESSOR_H
#define TRIANGLEACCESSOR_H

//! ------------------------------------------------------------------
//! Tool to get triangles from triangulation taking into account face
//! orientation and location
//! ------------------------------------------------------------------


//! --------------------------------------
//! replacement for opencascade::handle<>
//! --------------------------------------
#include "occhandle.h"

//! ----
//! OCC
//! ----
#include <TopoDS_Face.hxx>
#include <BRep_Tool.hxx>
#include <TopLoc_Location.hxx>
#include <TopAbs_Orientation.hxx>
#include <Poly_Triangulation.hxx>
#include <gp_Trsf.hxx>
#include <gp_Pnt.hxx>
#include <gp.hxx>

class TriangleAccessor
{
private:

    occHandle(Poly_Triangulation) myPoly;
    gp_Trsf myTrsf;
    int myNbTriangles;
    bool myInvert;

public:

    //! ------------
    //! constructor
    //! ------------
    TriangleAccessor (const TopoDS_Face& aFace);

    //! -------------------------------
    //! return the number of triangles
    //! -------------------------------
    int NbTriangles();

    //! -----------------------------------------
    //! get the i-th triangle and outward normal
    //! -----------------------------------------
    void GetTriangle (int iTri, gp_Vec &theNormal, gp_Pnt &thePnt1, gp_Pnt &thePnt2, gp_Pnt &thePnt3);
};

#endif // TRIANGLEACCESSOR_H
