#ifndef STLAPIWRITER_H
#define STLAPIWRITER_H

#include "occhandle.h"
#include <gp_Trsf.hxx>

class StlMesh_Mesh;
class TopoDS_Shape;
class TopoDS_Face;
class gp_Vec;
class gp_Pnt;
class Poly_Triangulation;

namespace STLAPIWriter
{

//! Tool to get triangles from triangulation taking into account face
//! orientation and location
class TriangleAccessor
{
public:

    TriangleAccessor (const TopoDS_Face& aFace);
    int NbTriangles () const { return myNbTriangles; }
    void GetTriangle (int iTri, gp_Vec &theNormal, gp_Pnt &thePnt1, gp_Pnt &thePnt2, gp_Pnt &thePnt3);

private:

    occHandle(Poly_Triangulation) myPoly;
    gp_Trsf myTrsf;
    int myNbTriangles;
    bool myInvert;
};

bool Write(const TopoDS_Shape& aShape, const char* aFileName);
bool WriteExtended(const TopoDS_Shape& aShape, const char *aFileName);

//! end namespace
}

#endif // STLAPIWRITER_H
