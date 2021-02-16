#ifndef GEOMTOOLSCLASS_H
#define GEOMTOOLSCLASS_H

//! ----
//! OCC
//! ----
#include <BRep_Tool.hxx>

#include <GeomAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Geom_Plane.hxx>

#include <GeomAdaptor_Curve.hxx>
#include <GeomAbs_CurveType.hxx>
#include <Geom_Curve.hxx>

#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax3.hxx>

#include <TopTools_ListOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Compound.hxx>

//! ----------------
//! custom includes
//! ----------------
#include "mydefines.h"
#include <geometrydatabase.h>

//! ----
//! C++
//! ----
#include <vector>

class TopoDS_Face;
class TopoDS_Edge;

namespace GeomToolsClass
{
// begin namespace

void getBoundingBox(const TopoDS_Shape &shape, double &L1, double &L2, double &L3);
void getFaceType(const TopoDS_Face &aFace, GeomAbs_SurfaceType &type);
void getCurveType(const TopoDS_Edge &anEdge, GeomAbs_CurveType &type);
bool getPlanarFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &centroid);
bool getCylindricalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P);
bool getSphericalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P);
bool getConicalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P);
bool getEdgeInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P);
bool getCircleInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P);
bool getEllipseInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P);
std::vector<double> calculateCentroid(geometryDataBase *gDB, const std::vector<GeometryTag> &vecLoc);

//! test if three points are aligned
bool isCollinear(double *A0, double *A1, double *A2);

//! rotation matrix for 3D->plane 2D
bool getRotationMatrix(double *P1, double *P2, double *P3, double *rotationMatrix);

//! get the center of a polygon in 3D space and its area
bool getFaceCenter(const std::vector<double*> &pointList,
                   double *rotationMatrix,
                   double *faceCenter, double &faceArea);

//! get the center of mass of a shape compound
gp_Pnt getCenterOfMass(const TopoDS_Compound &compound);

//! center of mass of a list of shapes
gp_Pnt getCenterOfMass(const TopTools_ListOfShape &aListOfShapes);

//! genus of a face
int faceGenus(const TopoDS_Face &aFace);

//! get the coefficients of a plane from three points4
std::tuple<double,double,double,double> getPlaneCoefficients(const std::vector<std::tuple<double,double,double>> &threePoints);

// end namespace
}

#endif // GEOMTOOLSCLASS_H
