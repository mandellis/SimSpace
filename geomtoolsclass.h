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

//! ---
//! Qt
//! ---
#include <QList>

class TopoDS_Face;
class TopoDS_Edge;

class GeomToolsClass
{
public:

    GeomToolsClass();
    ~GeomToolsClass();

    static void getBoundingBox(const TopoDS_Shape &shape, double &L1, double &L2, double &L3);
    static void getFaceType(const TopoDS_Face &aFace, GeomAbs_SurfaceType &type);
    static void getCurveType(const TopoDS_Edge &anEdge, GeomAbs_CurveType &type);
    static bool getPlanarFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &centroid);
    static bool getCylindricalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P);
    static bool getSphericalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P);
    static bool getConicalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P);
    static bool getEdgeInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P);
    static bool getCircleInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P);
    static bool getEllipseInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P);
    static QList<double> calculateCentroid(geometryDataBase *gDB, const std::vector<GeometryTag> &vecLoc);

    //! ---------------------------------
    //! test if three points are aligned
    //! ---------------------------------
    static bool isCollinear(const QList<double> &A0, const QList<double> &A1, const QList<double> &A2);
    static bool isCollinear(double *A0, double *A1, double *A2);
    static bool isCollinear(const QList<QList<double>> &threePoints);

    //! ---------------------------------
    //! rotation matrix for 3D->plane 2D
    //! ---------------------------------
    static bool getRotationMatrix(const QList<QList<double>> &threePoints, QList<double> &rotationMatrix);
    static bool getRotationMatrix(double *P1, double *P2, double *P3, double *rotationMatrix);

    //! -----------------------------------------------------
    //! get the center of a polygon in 3D space and its area
    //! -----------------------------------------------------
    static bool getFaceCenter(const QList<QList<double>> &pointList,
                              const QList<double> &rotationMatrix,
                              QList<double> &faceCenter, double &faceArea);

    static bool getFaceCenter(const std::vector<double*> &pointList,
                              double *rotationMatrix,
                              double *faceCenter, double &faceArea);

    //! --------------------------------
    //! function: linePlaneIntersection
    //! details:
    //! --------------------------------
    static bool linePlaneIntersection(const gp_Ax3 &thePlane, gp_Pnt &thePoint)
    {
        const occHandle(Geom_Plane) &transientPlane = new Geom_Plane(thePlane);
        double a,b,c,d;
        transientPlane->Coefficients(a,b,c,d);
        double x0 = thePlane.Location().X();
        double y0 = thePlane.Location().Y();
        double z0 = thePlane.Location().Z();
        double x1 = thePoint.X();
        double y1 = thePoint.Y();
        double z1 = thePoint.Z();
        double den = a*(x1-x0)+b*(y1-y0)+c*(z1-z0);
        //! no intersection
        if(den<gp::Resolution()) return false;

        double t = (a*x0+b*y0+c*z0+d)/den;
        thePoint.SetX(x0+(x1-x0)*t);
        thePoint.SetX(y0+(y1-y0)*t);
        thePoint.SetX(z0+(z1-z0)*t);
    }

    //! ------------------------------------------------------------------------
    //! function: planePolygonIntersection
    //! details:  perform the plane - polygon intersection
    //!           it scans the points of the polygon, defined through a list of
    //!           {x,y,z}. If the points are all at the same plane side no
    //!           intersection exists
    //! ------------------------------------------------------------------------
    static bool planePolygonIntersection(const gp_Ax3& aPlane,
                                         const QList<QList<double>> &polygon,
                                         QList<gp_Pnt> &intersectionPolySegment)
    {

        bool hasIntersection = false;

        //! --------------------------------------
        //! origin of the local coordinate system
        //! --------------------------------------
        const gp_Pnt &O = aPlane.Location();

        //! -------------
        //! plane normal
        //! -------------
        const gp_Ax1 &normal =  aPlane.Axis();
        const gp_Dir &normalDir = normal.Direction();

        //! ---------------------------
        //! first point of the polygon
        //! ---------------------------
        const QList<double> &firstPoint = polygon.first();

        gp_Pnt fP(firstPoint.at(0),firstPoint.at(1),firstPoint.at(2));
        gp_Vec fPO(O,fP);
        gp_Dir fPODir(fPO);

        double res_old = normalDir.Dot(fPODir);
        int n, NbPoints = polygon.length();
        QList<double> point;
        gp_Pnt P;
        double res;
        for(n=1; n<NbPoints; n++)
        {
            point = polygon.at(n);
            P.SetX(point.at(0));
            P.SetY(point.at(1));
            P.SetZ(point.at(2));
            gp_Vec PO(O,P);
            gp_Dir PODir(PO);
            res = normalDir.Dot(PODir);
            if(res_old*res>=0)
            {
                //! ---------------------------------------------------
                //! has intersection: calculate the intersection point
                //! ---------------------------------------------------
                hasIntersection = true;
                gp_Pnt intersectionPoint;
                GeomToolsClass::linePlaneIntersection(aPlane,intersectionPoint);
                intersectionPolySegment<<intersectionPoint;
            }
            res_old = res;
        }
        if(hasIntersection==true)
        {
            gp_Vec PO(O,P);
            gp_Dir PODir(PO);
            if(normalDir.Dot(PODir)*normalDir.Dot(PO)>=0)
            {
                //! ---------------------------------------------------
                //! has intersection: calculate the intersection point
                //! ---------------------------------------------------
                hasIntersection = true;
                gp_Pnt intersectionPoint;
                GeomToolsClass::linePlaneIntersection(aPlane,intersectionPoint);
                intersectionPolySegment<<intersectionPoint;
            }
        }
        //! --------------------------------------------------
        //! all the polygon points are at the same plane side
        //! no intersection exists
        //! --------------------------------------------------
        return hasIntersection;
    }

    //! -------------------------------------------
    //! get the center of mass of a shape compound
    //! -------------------------------------------
    static gp_Pnt getCenterOfMass(const TopoDS_Compound &compound);
    static gp_Pnt getCenterOfMass(const TopTools_ListOfShape &aListOfShapes);

    //! -----------------------
    //! calculate polygon area
    //! -----------------------
    static double polygon_area(const QList<QList<double>> &points);

    //static std::vector<double> getPolygonNormal(const std::vector<std::vector<double>> &vertices);
    static std::vector<double> getPlaneCoefficients(const std::vector<std::vector<double>> &threePoints);

    //! ----------
    //! faceGenus
    //! ----------
    static int faceGenus(const TopoDS_Face &aFace);
};

#endif // GEOMTOOLSCLASS_H
