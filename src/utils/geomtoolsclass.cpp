//! ----------------
//! custom includes
//! ----------------
#include "geomtoolsclass.h"
using namespace GeomToolsClass;

//! ----
//! OCC
//! ----
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Geom_Surface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <gp_Pnt.hxx>
#include <gp_Ax1.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Cone.hxx>
#include <gp_Lin.hxx>
#include <gp_Circ.hxx>
#include <gp_Elips.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Ellipse.hxx>
#include <gp_Parab.hxx>
#include <GProp_PGProps.hxx>
#include <BRep_Builder.hxx>

//! --------------------
//! function: faceGenus
//! details:
//! --------------------
int GeomToolsClass::faceGenus(const TopoDS_Face &aFace)
{
    if(aFace.IsNull()) return -1;
    int holes = 0;
    for(TopExp_Explorer anExp(aFace,TopAbs_WIRE);anExp.More();anExp.Next())
    {
        if(TopoDS::Wire(anExp.Current()).IsNull()) holes++;
    }
    return holes;
}

//! -------------------------
//! function: getBoundingBox
//! details:
//! -------------------------
void GeomToolsClass::getBoundingBox(const TopoDS_Shape &shape, double &L1, double &L2, double &L3)
{
    Bnd_Box boundingBox;
    BRepBndLib::Add(shape, boundingBox);
    double Xmin,Ymin,Zmin,Xmax,Ymax,Zmax;
    boundingBox.Get(Xmin,Ymin,Zmin,Xmax,Ymax,Zmax);
    L1 = fabs(Xmax-Xmin);
    L2 = fabs(Ymax-Ymin);
    L3 = fabs(Zmax-Zmin);
}

//! ----------------------
//! function: getFaceType
//! details:
//! ----------------------
void GeomToolsClass::getFaceType(const TopoDS_Face &aFace, GeomAbs_SurfaceType &type)
{
    const occHandle(Geom_Surface) &surface = BRep_Tool::Surface(aFace);
    GeomAdaptor_Surface adapt(surface);
    type = adapt.GetType();
}

//! ----------------------
//! function: getLineType
//! details:
//! ----------------------
void GeomToolsClass::getCurveType(const TopoDS_Edge &anEdge, GeomAbs_CurveType &type)
{
    double firstParam, lastParam;
    Q_UNUSED(firstParam)
    Q_UNUSED(lastParam)
    const occHandle(Geom_Curve) &curve = BRep_Tool::Curve(anEdge, firstParam, lastParam);
    GeomAdaptor_Curve adapt(curve);
    type = adapt.GetType();
}

//! -----------------------------------------------
//! function: getPlanarFaceInfo
//! details:  P is the centroid of the planar face
//! -----------------------------------------------
bool GeomToolsClass::getPlanarFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P)
{
    if(aFace.IsNull()) exit(1);
    const occHandle(Geom_Surface) &surface = BRep_Tool::Surface(aFace);
    GeomAdaptor_Surface adapt(surface);
    if(adapt.GetType()==GeomAbs_Plane)
    {
        occHandle(Geom_Plane) pln = occHandle(Geom_Plane)::DownCast(surface);
        theAxis = pln->Axis();
        GProp_GProps props;
        BRepGProp::SurfaceProperties(aFace,props,0.01);
        P = props.CentreOfMass();
        return true;
    }
    return false;
}

//! ---------------------------------
//! function: getCylindricalFaceInfo
//! details:
//! ---------------------------------
bool GeomToolsClass::getCylindricalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P)
{
    const occHandle(Geom_Surface) &surface = BRep_Tool::Surface(aFace);
    GeomAdaptor_Surface adapt(surface);
    if(adapt.GetType()==GeomAbs_Cylinder)
    {
        theAxis = adapt.Cylinder().Axis();
        P = theAxis.Location();
        return true;
    }
    return false;
}

//! -------------------------------
//! function: getSphericalFaceInfo
//! details:
//! -------------------------------
bool GeomToolsClass::getSphericalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P)
{
    const occHandle(Geom_Surface) &surface = BRep_Tool::Surface(aFace);
    GeomAdaptor_Surface adapt(surface);
    if(adapt.GetType()==GeomAbs_Sphere)
    {
        occHandle(Geom_SphericalSurface) sph = occHandle(Geom_SphericalSurface)::DownCast(surface);
        theAxis = sph->Axis();
        P = theAxis.Location();
        return true;
    }
    return false;
}

//! -----------------------------
//! function: getConicalFaceInfo
//! details:
//! -----------------------------
bool GeomToolsClass::getConicalFaceInfo(const TopoDS_Face &aFace, gp_Ax1 &theAxis, gp_Pnt &P)
{
    const occHandle(Geom_Surface) &surface = BRep_Tool::Surface(aFace);
    GeomAdaptor_Surface adapt(surface);
    if(adapt.GetType()==GeomAbs_Cone)
    {
        theAxis = adapt.Cone().Axis();
        P = theAxis.Location();
        return true;
    }
    return false;
}

//! ----------------------
//! function: getEdgeInfo
//! details:
//! ----------------------
bool GeomToolsClass::getEdgeInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P)
{
    double uin, uend;
    const occHandle(Geom_Curve) &curve = BRep_Tool::Curve(anEdge, uin, uend);
    GeomAdaptor_Curve adapt(curve);
    switch(adapt.GetType())
    {
    case GeomAbs_Line:
    {
        double u = 0.5*(uin+uend);
        P = curve->Value(u);
        theAxis.SetDirection(adapt.Line().Direction());
        theAxis.SetLocation(adapt.Line().Location());
        cerr<<"____get edge info P("<<P.X()<<", "<<P.Y()<<", "<<P.Z()<<")____"<<endl;
    }
        break;

    case GeomAbs_Ellipse:
    {
        const gp_Elips &ellipse = adapt.Ellipse();
        theAxis = ellipse.Axis();

        const gp_Pnt &F1 = ellipse.Focus1();
        const gp_Pnt &F2 = ellipse.Focus2();

        P.SetX(0.5*(F1.X()+F2.X()));
        P.SetY(0.5*(F1.Y()+F2.Y()));
        P.SetZ(0.5*(F1.Z()+F2.Z()));
    }
        break;

    case GeomAbs_Circle:
    {
        const gp_Circ &circle = adapt.Circle();
        P = circle.Position().Axis().Location();
        theAxis = circle.Axis();
    }
        break;

    case GeomAbs_Parabola:
    {
        const gp_Parab &parabola = adapt.Parabola();
        P = parabola.Location();
        theAxis = parabola.Axis();
    }
        break;

    default:
    {
        double first = curve->FirstParameter();
        double last = curve->LastParameter();
        double u = 0.5*(first+last);
        P = curve->Value(u);
    }
        break;
    }
    return true;
}


//! ----------------------------
//! function: getCircleEdgeInfo
//! details:
//! ----------------------------
bool GeomToolsClass::getCircleInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P)
{
    double uin, uend;
    const occHandle(Geom_Curve) &curve = BRep_Tool::Curve(anEdge, uin, uend);
    GeomAdaptor_Curve adapt(curve);
    if(adapt.GetType()==GeomAbs_Circle)
    {
        theAxis = adapt.Circle().Axis();
        P = theAxis.Location();
        return true;
    }
    return false;
}

//! -------------------------
//! function: getEllipseInfo
//! details:
//! -------------------------
bool GeomToolsClass::getEllipseInfo(const TopoDS_Edge &anEdge, gp_Ax1 &theAxis, gp_Pnt &P)
{
    double uin, uend;
    const occHandle(Geom_Curve) &curve = BRep_Tool::Curve(anEdge, uin, uend);
    GeomAdaptor_Curve adapt(curve);
    if(adapt.GetType()==GeomAbs_Ellipse)
    {
        const gp_Elips &ellipse = adapt.Ellipse();
        theAxis = ellipse.Axis();

        const gp_Pnt &F1 = ellipse.Focus1();
        const gp_Pnt &F2 = ellipse.Focus2();

        P.SetX(0.5*(F1.X()+F2.X()));
        P.SetY(0.5*(F1.Y()+F2.Y()));
        P.SetZ(0.5*(F1.Z()+F2.Z()));
        return true;
    }
    return false;
}

//! ----------------------------
//! function: calculateCentroid
//! details:
//! ----------------------------
std::vector<double> GeomToolsClass::calculateCentroid(geometryDataBase *gDB, const std::vector<GeometryTag> &vecLoc)
{
    if(gDB==NULL) return std::vector<double>(3,0);
    if(vecLoc.size()==0) return std::vector<double>(3,0);

    //! ----------------------------------------------
    //! compute the centroid of the selected entities
    //! ----------------------------------------------
    gp_Pnt CM;
    GProp_PGProps prop;
    BRep_Builder theBuilder;
    TopoDS_Compound theCompound;
    theBuilder.MakeCompound(theCompound);

    for(std::vector<GeometryTag>::const_iterator it = vecLoc.cbegin(); it!=vecLoc.cend(); ++it)
    {
        const GeometryTag &loc = *it;
        int parentShape = loc.parentShapeNr;
        int childShape = loc.subTopNr;
        TopAbs_ShapeEnum type = loc.subShapeType;
        TopoDS_Shape S;

        switch(type)
        {
        case TopAbs_COMPSOLID: S = gDB->MapOfBodyTopologyMap.value(parentShape).csolidMap.FindKey(childShape); break;
        case TopAbs_FACE: S = gDB->MapOfBodyTopologyMap.value(parentShape).faceMap.FindKey(childShape); break;
        case TopAbs_EDGE: S = gDB->MapOfBodyTopologyMap.value(parentShape).edgeMap.FindKey(childShape); break;
        case TopAbs_VERTEX: break; // intentionally not present
        }
        theBuilder.Add(theCompound,S);
    }

    BRepGProp::SurfaceProperties(theCompound,prop,Standard_False);
    CM = prop.CentreOfMass();
    return std::vector<double> { CM.X(), CM.Y(), CM.Z() };
}

//! ----------------------------
//! function: getRotationMatrix
//! details:
//! ----------------------------
bool GeomToolsClass::getRotationMatrix(double *P1, double *P2, double *P3, double *rotationMatrix)
{
    const double TINY = std::numeric_limits<double>::epsilon();
    double xP1 = P1[0], yP1 = P1[1], zP1 = P1[2];
    double xP2 = P2[0], yP2 = P2[1], zP2 = P2[2];
    double xP3 = P3[0], yP3 = P3[1], zP3 = P3[2];

    double a11,a12,a13,a21,a22,a23,a31,a32,a33;

    //! ----------------------
    //! local Z normal vector
    //!  i      j       k
    //! x2-x1   y2-y1   z2-z1
    //! x3-x1   y3-y1   z3-z1
    //! ----------------------
    a31 = (yP2-yP1)*(zP3-zP1)-(zP2-zP1)*(yP3-yP1);
    a32 = (zP2-zP1)*(xP3-xP1)-(xP2-xP1)*(zP3-zP1);
    a33 = (xP2-xP1)*(yP3-yP1)-(yP2-yP1)*(xP3-xP1);

    double Lz = sqrt(pow(a31,2)+pow(a32,2)+pow(a33,2));

    if(Lz>TINY)
    {
        a31 = a31/Lz;  a32 = a32/Lz;  a33 = a33/Lz;

        //! (2-1)   local X
        double Lx = sqrt(pow(xP2-xP1,2)+pow(yP2-yP1,2)+pow(zP2-zP1,2));

        if(Lx>TINY)
        {
            a11 = (xP2-xP1)/Lx;
            a12 = (yP2-yP1)/Lx;
            a13 = (zP2-zP1)/Lx;

            //! --------------------------------
            //! local Y = local Z cross local X
            //!  i      j       k
            //! a31     a32     a33
            //! a11     a12     a13
            //! --------------------------------
            a21 = a32*a13 - a33*a12;
            a22 = a33*a11 - a31*a13;
            a23 = a31*a12 - a32*a11;
            double Ly = sqrt(pow(a21,2)+pow(a22,2)+pow(a23,2));

            if(Ly>TINY)
            {
                a21 = a21/Ly; a22 = a22/Ly; a23 = a23/Ly;
            }
            else
            {
                cerr<<"---->Ly too small<----"<<endl;
                return false;
            }
        }
        else
        {
            cerr<<"---->Lx too small<----"<<endl;;
            return false;
        }
    }
    else
    {
        cerr<<"---->Lz too small<----"<<endl;
        return false;
    }

    rotationMatrix[0] = a11;
    rotationMatrix[1] = a12;
    rotationMatrix[2] = a13;
    rotationMatrix[3] = a21;
    rotationMatrix[4] = a22;
    rotationMatrix[5] = a23;
    rotationMatrix[6] = a31;
    rotationMatrix[7] = a32;
    rotationMatrix[8] = a33;
    return true;
}

//! ----------------------
//! function: isCollinear
//! details:  optimized
//! ----------------------
bool GeomToolsClass::isCollinear(double *A0, double *A1, double *A2)
{
    double x0 = A0[0], y0 = A0[1], z0 = A0[2];
    double x1 = A1[0], y1 = A1[1], z1 = A1[2];
    double x2 = A2[0], y2 = A2[1], z2 = A2[2];

    double L10 = sqrt(pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2));
    double L20 = sqrt(pow(x2-x0,2)+pow(y2-y0,2)+pow(z2-z0,2));
    double sp = (x1-x0)*(x2-x0)+(y1-y0)*(y2-y0)+(z1-z0)*(z2-z0);
    double r = fabs(sp/(L10*L20));
    if(r<std::numeric_limits<double>::epsilon()) return false;

    else return true;
}

/*
//! ------------------------
//! function: getFaceCenter
//! details:
//! ------------------------
bool GeomToolsClass::getFaceCenter(const QList<QList<double>> &pointList,
                                   const QList<double> &rotationMatrix,
                                   QList<double> &faceCenter, double &faceArea)
{
    const double TINY = std::numeric_limits<double>::epsilon();

    //! -------------------------
    //! the face rotation matrix
    //! -------------------------
    double a11 = rotationMatrix.at(0);
    double a12 = rotationMatrix.at(1);
    double a13 = rotationMatrix.at(2);
    double a21 = rotationMatrix.at(3);
    double a22 = rotationMatrix.at(4);
    double a23 = rotationMatrix.at(5);
    double a31 = rotationMatrix.at(6);
    double a32 = rotationMatrix.at(7);
    double a33 = rotationMatrix.at(8);

    QList<double> P1 = pointList.at(0);
    double xP1 = P1.at(0), yP1 = P1.at(1), zP1 = P1.at(2);

    //! ----------------------
    //! the transformed nodes
    //! ----------------------
    QList<QList<double>> transfPointList;

    //! ---------------------------------
    //! transform the points coordinates
    //! ---------------------------------
    for(int n=0; n<pointList.length(); n++)
    {
        QList<double> curPoint = pointList.at(n);
        QList<double> transPoint;

        double dxP = curPoint.at(0)-xP1;
        double dyP = curPoint.at(1)-yP1;
        double dzP = curPoint.at(2)-zP1;

        //! ---------------------------------
        //! from global to local coordinates
        //! ---------------------------------
        double XP = dxP*a11+dyP*a12+dzP*a13;
        double YP = dxP*a21+dyP*a22+dzP*a23;
        double ZP = dxP*a31+dyP*a32+dzP*a33;

        transPoint.push_back(XP);
        transPoint.push_back(YP);
        transPoint.push_back(ZP);
        transfPointList.append(transPoint);
    }

    //! ---------------------------------------------------------
    //! calculate the ->signed<- area and the center of the face
    //! ---------------------------------------------------------
    double Area = 0.0;
    double CX = 0, CY = 0;
    for(int i=0; i<pointList.length()-1; i++)
    {
        Area = Area + transfPointList.at(i).at(0)*transfPointList.at(i+1).at(1)-
                transfPointList.at(i+1).at(0)*transfPointList.at(i).at(1);
        //Area = Area + transfPointList[i][0]*transfPointList[i+1][1]-
        //        transfPointList[i+1][0]*transfPointList[i][1];

    }
    Area = Area*0.5;

    if(Area<-TINY || Area>TINY)
    {
        for(int i=0; i<pointList.length()-1; i++)
        {
            CX = CX + (transfPointList.at(i).at(0)+transfPointList.at(i+1).at(0))*
                    (transfPointList.at(i).at(0)*transfPointList.at(i+1).at(1)-transfPointList.at(i+1).at(0)*transfPointList.at(i).at(1));
            CY = CY + (transfPointList.at(i).at(1)+transfPointList.at(i+1).at(1))*
                    (transfPointList.at(i).at(0)*transfPointList.at(i+1).at(1)-transfPointList.at(i+1).at(0)*transfPointList.at(i).at(1));
        }
        CX = CX/(6.0*Area);
        CY = CY/(6.0*Area);

        //! -----------------------------------
        //! from local system to global system
        //! -----------------------------------
        double Cx = xP1 + CX*a11+CY*a21;
        double Cy = yP1 + CX*a12+CY*a22;
        double Cz = zP1 + CX*a13+CY*a23;
        faceCenter.push_back(Cx);
        faceCenter.push_back(Cy);
        faceCenter.push_back(Cz);
        faceArea = Area;
        return true;
    }
    else
    {
        faceArea = Area;
        return false;
    }
}
*/

//! ------------------------
//! function: getFaceCenter
//! details:
//! ------------------------
bool GeomToolsClass::getFaceCenter(const std::vector<double*> &pointList,
                                   double *rotationMatrix,
                                   double *faceCenter, double &faceArea)
{
    const double TINY = std::numeric_limits<double>::epsilon();

    //! the face rotation matrix
    double a11 = rotationMatrix[0];
    double a12 = rotationMatrix[1];
    double a13 = rotationMatrix[2];
    double a21 = rotationMatrix[3];
    double a22 = rotationMatrix[4];
    double a23 = rotationMatrix[5];
    double a31 = rotationMatrix[6];
    double a32 = rotationMatrix[7];
    double a33 = rotationMatrix[8];

    double xP1 = pointList[0][0];
    double yP1 = pointList[0][1];
    double zP1 = pointList[0][2];

    //! the transformed nodes
    struct Pt
    {
        double xt,yt,zt;
        Pt(double x=0, double y=0, double z=0):xt(x),yt(y),zt(z){;}
        Pt(const Pt &aP) { xt=aP.xt; yt=aP.yt; zt=aP.zt; }
        bool operator == (const Pt &aP) { if(xt==aP.xt && yt==aP.yt && zt==aP.zt) return true; return false; }
        Pt operator = (const Pt &aP) { xt=aP.xt; yt=aP.yt; zt=aP.zt; return *this; }
    };
    std::vector<Pt> transfPointList;

    //! transform the points coordinates
    for(int n=0; n<pointList.size(); n++)
    {
        double *curPoint = pointList[n];

        double dxP = curPoint[0]-xP1;
        double dyP = curPoint[1]-yP1;
        double dzP = curPoint[2]-zP1;

        //! from global to local coordinates
        double XP = dxP*a11+dyP*a12+dzP*a13;
        double YP = dxP*a21+dyP*a22+dzP*a23;
        double ZP = dxP*a31+dyP*a32+dzP*a33;

        transfPointList.push_back(Pt(XP,YP,ZP));
    }

    //! calculate the ->signed<- area and the center of the face
    double Area = 0.0;
    double CX = 0, CY = 0;
    for(int i=0; i<pointList.size()-1; i++)
    {
        Area = Area + transfPointList[i].xt*transfPointList[i+1].yt-transfPointList[i+1].xt*transfPointList[i].yt;
    }
    Area = Area*0.5;

    if(Area<-TINY || Area>TINY)
    {
        for(int i=0; i<pointList.size()-1; i++)
        {
            CX = CX + (transfPointList[i].xt+transfPointList[i+1].xt)*
                    (transfPointList[i].xt*transfPointList[i+1].yt-transfPointList[i+1].xt*transfPointList[i].yt);
            CY = CY + (transfPointList[i].yt+transfPointList[i+1].yt)*
                    (transfPointList[i].xt*transfPointList[i+1].yt-transfPointList[i+1].xt*transfPointList[i].yt);
        }
        CX /= (6.0*Area);
        CY /= (6.0*Area);

        //! from local system to global system
        double Cx = xP1 + CX*a11+CY*a21;
        double Cy = yP1 + CX*a12+CY*a22;
        double Cz = zP1 + CX*a13+CY*a23;
        faceCenter[0] = Cx; faceCenter[1] = Cy; faceCenter[2] = Cz;
        faceArea = Area;
        return true;
    }
    else
    {
        faceArea = Area;
        return false;
    }
}

//! -------------------------
//! fuction: getCenterOfMass
//! details:
//! -------------------------
gp_Pnt GeomToolsClass::getCenterOfMass(const TopoDS_Compound &compound)
{
    GProp_GProps props;
    BRepGProp::VolumeProperties(compound,props);
    gp_Pnt CM = props.CentreOfMass();
    return CM;
}

//! -------------------------
//! fuction: getCenterOfMass
//! details: overload
//! -------------------------
gp_Pnt GeomToolsClass::getCenterOfMass(const TopTools_ListOfShape &aListOfShapes)
{
    TopoDS_Builder theBuilder;
    TopoDS_Compound compound;
    theBuilder.MakeCompound(compound);
    for(TopTools_ListIteratorOfListOfShape it(aListOfShapes); it.More(); it.Next())
    {
        const TopoDS_Shape &curShape = it.Value();
        theBuilder.Add(compound,curShape);
    }
    GProp_GProps props;
    BRepGProp::VolumeProperties(compound,props);
    gp_Pnt CM = props.CentreOfMass();
    return CM;
}

/*
//! -------------------------------
//! function: getPlaneCoefficients
//! details:
//! -------------------------------
std::tuple<double,double,double,double> GeomToolsClass::getPlaneCoefficients(const std::vector<std::tuple<double,double,double>> &threePoints)
{
    const std::tuple<double,double,double> &first = threePoints[0];
    const std::tuple<double,double,double> &second = threePoints[1];
    const std::tuple<double,double,double> &third = threePoints[2];
    double x0 = std::get<0>(first); double y0 = std::get<1>(first); double z0 = std::get<2>(first);
    double x1 = std::get<0>(second); double y1 = std::get<1>(second); double z1 = std::get<2>(second);
    double x2 = std::get<0>(third); double y2 = std::get<1>(third); double z2 = std::get<2>(third);
    double a = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);
    double b = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
    double c = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);
    double d = -x0*a-y0*b-z0*c;
    return std::tuple<double,double,double,double> {a,b,c,d};
}
*/
