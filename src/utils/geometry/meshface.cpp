//! ----------------
//! custom includes
//! ----------------
#include "meshface.h"
#include "polygon.h"

//! ----
//! OCC
//! ----
#include <TColStd_PackedMapOfInteger.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>

//! ----
//! C++
//! ----
#include <vector>
#include <algorithm>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
meshFace::meshFace(const opencascade::handle<Ng_MeshVS_DataSourceFace> &theFaceMeshDS):
    myFaceMeshDS(theFaceMeshDS)
{
    if(myFaceMeshDS.IsNull()) return;
    if(myFaceMeshDS->GetAllElements().Extent()==0) return;
}

//! ----------------------
//! function: setGeometry
//! details:
//! ----------------------
void meshFace::setGeometry(const opencascade::handle<Ng_MeshVS_DataSourceFace> &theFaceMeshDS)
{
    myFaceMeshDS = theFaceMeshDS;
}

//! ---------------------------------------------------------
//! function: area
//! details:  supports first order triangles and quadrangles
//! ---------------------------------------------------------
double meshFace::area()
{
    //! ----------
    //! reset val
    //! ----------
    double val = 0;

    TColStd_PackedMapOfInteger emap = myFaceMeshDS->GetAllElements();

    MeshVS_EntityType atype;
    double buf[12];
    TColStd_Array1OfReal coords(*buf,1,12);
    int NbNodes;

    for(TColStd_MapIteratorOfPackedMapOfInteger it(emap); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        bool isFound = myFaceMeshDS->GetGeom(globalElementID,true,coords,NbNodes,atype);
        if(!isFound) continue;

        std::vector<polygon::Point> aPolygon;
        for(int i=0; i<NbNodes; i++)
        {
            int s = 3*i;
            aPolygon.push_back(polygon::Point(coords(s+1),coords(s+2),coords(s+3)));
        }

        double elementArea = polygon::area3D_Polygon(aPolygon);
        val += elementArea;
    }
    return val;
}

//! --------------------------
//! function: pointProjection
//! details:
//! --------------------------
bool meshFace::pointProjection(double *aPoint, double *theProjectedPoint, double eps)
{
    //cout<<"meshFace::pointProjection()->____function called for point("<<aPoint[0]<<", "<<
    //      aPoint[1]<<", "<<aPoint[2]<<")____"<<endl;

    std::vector<polygon::Point> validProjectedPoints;
    std::vector<double> distances;
    double curProjectedPoint[3];

    //! --------------------------
    //! iterate over the elements
    //! --------------------------
    MeshVS_EntityType atype;
    double buf[12];
    TColStd_Array1OfReal coords(*buf,1,12);
    int NbNodes;

    TColStd_PackedMapOfInteger emap = myFaceMeshDS->GetAllElements();

    for(TColStd_MapIteratorOfPackedMapOfInteger it(emap); it.More(); it.Next())
    {
        int globalElementID = it.Key();
        bool isFound = myFaceMeshDS->GetGeom(globalElementID,true,coords,NbNodes,atype);
        if(!isFound) continue;

        std::vector<polygon::Point> aPolygon;
        for(int i=0; i<NbNodes; i++)
        {
            int s=3*i;
            aPolygon.push_back(polygon::Point(coords(s+1),coords(s+2),coords(s+3)));
        }
        double a,b,c,d;
        polygon::planeCoefficients(aPolygon,a,b,c,d);

        double l = sqrt(a*a+b*b+c*c);
        a /= l; b /= l; c /= l; d /= l;

        //! ----------------------------------------------------
        //! projection of the point on the plane of the element
        //! ----------------------------------------------------
        double distance = (a*aPoint[0]+b*aPoint[1]+c*aPoint[2]+d)/sqrt(a*a+b*b+c*c);

        curProjectedPoint[0] = aPoint[0] - a*distance;
        curProjectedPoint[1] = aPoint[1] - b*distance;
        curProjectedPoint[2] = aPoint[2] - c*distance;

        polygon::Point PP(curProjectedPoint[0],curProjectedPoint[1],curProjectedPoint[2]);
        bool normalizeCrossProducts = true;
        double inOutTolerance = sqrt(polygon::area3D_Polygon(aPolygon));
        bool isInside = polygon::isPointInside(PP,aPolygon,normalizeCrossProducts,inOutTolerance);
        if(!isInside) continue;
        validProjectedPoints.push_back(PP);
        distances.push_back(fabs(distance));
    }
    if(validProjectedPoints.size()==0)
    {
        cerr<<"meshFace::pointProjection()->____failure in projection: list of points empty____"<<endl;
        //exit(10);
        return false;
    }
    //cout<<"meshFace::pointProjection()->____projection OK____"<<endl;
    int minElementIndex = std::min_element(distances.begin(),distances.end()) - distances.begin();
    theProjectedPoint[0] = validProjectedPoints[minElementIndex].x;
    theProjectedPoint[1] = validProjectedPoints[minElementIndex].y;
    theProjectedPoint[2] = validProjectedPoints[minElementIndex].z;

    if(distances.at(minElementIndex)>1.0)
    {
        cerr<<"meshFace::pointProjection()->____failure in projection: projection point too far____"<<endl;
        //exit(11);
        return false;
    }
    return true;
}

//! ------------------------
//! function: randomPointOn
//! details:
//! ------------------------
#include <random>
void meshFace::randomPointOn(double *point)
{
    //! ----------------------------------------
    //! sample a random position in (u,v) space
    //! ----------------------------------------
    std::random_device rd;
    std::mt19937 e(rd());
    std::uniform_int_distribution<> dist(1,myFaceMeshDS->GetAllElements().Extent());
    int sampledElementNr = dist(e);

    int globalNodeID = myFaceMeshDS->myNodesMap.FindKey(sampledElementNr);
    int NbNodes;
    double buf[24];
    MeshVS_EntityType aType;
    TColStd_Array1OfReal x(*buf,1,24);
    myFaceMeshDS->GetGeom(globalNodeID,true,x,NbNodes,aType);
    double x1,y1,z1,x2,y2,z2,x3,y3,z3;
    x1 = x(1); y1 = x(2); z1 = x(3);
    x2 = x(4); y2 = x(5); z2 = x(6);
    x3 = x(7); y3 = x(8); z3 = x(9);

    //! ---------------------------------
    //! transform into a planar triangle
    //! ---------------------------------
    double exx = x2-x1;
    double exy = y2-y1;
    double exz = z2-z1;
    double lx = sqrt(exx*exx+exy*exy+exz*exz);
    exx /= lx; exy /= lx; exz /= lx;

    double a = x3-x1;
    double b = y3-y1;
    double c = z3-z1;

    double ezx = exy*c-exz*b;
    double ezy = exz*a-exx*c;
    double ezz = exx*b-exy*a;
    double lz = sqrt(ezx*ezx+ezy*ezy+ezz*ezz);
    ezx /= lz; ezy /= lz; ezz /= lz;

    double eyx = ezy*exz-ezz*exy;
    double eyy = ezz*exx-ezx*exy;
    double eyz = ezx*exy-ezy*exx;
    double ly = sqrt(eyx*eyx+eyy*eyy+eyz*eyz);
    eyx /= ly; eyy /= ly; eyz /= ly;

    double X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3;

    X1 = 0; Y1 = 0; Z1 = 0;
    X2 = lx; Y2 = 0; Z2 = 0;
    X3 = (x3-x1)*exx+(y3-y1)*eyx+(z3-z1)*ezx;
    Y3 = (x3-x1)*exz+(y3-y1)*eyz+(z3-z1)*ezz;
    Z3 = 0;

    //! ------------------------------------
    //! a random (u,v) in range [0,1]x[0,1]
    //! ------------------------------------
    std::mt19937 e1(rd());
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> distu(0, 1);
    std::uniform_real_distribution<> distv(0, 1);
    double u = distu(e1);
    double v = distv(e2);

    //! ---------------------------------
    //! get the random point coordinates
    //! ---------------------------------
    double XS = X1 + (X2-X1)*u+(X3-X1)*v;
    double YS = Y1 + (Y2-Y1)*u+(Y3-Y1)*v;

    //! ---------------
    //! transform back
    //! ---------------
    point[1] = x1+XS*exx+YS*eyx;
    point[2] = y1+XS*exy+YS*eyy;
    point[3] = z1+XS*exz+YS*eyz;
}
