#include <algorithm>
#include <limits>
#include "polygon.h"
using namespace polygon;

//! -----------------------
//! function: triangleArea
//! details:
//! -----------------------
double polygon::triangleArea(const polygon::Point &P0, const polygon::Point &P1, const polygon::Point &P2)
{
    double x1 = P0.x; double y1 = P0.y; double z1 = P0.z;
    double x2 = P1.x; double y2 = P1.y; double z2 = P1.z;
    double x3 = P2.x; double y3 = P2.y; double z3 = P2.z;
    double L21x = x2-x1; double L21y = y2-y1; double L21z = z2-z1;
    double L31x = x3-x1; double L31y = y3-y1; double L31z = z3-z1;
    double area = 0.5*(sqrt(pow(L21y*L31z-L21z*L31y,2)+
                            pow(L21z*L31x-L21x*L31z,2)+
                            pow(L21x*L31y-L21y*L31x,2)));
    return area;
}
//! -----------------------------------------------------------------------------------------------------
//! function: isPointInside
//! details:  check if a point is inside a polygon in 3D space
//!
//!           ALGO in:  https://math.stackexchange.com/questions/914795/
//!           considering-a-convex-polygon-lying-on-a-plane-in-3d-space-how-can-i-know-if-a-p
//!
//!           The following solution is elegant and efficient, and relatively simple to implement.
//!           First the idea, then the algorithm. Consider a point x and an n-sided convex polygon P
//!           with ordered boundary vertices (p1,p2,...pn), all lying in a plane in 3D space.
//!           We want to know whether x lies inside P, outside of it, or on its boundary.
//!           Notice that triangles can be constructed using x and any pair of adjacent pi and pi+1.
//!           If x lies inside the polygon, then the normals of all triangles (x,pi,pi+1) point in
//!           the same direction.
//!           Now the algorithm. First, calculate the cross-product N=(p1−x)×(p2−x); we will consider
//!           N to be the normal direction of P, and compare all other normals against it.
//!           Next, we compute all remaining (pi−x)×(pi+1−x). Each of these cross-products either points
//!           in the same direction as N, or in the direction of −N. To find out in which way they point,
//!           we take the dot product of (pi−x)×(pi+1−x) with N. If every dot product ((pi−x)×(pi+1−x))⋅N
//!           is positive, then x lies inside P; if at least one dot product is negative, then x lies
//!           outside of P; if a dot product is zero for some pi,pi+1, then x lies on the boundary of
//!           the polygon; in particular, it lies on the edge connecting pi to pi+1.
//!
//! -----------------------------------------------------------------------------------------------------
bool polygon::isPointInside(const polygon::Point &P,
                            const std::vector<polygon::Point> &V,
                            bool normalizeCrossProducts,
                            double tolerance)
{
    Point P0 = V[0];
    Point P1 = V[1];

    //! -----------------------------------------------------
    //! cross product (P0-P)x(P1-P)
    //! this is the normal N to the polygon (not normalized)
    //! -----------------------------------------------------
    double lx = P0.x-P.x; double ly = P0.y-P.y; double lz = P0.z-P.z;
    double mx = P1.x-P.x; double my = P1.y-P.y; double mz = P1.z-P.z;
    double nx = ly*mz-lz*my; double ny = lz*mx-lx*mz; double nz = lx*my-ly*mx;

    //! ----------
    //! normalize
    //! ----------
    if(normalizeCrossProducts)
    {
        double L = sqrt(nx*nx+ny*ny+nz*nz);
        lx /= L; ly /= L; lz /= L;
    }

    //! -------------------------------------------------------------------------
    //! work on the other polygon triangles. If NbPoints = 4 these triangles are
    //! built with (P,1,2) (P,2,3) (P3,0)
    //! -------------------------------------------------------------------------
    size_t NbPoints = V.size();
    for(int i=1; i<NbPoints; i++)
    {
        Point Pi = V[i%NbPoints];
        Point Pii = V[(i+1)%NbPoints];

        //! ----------------------------------------------
        //! cross product (Pi-P)x(Pii-P) (not normalized)
        //! ----------------------------------------------
        double llx = Pi.x-P.x; double lly = Pi.y-P.y; double llz = Pi.z-P.z;
        double mmx = Pii.x-P.x; double mmy = Pii.y-P.y; double mmz = Pii.z-P.z;
        double nnx = lly*mmz-llz*mmy; double nny = llz*mmx-llx*mmz; double nnz = llx*mmy-lly*mmx;

        if(normalizeCrossProducts)
        {
            double LL = sqrt(nnx*nnx+nny*nny+nnz*nnz);
            llx /= LL; lly /= LL; llz /= LL;
        }

        //! -----------------------------------------
        //! dot product between (Pi-P)x(Pii-P) and N
        //! -----------------------------------------
        double dot = nx*nnx+ny*nny+nz*nnz;

        //! ---------------------------------------
        //! "a bit negative" means "a bit outside"
        //! ---------------------------------------
        if(dot<-tolerance) return false;
    }
    return true;
}

//! ---------------------------
//! function: getPolygonCenter
//! details:
//! ---------------------------
bool polygon::getPolygonCenter(const std::vector<polygon::Point> &aPolygon,
                               polygon::Point &faceCenter, double &faceArea)
{
    //! --------------------------
    //! the polygon is a triangle
    //! --------------------------
    if(aPolygon.size()==3)
    {
        faceCenter.x = (aPolygon[0].x+aPolygon[1].x+aPolygon[2].x)/3.0;
        faceCenter.y = (aPolygon[0].y+aPolygon[1].y+aPolygon[2].y)/3.0;
        faceCenter.z = (aPolygon[0].z+aPolygon[1].z+aPolygon[2].z)/3.0;
        faceArea = polygon::triangleArea(aPolygon[0],aPolygon[1],aPolygon[2]);
        return true;
    }

    //! -------------------------
    //! the face rotation matrix
    //! -------------------------
    std::vector<double> rotationMatrix;
    bool isDone = getRotationMatrix(aPolygon,rotationMatrix);
    if(!isDone)
    {
        cerr<<"polygon::getPolygonCenter()->____error in retrieving the rotation matrix____"<<endl;
        return false;
    }

    double a11 = rotationMatrix[0];
    double a12 = rotationMatrix[1];
    double a13 = rotationMatrix[2];
    double a21 = rotationMatrix[3];
    double a22 = rotationMatrix[4];
    double a23 = rotationMatrix[5];
    double a31 = rotationMatrix[6];
    double a32 = rotationMatrix[7];
    double a33 = rotationMatrix[8];

    const polygon::Point &P1 = aPolygon[0];
    double xP1 = P1.x; double yP1 = P1.y; double zP1 = P1.z;

    //! ----------------------
    //! the transformed nodes
    //! ----------------------
    std::vector<polygon::Point> transfPointList;

    //! ---------------------------------
    //! transform the points coordinates
    //! ---------------------------------
    size_t NbPoints = aPolygon.size();
    for(int n=0; n<NbPoints; n++)
    {
        const polygon::Point &curPoint = aPolygon[n];
        polygon::Point transPoint;

        double dxP = curPoint.x-xP1;
        double dyP = curPoint.y-yP1;
        double dzP = curPoint.z-zP1;

        //! ---------------------------------
        //! from global to local coordinates
        //! ---------------------------------
        double XP = dxP*a11+dyP*a12+dzP*a13;
        double YP = dxP*a21+dyP*a22+dzP*a23;
        double ZP = dxP*a31+dyP*a32+dzP*a33;

        transPoint.x = XP;
        transPoint.y = YP;
        transPoint.z = ZP;
        transfPointList.push_back(transPoint);
    }

    //! -----------------------------------------------------
    //! calculate the signed area and the center of the face
    //! -----------------------------------------------------
    double Area = 0.0, CX = 0, CY = 0;
    for(int i=0; i<NbPoints-1; i++)
    {
        const polygon::Point &Pi = transfPointList[i];
        const polygon::Point &Pii = transfPointList[i+1];
        Area += Pi.x*Pii.y-Pii.x*Pi.y;
    }
    Area = Area*0.5;
    if(fabs(Area)<1e-16)
    {
        faceArea = Area;
        return false;
    }

    for(int i=0; i<NbPoints-1; i++)
    {
        const polygon::Point &Pi = transfPointList[i];
        const polygon::Point &Pii = transfPointList[i+1];

        CX = CX + (Pi.x+Pii.x)*(Pi.x*Pii.y-Pii.x*Pi.y);
        CY = CY + (Pi.y+Pii.y)*(Pi.x*Pii.y-Pii.x*Pi.y);
    }
    CX = CX/(6*Area);
    CY = CY/(6*Area);

    //! -----------------------------------
    //! from local system to global system
    //! -----------------------------------
    double Cx = xP1 + CX*a11+CY*a21;
    double Cy = yP1 + CX*a12+CY*a22;
    double Cz = zP1 + CX*a13+CY*a23;
    faceCenter.x = Cx;
    faceCenter.y = Cy;
    faceCenter.z = Cz;
    faceArea = Area;
    return true;
}

//! ---------------------------------
//! function: polygonIntersectsPlane
//! details:
//! ---------------------------------
bool polygon::testPolygonPlaneIntersection(const std::vector<polygon::Point> &aPolygon, double a, double b, double c, double d/*, int &NbPoints*/)
{
    bool intersection = false;

    double l0 = sqrt(a*a+b*b+c*c);
    a /= l0; b /= l0; c /= l0; d /= l0;

    double ap,bp,cp,dp;
    polygon::planeCoefficients(aPolygon,ap,bp,cp,dp);
    double l = sqrt(ap*ap+bp*bp+cp*cp);
    ap /= l; bp /= l; cp /= l; dp /= l;

    double dotplane = a*ap+b*bp+c*cp;
    if(fabs(dotplane)>=1.0)
    {
        if(d != dp) return false;   //! the polygon plane is parallel to the input plane
        else return true;           //! the polygon lays on the plane
    }

    //! ---------------------
    //! point plane distance
    //! ---------------------
    const polygon::Point &aP = aPolygon[0];
    double d_old = (a*aP.x+b*aP.y+c*aP.z+d);
    for(int i=1; i<aPolygon.size(); i++)
    {
        const polygon::Point &aP = aPolygon[i];
        double d_new = (a*aP.x+b*aP.y+c*aP.z+d);

        //! -----
        //! sign
        //! -----
        if(d_new*d_old<=0)
        {
            intersection = true;
            break;
        }
    }
    return intersection;
}


//! ----------------------------
//! function: getRotationMatrix
//! details:
//! ----------------------------
bool polygon::getRotationMatrix(const std::vector<polygon::Point> &aPolygon,
                                std::vector<double> &rotationMatrix)
{
    const polygon::Point &P1 = aPolygon[0];
    const polygon::Point &P2 = aPolygon[1];

    double xP1 = P1.x; double yP1 = P1.y; double zP1 = P1.z;
    double xP2 = P2.x; double yP2 = P2.y; double zP2 = P2.z;

    //! -------------------------------------
    //! search for a third not aligned point
    //! -------------------------------------
    bool found = false;
    polygon::Point P3;
    int NbPoints = int(aPolygon.size());
    for(int n=2; n<NbPoints; n++)
    {
        P3 = aPolygon[n];
        bool aligned = areCollinear(P1,P2,P3);
        if(aligned == true) continue;
        else
        {
            found = true;
            break;
        }
    }
    if(!found) return false;

    double xP3 = P3.x; double yP3 = P3.y; double zP3 = P3.z;
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
    if(Lz<1e-16)
    {
        cout<<"polygon::getRotationMatrix()->____returning false because of Lz____"<<endl;
        return false;
    }
    a31 = a31/Lz;  a32 = a32/Lz;  a33 = a33/Lz;

    //! ----------------
    //! (2-1)   local X
    //! ----------------
    double Lx = sqrt(pow(xP2-xP1,2)+pow(yP2-yP1,2)+pow(zP2-zP1,2));

    if(Lx<1e-16)
    {
        cout<<"polygon::getRotationMatrix()->____returning false because of Lx____"<<endl;
        return false;
    }

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

    if(Ly<1e-16)
    {
        cout<<"polygon::getRotationMatrix()->____returning false because of Ly____"<<endl;
        return false;
    }

    rotationMatrix.push_back(a11);
    rotationMatrix.push_back(a12);
    rotationMatrix.push_back(a13);
    rotationMatrix.push_back(a21);
    rotationMatrix.push_back(a22);
    rotationMatrix.push_back(a23);
    rotationMatrix.push_back(a31);
    rotationMatrix.push_back(a32);
    rotationMatrix.push_back(a33);
    return true;
}


//! -----------------------------
//! function: pointPlaneDistance
//! details:  signed distance
//! -----------------------------
double polygon::pointPlaneDistance(const polygon::Point &P, double a, double b,double c, double d)
{
    return (a*P.x+b*P.y+c*P.z+d)/sqrt(a*a+b*b+c*c);
}


//! -------------------------------------------------------------------
//! function: planeCoefficients
//! details:  input => a polygon; output => polygon plane coefficients
//! -------------------------------------------------------------------
bool polygon::planeCoefficients(const std::vector<polygon::Point> &points, double &a, double &b, double &c, double &d)
{
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    size_t NbPoints = points.size();

    //! -----------------
    //! first two points
    //! -----------------
    Point P0 = points[0];
    Point P1 = points[1];
    x0 = P0.x; y0 = P0.y; z0 = P0.z;
    x1 = P1.x; y1 = P1.y; z1 = P1.z;

    //! -------------------------------------
    //! third point: should not be collinear
    //! -------------------------------------
    for(int i=2; i<NbPoints; i++)
    {
        Point P2 = points[i];
        if(!areCollinear(P0,P1,P2))
        {
            x2 = P2.x; y2 = P2.y; z2 = P2.z;
            a = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);
            b = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
            c = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);
            d = -x0*a-y0*b-z0*c;
            return true;
        }
    }
    return false;
}


//! ----------------------------------------
//! function: areCollinear
//! details:  uses the area of the triangle
//! ----------------------------------------
bool polygon::areCollinear(const polygon::Point &P0, const polygon::Point &P1, const polygon::Point &P2)
{
    double x1 = P0.x; double y1 = P0.y; double z1 = P0.z;
    double x2 = P1.x; double y2 = P1.y; double z2 = P1.z;
    double x3 = P2.x; double y3 = P2.y; double z3 = P2.z;

    double L21x = x2-x1; double L21y = y2-y1; double L21z = z2-z1;
    double L31x = x3-x1; double L31y = y3-y1; double L31z = z3-z1;

    double area = sqrt(pow(L21y*L31z-L21z*L31y,2)+
                       pow(L21z*L31x-L21x*L31z,2)+
                       pow(L21x*L31y-L21y*L31x,2));

    if(area<=1e-20) return true; return false;
}

//! ----------------------------------------------------------------------------------------
//! https://stackoverflow.com/questions/33893856/area-of-n-sided-planar-polygon-in-3d-space
//!        |       [ (x_(i) - x0) ]   [ (x_((i+1) % N) - x0) ] |
//! area = | sum_i [ (y_(i) - y0) ] x [ (y_((i+1) % N) - y0) ] | / 2
//!        |       [ (z_(i) - z0) ]   [ (z_((i+1) % N) - z0) ] |
//! ----------------------------------------------------------------------------------------
double polygon::area3D_Polygon(const std::vector<polygon::Point> &V)
{
    std::size_t NbPoints = V.size();
    const Point &P0 = V[0];
    double S = 0, Scompx = 0, Scompy = 0, Scompz = 0;
    for(int i=0; i<NbPoints; i++)
    {
        //! -----------------------------------------------------------------------
        //!   | i           j           k        |       | i       j       k     |
        //!   | (x_1-x0)    (y_1-y0)    (z_1-z0) |  =>   | l10_x   l10_y   l10_z |
        //!   | (x_2-x0)    (y_2-y0)    (z_2-z0) |       | l20_x   l20_y   l20_z |
        //! -----------------------------------------------------------------------
        const Point &P1 = V[i];
        const Point &P2 = V[(i+1)%NbPoints];
        double l10_x = P1.x-P0.x; double l10_y = P1.y-P0.y; double l10_z = P1.z-P0.z;
        double l20_x = P2.x-P0.x; double l20_y = P2.y-P0.y; double l20_z = P2.z-P0.z;
        double compx = l10_y*l20_z - l10_z*l20_y;
        double compy = l10_z*l20_x - l10_x*l20_z;
        double compz = l10_x*l20_y - l10_y*l20_x;
        Scompx += compx; Scompy += compy; Scompz += compz;
    }
    S = sqrt(Scompx*Scompx+Scompy*Scompy+Scompz*Scompz);
    return fabs(S);
}

//! ----------------------------------------------------------
//! function: getNormal
//! details:  return (0,0,0) if the normal cannot be computed
//! ----------------------------------------------------------
std::vector<double> polygon::getNormal(const std::vector<polygon::Point> &points)
{
    std::vector<double> zeroNormal{0,0,0};  // null normal
    size_t NbPoints = points.size();
    if(NbPoints<3) return zeroNormal;
    int start = 1;
    const Point &P0 = points[start%NbPoints];
    const Point &P1 = points[(start+1)%NbPoints];
    Point P2;

    for(int i=2; i<NbPoints; i++)
    {
        //P2 = points[i];
        P2 = points[(start+i)%NbPoints];
        if(areCollinear(P0,P1,P2)) continue;
        double x0 = P1.x-P0.x; double y0 = P1.y-P0.y; double z0 = P1.z-P0.z;
        double x1 = P2.x-P0.x; double y1 = P2.y-P0.y; double z1 = P2.z-P0.z;
        double n1 = y0*z1-z0*y1;
        double n2 = z0*x1-x0*z1;
        double n3 = x0*y1-y0*x1;
        double L = sqrt(n1*n1+n2*n2+n3*n3);
        if(L<=std::numeric_limits<double>::epsilon())
        {
            cerr<<"polygon::getNormal()->____normal too small. Returning a zero normal____"<<endl;
            return zeroNormal;
        }
        std::vector<double> normal {n1/L,n2/L,n3/L};
        return normal;
    }
    cerr<<"polygon::getNormal()->____returning a zero normal____"<<endl;
    return zeroNormal;
}

//! ---------------------------
//! function: triangle_area
//! details:  counterclockwise
//! ---------------------------
double polygon::triangle_area (double xa, double ya, double xb, double yb, double xc, double yc)
{
    return 0.5*((xb-xa)*(yc-ya)-(xc-xa)*(yb-ya));
}

//! ------------------------------
//! function: polygon_triangulate
//! details:  in 3D space
//! ------------------------------
bool polygon::polygon_triangulate(const std::vector<Point> &aPolygon, std::vector<std::vector<Point>> &triangles)
{
    int n = int(aPolygon.size());    //! number of points
    if(n<3) return false;

    //! ----------------------------------------
    //! perform local coordinate transformation
    //! 1. select three not collinear points
    //! ----------------------------------------
    const polygon::Point &A = aPolygon[0];
    const polygon::Point &B = aPolygon[1];
    polygon::Point C;
    bool foundNotCollinear = false;
    for(int i=2; i<n; i++)
    {
        C = aPolygon[i];
        if(areCollinear(A,B,C)) continue;
        foundNotCollinear = true;
        break;
    }
    if(foundNotCollinear==false) return false;

    //cout<<"\n____COSINES____"<<endl;
    //! ------------
    //! 2. eX local
    //! ------------
    double eXx = B.x-A.x;
    double eXy = B.y-A.y;
    double eXz = B.z-A.z;
    double lX = sqrt(eXx*eXx+eXy*eXy+eXz*eXz);
    eXx /= lX; eXy /= lX; eXz /= lX;
    //cout<<"____ey ("<<eXx<<", "<<eXy<<", "<<eXz<<")____"<<endl;

    //! ----
    //! C-A
    //! ----
    double xCA = C.x-A.x;
    double yCA = C.y-A.y;
    double zCA = C.z-A.z;

    //! -------------------
    //! 3. eZ local
    //! i       j      k
    //! eXx     eXy    eXz
    //! xCA     yCA    zCA
    //! -------------------
    double eZx = eXy*zCA-eXz*yCA;
    double eZy = eXz*xCA-eXx*zCA;
    double eZz = eXx*yCA-eXy*xCA;
    double lZ = sqrt(eZx*eZx+eZy*eZy+eZz*eZz);
    eZx /= lZ; eZy /= lZ; eZz /= lZ;
    //cout<<"____ez ("<<eZx<<", "<<eZy<<", "<<eZz<<")____"<<endl;

    //! -------------------
    //! 3. eY local
    //! i       j       k
    //! eZx     eZy     eZz
    //! eXx     eXy     eXz
    //! --------------------
    double eYx = eZy*eXz-eZz*eXy;
    double eYy = eZz*eXx-eZx*eXz;
    double eYz = eZx*eXy-eZy*eXx;
    double lY = sqrt(eYx*eYx+eYy*eYy+eYz*eYz);
    eYx /= lY; eYz /= lY; eZy /= lY;
    //cout<<"____ey ("<<eYx<<", "<<eYy<<", "<<eYz<<")____"<<endl;

    //! ---------------------------------------------
    //! definition of the polygon on the local plane
    //! ---------------------------------------------
    //cout<<"\n____POLYGON ON LOCAL PLANE____"<<endl;
    std::vector<double> x,y;
    for(int i=0; i<n; i++)
    {
        const polygon::Point &aP = aPolygon[i];

        //! ---------------------------------
        //! from global to local coordinates
        //! ---------------------------------
        double dxP = aP.x-A.x;
        double dyP = aP.y-A.y;
        double dzP = aP.z-A.z;

        //! ------------------
        //! transformed point
        //! ------------------
        double XP = dxP*eXx+dyP*eXy+dzP*eXz;
        double YP = dxP*eYx+dyP*eYy+dzP*eYz;
        //double ZP = dxP*eZx+dyP*eZy+dzP*eZz;        // should be zero
        x.push_back(XP);
        y.push_back(YP);
        //cout<<"____("<<XP<<", "<<YP<<", "<<ZP<<")____"<<endl;
    }

    //! -------------------------------
    //! triangulate the planar polygon
    //! -------------------------------
    std::vector<int> trianglesDefinitions;
    bool isDone = polygon::polygon_triangulate_planar(x,y,trianglesDefinitions);
    if(!isDone) return false;

    //! ---------------------------------------
    //! return to the global coordinate system
    //! ---------------------------------------
    //cout<<"____"<<trianglesDefinitions.size()<<"____"<<endl;
    //size_t Nt = trianglesDefinitions.size()/3;
    //for(int i=0; i<Nt; i++)
    //{
    //    cout<<"____triangle: "<<i<<"____"<<endl;
    //    int s=3*i;
    //    cout<<"____"<<trianglesDefinitions[s]<<", "<<trianglesDefinitions[s+1]<<", "<<trianglesDefinitions[s+2]<<"____"<<endl;
    //}

    size_t NbTriangles = trianglesDefinitions.size()/3;
    for(int n=0; n<NbTriangles; n++)
    {
        int s = n*3;
        int i0 = trianglesDefinitions[s+0];
        double x0 = x[i0];
        double y0 = y[i0];
        double x0_global = A.x+eXx*x0+eYx*y0;
        double y0_global = A.y+eXy*x0+eYy*y0;
        double z0_global = A.z+eXz*x0+eYz*y0;

        int i1 = trianglesDefinitions[s+1];
        double x1 = x[i1];
        double y1 = y[i1];
        double x1_global = A.x+eXx*x1+eYx*y1;
        double y1_global = A.y+eXy*x1+eYy*y1;
        double z1_global = A.z+eXz*x1+eYz*y1;

        int i2 = trianglesDefinitions[s+2];
        double x2 = x[i2];
        double y2 = y[i2];
        double x2_global = A.x+eXx*x2+eYx*y2;
        double y2_global = A.y+eXy*x2+eYy*y2;
        double z2_global = A.z+eXz*x2+eYz*y2;

        std::vector<polygon::Point> a3DTriangle
        {
            polygon::Point(x0_global,y0_global,z0_global),
                    polygon::Point(x1_global,y1_global,z1_global),
                    polygon::Point(x2_global,y2_global,z2_global)
        };

        //cout<<"____("<<i0<<", "<<i1<<", "<<i2<<")____"<<endl;
        triangles.push_back(a3DTriangle);
    }
    return true;
}

//! -----------------------------------------------------------------------------
//! function: polygon_triangulate_planar
//! details:  Input, int N, the number of vertices.
//!           Input, double X[N], Y[N], the coordinates of each vertex.
//!           Output, int TRIANGLES[3*(N-2)], the triangles of the triangulation
//! -----------------------------------------------------------------------------
bool polygon::polygon_triangulate_planar(const std::vector<double> &x,
                                         const std::vector<double> &y,
                                         std::vector<int> &triangles)
{
    //cout<<"polygon::polygon_triangulate_planar()->____function called____"<<endl;
    size_t n = x.size();
    std::vector<bool> ear(n);
    std::vector<int> prev_node(n), next_node(n);
    triangles = std::vector<int>(3*(n-2));
    //triangles.reserve(3*(n-2));

    int i0, i1, i2, i3, i4;
    int triangle_num;

    //! -------------------------------------------------------------
    //! prev_node and next_node point to the previous and next nodes
    //! -------------------------------------------------------------
    int i = 0;
    prev_node[i] = n - 1;
    next_node[i] = i + 1;

    for (i=1; i<n-1; i++)
    {
        prev_node[i] = i-1;
        next_node[i] = i+1;
    }

    i = n-1;
    prev_node[i] = i-1;
    next_node[i] = 0;

    //! -----------------------------------------------------------------------
    //! EAR indicates whether the node and its immediate neighbors form an ear
    //! that can be sliced off immediately
    //! -----------------------------------------------------------------------
    for (i=0; i<n; i++)
    {
        ear[i] = polygon::diagonal(prev_node[i], next_node[i], prev_node, next_node, x, y);
    }
    triangle_num = 0;
    i2 = 0;

    int tryn=0;
    while (triangle_num<n-3)
    {
        tryn++;
        if(tryn>n-1)    //! avoid entering in an endless loop
        {
            ear.clear();
            next_node.clear();
            prev_node.clear();
            return false;
        }

        //!  If I2 is an ear, gather information necessary to carry out
        //!  the slicing operation and subsequent "healing"
        if (ear[i2])
        {
            i3 = next_node[i2];
            i4 = next_node[i3];
            i1 = prev_node[i2];
            i0 = prev_node[i1];

            //! make vertex I2 disappear
            next_node[i1] = i3;
            prev_node[i3] = i1;

            //! update the earity of I1 and I3, because I2 disappeared.
            ear[i1] = polygon::diagonal(i0, i3, prev_node, next_node, x, y);
            ear[i3] = polygon::diagonal(i1, i4, prev_node, next_node, x, y);

            //! add the diagonal [I3, I1, I2] to the list.
            triangles[0+triangle_num*3] = i3;
            triangles[1+triangle_num*3] = i1;
            triangles[2+triangle_num*3] = i2;
            triangle_num = triangle_num + 1;
        }
        //!  try the next vertex
        i2 = next_node[i2];
    }

    //!  The last triangle is formed from the three remaining vertices.
    i3 = next_node[i2];
    i1 = prev_node[i2];

    triangles[0+triangle_num*3] = i3;
    triangles[1+triangle_num*3] = i1;
    triangles[2+triangle_num*3] = i2;
    triangle_num = triangle_num + 1;

    return true;
}

//! -------------------
//! function: diagonal
//! details:
//! -------------------
bool polygon::diagonal (int im1, int ip1, const std::vector<int> &prev_node, const std::vector<int> &next_node, const std::vector<double> &x, const std::vector<double> &y )
{
    bool value1 = in_cone (im1, ip1, prev_node, next_node, x, y);
    bool value2 = in_cone (ip1, im1, prev_node, next_node, x, y);
    bool value3 = diagonalie (im1,ip1, next_node, x, y);
    bool value = (value1 && value2 && value3);
    return value;
}

//! ---------------------------------------------------------------------------------------
//! function: in_cone
//! details:  IN_CONE is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal
//!           Parameters:
//!
//!           Input, int IM1, IP1, the indices of two vertices.
//!           Input, int N, the number of vertices.
//!           Input, int PREV_NODE[N], the previous neighbor of each vertex
//!           Input, int NEXT_NODE[N], the next neighbor of each vertex
//!           Input, double X[N], Y[N], the coordinates of each vertex
//! ---------------------------------------------------------------------------------------
bool polygon::in_cone(int im1, int ip1, const std::vector<int> &prev_node, const std::vector<int> &next_node, const std::vector<double> &x, const std::vector<double> &y)
{
    bool value;
    int im2 = prev_node[im1];
    int i = next_node[im1];
    double t1 = polygon::triangle_area(x[im1], y[im1], x[i], y[i], x[im2], y[im2]);
    if ( 0.0 <= t1 )
    {
        double t2 = polygon::triangle_area(x[im1], y[im1], x[ip1], y[ip1], x[im2], y[im2]);
        double t3 = polygon::triangle_area(x[ip1], y[ip1], x[im1], y[im1], x[i], y[i]);
        value = ( ( 0.0 < t2 ) && ( 0.0 < t3 ) );
    }
    else
    {
        double t4 = polygon::triangle_area (x[im1], y[im1], x[ip1], y[ip1], x[i], y[i]);
        double t5 = polygon::triangle_area (x[ip1], y[ip1], x[im1], y[im1], x[im2], y[im2]);
        value = ! ((0.0 <= t4) && (0.0 <= t5));
    }
    return value;
}

//! ---------------------
//! function: diagonalie
//! details:
//! ---------------------
bool polygon::diagonalie (int im1, int ip1, const std::vector<int> &next_node, const std::vector<double> &x, const std::vector<double> &y)
{
    int first = im1;
    int j = first;
    int jp1 = next_node[first];

    bool value = true;
    bool value2;

    //! For each edge VERTEX(J):VERTEX(JP1) of the polygon:
    while (1)
    {
        //! Skip any edge that includes vertex IM1 or IP1
        if (j == im1 || j == ip1 || jp1 == im1 || jp1 == ip1)
        {
            //! do nothing
        }
        else
        {
            value2 = polygon::intersect (x[im1], y[im1], x[ip1], y[ip1], x[j], y[j], x[jp1], y[jp1]);
            if(value2)
            {
                value = false;
                break;
            }
        }
        j = jp1;
        jp1 = next_node[j];
        if (j == first) break;
    }
    return value;
}

//!---------------------------------------------------------------------
//! function: intersect
//! details:  INTERSECT is true if lines VA:VB and VC:VD intersect
//!           Parameters:
//!           Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y
//!           coordinates of the four vertices
//!---------------------------------------------------------------------
bool polygon::intersect(double xa, double ya, double xb, double yb, double xc,
                        double yc, double xd, double yd)
{
    bool value;
    if(polygon::intersect_prop(xa, ya, xb, yb, xc, yc, xd, yd)) value = true;
    else if(polygon::between(xa, ya, xb, yb, xc, yc)) value = true;
    else if(polygon::between(xa, ya, xb, yb, xd, yd)) value = true;
    else if(polygon::between(xc, yc, xd, yd, xa, ya)) value = true;
    else if(polygon::between(xc, yc, xd, yd, xb, yb)) value = true;
    else value = false;
    return value;
}

//! --------------------
//! function: collinear
//! details:
//! --------------------
bool polygon::collinear(double xa, double ya, double xb, double yb, double xc, double yc)
{
    const double epsZero = std::numeric_limits<double>::epsilon();
    bool value;

    double area = triangle_area ( xa, ya, xb, yb, xc, yc );
    double side_ab_sq = pow(xa-xb,2) + pow(ya-yb,2);
    double side_bc_sq = pow(xb-xc,2) + pow(yb-yc,2);
    double side_ca_sq = pow(xc-xa,2) + pow(yc-ya,2);
    double side_max_sq = std::max<double>(side_ab_sq, std::max<double>(side_bc_sq, side_ca_sq));

    if (side_max_sq <= epsZero) value = true;
    else if (2.0*fabs(area) <= epsZero * side_max_sq) value = true;
    else value = false;
    return value;
}

//! --------------------------------------------------------------------------------
//! function: between
//! details:  return TRUE if the vertex C is between vertices A and B.
//!           The points must be (numerically) collinear.
//!           Given that condition, we take the greater of XA - XB and YA - YB
//!           as a "scale" and check where C's value lies.
//!           Parameters:
//!           Input, double XA, YA, XB, YB, XC, YC, the coordinates of the vertices
//! --------------------------------------------------------------------------------
bool polygon::between (double xa, double ya, double xb, double yb, double xc, double yc)
{
    bool value;
    double xmax, xmin, ymax, ymin;
    if (!polygon::collinear(xa, ya, xb, yb, xc, yc)) value = false;
    else if (fabs(ya-yb) < fabs(xa-xb))
    {
        xmax = std::max(xa, xb);
        xmin = std::min(xa, xb);
        value = (xmin <= xc && xc <= xmax);
    }
    else
    {
        ymax = std::max<double>(ya, yb);
        ymin = std::min(ya, yb);
        value = (ymin <= yc && yc <= ymax);
    }
    return value;
}

//!--------------------------------------------------------------------------------------
//! function: intesect_prop
//! details:  INTERSECT_PROP is TRUE if lines VA:VB and VC:VD have a proper intersection
//!           Parameters:
//!           Input, double XA, YA, XB, YB, XC, YC, XD, YD, the X and Y
//!           coordinates of the four vertices.
//!           Output, bool INTERSECT_PROP, the result of the test
//!--------------------------------------------------------------------------------------
bool polygon::intersect_prop(double xa, double ya, double xb, double yb, double xc, double yc, double xd, double yd)
{
    bool value;
    if(polygon::collinear(xa, ya, xb, yb, xc, yc)) value = false;
    else if(polygon::collinear(xa, ya, xb, yb, xd, yd)) value = false;
    else if(polygon::collinear(xc, yc, xd, yd, xa, ya)) value = false;
    else if(polygon::collinear(xc, yc, xd, yd, xb, yb)) value = false;
    else
    {
        double t1 = polygon::triangle_area(xa, ya, xb, yb, xc, yc);
        double t2 = polygon::triangle_area(xa, ya, xb, yb, xd, yd);
        double t3 = polygon::triangle_area(xc, yc, xd, yd, xa, ya);
        double t4 = polygon::triangle_area(xc, yc, xd, yd, xb, yb);
        bool value1 = (0.0 < t1);
        bool value2 = (0.0 < t2);
        bool value3 = (0.0 < t3);
        bool value4 = (0.0 < t4);
        value = (XOR(value1, value2)) && (XOR(value3, value4));
    }
    return value;
}
