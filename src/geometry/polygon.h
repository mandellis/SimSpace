#ifndef POLYGON_H
#define POLYGON_H

#define SIGN(x) ((0 < x) - (x < 0))

#include <vector>
#include <iostream>
using namespace std;

class polygon
{

public:

    //! ------------------------
    //! a Point in the 2D space
    //! ------------------------
    struct Point2D
    {
        double x,y;
        Point2D(double x0 = 0, double y0 = 0):x(x0),y(y0){;}
        Point2D(const Point2D &rhs){ x=rhs.x; y = rhs.y; }
        bool operator == (const Point2D &rhs) { if(x==rhs.x && y==rhs.y) return true; return false; }
        Point2D operator = (const Point2D &rhs) {x=rhs.x; y = rhs.y; return *this; }
    };

    //! ------------------------
    //! a Point in the 3D space
    //! ------------------------
    struct Point
    {
        double x, y, z;
        Point(double x0=0, double y0=0, double z0=0):x(x0),y(y0),z(z0){;}
        Point(const Point &other) {x = other.x; y = other.y; z = other.z; }
        Point operator = (const Point &other) { x = other.x; y = other.y; z = other.z; return *this; }
        bool operator == (const Point &other) { if(x == other.x && y == other.y && z == other.z) return true; return false; }
    };

    //! -------------------------
    //! function: triangleArea
    //! details:  absolute value
    //! -------------------------
    static double triangleArea(const Point &P0, const Point &P1, const Point &P2)
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

    //! ----------------------------------------
    //! function: areCollinear
    //! details:  uses the area of the triangle
    //! ----------------------------------------
    static bool areCollinear(const Point &P0, const Point &P1, const Point &P2)
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
    //!
    //!        |       [ (x_(i) - x0) ]   [ (x_((i+1) % N) - x0) ] |
    //! area = | sum_i [ (y_(i) - y0) ] x [ (y_((i+1) % N) - y0) ] | / 2
    //!        |       [ (z_(i) - z0) ]   [ (z_((i+1) % N) - z0) ] |
    //! ----------------------------------------------------------------------------------------
    static double area3D_Polygon(const std::vector<Point> &V)
    {
        std::size_t NbPoints = V.size();
        Point P0 = V.at(0);
        double S = 0;
        for(int i=0; i<NbPoints; i++)
        {
            Point P1 = V.at(i);
            Point P2 = V.at((i+1)%NbPoints);
            S = S + triangleArea(P0,P1,P2);
        }
        return fabs(S);
    }

    //! ----------------------------------------------------------
    //! function: getNormal
    //! details:  return (0,0,0) if the normal cannot be computed
    //! ----------------------------------------------------------
    static std::vector<double> getNormal(const std::vector<Point> &points)
    {
        std::vector<double> zeroNormal{0,0,0};
        size_t NbPoints = points.size();
        if(NbPoints<3)
        {
            cout<<"polygon::getNormal()->____wrong number of nodes____"<<endl;
            return zeroNormal;
        }
        int start = 1;
        Point P0 = points[start%NbPoints];
        Point P1 = points[(start+1)%NbPoints];
        Point P2;

        for(int i=2; i<NbPoints; i++)
        {
            //P2 = points[i];
            P2 = points[(start+i)%NbPoints];
            if(areCollinear(P0,P1,P2)) continue;
            double x0 = P1.x-P0.x; double y0 = P1.y-P0.y; double z0 = P1.z-P0.z;
            double x1 = P2.x-P0.x; double y1 = P2.y-P0.y; double z1 = P2.z-P0.z;

            double n1= y0*z1-z0*y1; double n2 = z0*x1-x0*z1; double n3 = x0*y1-y0*x1;
            double L = sqrt(n1*n1+n2*n2+n3*n3);
            if(L<1e-10)
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

    //! ----------------------------------------
    //! function: isPointBetween
    //! details:  check if C is between A and B
    //! ----------------------------------------
    static bool isPointBetween(const polygon::Point &A, const polygon::Point &B, const polygon::Point &C, const double tolerance=1e-6)
    {
        //! -----------------
        //! check alignement
        //! -----------------
        double ABx = A.x-B.x; double ABy = A.y-B.y; double ABz = A.z-B.y;
        double ACx = A.x-C.x; double ACy = A.y-C.y; double ACz = A.z-C.y;

        double llx = ABy*ACz-ACy*ABz;
        double lly = ABz*ACy-ABy*ACz;
        double llz = ABx*ACy-ABy*ACx;

        double L = sqrt(llx*llx+lly*lly+llz*llz);
        if(L>tolerance) return false;                    //! the points are not aligned

        double KAC = ABx*ACx+ABy*ACy+ABz*ACz;
        double KAB = ABx*ABx+ABy*ABy+ABz*ABz;

        //! ------
        //! rules
        //! ------
        //! if(KAC<0) return false;
        //! if(KAC>KAB) return false;
        //! if(KAC==0) return true;
        //! if(KAC==KAB) return true;
        //! if(KAC>0 && KAC<KAB) return true;

        if(KAC<tolerance) return false;
        if(KAC-KAB>tolerance) return false;
        if(fabs(KAC-tolerance)<=0) return true;
        if(fabs(KAC-KAB)<0) return true;
        if((KAC-tolerance)>=0 && fabs(KAB-KAC)>=tolerance) return true;
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
    static bool isPointInside(const Point &P,
                              const std::vector<Point> &V,
                              bool normalizeCrossProducts,
                              double inOutTolerance=1e-5)
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
            if(dot<-inOutTolerance) return false;
        }
        return true;
    }

    //! -------------------------------------------------------------------
    //! function: planeCoefficients
    //! details:  input => a polygon; output => polygon plane coefficients
    //! -------------------------------------------------------------------
    static bool planeCoefficients(const std::vector<Point> &points, double &a, double &b, double &c, double &d)
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

    //! -----------------------------
    //! function: pointPlaneDistance
    //! details:  signed distance
    //! -----------------------------
    static double pointPlaneDistance(const polygon::Point &P, double a, double b,double c, double d)
    {
        return (a*P.x+b*P.y+c*P.z+d)/sqrt(a*a+b*b+c*c);
    }

    //! --------------------------------
    //! function: isPointCloseToPolygon
    //! details:
    //! --------------------------------
    static bool isPointCloseToPolygon(polygon::Point P, std::vector<polygon::Point> aPolygon, double proximity, int mode = 0)
    {
        switch(mode)
        {
        case 0:
        {
            //! ------------------------------------------------------
            //! all polygon points should be close to the input point
            //! ------------------------------------------------------
            size_t NbPoints = aPolygon.size();
            for(int n=0; n<NbPoints; n++)
            {
                const polygon::Point aP = aPolygon[n];
                double d = sqrt(pow(aP.x-P.x,2)+pow(aP.y-P.y,2)+pow(aP.z-P.z,2));
                if(d>proximity) return false;
            }
            return true;
        }
            break;

        case 1:
        {
            //! -----------------------------------------------------
            //! at least one point of the polygon should be close to
            //! the input point for returning true
            //! -----------------------------------------------------
            size_t NbPoints = aPolygon.size();
            for(int n=0; n<NbPoints; n++)
            {
                const polygon::Point &aP = aPolygon.at(n);
                double d = sqrt(pow(aP.x-P.x,2)+pow(aP.y-P.y,2)+pow(aP.z-P.z,2));
                if(d<=proximity) return true;
            }
            return false;
        }
            break;

        case 2:
        {
            //! --------------------------------------------------------------------
            //! use the polygon center for evaluating the point to polygon distance
            //! --------------------------------------------------------------------
            polygon::Point C;
            double area = 0;
            bool isDone = polygon::getPolygonCenter(aPolygon,C,area);
            if(!isDone) return false;
            double d = sqrt(pow(C.x-P.x,2)+pow(C.y-P.y,2)+pow(C.z-P.z,2));
            if(d>proximity) return false;
            return true;
        }
            break;
        }
    }

    //! ----------------------------
    //! function: getRotationMatrix
    //! details:
    //! ----------------------------
    static bool getRotationMatrix(const std::vector<polygon::Point> &aPolygon,
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
            //P3 = aPolygon.at(n);
            P3 = aPolygon[n];

            bool aligned = areCollinear(P1,P2,P3);
            if(aligned == true)
            {
                //cout<<"getRotationMatrix()->____third not found____"<<endl;
            }
            else
            {
                //cout<<"getRotationMatrix()->____third found____"<<endl;
                found = true;
                break;
            }
        }
        if(!found) return false;

        //cout<<"getRotationMatrix()->____("<<P1.x<<", "<<P1.y<<", "<<P1.z<<")____"<<endl;
        //cout<<"getRotationMatrix()->____("<<P2.x<<", "<<P2.y<<", "<<P2.z<<")____"<<endl;
        //cout<<"getRotationMatrix()->____("<<P3.x<<", "<<P3.y<<", "<<P3.z<<")____"<<endl;

        if(!found)
        {
            //cout<<"getRotationMatrix()->____collinear____"<<endl;
            return false;
        }
        //cout<<"getRotationMatrix()->____CAN calculate a plane____"<<endl;

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

        /*
        rotationMatrix[0]=a11;
        rotationMatrix[1]=a12;
        rotationMatrix[2]=a13;
        rotationMatrix[3]=a21;
        rotationMatrix[4]=a22;
        rotationMatrix[5]=a23;
        rotationMatrix[6]=a31;
        rotationMatrix[7]=a32;
        rotationMatrix[8]=a33;
        */

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


    //! ---------------------------
    //! function: getPolygonCenter
    //! details:
    //! ---------------------------
    static bool getPolygonCenter(const std::vector<polygon::Point> &aPolygon,
                                 polygon::Point &faceCenter, double &faceArea)
    {
        //! ------------------------------------------
        //! the polygon is actually a triangle
        //! direct (faster) calculation of the center
        //! ------------------------------------------
        if(aPolygon.size()==3)
        {
            faceCenter.x = (aPolygon.at(0).x+aPolygon.at(1).x+aPolygon.at(2).x)/3.0;
            faceCenter.y = (aPolygon.at(0).y+aPolygon.at(1).y+aPolygon.at(2).y)/3.0;
            faceCenter.z = (aPolygon.at(0).z+aPolygon.at(1).z+aPolygon.at(2).z)/3.0;
            faceArea = polygon::triangleArea(aPolygon.at(0),aPolygon.at(1),aPolygon.at(2));
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
            exit(100);
            return false;
        }

        double a11 = rotationMatrix.at(0);
        double a12 = rotationMatrix.at(1);
        double a13 = rotationMatrix.at(2);
        double a21 = rotationMatrix.at(3);
        double a22 = rotationMatrix.at(4);
        double a23 = rotationMatrix.at(5);
        double a31 = rotationMatrix.at(6);
        double a32 = rotationMatrix.at(7);
        double a33 = rotationMatrix.at(8);

        const polygon::Point P1 = aPolygon.at(0);
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
            const polygon::Point &curPoint = aPolygon.at(n);
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
            const polygon::Point &Pi = transfPointList.at(i);
            const polygon::Point &Pii = transfPointList.at(i+1);

            //const polygon::Point &Pi = transfPointList[i];
            //const polygon::Point &Pii = transfPointList[i+1];

            //Area += Pi.x*Pii.y-Pii.x*Pi.x;    // check code
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
            const polygon::Point &Pi = transfPointList[0];
            const polygon::Point &Pii = transfPointList[1];

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
        //cout<<"____polygon center("<<Cx<<", "<<Cy<<", "<<Cz<<")____"<<endl;
        faceArea = Area;
        return true;
    }

    //! ---------------------------------
    //! function: polygonIntersectsPlane
    //! details:
    //! ---------------------------------
    static bool testPolygonPlaneIntersection(const std::vector<polygon::Point> &aPolygon, double a, double b, double c, double d/*, int &NbPoints*/)
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
            if(d != dp)
            {
                //! -------------------------------------------------
                //! the polygon plane is parallel to the input plane
                //! -------------------------------------------------
                return false;
            }
            else
            {
                //! ------------------------------
                //! the polygon lays on the plane
                //! ------------------------------
                return true;
            }
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
};

#endif // POLYGON_H
