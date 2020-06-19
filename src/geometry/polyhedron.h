#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include <polygon_triangulate.h>
#include <polygon.h>
#include <vector>
using namespace std;

struct polyhedron
{
    std::vector<std::vector<polygon::Point>> faces;

    //! -------------
    //! constructors
    //! -------------
    polyhedron() { ; }
    polyhedron(std::vector<std::vector<polygon::Point>> faces)
    {
        int NbFaces = int(faces.size());
        for(int i=0; i<NbFaces; i++)
        {
            faces.push_back(faces.at(i));
        }
    }

    //! -----------------
    //! copy constructor
    //! -----------------
    polyhedron(const polyhedron &other)
    {
        int NbFaces = int(other.faces.size());
        for(int i=0; i<NbFaces; i++)
        {
            const std::vector<polygon::Point> &aFace = other.faces.at(i);
            faces.push_back(aFace);
        }
    }

    //! ---------------------------
    //! getPolyhedronTriangulation
    //! ---------------------------
    std::vector<std::vector<polygon::Point>> getPolyhedronTriangulation()
    {
        //! -----------------------------
        //! triangulation of the polygon
        //! -----------------------------
        std::vector<std::vector<polygon::Point>> polygonTriangulation;

        //! -----------------------------------------
        //! scan the faces and triangulate each face
        //! -----------------------------------------
        int NbFaces = int(faces.size());
        for(int i=0; i<NbFaces; i++)
        {
            std::vector<polygon::Point> curFace = faces[i];
            int NbPoints = int(curFace.size());

            //! ------------------------------------
            //! get the rotation matrix of the face
            //! ------------------------------------
            std::vector<double> rm;
            polygon::getRotationMatrix(curFace,rm);
            double a11 = rm[0], a12 = rm[1], a13 = rm[2];
            double a21 = rm[3], a22 = rm[4], a23 = rm[5];
            double a31 = rm[6], a32 = rm[7], a33 = rm[8];

            //! ------------------------
            //! first point of the face
            //! ------------------------
            polygon::Point P1 = curFace[0];

            //! ----------------------------------
            //! transform the points of the faces
            //! from global to local coordinates
            //! ----------------------------------
            std::vector<polygon::Point> face_t;
            for(int j=0; j<NbPoints; j++)
            {
                polygon::Point curPoint = curFace.at(j);

                double dxP = curPoint.x-P1.x;
                double dyP = curPoint.y-P1.y;
                double dzP = curPoint.z-P1.z;

                double XP = dxP*a11+dyP*a12+dzP*a13;
                double YP = dxP*a21+dyP*a22+dzP*a23;
                double ZP = dxP*a31+dyP*a32+dzP*a33;

                polygon::Point Pt(XP,YP,ZP);
                face_t.push_back(Pt);
            }

            //! ---------------------
            //! triangulate the face
            //! ---------------------
            double *x = new double[NbPoints];
            double *y = new double[NbPoints];
            for(int j=0; j<NbPoints; j++)
            {
                polygon::Point aP = face_t.at(j);
                x[j] = aP.x; y[j] = aP.y;
            }
            int *triangles = new int[3*(NbPoints-2)];
            polygon_triangulate(NbPoints,x,y,triangles);

            //! ---------------
            //! transform back
            //! ---------------
            int NbTriads = NbPoints-2;
            for(int j=0; j<NbTriads; j++)
            {
                std::vector<polygon::Point> aTriangle;
                int s = 3*j;
                for(int n=0; n<3; n++)
                {
                    int i=triangles[s+n];
                    double X = x[i];
                    double Y = y[i];
                    double xP = P1.x+a11*X+a21*Y;
                    double yP = P1.y+a12*X+a22*Y;
                    double zP = P1.z+a13*X+a23*Y;
                    polygon::Point P(xP,yP,zP);
                    aTriangle.push_back(P);
                    polygonTriangulation.push_back(aTriangle);
                }
            }
            //polygonTriangulation.push_back(aTriangle);
        }
        return polygonTriangulation;
    }

    static double det3(double a11, double a12, double a13,
                       double a21, double a22, double a23,
                       double a31, double a32, double a33)
    {
        double val = a11*(a22*a33-a23*a32)+
                a12*(a23*a31-a21*a33)+
                a13*(a21*a32-a22*a31);
        return val;
    }

    //! ----------------------
    //! volume of tetrahedron
    //! ----------------------
    static double tetrahedronVolume(std::vector<polygon::Point> points)
    {
        double x1 = points[0].x; double y1 = points[0].y; double z1 = points[0].z;
        double x2 = points[1].x; double y2 = points[1].y; double z2 = points[1].z;
        double x3 = points[2].x; double y3 = points[2].y; double z3 = points[2].z;
        double x4 = points[3].x; double y4 = points[3].y; double z4 = points[3].z;

        //! --------------
        //! x1  y1  z1  1
        //! x2  y2  z2  1
        //! x3  y3  z3  1
        //! x4  y4  z4  1
        //! --------------
        double vol = x1*det3(y2,z2,1,y3,z3,1,y4,z4,1)-
                y1*det3(x2,z2,1,x3,z3,1,x4,z4,1)+
                z1*det3(x2,y2,1,x3,y3,1,x4,y4,1)-
                1*det3(x2,y2,z2,x3,y3,z3,x4,y4,z4);
        return vol;
    }

    //! -------
    //! volume
    //! -------
    double volume()
    {
        double volume = 0;
        std::vector<std::vector<polygon::Point>> triangulation = getPolyhedronTriangulation();
        int NbTriangles = int(triangulation.size());
        for(int i=0; i<NbTriangles; i++)
        {
            std::vector<polygon::Point> aTriangle = triangulation.at(i);
            polygon::Point O(0,0,0);
            aTriangle.push_back(O);
            double dV = tetrahedronVolume(aTriangle);
            volume = volume + dV;
        }
        return volume;
    }
};

#endif // POLYHEDRON_H
