#include "polyhedron.h"
#include "limits"
using namespace polyhedron;

//! do not bring away this macros from .cpp
#define det3(a11,a12,a13,a21,a22,a23,a31,a32,a33) (a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31))
#define EPS 10.0*std::numeric_limits<double>::epsilon()

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
polyhedron::cell::cell(const std::vector<face> &faces):
    myFaces(faces)
{
    //! ---------------------------------
    //! compute polyhedron triangulation
    //! ---------------------------------
    for(int faceNr=0; faceNr<myFaces.size(); faceNr++)
    {
        //! triangulation of a face
        const std::vector<polygon::Point> &aFace = this->getFace(faceNr);
        std::vector<std::vector<polygon::Point>> aFaceTriangulation;
        polygon::polygon_triangulate(aFace,aFaceTriangulation);

        //! pile up the triangles into the overall triangulation
        for(int i=0; i<aFaceTriangulation.size(); i++) myTriangulation.push_back(aFaceTriangulation[i]);
    }
}

//! ----------------------
//! function: getCentroid
//! details:
//! ----------------------
bool polyhedron::cell::getCentroid(double &xc, double &yc, double &zc) const
{
    //! ----------------------------------------
    //! calculate the center of the element
    //! working on the polyhedron triangulation
    //! ----------------------------------------
    double Area_i = 0.0, AreaTot = 0.0;
    double AiRi_x = 0.0, AiRi_y = 0.0, AiRi_z = 0.0;

    size_t NbTriangles = myTriangulation.size();
    for(int j=0; j<NbTriangles; j++)
    {
        const face &aTriangle = myTriangulation[j];

        //! the three points defining the triangle
        double xAi = aTriangle[0].x; double yAi = aTriangle[0].y; double zAi = aTriangle[0].z;
        double xBi = aTriangle[1].x; double yBi = aTriangle[1].y; double zBi = aTriangle[1].z;
        double xCi = aTriangle[2].x; double yCi = aTriangle[2].y; double zCi = aTriangle[2].z;

        std::vector<polygon::Point> trig { polygon::Point(xAi,yAi,zAi),polygon::Point(xBi,yBi,zBi),polygon::Point(xCi,yCi,zCi) };
        Area_i = polygon::area3D_Polygon(trig);
        AreaTot += Area_i;

        double Ri_x = (xAi+xBi+xCi)/3.0;
        double Ri_y = (yAi+yBi+yCi)/3.0;
        double Ri_z = (zAi+zBi+zCi)/3.0;

        AiRi_x += Ri_x * Area_i;
        AiRi_y += Ri_y * Area_i;
        AiRi_z += Ri_z * Area_i;
    }
    if(AreaTot>EPS)
    {
        xc = AiRi_x/AreaTot;
        yc = AiRi_y/AreaTot;
        zc = AiRi_z/AreaTot;
        return true;
    }
    else return false;
}

//! --------------------
//! function: getVolume
//! details:
//! --------------------
double polyhedron::cell::getVolume() const
{
    double vol = 0;
    polygon::Point P;
    this->getCentroid(P.x,P.y,P.z);
    for(int i=0; i<myTriangulation.size(); i++)
    {
        std::vector<polygon::Point> aTet = myTriangulation[i];
        aTet.push_back(P);
        double x1 = aTet[0].x; double y1 = aTet[0].y; double z1 = aTet[0].z;
        double x2 = aTet[1].x; double y2 = aTet[1].y; double z2 = aTet[1].z;
        double x3 = aTet[2].x; double y3 = aTet[2].y; double z3 = aTet[2].z;
        double x4 = aTet[3].x; double y4 = aTet[3].y; double z4 = aTet[3].z;

        //! --------------
        //! x1  y1  z1  1
        //! x2  y2  z2  1
        //! x3  y3  z3  1
        //! x4  y4  z4  1
        //! --------------
        double volT = (1/6.0)*(
                    +x1*det3(y2, z2, 1, y3, z3, 1, y4, z4, 1)
                    -y1*det3(x2, z2, 1, x3, z3, 1, x4, z4, 1)
                    +z1*det3(x2, y2, 1, x3, y3, 1, x4, y4, 1)
                    -1.0*det3(x2, y2, z2, x3, y3, z3, x4, y4, z4));
        vol += volT;
    }
    return vol;
}
