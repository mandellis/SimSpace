#ifndef MESHELEMENTBYCOORDS_H
#define MESHELEMENTBYCOORDS_H

//! ------------
//! a11 a12 a13
//! a21 a22 a23
//! a31 a32 a33
//! ------------
#define DET3(a11,a12,a13,a21,a22,a23,a31,a32,a33) a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31)

//! ---
//! Qt
//! ---
#include <QList>

//! ----------------
//! custom includes
//! ----------------
#include <mesh.h>
#include <elementtypes.h>
#include <polygon.h>

struct meshElementByCoords
{
    int ID;
    ElemType type;
    QList<mesh::meshPoint> pointList;
    
    //! --------------------
    //! default constructor
    //! --------------------
    meshElementByCoords()
    {
        type = NULL_ELEMENT;
        pointList = QList<mesh::meshPoint>();
    }

    //! ----------------------------------
    //! compute jacobian of a tetrahedron
    //! ----------------------------------
    double jacobianTet()
    {
        //! ---------------
        //! 1   1   1   1
        //! x1  x2  x3  x4
        //! y1  y2  y3  y4
        //! z1  z2  z3  z4
        //! ---------------
        double x1 = pointList.at(0).x; double y1 = pointList.at(0).y; double z1 = pointList.at(0).z;
        double x2 = pointList.at(1).x; double y2 = pointList.at(1).y; double z2 = pointList.at(1).z;
        double x3 = pointList.at(2).x; double y3 = pointList.at(2).y; double z3 = pointList.at(2).z;
        double x4 = pointList.at(3).x; double y4 = pointList.at(3).y; double z4 = pointList.at(3).z;
        double J =  DET3(x2,x3,x4,y2,y3,y4,z2,z3,z4);
        J -= DET3(x1,x3,x4,y1,y3,y4,z1,z3,z4);
        J += DET3(x1,x2,x4,y1,y2,y4,z1,z2,z4);
        J -= DET3(x1,x2,x3,y1,y2,y3,z1,z2,z3);
        return J;
    }

    //! ------------------------------------
    //! function: checkPyramidBasePlanarity
    //! details:  returns is base flat?
    //! ------------------------------------
    bool checkPyramidBasePlanarity()
    {
        double a,b,c,d;
        const mesh::meshPoint &A0 = pointList.at(1);
        const mesh::meshPoint &A1 = pointList.at(2);
        const mesh::meshPoint &A2 = pointList.at(3);
        const mesh::meshPoint &A3 = pointList.at(4);
        std::vector<polygon::Point> meshTriangle;
        meshTriangle.push_back(polygon::Point(A0.x,A0.y,A0.z));
        meshTriangle.push_back(polygon::Point(A1.x,A1.y,A1.z));
        meshTriangle.push_back(polygon::Point(A2.x,A2.y,A2.z));
        polygon::planeCoefficients(meshTriangle,a,b,c,d);

        polygon::Point P(A3.x,A3.y,A3.z);
        if(a*P.x+b*P.y+c*P.z+d==0.0) return true;

        //! ------------------------------------------------
        //! the nodes are not coplanar
        //! correct the position of A3 using the projection
        //! ------------------------------------------------
        //cout<<"____the base points are not coplanar: ";

        double t = polygon::pointPlaneDistance(P,a,b,c,d);
        double l = sqrt(a*a+b*b+c*c);
        a /= l;  b /= l; c /= l;
        double xP = P.x-a*t;
        double yP = P.y-b*t;
        double zP = P.z-c*t;

        //cout<<std::scientific;
        //cout<<" d = "<<t<<" correcting point from ("<<P.x<<", "<<P.y<<", "<<P.z<<") to ("<<xP<<", "<<yP<<", "<<zP<<")____"<<endl;
        pointList[4].x = xP;
        pointList[4].y = yP;
        pointList[4].z = zP;
        return false;
    }

    //! ------------------------
    //! checkPyramidBaseWarping
    //! ------------------------
    bool checkPyramidBaseWarping()
    {
        const mesh::meshPoint & A0 = pointList.at(1);
        const mesh::meshPoint & A1 = pointList.at(2);
        const mesh::meshPoint & A2 = pointList.at(3);
        const mesh::meshPoint & A3 = pointList.at(4);

        double Area1 = polygon::triangleArea(polygon::Point(A0.x,A0.y,A0.z),
                                             polygon::Point(A1.x,A1.y,A1.z),
                                             polygon::Point(A2.x,A2.y,A2.z));

        double Area2 = polygon::triangleArea(polygon::Point(A0.x,A0.y,A0.z),
                                             polygon::Point(A2.x,A2.y,A2.z),
                                             polygon::Point(A3.x,A3.y,A3.z));

        if(Area1*Area2>0.0) return false;
        else
        {
            cout<<"____invalid pyramid base____"<<endl;
            return true;
        }
    }

    //! --------------------------------------------------------------
    //! constructor of a prismatic element
    //! face1 {1,2,3,4} is the base face
    //! face2 {1',2',3',4'} is the inflated face
    //! the element construction takes into account the locked nodes.
    //! Degenerated elements are created for some combination of
    //! locked nodes
    //!
    //!              normal
    //!       2'/|--/\----/| 3'
    //!        / |  ||   / |      "Extruded" meshElement2D
    //!     1'/__|______/ 4'
    //!       |  |      |  |
    //!       |  |      |  |
    //!       | 2/------|--/ 3
    //!       | /       | /       "Base" meshElement2D
    //!     1 |/________|/ 4
    //!            ||
    //!            \/ normal
    //! --------------------------------------------------------------
    meshElementByCoords(const QList<mesh::meshPoint> &face1, const QList<mesh::meshPoint> &face2)
    {
        //! --------------------
        //! add the mesh points
        //! --------------------
        QList<int> listOfIndexToDiscard;
        for(int i=0; i<face1.length(); i++)
        {
            //! ---------------------------------------
            //! points of the face 1 and of the face 2
            //! ---------------------------------------
            mesh::meshPoint pointFace1 = face1.at(i);
            mesh::meshPoint pointFace2 = face2.at(i);

            //!cout<<"____face 1("<<pointFace1.x<<", "<<pointFace1.y<<", "<<pointFace1.z<<")____"<<endl;
            //!cout<<"____face 2("<<pointFace2.x<<", "<<pointFace2.y<<", "<<pointFace2.z<<")____"<<endl;

            if(pointFace1==pointFace2)
            {
                listOfIndexToDiscard<<i;
                //cout<<"____found duplicated mesh point____"<<endl;
            }

            pointFace1.ID = i+1;
            pointList<<pointFace1;
        }
        int NbPoints = pointList.length();
        int k=0;
        for(int i=0; i<face2.length(); i++)
        {
            if(!listOfIndexToDiscard.contains(i))
            {
                k++;
                mesh::meshPoint pointFace2 = face2.at(i);
                //cout<<"____("<<pointFace2.x<<", "<<pointFace2.y<<", "<<pointFace2.z<<")____"<<endl;
                pointFace2.ID = NbPoints+k;
                pointList<<pointFace2;
            }
        }

        //! ---------------------
        //! set the element type
        //! build the element
        //! ---------------------
        NbPoints = pointList.length();

        //! -----------------------------------
        //! the base face is a linear triangle
        //! -----------------------------------
        if(face1.length()==3)
        {
            switch(NbPoints)
            {
            case 3:
            {
                type = NULL_ELEMENT;
            }
                break;

            case 4:
            {
                type = TET;
            }
                break;

            case 5:
            {
                type = PYRAM;
                switch (listOfIndexToDiscard.at(0))
                {
                case 0:
                {
                    QList<mesh::meshPoint> pl; pl.append(pointList);
                    pointList.clear();
                    pointList<<pl.at(0)<<pl.at(3)<<pl.at(4)<<pl.at(2)<<pl.at(1);
                }
                    break;
                case 1:
                {
                    QList<mesh::meshPoint> pl; pl.append(pointList);
                    pointList.clear();
                    pointList<<pl.at(1)<<pl.at(4)<<pl.at(3)<<pl.at(0)<<pl.at(2);
                    }
                    break;
                case 2:
                {
                    QList<mesh::meshPoint> pl; pl.append(pointList);
                    pointList.clear();
                    pointList<<pl.at(2)<<pl.at(3)<<pl.at(4)<<pl.at(1)<<pl.at(0);
                }
                    break;
                }
            }
                break;

            case 6:
            {
                type = PRISM;
                QList<mesh::meshPoint> pl; pl.append(pointList);
                pointList.clear();
                pointList<<pl.at(2)<<pl.at(1)<<pl.at(0)<<pl.at(5)<<pl.at(4)<<pl.at(3);      //old working with calculix
            }
                break;
                //! --------------------------------------
                //! other elements (2nd order); to do ...
                //! --------------------------------------
            }
        }
    }

    //! ----------------------
    //! operator = assignment
    //! ----------------------
    meshElementByCoords operator = (const meshElementByCoords &rhs)
    {
        ID = rhs.ID;
        type = rhs.type;
        for(int n=0; n<rhs.pointList.length(); n++) pointList<<rhs.pointList.at(n);
        return *this;
    }

    //! ---------------------
    //! operator == identity
    //! "strong" identity
    //! ---------------------
    bool operator == (const meshElementByCoords &rhs)
    {
        if(ID != rhs.ID) return false;
        if(type != rhs.type) return false;
        //! to be implemented; to do...
    }

    //! --------------------
    //! function: getPoints
    //! details:
    //! --------------------
    std::vector<mesh::meshPoint> getPoints() const
    {
        return pointList.toVector().toStdVector();
    }
};

#endif // MESHELEMENTBYCOORDS_H
