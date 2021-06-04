#include "tetsplitter.h"

//! -----------------------
//! function: pointBetween
//! details:  helper
//! -----------------------
mesh::meshPoint pointBetween(const mesh::meshPoint &A, const mesh::meshPoint &B)
{
    double x = (A.x+B.x)/2;
    double y = (A.y+B.y)/2;
    double z = (A.z+B.z)/2;
    return mesh::meshPoint(x,y,z);
}

//! ------------------------------------------------
//! function: splitTet
//! details:  split a tetrahedron in four hexahedra
//! ------------------------------------------------
std::vector<meshElementByCoords> tetSplitter::splitTet(const meshElementByCoords &aTet)
{
    //! ----------------------
    //! define points 1,2,3,4
    //! ----------------------
    mesh::meshPoint P1 = aTet.pointList[0];
    mesh::meshPoint P2 = aTet.pointList[1];
    mesh::meshPoint P3 = aTet.pointList[2];
    mesh::meshPoint P4 = aTet.pointList[3];

    double xC,yC,zC;

    //! ----------------------------------------------
    //! define the center of the faces {1,2,3} => P14
    //! ----------------------------------------------
    xC = (P1.x+P2.x+P3.x)/3;
    yC = (P1.y+P2.y+P3.y)/3;
    zC = (P1.z+P2.z+P3.z)/3;
    mesh::meshPoint P14(xC,yC,zC);

    //! ----------------------------------------------
    //! define the center of the faces {2,3,4} => P11
    //! ----------------------------------------------
    xC = (P2.x+P3.x+P4.x)/3;
    yC = (P2.y+P3.y+P4.y)/3;
    zC = (P2.z+P3.z+P4.z)/3;
    mesh::meshPoint P11(xC,yC,zC);

    //! ----------------------------------------------
    //! define the center of the faces {1,3,4} => P13
    //! ----------------------------------------------
    xC = (P1.x+P3.x+P4.x)/3;
    yC = (P1.y+P3.y+P4.y)/3;
    zC = (P1.z+P3.z+P4.z)/3;
    mesh::meshPoint P13(xC,yC,zC);

    //! ----------------------------------------------
    //! define the center of the faces {1,2,4} => P12
    //! ----------------------------------------------
    xC = (P1.x+P2.x+P4.x)/3;
    yC = (P1.y+P2.y+P4.y)/3;
    zC = (P1.z+P2.z+P4.z)/3;
    mesh::meshPoint P12(xC,yC,zC);

    //! --------------------------
    //! center of the tetrahedron
    //! --------------------------
    xC = (P1.x+P2.x+P3.x+P4.x)/4;
    yC = (P1.y+P2.y+P3.y+P4.y)/4;
    zC = (P1.z+P2.z+P3.z+P4.z)/4;
    mesh::meshPoint C(xC,yC,zC);

    //! ---------------
    //! points between
    //! ---------------
    mesh::meshPoint P5 = pointBetween(P1,P2);
    mesh::meshPoint P6 = pointBetween(P2,P3);
    mesh::meshPoint P7 = pointBetween(P2,P4);
    mesh::meshPoint P8 = pointBetween(P3,P4);
    mesh::meshPoint P9 = pointBetween(P1,P4);
    mesh::meshPoint P10 = pointBetween(P1,P3);

    //! -------------------
    //! build the elements
    //! -------------------
    meshElementByCoords element1, element2, element3, element4;
    element1.type = HEXA; element2.type = HEXA; element3.type = HEXA; element4.type = HEXA;
    element1.pointList<<P1<<P5<<P14<<P10<<P9<<P12<<C<<P13;
    element2.pointList<<P10<<P14<<P6<<P3<<P13<<C<<P11<<P8;
    element3.pointList<<P5<<P2<<P6<<P14<<P12<<P7<<P11<<C;
    element4.pointList<<P9<<P12<<C<<P13<<P4<<P7<<P11<<P8;

    std::vector<meshElementByCoords> decomposition;
    decomposition.push_back(element1);
    decomposition.push_back(element2);
    decomposition.push_back(element3);
    decomposition.push_back(element4);
    return decomposition;
}

//! -------------------
//! function: splitTri
//! details:
//! -------------------
std::vector<meshElementByCoords> tetSplitter::splitTri(const meshElementByCoords &aTri)
{
    //! --------------------
    //! define points 1,2,3
    //! --------------------
    mesh::meshPoint P1 = aTri.pointList[0];
    mesh::meshPoint P2 = aTri.pointList[1];
    mesh::meshPoint P3 = aTri.pointList[2];

    mesh::meshPoint P4 = pointBetween(P1,P2);
    mesh::meshPoint P5 = pointBetween(P2,P3);
    mesh::meshPoint P6 = pointBetween(P3,P1);

    mesh::meshPoint C((P1.x+P2.x+P3.x)/3.0,(P1.y+P2.y+P3.y)/3.0,(P1.z+P2.z+P3.z)/3.0);
    meshElementByCoords element1, element2, element3;
    element1.type = QUAD; element2.type = QUAD; element3.type = QUAD;
    element1.pointList<<P1<<P4<<C<<P6;
    element2.pointList<<P4<<P2<<P5<<C;
    element3.pointList<<C<<P5<<P3<<P6;

    std::vector<meshElementByCoords> vecElements;
    vecElements.push_back(element1);
    vecElements.push_back(element2);
    vecElements.push_back(element3);
    return vecElements;
}
