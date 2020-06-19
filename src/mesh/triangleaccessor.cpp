#include "triangleaccessor.h"

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TriangleAccessor::TriangleAccessor (const TopoDS_Face& aFace)
{
    TopLoc_Location aLoc;
    myPoly = BRep_Tool::Triangulation (aFace, aLoc);
    myTrsf = aLoc.Transformation();
    myNbTriangles = (myPoly.IsNull() ? 0 : myPoly->Triangles().Length());
    myInvert = (aFace.Orientation() == TopAbs_REVERSED);
    if(myTrsf.IsNegative()) myInvert = ! myInvert;
}

//! -------------------------------------------------------
//! function: getTriangle
//! details:  get the i-th triangle and the outward normal
//! -------------------------------------------------------
void TriangleAccessor::GetTriangle (int iTri, gp_Vec &theNormal, gp_Pnt &thePnt1, gp_Pnt &thePnt2, gp_Pnt &thePnt3)
{
    //! -----------------------
    //! get positions of nodes
    //! -----------------------
    int iNode1, iNode2, iNode3;
    myPoly->Triangles()(iTri).Get (iNode1, iNode2, iNode3);
    thePnt1 = myPoly->Nodes()(iNode1);
    thePnt2 = myPoly->Nodes()(myInvert ? iNode3 : iNode2);
    thePnt3 = myPoly->Nodes()(myInvert ? iNode2 : iNode3);

    //! ------------------------------------
    //! apply transormation if not identity
    //! ------------------------------------
    if (myTrsf.Form() != gp_Identity)
    {
        thePnt1.Transform (myTrsf);
        thePnt2.Transform (myTrsf);
        thePnt3.Transform (myTrsf);
    }

    //! ---------------------
    //! calculate the normal
    //! ---------------------
    theNormal = (thePnt2.XYZ() - thePnt1.XYZ()) ^ (thePnt3.XYZ() - thePnt1.XYZ());
    double aNorm = theNormal.Magnitude();
    if (aNorm > gp::Resolution()) theNormal /= aNorm;
}

//! ----------------------
//! function: NbTriangles
//! details:
//! ----------------------
int TriangleAccessor::NbTriangles()
{
    return myNbTriangles;
}
