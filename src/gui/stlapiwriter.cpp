#include "stlapiwriter.h"
using namespace STLAPIWriter;

#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp_Explorer.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

STLAPIWriter::TriangleAccessor::TriangleAccessor(const TopoDS_Face& aFace)
{
    TopLoc_Location aLoc;
    myPoly = BRep_Tool::Triangulation (aFace, aLoc);
    myTrsf = aLoc.Transformation();
    myNbTriangles = (myPoly.IsNull() ? 0 : myPoly->Triangles().Length());
    myInvert = (aFace.Orientation() == TopAbs_REVERSED);
    if (myTrsf.IsNegative()) myInvert = ! myInvert;
}

void STLAPIWriter::TriangleAccessor::GetTriangle(int iTri, gp_Vec &theNormal, gp_Pnt &thePnt1, gp_Pnt &thePnt2, gp_Pnt &thePnt3)
{
    int iNode1, iNode2, iNode3;
    myPoly->Triangles()(iTri).Get (iNode1, iNode2, iNode3);
    thePnt1 = myPoly->Nodes()(iNode1);
    thePnt2 = myPoly->Nodes()(myInvert ? iNode3 : iNode2);
    thePnt3 = myPoly->Nodes()(myInvert ? iNode2 : iNode3);

    //! apply transormation if not identity
    if (myTrsf.Form() != gp_Identity)
    {
        thePnt1.Transform (myTrsf);
        thePnt2.Transform (myTrsf);
        thePnt3.Transform (myTrsf);
    }

    //! calculate normal
    theNormal = (thePnt2.XYZ() - thePnt1.XYZ()) ^ (thePnt3.XYZ() - thePnt1.XYZ());
    Standard_Real aNorm = theNormal.Magnitude();
    if (aNorm > gp::Resolution()) theNormal /= aNorm;
}

//! ----------------
//! function: Write
//! details:
//! ----------------
bool STLAPIWriter::Write(const TopoDS_Shape& theShape, const char *theFileName)
{
    FILE* aFile = fopen(theFileName, "wb");
    if (aFile==NULL) return false;
    fprintf (aFile, "solid shape, STL ascii file, created with Open CASCADE Technology\n");
    for (TopExp_Explorer exp (theShape, TopAbs_FACE); exp.More(); exp.Next())
    {
        TriangleAccessor aTool (TopoDS::Face (exp.Current()));
        for (int iTri = 1; iTri <= aTool.NbTriangles(); iTri++)
        {
            gp_Vec aNorm;
            gp_Pnt aPnt1, aPnt2, aPnt3;
            aTool.GetTriangle (iTri, aNorm, aPnt1, aPnt2, aPnt3);
            fprintf (aFile,
                     " facet normal %.6e %.6e %.6e\n"
                     "   outer loop\n"
                     "     vertex %.6e %.6e %.6e\n"
                     "     vertex %.6e %.6e %.6e\n"
                     "     vertex %.6e %.6e %.6e\n"
                     "   endloop\n"
                     " endfacet\n",
                     aNorm.X(), aNorm.Y(), aNorm.Z(),
                     aPnt1.X(), aPnt1.Y(), aPnt1.Z(),
                     aPnt2.X(), aPnt2.Y(), aPnt2.Z(),
                     aPnt3.X(), aPnt3.Y(), aPnt3.Z());
        }
    }
    fprintf (aFile, "endsolid shape\n");
    fclose (aFile);
    return true;
}

//! ------------------------
//! function: WriteExtended
//! details:
//! ------------------------
bool STLAPIWriter::WriteExtended(const TopoDS_Shape& theShape, const char* theFileName)
{
    FILE* aFile = fopen(theFileName, "w");
    if (aFile==NULL) return false;
    fprintf (aFile, "solid shape, STL ascii file, each triangle has an int tag\n");
    TopTools_IndexedMapOfShape M;
    TopExp::MapShapes(theShape,TopAbs_FACE,M);
    for(TopExp_Explorer exp (theShape, TopAbs_FACE); exp.More(); exp.Next())
    {
        int tag = M.FindIndex(exp.Current());
        TriangleAccessor aTool (TopoDS::Face (exp.Current()));
        for(int iTri = 1; iTri <= aTool.NbTriangles(); iTri++)
        {
            gp_Vec aNorm;
            gp_Pnt aPnt1, aPnt2, aPnt3;
            aTool.GetTriangle (iTri, aNorm, aPnt1, aPnt2, aPnt3);

            fprintf (aFile,
                     " facet normal %.6e %.6e %.6e\n"
                     "   outer loop\n"
                     "     vertex %.6e %.6e %.6e\n"
                     "     vertex %.6e %.6e %.6e\n"
                     "     vertex %.6e %.6e %.6e\n"
                     "   endloop\n"
                     "   %d\n"
                     " endfacet\n",
                     aNorm.X(), aNorm.Y(), aNorm.Z(),
                     aPnt1.X(), aPnt1.Y(), aPnt1.Z(),
                     aPnt2.X(), aPnt2.Y(), aPnt2.Z(),
                     aPnt3.X(), aPnt3.Y(), aPnt3.Z(),
                     tag);
        }
    }
    fprintf (aFile, "endsolid shape\n");
    fclose (aFile);
    return true;
}
