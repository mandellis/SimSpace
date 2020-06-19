#include "geometryhealing.h"

//! OCC
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <BRep_Tool.hxx>

#include <ShapeFix_Shape.hxx>

#include <ShapeBuild_ReShape.hxx>
#include <ShapeFix_Face.hxx>

#include <Precision.hxx>

//! C++
#include <stdio.h>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
geometryHealing::geometryHealing(QObject *parent): QObject(parent)
{
    ;
}


geometryHealing::geometryHealing(const TopoDS_Shape &aShape, QObject *parent): QObject(parent), myTopoDS_Shape(aShape)
{
    ;
}

//! --------------------------------
//! function: removeDegenerateEdges
//! details:
//! --------------------------------
void geometryHealing::removeDegenerateEdges(TopoDS_Shape &shape)
{
    TopExp_Explorer exp1;
    occHandle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
    rebuild->Apply(shape);
    int count = 0;
    for (exp1.Init(shape, TopAbs_EDGE); exp1.More(); exp1.Next())
    {
        TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
        if (BRep_Tool::Degenerated(edge)) rebuild->Remove(edge);
        count++;
    }
    cout<<"geometryHealing::removeDegenerateEdges()->____degenerate edges: "<<count<<"____"<<endl;
    rebuild->Apply(shape);
}

//! -------------------
//! function: fixFaces
//! details:
//! -------------------
void geometryHealing::fixFaces(TopoDS_Shape &shape)
{
    occHandle(ShapeBuild_ReShape) rebuild = new ShapeBuild_ReShape;
    TopExp_Explorer exp0;
    for (exp0.Init (shape, TopAbs_FACE); exp0.More(); exp0.Next())
    {
        TopoDS_Face face = TopoDS::Face(exp0.Current());

        occHandle(ShapeFix_Face) sff = new ShapeFix_Face (face);
        sff->FixAddNaturalBoundMode() = Standard_True;
        sff->FixSmallAreaWireMode() = Standard_True;
        sff->Perform();

        if(sff->Status(ShapeExtend_DONE1) ||
                sff->Status(ShapeExtend_DONE2) ||
                sff->Status(ShapeExtend_DONE3) ||
                sff->Status(ShapeExtend_DONE4) ||
                sff->Status(ShapeExtend_DONE5))
        {
            cout << "repaired face ";
            if(sff->Status(ShapeExtend_DONE1))
                cout << "(some wires are fixed)" <<endl;
            else if(sff->Status(ShapeExtend_DONE2))
                cout << "(orientation of wires fixed)" <<endl;
            else if(sff->Status(ShapeExtend_DONE3))
                cout << "(missing seam added)" <<endl;
            else if(sff->Status(ShapeExtend_DONE4))
                cout << "(small area wire removed)" <<endl;
            else if(sff->Status(ShapeExtend_DONE5))
                cout << "(natural bounds added)" <<endl;
            TopoDS_Face newface = sff->Face();

            rebuild->Replace(face, newface);
        }
    }
    shape = rebuild->Apply(shape);
}

//! ------------------
//! function: perform
//! details:
//! ------------------
#include <ShapeFix.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <ShapeFix_FixSmallFace.hxx>

//! -----------------------------------------------------------
//! function: perform
//! details:  use "precision" to set "tolerance": no words ...
//! -----------------------------------------------------------
void geometryHealing::perform(double Prec)
{
    TopoDS_Shape shape = myTopoDS_Shape;

    double maxTol = 10*Prec;
    double minTol = Prec;

    cout<<"/--------------------------------------------/"<<endl;
    cout<<"____min edge length: "<<getMinimumEdge(shape)<<"____"<<endl;
    cout<<"____precision: "<<Prec<<"____"<<endl;
    cout<<"____maxTol: "<<maxTol<<"____"<<endl;
    cout<<"____minTol: "<<minTol<<"____"<<endl;

    /*
    //! to be studied - at the moment it has no effect
    occHandle(ShapeFix_Shape) sfs = new ShapeFix_Shape;
    sfs->Init(shape);
    sfs->SetPrecision (Prec);
    sfs->SetMaxTolerance (maxTol);
    sfs->SetMinTolerance (minTol);

    sfs->FixWireTool()->ModifyTopologyMode() = true;
    sfs->FixWireTool()->ModifyGeometryMode() = true;    //! true by default

    sfs->FixWireTool()->FixReorderMode() = true;
    sfs->FixWireTool()->FixReorder();

    sfs->FixWireTool()->FixSelfIntersectionMode() = true;
    sfs->FixWireTool()->FixSelfIntersection();

    //sfs->FixFaceTool()->FixOrientationMode() = true;
    //sfs->FixFaceTool()->FixOrientation();

    sfs->FixWireTool()->FixSmallMode() = true;
    sfs->FixWireTool()->FixSmall(true,Prec);

    sfs->FixWireTool()->FixDegeneratedMode() = true;
    sfs->FixWireTool()->FixDegenerated();

    shape = sfs->Shape();
    */

    //! ------------------
    //! work on wireframe
    //! ------------------
    occHandle(ShapeFix_Wireframe) sfwf = new ShapeFix_Wireframe(shape);

    sfwf->SetPrecision(Prec);
    sfwf->SetMinTolerance(minTol);
    sfwf->SetMaxTolerance(maxTol);

    sfwf->ModeDropSmallEdges()= Standard_True;
    sfwf->FixSmallEdges();
    sfwf->FixWireGaps();

    shape = sfwf->Shape();

    //! --------------
    //! work on faces
    //! ---------------
    occHandle(ShapeFix_FixSmallFace) sff = new ShapeFix_FixSmallFace();
    sff->Init(shape);
    //sff->SetPrecision(pow(Prec,2));
    sff->SetPrecision(Prec);
    sff->SetMaxTolerance(maxTol);

    /*shape = */sff->FixSpotFace();
    /*shape = */sff->FixStripFace();

    sff->Perform();
    shape = sff->Shape();

    //! -------------------
    //! work on wirefreame
    //! -------------------
    sfwf = new ShapeFix_Wireframe(shape);

    sfwf->SetPrecision(Prec);
    sfwf->SetMaxTolerance(maxTol);

    sfwf->FixWireGaps();

    sfwf->ModeDropSmallEdges()= Standard_True;
    sfwf->FixSmallEdges();

    shape = sfwf->Shape();

    cout<<"____min edge length: "<<getMinimumEdge(shape)<<"____"<<endl;
    cout<<"/--------------------------------------------/"<<endl;
    myTopoDS_Shape = shape;
}


double geometryHealing::getMinimumEdge(const TopoDS_Shape &aShape)
{
    GProp_GProps props;
    double minlength = 1e99;
    for(TopExp_Explorer exp(aShape,TopAbs_EDGE);exp.More();exp.Next())
    {
        if(BRep_Tool::Degenerated(TopoDS::Edge(exp.Current()))) continue;
        BRepGProp::LinearProperties(exp.Current(),props);
        if(props.Mass()<minlength)minlength=props.Mass();
    }
    return minlength;
}
