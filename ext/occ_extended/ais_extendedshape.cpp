//! ----------------
//! custom includes
//! ----------------
#include "ais_extendedshape.h"

//! ----
//! OCC
//! ----
#include <StdPrs_HLRPolyShape.hxx>
#include <StdPrs_ShadedShape.hxx>
#include <StdPrs_WFShape.hxx>
#include <Prs3d_LineAspect.hxx>
#include <GProp_GProps.hxx>
#include <GProp_PrincipalProps.hxx>
#include <BRepGProp.hxx>
#include <gp_Pnt.hxx>

IMPLEMENT_STANDARD_HANDLE(AIS_ExtendedShape,AIS_ColoredShape) // implements the DownCast() method of the Handle class
IMPLEMENT_STANDARD_RTTIEXT(AIS_ExtendedShape,AIS_ColoredShape)

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
AIS_ExtendedShape::AIS_ExtendedShape(const TopoDS_Shape &aTopoDSShape):AIS_ColoredShape(aTopoDSShape)
{
    myTopoDSShape = aTopoDSShape;
    myName = "";
    myIndex = -1;
}

//! ---------------------------
//! function: Compute
//! details:  hidden line mode
//! ---------------------------
void AIS_ExtendedShape::Compute(const Handle_Prs3d_Projector& aProjector, const Handle_Prs3d_Presentation& aPresentation)
{
    // hidden line mode
    myDrawer->EnableDrawHiddenLine();
    occHandle(Prs3d_LineAspect) lineAspectHidden = new Prs3d_LineAspect(Quantity_NOC_BLACK,Aspect_TOL_DASH,2);
    myDrawer->SetHiddenLineAspect(lineAspectHidden);
    StdPrs_HLRPolyShape::Add(aPresentation,myTopoDSShape,myDrawer,aProjector);
}

//! ------------------------------------
//! function: Compute method
//! details:  wireframe and shaded mode
//! ------------------------------------
void AIS_ExtendedShape::Compute(const Handle_PrsMgr_PresentationManager3d &/*aPresentationManager*/,
                               const Handle_Prs3d_Presentation& aPresentation,
                               const Standard_Integer aMode)
{
    occHandle(Prs3d_LineAspect) lineAspect1;
    occHandle(Prs3d_LineAspect) lineAspect2;

    switch(aMode)
    {
    case 0:
        lineAspect1 = new Prs3d_LineAspect(Quantity_NOC_BLACK,Aspect_TOL_SOLID,2);
        lineAspect2 = new Prs3d_LineAspect(Quantity_NOC_RED,Aspect_TOL_SOLID,2);
        myDrawer->SetUnFreeBoundaryAspect(lineAspect1);
        myDrawer->SetFreeBoundaryAspect(lineAspect2);
        StdPrs_WFShape::Add(aPresentation,myTopoDSShape,myDrawer);
        break;

    case 1:
        StdPrs_ShadedShape::Add(aPresentation,myTopoDSShape,myDrawer);
        break;
    }
}
