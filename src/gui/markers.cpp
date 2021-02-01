//! custom includes
#include "markers.h"
#include "src/utils/myconstant.h"

//! OCC
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <Prs3d_IsoAspect.hxx>

#include <TopoDS_Solid.hxx>
#include <TopoDS_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>

#include <AIS_Shape.hxx>
#include <AIS_DisplayMode.hxx>
#include <AIS_Trihedron.hxx>
#include <gp_Trsf.hxx>

#include <Geom_Axis2Placement.hxx>
#include <Prs3d_DatumAspect.hxx>
#include <Prs3d_ArrowAspect.hxx>
#include <Prs3d_ShadingAspect.hxx>
#include <Prs3d_LineAspect.hxx>
#include <Aspect_WidthOfLine.hxx>

//! Qt
#include <QVector>
#include <Quantity_NameOfColor.hxx>

//! -----------------------------------
//! function: buildCustomTrihedron
//! details:  build a custom trihedron
//! -----------------------------------
AIS_CustomTrihedron_handle_reg markers::buildCustomTrihedron(QVector<double> Origin,
                                                             QVector<QVector<double>> directionalData,
                                                             double axisLength)
{
    cout<<"markers::buildTrihedron()->____function called____"<<endl;
    gp_Pnt P(Origin.at(0),Origin.at(1),Origin.at(2));

    QVector<double> xAxisData = directionalData.at(0);
    QVector<double> zAxisData = directionalData.at(2);

    gp_Dir N(zAxisData.at(0),zAxisData.at(1),zAxisData.at(2));
    gp_Dir V(xAxisData.at(0),xAxisData.at(1),xAxisData.at(2));
    gp_Ax2 theCoordinateSystem(P,N,V);

    //! ------------------------
    //! dimensions of the arrow
    //! ------------------------
    double D = axisLength/4.0;
    double H_cylinder = 4*D;
    double H_cone = 1.75*D;
    double R_cylinder = D/2.0;
    double R_cone = 1.25*R_cylinder;

    //! -----------------------
    //! build the cylinder "Z"
    //! -----------------------
    gp_Ax2 axes(P,N,V);
    TopoDS_Solid cylinderZ = BRepPrimAPI_MakeCylinder(axes, R_cylinder, H_cylinder).Solid();
    gp_Pnt point;
    gp_Vec vec(axes.Direction());
    vec.Multiply(H_cylinder);
    point.Translate(vec);
    gp_Ax2 axes_copy = axes;
    axes_copy.Translate(vec);

    //! -------------------
    //! build the cone "Z"
    //! -------------------
    TopoDS_Solid coneZ = BRepPrimAPI_MakeCone(axes_copy, 1.5*R_cone, 0.0, H_cone).Solid();
    TopoDS_Compound trihedron;
    TopoDS_Builder theBuilder;

    //! ------------------------------------
    //! build the cone and the cylinder "X"
    //! ------------------------------------
    gp_Trsf rot1;
    rot1.SetRotation(gp_Ax1(P,theCoordinateSystem.YDirection()),PIGREEK_HALVED);
    TopoDS_Shape cylinderX = BRepBuilderAPI_Transform (cylinderZ,rot1,true).ModifiedShape(cylinderZ);
    TopoDS_Shape coneX = BRepBuilderAPI_Transform (coneZ,rot1,true).ModifiedShape(coneZ);

    //! ------------------------------------
    //! build the cone and the cylinder "Y"
    //! ------------------------------------
    gp_Trsf rot2;
    rot2.SetRotation(gp_Ax1(P,theCoordinateSystem.Direction()),PIGREEK_HALVED);

    const TopoDS_Solid cylinderY = TopoDS::Solid(BRepBuilderAPI_Transform (cylinderX,rot2,true).ModifiedShape(cylinderX));
    const TopoDS_Solid coneY = TopoDS::Solid(BRepBuilderAPI_Transform (coneX,rot2,true).ModifiedShape(coneX));

    //! ---------------
    //! build a sphere
    //! ---------------
    TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(P,R_cylinder);

    theBuilder.MakeCompound(trihedron);
    theBuilder.Add(trihedron,cylinderZ);
    theBuilder.Add(trihedron,coneZ);
    theBuilder.Add(trihedron,cylinderX);
    theBuilder.Add(trihedron,coneX);
    theBuilder.Add(trihedron,cylinderY);
    theBuilder.Add(trihedron,coneY);
    theBuilder.Add(trihedron,sphere);

    Handle_AIS_CustomTrihedron AIS_trihedron = new AIS_CustomTrihedron(trihedron);
    TopTools_IndexedMapOfShape sMap;
    TopExp::MapShapes(trihedron,TopAbs_SOLID,sMap);

    AIS_trihedron->SetCustomColor(sMap.FindKey(1),Quantity_NOC_BLUE1);
    AIS_trihedron->SetCustomColor(sMap.FindKey(2),Quantity_NOC_BLUE1);
    AIS_trihedron->SetCustomColor(sMap.FindKey(3),Quantity_NOC_RED);
    AIS_trihedron->SetCustomColor(sMap.FindKey(4),Quantity_NOC_RED);
    AIS_trihedron->SetCustomColor(sMap.FindKey(5),Quantity_NOC_GREEN);
    AIS_trihedron->SetCustomColor(sMap.FindKey(6),Quantity_NOC_GREEN);
    AIS_trihedron->SetCustomColor(sMap.FindKey(7),Quantity_NOC_WHITE);
    AIS_trihedron->Attributes()->SetFaceBoundaryDraw(false);
    AIS_trihedron->SetDisplayMode(AIS_Shaded);

    return AIS_trihedron;
}

//! ---------------------------------
//! function: buildDoubleArrowMarker
//! details:
//! ---------------------------------
AIS_DoubleArrowMarker_handle_reg markers::buildDoubleArrowMarker(const gp_Pnt &P, const gp_Dir &V, double D)
{
    //! ------------------------
    //! dimensions of the arrow
    //! ------------------------
    double H_cylinder = 4*D;
    double H_cone = 1.75*D;
    double R_cylinder = D/2;
    double R_cone = 1.25*R_cylinder;

    //! ----------------
    //! shift the point
    //! ----------------
    gp_Vec s(V);
    s.Multiply(-(H_cylinder+H_cone));
    gp_Pnt P_copy = P;
    P_copy.Translate(s);

    //! ------------
    //! first arrow
    //! ------------
    gp_Ax2 axes(P_copy, V);
    TopoDS_Solid cylinder = BRepPrimAPI_MakeCylinder(axes, R_cylinder, H_cylinder).Solid();
    gp_Pnt point;
    gp_Vec vec(V);
    vec.Multiply(H_cylinder);
    point.Translate(vec);
    axes.Translate(vec);
    TopoDS_Solid cone = BRepPrimAPI_MakeCone(axes, 1.5*R_cone, 0.0, H_cone).Solid();

    TopoDS_Compound arrow;
    TopoDS_Builder theBuilder;
    theBuilder.MakeCompound(arrow);
    theBuilder.Add(arrow,cylinder);
    theBuilder.Add(arrow,cone);

    //! -------------
    //! second arrow
    //! -------------
    gp_Ax2 axes1(P_copy, V);
    vec.Multiply((1/H_cylinder)*(H_cylinder+H_cone));
    axes1.Translate(vec);

    gp_Trsf trsf;
    trsf.SetMirror(axes1);
    BRepBuilderAPI_Transform t(arrow,trsf,true);
    theBuilder.Add(arrow,t.ModifiedShape(arrow));

    AIS_DoubleArrowMarker_handle_reg AIS_DoubleArrow = new AIS_DoubleArrowMarker(arrow);
    AIS_DoubleArrow->Attributes()->SetFaceBoundaryDraw(false);
    AIS_DoubleArrow->SetColor(Quantity_NOC_RED);
    AIS_DoubleArrow->SetDisplayMode(AIS_Shaded);
    return AIS_DoubleArrow;
}

//! ---------------------
//! function: buildArrow
//! details:
//! ---------------------
AIS_ArrowMarker_handle_reg markers::buildArrowMarker(const gp_Pnt &P, const gp_Dir &V, double D)
{
    //! ------------------------
    //! dimensions of the arrow
    //! ------------------------
    double H_cylinder = 4*D;
    double H_cone = 1.75*D;
    double R_cylinder = D/2;
    double R_cone = 1.25*R_cylinder;

    gp_Ax2 axes(P, V);
    TopoDS_Solid cylinder = BRepPrimAPI_MakeCylinder(axes, R_cylinder, H_cylinder).Solid();
    gp_Pnt point;
    gp_Vec vec(V);
    vec.Multiply(H_cylinder);
    point.Translate(vec);
    axes.Translate(vec);
    TopoDS_Solid cone = BRepPrimAPI_MakeCone(axes, 1.5*R_cone, 0.0, H_cone).Solid();
    TopoDS_Compound arrow;
    TopoDS_Builder theBuilder;
    theBuilder.MakeCompound(arrow);
    theBuilder.Add(arrow,cylinder);
    theBuilder.Add(arrow,cone);
    AIS_ArrowMarker_handle_reg AIS_arrow = new AIS_ArrowMarker(arrow);
    AIS_arrow->SetDisplayMode(AIS_Shaded);
    AIS_arrow->Attributes()->SetFaceBoundaryDraw(false);
    AIS_arrow->SetColor(Quantity_NOC_RED);
    return AIS_arrow;
}

//! -----------------------------------
//! function: buildSphereMarker
//! details:  build a spherical marker
//! -----------------------------------
AIS_SphereMarker_handle_reg markers::buildSphereMarker(gp_Pnt C, double radius)
{
    TopoDS_Shape sphere = BRepPrimAPI_MakeSphere(C,radius);
    AIS_SphereMarker_handle_reg AIS_sphere = new AIS_SphereMarker(sphere);
    AIS_sphere->SetColor(Quantity_NOC_BLUE1);

    //! ---------------------
    //! do not show isolines
    //! ---------------------
    occHandle(Prs3d_IsoAspect) isoLineAspect = new Prs3d_IsoAspect(Quantity_NOC_GRAY,Aspect_TOL_SOLID,1.0,1);
    isoLineAspect->SetNumber(0);
    AIS_sphere->Attributes()->SetUIsoAspect(isoLineAspect);
    AIS_sphere->Attributes()->SetVIsoAspect(isoLineAspect);

    AIS_sphere->Attributes()->SetFaceBoundaryDraw(false);
    AIS_sphere->SetDisplayMode(AIS_Shaded);
    return AIS_sphere;
}

//! -------------------------
//! function: buildTrihedron
//! details:
//! -------------------------
AIS_Trihedron_handle_reg markers::buildTrihedron(QVector<double> Origin, QVector<QVector<double>> directionalData, int axisLength)
{
    gp_Pnt P(Origin.at(0),Origin.at(1),Origin.at(2));
    QVector<double> xAxisData = directionalData.at(0);
    QVector<double> zAxisData = directionalData.at(2);

    gp_Dir N(zAxisData.at(0),zAxisData.at(1),zAxisData.at(2));
    gp_Dir V(xAxisData.at(0),xAxisData.at(1),xAxisData.at(2));
    gp_Ax2 theCoordinateSystem(P,N,V);

    occHandle(Geom_Axis2Placement) thePlacement =  new Geom_Axis2Placement(theCoordinateSystem);

    AIS_Trihedron_handle_reg theTrihedron = new AIS_Trihedron(thePlacement);
    theTrihedron->SetTextColor(Quantity_NOC_BLACK);

    theTrihedron->Attributes()->DatumAspect()->SetAxisLength(axisLength,axisLength,axisLength);
    theTrihedron->Attributes()->DatumAspect()->FirstAxisAspect()->SetColor(Quantity_NOC_RED);
    theTrihedron->Attributes()->DatumAspect()->SecondAxisAspect()->SetColor(Quantity_NOC_GREEN);
    theTrihedron->Attributes()->DatumAspect()->ThirdAxisAspect()->SetColor(Quantity_NOC_BLUE1);
    theTrihedron->Attributes()->DatumAspect()->FirstAxisAspect()->SetWidth(Aspect_WOL_VERYTHICK);
    theTrihedron->Attributes()->DatumAspect()->SecondAxisAspect()->SetWidth(Aspect_WOL_VERYTHICK);
    theTrihedron->Attributes()->DatumAspect()->ThirdAxisAspect()->SetWidth(Aspect_WOL_VERYTHICK);

    //! -------------
    //! arrow aspect
    //! -------------
    const occHandle(Prs3d_ArrowAspect) &arrowAspect = theTrihedron->Attributes()->ArrowAspect();
    arrowAspect->SetColor(Quantity_NOC_BLACK);
    theTrihedron->Attributes()->SetArrowAspect(arrowAspect);

    occHandle(Graphic3d_AspectFillArea3d) g = new Graphic3d_AspectFillArea3d();
    Graphic3d_MaterialAspect ma(Graphic3d_NOM_BRASS);
    g->SetBackMaterial(ma);
    g->SetFrontMaterial(ma);

    occHandle(Prs3d_ShadingAspect) shadingAspect = new Prs3d_ShadingAspect(g);
    theTrihedron->Attributes()->SetShadingAspect(shadingAspect);
    return theTrihedron;
}

//! ---------------------------
//! function: buildCurvedArrow
//! details:
//! ---------------------------
AIS_CurvedArrowMarker_handle_reg markers::buildCurvedArrow(const gp_Ax2 &axes, double Rin, bool invert)
{
    double Rext = 1.2*Rin;
    double H = Rext-Rin;

    const TopoDS_Solid &cylinderIn = BRepPrimAPI_MakeCylinder(axes,Rin,H,1.5*PIGREEK);
    const TopoDS_Solid &cylinderOut = BRepPrimAPI_MakeCylinder(axes,Rext,H,1.5*PIGREEK);

    TopTools_ListOfShape s1,s2;
    s1.Append(cylinderOut);
    s2.Append(cylinderIn);
    TopoDS_Shape ring = BRepAlgoAPI_Cut(cylinderOut,cylinderIn).Shape();
    gp_Trsf aTrans;
    gp_Vec trans(axes.Direction());
    trans.Multiply(-H/2.0);
    aTrans.SetTranslation(trans);
    ring = BRepBuilderAPI_Transform(ring,aTrans,false).ModifiedShape(ring);

    gp_Ax2 coneAxes = axes;
    gp_Dir Ydir = axes.YDirection();
    gp_Vec PO(Ydir);
    double t = -(Rext+Rin)/2.0;
    PO.Multiply(t);
    coneAxes.Translate(PO);

    gp_Trsf aRot;
    gp_Pnt P;
    P.Translate(PO);
    gp_Ax1 ax1(P,Ydir);
    aRot.SetRotation(ax1,PIGREEK_HALVED);
    coneAxes.Transform(aRot);

    TopoDS_Solid cone = BRepPrimAPI_MakeCone(coneAxes, 1.5*H, 0.0, 2*H).Solid();
    TopoDS_Builder theBuilder;
    TopoDS_Compound curvedArrow;
    theBuilder.MakeCompound(curvedArrow);
    theBuilder.Add(curvedArrow,ring);
    theBuilder.Add(curvedArrow,cone);

    //! ------------------
    //! reverse the arrow
    //! ------------------
    if(invert==true)
    {
        //gp_Trsf aTrans;
        //aTrans.SetMirror(coneAxes.Location());
        //curvedArrow = BRepBuilderAPI_Transform(curvedArrow,aTrans,false).ModifiedShape(curvedArrow);
    }

    AIS_CurvedArrowMarker_handle_reg marker = new AIS_CurvedArrowMarker(curvedArrow);
    marker->SetDisplayMode(AIS_Shaded);
    marker->Attributes()->SetFaceBoundaryDraw(false);
    marker->SetColor(Quantity_NOC_RED);
    return marker;
}
