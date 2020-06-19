#include "wbtrihedron.h"

//! --------------------------
//! for drawing the trihedron
//! --------------------------
#include <Prs3d_ToolCylinder.hxx>
#include <Prs3d_ToolDisk.hxx>
#include <Prs3d_ToolSphere.hxx>
#include <gp_Dir.hxx>
#include <gp_Trsf.hxx>
#include <Graphic3d_Group.hxx>
#include <Graphic3d_ArrayOfTriangles.hxx>
#include <Prs3d_ShadingAspect.hxx>


#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
WBTrihedron::WBTrihedron(double size):mySize(size)
{
    //this->setSelectionMode(selectionMode_singleAxisSelect);
    //this->setSelectionMode(selectionMode_wholeSelect);
    this->setGlobalSelMode(1);
}

//! ------------------
//! function: Compute
//! details:
//! ------------------
void WBTrihedron::Compute(const occHandle(PrsMgr_PresentationManager3d) &/*thePresentationManager*/,
                          const occHandle(Prs3d_Presentation) &thePresentation,
                          const Standard_Integer /*theMode*/)
{
    cout<<"void WBTrihedron::Compute()->____function called____"<<endl;
    const double pi = 3.14159265359;

    //! --------------------------------------------------------
    //! the base dimension (size of the box) should be computed
    //! on the basis of the windows size
    //! --------------------------------------------------------
    double baseDim = 40;
    int NbSlices = 36;
    int NbStacks = 10;

    double cylinderRadius = baseDim/12.0;
    double sphereRadius = cylinderRadius*1.5;
    double height = 1.5*baseDim;
    double tipRadius = cylinderRadius*1.65;
    double tipHeight = tipRadius*3.00;

    //! -------------------
    //! identity transform
    //! -------------------
    gp_Trsf aTrsf;

    occHandle(Prs3d_ShadingAspect) aShadingAspect = new Prs3d_ShadingAspect();
    aShadingAspect->SetColor(Quantity_NOC_WHITE);

    //! -----------
    //! the sphere
    //! -----------
    occHandle(Graphic3d_ArrayOfTriangles) sphere = Prs3d_ToolSphere::Create(sphereRadius,NbSlices,NbStacks,aTrsf);
    occHandle(Graphic3d_Group) aGroup0 = thePresentation->NewGroup();
    aGroup0->AddPrimitiveArray(sphere);
    aGroup0->SetGroupPrimitivesAspect(aShadingAspect->Aspect());

    cout<<"void WBTrihedron::Compute()->____sphere done____"<<endl;

    //! ----------------------
    //! axis Z - the cylinder
    //! ----------------------
    occHandle(Graphic3d_ArrayOfTriangles) AxisZ;
    Prs3d_ToolCylinder toolCylinderZ(cylinderRadius,cylinderRadius,height,NbSlices,NbStacks);
    toolCylinderZ.FillArray(AxisZ,polyTriangulationAxisZ,aTrsf);

    cout<<"void WBTrihedron::Compute()->____cylinder Z done____"<<endl;

    //! -----------------
    //! axis Z - the tip
    //! -----------------
    occHandle(Graphic3d_ArrayOfTriangles) AxisZ_cone;
    aTrsf.SetTranslation(gp_Vec(0,0,height));
    Prs3d_ToolCylinder toolConeZ(tipRadius,0,tipHeight,NbSlices,NbStacks);
    toolConeZ.FillArray(AxisZ_cone,polyTriangulationTipAxisZ,aTrsf);

    cout<<"void WBTrihedron::Compute()->____cylinder tip Z done____"<<endl;

    //! -------------------------
    //! axis Z - the hollow disk
    //! -------------------------
    occHandle(Graphic3d_ArrayOfTriangles) hollowDiskZ;
    Prs3d_ToolDisk toolDiskZ(cylinderRadius,tipRadius,NbSlices,NbStacks);
    toolDiskZ.FillArray(hollowDiskZ,polyTriangulationTipHollowDiskZ,aTrsf);

    cout<<"void WBTrihedron::Compute()->____hollow disk Z done____"<<endl;

    occHandle(Graphic3d_Group) GroupZ = thePresentation->NewGroup();
    GroupZ->AddPrimitiveArray (AxisZ);
    GroupZ->AddPrimitiveArray (AxisZ_cone);
    GroupZ->AddPrimitiveArray (hollowDiskZ);

    occHandle(Prs3d_ShadingAspect) aShadingAspectZ = new Prs3d_ShadingAspect();
    aShadingAspectZ->SetColor(Quantity_NOC_BLUE1);
    GroupZ->SetGroupPrimitivesAspect(aShadingAspectZ->Aspect());

    cout<<"void WBTrihedron::Compute()->____axis Z array of primitives done____"<<endl;

    //! ----------------------
    //! axis Y - the cylinder
    //! ----------------------
    gp_Trsf aTrsf1;
    aTrsf1.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(1,0,0)),-pi/2);

    occHandle(Graphic3d_ArrayOfTriangles) AxisY;
    Prs3d_ToolCylinder toolCilynderY(cylinderRadius,cylinderRadius,height,NbSlices,NbStacks);
    toolCilynderY.FillArray(AxisY,polyTriangulationAxisY,aTrsf1);

    cout<<"void WBTrihedron::Compute()->____cylinder Y done____"<<endl;

    //! -----------------
    //! axis Y - the tip
    //! -----------------
    gp_Trsf t;
    t.SetTranslation(gp_Vec(0,0,height));
    aTrsf1.Multiply(t);

    occHandle(Graphic3d_ArrayOfTriangles) AxisY_cone;
    Prs3d_ToolCylinder toolConeY(tipRadius,0,tipHeight,NbSlices,NbStacks);
    toolConeY.FillArray(AxisY_cone,polyTriangulationTipAxisY,aTrsf1);

    cout<<"void WBTrihedron::Compute()->____cylinder tip Y done____"<<endl;

    //! -------------------------
    //! axis Y - the hollow disk
    //! -------------------------
    occHandle(Graphic3d_ArrayOfTriangles) hollowDiskY;
    Prs3d_ToolDisk toolDiskY(cylinderRadius,tipRadius,NbSlices,NbStacks);
    toolDiskY.FillArray(hollowDiskY,polyTriangulationTipHollowDiskY,aTrsf1);

    cout<<"void WBTrihedron::Compute()->____hollow disk Y done____"<<endl;

    occHandle(Graphic3d_Group) aGroupY = thePresentation->NewGroup();
    aGroupY->AddPrimitiveArray (AxisY);
    aGroupY->AddPrimitiveArray (AxisY_cone);
    aGroupY->AddPrimitiveArray (hollowDiskY);

    occHandle(Prs3d_ShadingAspect) aShadingAspectY = new Prs3d_ShadingAspect();
    aShadingAspectY->SetColor(Quantity_NOC_GREEN);
    aGroupY->SetGroupPrimitivesAspect(aShadingAspectY->Aspect());

    cout<<"void WBTrihedron::Compute()->____axis Y array of primitives done____"<<endl;

    //! ----------------------
    //! axis X - the cylinder
    //! ----------------------
    gp_Trsf aTrsf2;
    aTrsf2.SetRotation(gp_Ax1(gp_Pnt(0,0,0),gp_Dir(0,1,0)),pi/2);

    occHandle(Graphic3d_ArrayOfTriangles) AxisX;
    Prs3d_ToolCylinder toolCylinderX(cylinderRadius,cylinderRadius,height,NbSlices,NbStacks);
    toolCylinderX.FillArray(AxisX,polyTriangulationAxisX,aTrsf2);

    cout<<"void WBTrihedron::Compute()->____cylinder X done____"<<endl;

    //! -----------------
    //! axis X - the tip
    //! -----------------
    gp_Trsf t2;
    t2.SetTranslation(gp_Vec(0,0,height));
    aTrsf2.Multiply(t2);

    occHandle(Graphic3d_ArrayOfTriangles) AxisX_cone;
    Prs3d_ToolCylinder toolConeX(tipRadius,0,tipHeight,NbSlices,NbStacks);
    toolConeX.FillArray(AxisX_cone,polyTriangulationTipAxisX,aTrsf2);

    cout<<"void WBTrihedron::Compute()->____cylinder tip X done____"<<endl;

    //! -------------------------
    //! axis X - the hollow disk
    //! -------------------------
    occHandle(Graphic3d_ArrayOfTriangles) hollowDiskX;
    Prs3d_ToolDisk toolDiskX(cylinderRadius,tipRadius,NbSlices,NbStacks);
    toolDiskX.FillArray(hollowDiskX,polyTriangulationTipHollowDiskX,aTrsf2);

    cout<<"void WBTrihedron::Compute()->____hollow disk X done____"<<endl;

    occHandle(Graphic3d_Group) aGroupX = thePresentation->NewGroup();
    aGroupX->AddPrimitiveArray (AxisX);
    aGroupX->AddPrimitiveArray (AxisX_cone);
    aGroupX->AddPrimitiveArray (hollowDiskX);

    occHandle(Prs3d_ShadingAspect) aShadingAspectX = new Prs3d_ShadingAspect();
    aShadingAspectX->SetColor(Quantity_NOC_RED);
    aGroupX->SetGroupPrimitivesAspect(aShadingAspectX->Aspect());

    cout<<"void WBTrihedron::Compute()->____axis X array of primitives done____"<<endl;
}

//! ---------------------------
//! function: ComputeSelection
//! details:
//! ---------------------------
#include <SelectMgr_EntityOwner.hxx>
#include <Select3D_SensitiveTriangulation.hxx>
#include <Standard_PrimitiveTypes.hxx>
void WBTrihedron::ComputeSelection (const occHandle(SelectMgr_Selection)& theSel,
                                    const Standard_Integer aMode)
{
    cout<<"WBTrihedron::ComputeSelection()->____function called____"<<endl;
    bool isInterior = true;
    /*
    switch(aMode)
    {
    case 0:
    {
        cout<<"void WBTrihedron::ComputeSelection()->____mode whole trihedron____"<<endl;

        //! ------
        //! owner
        //! ------
        TopLoc_Location zLoc, yLoc, xLoc;
        occHandle(SelectMgr_EntityOwner) anOwner = new SelectMgr_EntityOwner();

        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesZ =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationAxisZ,zLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesZ1 =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationTipAxisZ,zLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesZ2 =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationTipHollowDiskZ,zLoc,isInterior);

        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesY =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationAxisY,yLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesY1 =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationTipAxisY,yLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesY2 =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationTipHollowDiskY,yLoc,isInterior);

        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesX =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationAxisY,xLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesX1 =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationTipAxisY,xLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesX2 =
                new Select3D_SensitiveTriangulation(anOwner,polyTriangulationTipHollowDiskY,xLoc,isInterior);


        theSel->Add(arrayOfTrianglesZ);
        theSel->Add(arrayOfTrianglesZ1);
        theSel->Add(arrayOfTrianglesZ2);
        theSel->Add(arrayOfTrianglesY);
        theSel->Add(arrayOfTrianglesY1);
        theSel->Add(arrayOfTrianglesY2);
        theSel->Add(arrayOfTrianglesX);
        theSel->Add(arrayOfTrianglesX1);
        theSel->Add(arrayOfTrianglesX2);

    }
        break;

    case 1:
    {
    */
        cout<<"void WBTrihedron::ComputeSelection()->____mode single axis____"<<endl;

        //! --------------------
        //! owner of the Z axis
        //! --------------------
        TopLoc_Location zLoc;
        occHandle(SelectMgr_EntityOwner) anOwner_Z_axis = new SelectMgr_EntityOwner();
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesZ =
                new Select3D_SensitiveTriangulation(anOwner_Z_axis,polyTriangulationAxisZ,zLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesZ1 =
                new Select3D_SensitiveTriangulation(anOwner_Z_axis,polyTriangulationTipAxisZ,zLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesZ2 =
                new Select3D_SensitiveTriangulation(anOwner_Z_axis,polyTriangulationTipHollowDiskZ,zLoc,isInterior);
        theSel->Add(arrayOfTrianglesZ);
        theSel->Add(arrayOfTrianglesZ1);
        theSel->Add(arrayOfTrianglesZ2);

        //! --------------------
        //! owner of the Y axis
        //! --------------------
        TopLoc_Location yLoc;
        occHandle(SelectMgr_EntityOwner) anOwner_Y_axis = new SelectMgr_EntityOwner();
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesY =
                new Select3D_SensitiveTriangulation(anOwner_Y_axis,polyTriangulationAxisY,yLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesY1 =
                new Select3D_SensitiveTriangulation(anOwner_Y_axis,polyTriangulationTipAxisZ,yLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesY2 =
                new Select3D_SensitiveTriangulation(anOwner_Y_axis,polyTriangulationTipHollowDiskY,yLoc,isInterior);
        theSel->Add(arrayOfTrianglesY);
        theSel->Add(arrayOfTrianglesY1);
        theSel->Add(arrayOfTrianglesY2);

        //! --------------------
        //! owner of the X axis
        //! --------------------
        TopLoc_Location xLoc;
        occHandle(SelectMgr_EntityOwner) anOwner_X_axis = new SelectMgr_EntityOwner();
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesX =
                new Select3D_SensitiveTriangulation(anOwner_X_axis,polyTriangulationAxisX,xLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesX1 =
                new Select3D_SensitiveTriangulation(anOwner_X_axis,polyTriangulationTipAxisX,xLoc,isInterior);
        occHandle(Select3D_SensitiveTriangulation) arrayOfTrianglesX2 =
                new Select3D_SensitiveTriangulation(anOwner_X_axis,polyTriangulationTipHollowDiskX,xLoc,isInterior);
        theSel->Add(arrayOfTrianglesX);
        theSel->Add(arrayOfTrianglesX1);
        theSel->Add(arrayOfTrianglesX2);
        /*
    }
        break;
    }
    */
}

/*
//! ---------------------------
//! function: setSelectionMode
//! details:
//! ---------------------------
void WBTrihedron::setSelectionMode(selectionMode aSelectionMode)
{
    this->setGlobalSelMode((int)aSelectionMode);
}
*/
