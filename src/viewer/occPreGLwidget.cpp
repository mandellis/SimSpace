//#define USE_FACE_DS_BUILDER

//! ----------------
//! custom includes
//! ----------------
#include <occPreGLwidget.h>
#include "simulationdatabase.h"
#include "mydefines.h"
#include "stlapiwriter.h"

#include "meshtools.h"
#include "exportingtools.h"
#include "stepimporter.h"
#include "ais_colorscaleextended.h"
#include "qbackgroundevent.h"
#include "ccout.h"
#include "contextmenubuilder.h"
#include "simulationmanager.h"
#include "tools.h"
#include "resultstoolbar.h"
#include "ais_meshsegmentmarker.h"
#include "modelloader.h"

#include <ng_meshvs_datasource1d.h>
#include <ng_meshvs_datasource3d.h>
#include <ng_meshvs_datasource2d.h>

#include <curvature.h>

//! ---
//! Qt
//! ---
#include <QActionGroup>
#include <QMessageBox>
#include <QMenu>
#include <QDebug>
#include <QFileDialog>
#include <QKeyEvent>

//! ----
//! OCC
//! ----
#include <AIS_ColoredShape.hxx>
#include <AIS_ListOfInteractive.hxx>
#include <AIS_ListIteratorOfListOfInteractive.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_CompSolid.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopoDS_Builder.hxx>
#include <TopExp_Explorer.hxx>
#include <MeshVS_DataSource.hxx>
#include <IGESControl_Reader.hxx>
#include <STEPControl_Reader.hxx>
#include <XSControl_WorkSession.hxx>
#include <Interface_InterfaceModel.hxx>
#include <StepRepr_Representation.hxx>
#include <XSControl_TransferReader.hxx>
#include <Interface_Static.hxx>
#include <TCollection_HAsciiString.hxx>
#include <TCollection_AsciiString.hxx>
#include <TColStd_ListIteratorOfListOfAsciiString.hxx>
#include <Prs3d_IsoAspect.hxx>
#include <MeshVS_MeshPrsBuilder.hxx>
#include <MeshVS_Drawer.hxx>
#include <MeshVS_DrawerAttribute.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <BRepExtrema_DistanceSS.hxx>
#include <BRepExtrema_SeqOfSolution.hxx>
#include <BRepExtrema_SolutionElem.hxx>
#include <BRepBndLib.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <NCollection_Sequence.hxx>
#include <BRepExtrema_SeqOfSolution.hxx>
#include <BRepExtrema_SupportType.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Geom_Surface.hxx>
#include <Extrema_ExtSS.hxx>
#include <BRepTools.hxx>
#include <gp_Pnt2d.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_HSurface.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TColStd_MapIteratorOfPackedMapOfInteger.hxx>
#include <Poly_Triangulation.hxx>
#include <BRepMesh_ShapeTool.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRep_Builder.hxx>
#include <Prs3d_PointAspect.hxx>
#include <Graphic3d_TransformPers.hxx>
#include <Quantity_NameOfColor.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <CPnts_AbscissaPoint.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GeomAdaptor_Curve.hxx>

using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
occPreGLWidget::occPreGLWidget(QWidget *parent):occGLWidget(parent),
    myCurWorkingMode(curWorkingMode_none),
    isMeshViewVolume(false)
{
    cout<<"occPreGLWidget::occPreGLWidget()->____CONSTRUCTOR CALLED____"<<endl;
    //this->setFocusPolicy(Qt::StrongFocus);
    //this->setFocusPolicy(Qt::ClickFocus);

    //if(occContext.IsNull()) this->init();
    //occMeshContext = new AIS_InteractiveContext(occViewer);

    this->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(this,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(ShowContextMenu1(const QPoint&)));
    connect(this,SIGNAL(selectionChanged()),this,SLOT(selectionProperties()));
    connect(this,SIGNAL(highlightmeshface()),this,SLOT(showFaceMesh()));
    connect(this,SIGNAL(highlightmeshedge()),this,SLOT(showEdgeMesh()));

    connect(this,SIGNAL(selectionChanged()),this,SLOT(printTopologyNumber()));    //! print the topology numbers

    //! -----------------------------
    //! actions for the context menu
    //! -----------------------------
    actionGenerateMeshOnSelectedBodies = new QAction("Generate mesh on selected bodies",this);
    actionGenerateMeshOnSelectedBodies->setIcon(QIcon(":/icons/icon_volume mesh.png"));
    actionGenerateMeshOnSelectedBodies->setData(16);

    actionPreviewSurfaceOnSelectedBodies = new QAction("Preview mesh on selected bodies",this);
    actionPreviewSurfaceOnSelectedBodies->setIcon(QIcon(":/icons/icon_surface mesh.png"));
    actionPreviewSurfaceOnSelectedBodies->setData(17);

    actionClearGeneratedDataOnSelectedBodies = new QAction("Clear generated data on selected bodies",this);
    actionClearGeneratedDataOnSelectedBodies->setIcon(QIcon(":/icons/icon_clear data.png"));
    actionClearGeneratedDataOnSelectedBodies->setData(18);
}

//! --------------------------------------------------------------
//! function: buildMeshToolsContextMenu
//! details:  build a piece of context menu containing mesh tools
//! --------------------------------------------------------------
void occPreGLWidget::buildMeshToolsContextMenu(QMenu *aContextMenu)
{
    QMenu *meshToolsMenu = aContextMenu->addMenu("Mesh tools");
    aContextMenu->addSeparator();

    //! -------------------------------------
    //! show the mesh for the selected faces
    //! -------------------------------------
    QAction *actionShowFaceMesh = meshToolsMenu->addAction("Show face mesh");
    actionShowFaceMesh->setIcon(QIcon(":/icons/icon_mesh face sizing.png"));
    actionShowFaceMesh->setData(300);

    //! ----------------------------------------------
    //! experimental - trigger the custom face mesher
    //! ----------------------------------------------
    QAction *actionMeshThisFace = meshToolsMenu->addAction("Mesh this face (experimental)");
    actionMeshThisFace->setIcon(QIcon(":/icons/icon_surface mesh.png"));
    actionMeshThisFace->setData(306);

    QAction *actionShowEdgeMesh = meshToolsMenu->addAction("Show edge mesh");
    actionShowEdgeMesh->setIcon(QIcon(":/icons/icon_edge mesh.png"));
    actionShowEdgeMesh->setData(301);

    QAction *actionShowHealingElements = meshToolsMenu->addAction("Show healing elements");
    actionShowHealingElements->setIcon(QIcon(":/icons/icon_healing patch.png"));
    actionShowHealingElements->setData(304);

    meshToolsMenu->addSeparator();

    QAction *actionShowElementsByType = meshToolsMenu->addAction("Show elements by type");
    actionShowElementsByType->setIcon(QIcon(":/icons/icon_show elements by type.png"));

    QMenu *menuShowElementByType = new QMenu(this);
    QAction *actionShowTets = menuShowElementByType->addAction("Show tetrahedra");
    actionShowTets->setIcon(QIcon(":/icons/icon_tetrahedron.png"));
    actionShowTets->setData(307);

    QAction *actionShowPyramids = menuShowElementByType->addAction("Show pyramids");
    actionShowPyramids->setIcon(QIcon(":/icons/icon_pyramid.png"));
    actionShowPyramids->setData(308);

    QAction *actionShowPrisms = menuShowElementByType->addAction("Show prisms");
    actionShowPrisms->setIcon(QIcon(":/icons/icon_prism.png"));
    actionShowPrisms->setData(309);

    QAction *actionShowHexa = menuShowElementByType->addAction("Show hexahedra");
    actionShowHexa->setIcon(QIcon(":/icons/icon_hexahedron.png"));
    actionShowHexa->setData(310);

    actionShowElementsByType->setMenu(menuShowElementByType);

    meshToolsMenu->addSeparator();

    QAction *actionExportSTL = meshToolsMenu->addAction("Export as STL");
    actionExportSTL->setIcon(QIcon(":/icons/icon_STL format.png"));
    actionExportSTL->setData(302);

    QAction *actionExportCloud = meshToolsMenu->addAction("Export cloud");
    actionExportCloud->setIcon(QIcon(":/icons/icon_cloud points"));
    actionExportCloud->setData(303);

    meshToolsMenu->addSeparator();

    QAction *actionDisplayCurvature = meshToolsMenu->addAction("Display curvature");
    actionDisplayCurvature->setIcon(QIcon(":/icons/icon_curvature.png"));
    actionDisplayCurvature->setData(305);
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
occPreGLWidget::~occPreGLWidget()
{
    cout<<"occPreGLWidget::~occPreGLWidget()->____DESTRUCTOR CALLED____"<<endl;
}

//! ----------------------------------
//! function: createInteractiveShapes
//! details:
//! ----------------------------------
void occPreGLWidget::createInteractiveShapes()
{
    //cout<"occPreGLWidget::createInteractiveShapes()->____function called____"<<endl;

    //! --------------------------
    //! number of 1D, 2D, 3D body
    //! --------------------------
    int N3D = myDS2->N3D();
    int N2D = myDS2->N2D();
    int N1D = myDS2->N1D();

    //! -----------------------
    //! total number of bodies
    //! -----------------------
    int Nt = N3D+N2D+N1D;

    if(Nt>0)
    {
        //! -------------------------------------------------
        //! retrieve the geometry from the geometry database
        //! -------------------------------------------------
        for(QMap<int,TopoDS_Shape>::iterator it = myDS2->bodyMap.begin(); it!=myDS2->bodyMap.end(); ++it)
        {
            //! the main shape of the body
            const TopoDS_Shape &aShape = it.value();

            //! body index
            int index = it.key();

            //! the interactive object
            occHandle(AIS_ExtendedShape) anAISShape = new AIS_ExtendedShape(aShape);

            //! sets the index of the shape - the AIS Shape and the underlying shape in the geometry data structure have the same index
            anAISShape->setIndex(index);

            //! set the name of the interactive shape
            anAISShape->setName(myDS2->MapOfBodyNames.value(index).toStdString().c_str());

            //! -------------------------------------
            //! fill the array of interactive shapes
            //! -------------------------------------
            myMapOfInteractiveShapes.insert(index, anAISShape);
            //std::pair<int,occHandle(AIS_ExtendedShape)> apair;
            //apair.first = index;
            //apair.second = anAISShape;
            //myMapOfInteractiveShapes.insert(apair);

            //! ------------------------------------------------------------
            //! fill the array of mesh interactive objects with null values
            //! ------------------------------------------------------------
            //std::pair<int,occHandle(MeshVS_Mesh)> apair_;
            //apair_.first = index;
            //apair_.second = occHandle(MeshVS_Mesh)();
            //myMapOfInteractiveMeshes.insert(apair_);
            myMapOfInteractiveMeshes.insert(index,occHandle(MeshVS_Mesh)());
        }
    }
}

//! ---------------------
//! function: displayCAD
//! details:
//! ---------------------
void occPreGLWidget::displayCAD(bool onlyLoad)
{
    cout<<"occPreGLWidget::displayCAD()->____function called____"<<endl;

    //! ------------------------------
    //! Sets the current working mode
    //! ------------------------------
    myCurWorkingMode = curWorkingMode_onModel;

    //! ------------------------------------------------------------------------------------
    //! When a CAD model is loaded and displayed, the viewer must be re-initialized.
    //! At the CAD loading time, the MainWindow resets the status of its toolbar buttons,
    //! but it cannot reinitialize the occGLWidget myCurSelectionMode, myCurAction3D,
    //! myCurGlobalSelectionMode because they are protected members: unsetViewOperations(),
    //! and unsetSelectionModes() are provided
    //! ------------------------------------------------------------------------------------
    this->unsetViewOperations();
    this->unsetSelectionModes();

    //! Reset the status of the buttons within the "Cursor <sub> menu"
    //! (this widget context menu)
    QList<QAction *>listOfActions = myCursorModeMenu->actions();
    for(int i=0;i<listOfActions.size();i++) listOfActions.value(i)->setChecked(false);

    //! This closes all contexts - returns to the NEUTRAL POINT
    occContext->CloseAllContexts(true);

    //! Set the number of the context - now "0", NEUTRAL POINT
    myLocalCtxNumber = occContext->IndexOfCurrentLocal();

    //! Not interested in selection in NEUTRAL POINT
    //occContext->SetAutoActivateSelection(Standard_False);

    //! Set the view mode
    this->setShadedExteriorAndEdgesView();

    //! Hide UV isolines
    occHandle(Prs3d_IsoAspect) isoLineAspect = new Prs3d_IsoAspect(Quantity_NOC_GRAY,Aspect_TOL_SOLID,1.0,1);
    isoLineAspect->SetNumber(0);

    //! -------------------
    //! display the shapes
    //! -------------------
    myCurDisplayMode = CurDisplayMode_ShadedExteriorAndEdges;
    occHandle(Graphic3d_AspectLine3d) la = new Graphic3d_AspectLine3d(static_cast<Quantity_Color>(Quantity_NOC_BLACK),Aspect_TOL_DASH ,2.0);
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin();
        it!=myMapOfInteractiveShapes.end(); ++it)

    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin();
    //    it!=myMapOfInteractiveShapes.end(); ++it)
    {
        //! the current AIS shape
        const occHandle(AIS_ExtendedShape) &anAISShape=occHandle(AIS_ExtendedShape)::DownCast(it.value());

        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //const occHandle(AIS_ExtendedShape) &anAISShape=occHandle(AIS_ExtendedShape)::DownCast(apair.second);

        //! set the visibility flag of the shape
        anAISShape->setShapeVisibility(Standard_True);

        //! -----------------------------------------------------------------
        //! set the color of the shape: use the interactive context function
        //! -----------------------------------------------------------------
        int k = it.key();
        //int k = apair.first;
        occContext->SetColor(anAISShape,myShapeColor.getColor(k),false);

        //! --------------------------------------------------------------------
        //! unset isolines: use the interactive object function: "deep" feature
        //! --------------------------------------------------------------------
        anAISShape->Attributes()->SetUIsoAspect(isoLineAspect);
        anAISShape->Attributes()->SetVIsoAspect(isoLineAspect);

        //! -------------
        //! transparency
        //! -------------
        //anAISShape->SetTransparency(TRANSPARENCY_IN_WORKING_MODE_MODEL);

        //! -----------------------------------------------------
        //! The selection or highlight is on the whole shape
        //! It can be: Aspect_TOHM_BOUNDBOX or Aspect_TOHM_COLOR
        //! -----------------------------------------------------
        //Aspect_TypeOfHighlightMethod typeOfHighlightMethod = Aspect_TOHM_COLOR;
        //anAISShape->Attributes()->SetMethod(typeOfHighlightMethod);

        if(onlyLoad)
        {

            //! only load into the context
            occContext->Load(anAISShape,-1,false);
        }
        else
        {
            //! call context Display()
            occContext->Display(anAISShape,false);

            double angle = anAISShape->Attributes()->DeviationAngle();
            double deviation = anAISShape->Attributes()->DeviationCoefficient();

            cout<<"@ Shape nr."<<k<<endl;
            cout<<"@ deviation angle: "<<angle<<endl;
            cout<<"@ deviation coefficient: "<<deviation<<endl;
        }
        cout<<"___AIS shape nr: "<<anAISShape->index()<<" has visibility "<<(anAISShape->isVisible()? "VISIBLE":"HIDDEN")<<"____"<<endl;
    }

    //! do not highlight a selected object
    occContext->SetToHilightSelected(false);

    //! fit the view
    occView->ZFitAll();
    occView->FitAll();

    //! open a local context for model access (selection, ...)
    myLocalCtxNumber=occContext->OpenLocalContext();

    //! set the selection style
    //setSelectionStyle(Quantity_NOC_GREEN,float(0.1));

    //! initialize the center of rotation using the camera target
    Standard_Real COV_Xin = occView->Camera()->Center().X();
    Standard_Real COV_Yin = occView->Camera()->Center().Y();
    Standard_Real COV_Zin = occView->Camera()->Center().Z();

    myCOR.SetX(COV_Xin);
    myCOR.SetY(COV_Yin);
    myCOR.SetX(COV_Zin);

    //! diagnostic
    //this->printSummary();
}

//! ----------------------------
//! function: getTopologyNumber
//! details:
//! ----------------------------
void occPreGLWidget::getTopologyNumber(const TopoDS_Shape &theParentShape,
                                       const TopoDS_Shape &theSubShape,
                                       int &parentShapeIndex,
                                       int &subShapeIndex)
{
    parentShapeIndex = myDS2->bodyMap.key(theParentShape);
    switch(theSubShape.ShapeType())
    {
    case(TopAbs_FACE):
        subShapeIndex = myDS2->MapOfBodyTopologyMap.value(parentShapeIndex).faceMap.FindIndex(theSubShape);
        break;
    case(TopAbs_EDGE):
        subShapeIndex = myDS2->MapOfBodyTopologyMap.value(parentShapeIndex).edgeMap.FindIndex(theSubShape);
        break;
    case(TopAbs_VERTEX):
        subShapeIndex = myDS2->MapOfBodyTopologyMap.value(parentShapeIndex).vertexMap.FindIndex(theSubShape);
        break;
    }
}

//! ----------------------------------------------------------
//! function: getNormal
//! details:  normal at the "center" of the (U, V) projection
//! ----------------------------------------------------------
static void getNormal(const TopoDS_Face& aFace, gp_Dir& aDNS)
{
    gp_Pnt aPoint;
    gp_Vec aD1U, aD1V;
    const occHandle(Geom_Surface) &aSurface = BRep_Tool::Surface(aFace);

    Standard_Real U1, U2, V1, V2;
    aSurface->Bounds(U1, U2, V1, V2);
    Standard_Real U = 0.5*(U1+U2);
    Standard_Real V = 0.5*(V1+V2);

    //! compute the first derivative (gradient)
    aSurface->D1(U, V, aPoint, aD1U, aD1V);

    //! directions
    gp_Dir aDD1U(aD1U);
    gp_Dir aDD1V(aD1V);
    aDNS=aDD1U^aDD1V;
    if(aFace.Orientation()==TopAbs_REVERSED)
    {
        aDNS.Reverse();
    }
}

//! --------------------------------------
//! function: normal to a face in a point
//! details:
//! --------------------------------------
static void getNormal(Standard_Real U, Standard_Real V,const TopoDS_Face& aFace, gp_Dir& aDNS)
{
    gp_Pnt aPoint;
    gp_Vec aD1U, aD1V;

    const occHandle(Geom_Surface) &aSurface = BRep_Tool::Surface(aFace);
    // computes the first derivative (gradient)
    aSurface->D1(U, V, aPoint, aD1U, aD1V);
    // directions
    gp_Dir aDD1U(aD1U);
    gp_Dir aDD1V(aD1V);
    aDNS=aDD1U^aDD1V;
    if(aFace.Orientation()==TopAbs_REVERSED)
    {
        aDNS.Reverse();
    }
}

//! ---------------------------------------------------------------
//! function: angleInRange
//! details:  check if the angle between normals is whitin a range
//!           the reverse() method has been introduced because the
//!           this method is used for contact search purposes
//! ---------------------------------------------------------------
static bool angleInRange(const gp_Dir &dir1, gp_Dir dir2)
{
    //! efficaci ragionamenti serali
    dir2.Reverse();
    //! fine degli efficaci ragionamenti serali

    const double limitAngle = 15;
    Standard_Real theDeviationAngle = dir1.Angle(dir2)*180/M_PI;
    //! check normals
    bool rv = (theDeviationAngle<=limitAngle) && (-limitAngle<=theDeviationAngle)? true: false;
    return rv;
}

//! ------------------------------------------------
//! function - static: utility; (x, y, z) to (U, V)
//! details:
//! ------------------------------------------------
static gp_Pnt2d FaceParameters(const TopoDS_Face &face,const gp_Pnt &pt)
{
    //! get face as surface
    const occHandle(Geom_Surface) &surface = BRep_Tool::Surface(face);
    //! create shape analysis object
    ShapeAnalysis_Surface sas(surface);
    //! get UV of point on surface
    gp_Pnt2d uv = sas.ValueOfUV(pt, 0.01);
    //! return parameters of point on face
    return uv;
}

//! ------------------------------------------------------------
//! function - static: utility; check it two faces are parallel
//! details:
//! ------------------------------------------------------------
static bool areParallel(const TopoDS_Shape &theFace1, const TopoDS_Shape &theFace2)
{
    TopoDS_Face face1 = TopoDS::Face(theFace1);
    TopoDS_Face face2 = TopoDS::Face(theFace2);
    BRepAdaptor_Surface Surf1(face1);
    BRepAdaptor_Surface Surf2(face2);

    occHandle(BRepAdaptor_HSurface) HS1 = new BRepAdaptor_HSurface(Surf1);
    occHandle(BRepAdaptor_HSurface) HS2 = new BRepAdaptor_HSurface(Surf2);
    Standard_Real Tol1 = Min(BRep_Tool::Tolerance(face1), Precision::Confusion());
    Tol1 = Min(Surf1.UResolution(Tol1), Surf1.VResolution(Tol1));
    Tol1 = Max(Tol1, Precision::PConfusion());

    Standard_Real U1_1, U2_1, V1_1, V2_1;
    Standard_Real U1_2, U2_2, V1_2, V2_2;
    BRepTools::UVBounds(face1, U1_1, U2_1, V1_1, V2_1);
    BRepTools::UVBounds(face2, U1_2, U2_2, V1_2, V2_2);
    Extrema_ExtSS myExtSS;
    myExtSS.Initialize(Surf2,U1_2,U2_2,V1_2,V2_2,Tol1);
    myExtSS.Perform(Surf1,U1_1,U2_1,V1_1,V2_1,Tol1);
    return(myExtSS.IsParallel());
}

//! ----------------
//! function: reset
//! details:
//! ----------------
void occPreGLWidget::reset()
{
    //! since there will be no body on the screen the item
    //! "pick point" in the context sub-menu must disappear
    myCursorModeMenu->removeAction(selectPickPointCoordinates);

    QList<QAction *>listOfActions = myCursorModeMenu->actions();

    //cout<<"____actual number of actions in the submenu: "<<listOfActions.size()<<endl;
    for(int i=0;i<listOfActions.size();i++)
    {
        listOfActions.value(i)->setChecked(false);
    }

    //! unset selection modes (internal OCC engine)
    unsetSelectionModes();

    //! unset the view operations (internal widget status)
    unsetViewOperations();

    //! close all the contexts and clear the scene
    occContext->CloseAllContexts(false);
    myLocalCtxNumber = occContext->IndexOfCurrentLocal();
    occContext->RemoveAll(true);

    //! Mesh context: closes all the contexts and clear the mesh scene
    occMeshContext->CloseAllContexts(false);
    occMeshContext->RemoveAll(true);

    //! clear the map of shapes
    myMapOfInteractiveShapes.clear();

    //! clear the map of meshes
    myMapOfInteractiveMeshes.clear();
}

//! ------------------------------
//! function: printTopologyNumber
//! details:
//! ------------------------------
void occPreGLWidget::printTopologyNumber()
{
    int index;
    for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
    {
        const TopoDS_Shape &theSelectedShape = occContext->SelectedShape();

        const occHandle(AIS_ExtendedShape) &theSelectedAIS_Shape=
                occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());

        int AIS_shapeIndex = theSelectedAIS_Shape->index();

        switch(theSelectedShape.ShapeType())
        {
        case TopAbs_SOLID:
            //! if the data structure contains only solids, when selecting
            //! a solid its index must be equal to the index of the AIS_ExtendedShape
            //! built with it
            index= myDS2->bodyMap.key(theSelectedShape);
            cout<<"AIS_Shape index: "<<AIS_shapeIndex<<" - TopoDS_Solid index: "<<index<<endl;
            break;

        case TopAbs_FACE:
            index= (myDS2->MapOfBodyTopologyMap.value(AIS_shapeIndex)).faceMap.FindIndex(theSelectedShape);
            cout<<"AIS_Shape index: "<<AIS_shapeIndex<<" - TopoDS_Face index: "<<index<<endl;
            break;

        case TopAbs_EDGE:
            index= (myDS2->MapOfBodyTopologyMap.value(AIS_shapeIndex)).edgeMap.FindIndex(theSelectedShape);
            cout<<"AIS_Shape index: "<<AIS_shapeIndex<<" - TopoDS_Edge index: "<<index<<endl;
            break;

        case TopAbs_VERTEX:
            index= (myDS2->MapOfBodyTopologyMap.value(AIS_shapeIndex)).vertexMap.FindIndex(theSelectedShape);
            cout<<"AIS_Shape index: "<<AIS_shapeIndex<<" - TopoDS_Vertex index: "<<index<<endl;
            break;
        }
    }
}

//! -------------------------------
//! function: identifyTheSelection
//! details:
//! -------------------------------
void occPreGLWidget::identifyTheSelection(int &vertexIndex, int &edgeIndex, int &faceIndex, int &solidIndex)
{
    for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
    {
        //! This is the selected topology
        const TopoDS_Shape &theSelectedShape = occContext->SelectedShape();

        //! This is the AIS_Shape providing the selected topology
        const occHandle(AIS_ExtendedShape) &theSelectedAIS_Shape=
                occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());

        //! This is the index of the AIS_ExtendedShape.
        //! Note: at the loading time the AIS_ExtendedShapes are built
        //! using the information contained in "myDS" (instance of the
        //! simulationDataStructure class). The index of the AIS_ExtendedShape
        //! is set as the same of the element in the map of 3D (2D, or 1D) bodies
        int index = theSelectedAIS_Shape->index();

        switch(theSelectedShape.ShapeType())
        {
        case TopAbs_SOLID:
            //! if the data structure contains only solids, when selecting
            //! a solid its index must be equal to the index of the AIS_ExtendedShape
            //! built with it
            solidIndex= myDS2->bodyMap.key(theSelectedShape);
            //cout<<"AIS_Shape index: "<<index<<" - TopoDS_Solid index: "<<solidIndex<<endl;
            break;

        case TopAbs_FACE:
            faceIndex= (myDS2->MapOfBodyTopologyMap.value(index)).faceMap.FindIndex(theSelectedShape);
            cout<<"AIS_SHape index: "<<index<<" - TopoDS_Face index: "<<faceIndex<<endl;
            break;

        case TopAbs_EDGE:
            edgeIndex= (myDS2->MapOfBodyTopologyMap.value(index)).edgeMap.FindIndex(theSelectedShape);
            //cout<<"AIS_Shape index: "<<index<<" - TopoDS_Edge index: "<<edgeIndex<<endl;
            break;

        case TopAbs_VERTEX:
            vertexIndex= (myDS2->MapOfBodyTopologyMap.value(index)).vertexMap.FindIndex(theSelectedShape);
            //cout<<"AIS_Shape index: "<<index<<" - TopoDS_Vertex index: "<<vertexIndex<<endl;
            break;
        }
    }
}

//! ---------------------------
//! function: displayAllMeshes
//! details:
//! ---------------------------
void occPreGLWidget::displayAllMeshes(bool meshNodesVisible, Graphic3d_NameOfMaterial theMaterial, bool updateViewers)
{
    this->hideAllMeshes();

    Graphic3d_MaterialAspect myAspect(theMaterial);
    MeshVS_DisplayModeFlags theMeshDisplayMode;

    //! ----------------------------------
    //! establish how to display the mesh
    //! ----------------------------------
    if(myCurDisplayMode==CurDisplayMode_Wireframe) theMeshDisplayMode = MeshVS_DMF_WireFrame;
    else theMeshDisplayMode = MeshVS_DMF_Shading;

    //! -----------------------------------------
    //! scan all the shapes, using signature "0"
    //! (only "geometry" entities)
    //! -----------------------------------------
    AIS_ListOfInteractive aList;
    occContext->ObjectsInside(aList,AIS_KOI_Shape,0);
    if(aList.IsEmpty()) return;

    //! ---------------------
    //! non-prismatic meshes
    //! ---------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it= myMapOfInteractiveMeshes.begin(); it!= myMapOfInteractiveMeshes.end(); ++it)
    {
        const occHandle(AIS_InteractiveObject) &curInteractiveMesh = it.value();
        if(curInteractiveMesh.IsNull()) continue;

        const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(myMapOfInteractiveShapes.value(it.key()));
        int bodyIndex = theAISShape->index();

        //! ----------------------------------------------------------------------------
        //! display the surface mesh if:
        //! 1. the shape is visible
        //! 2. the shape mesh has not prismatic faces: in case the shape is visible and
        //!    HAS prismatic faces, the prismatic mesh will be shown by the code [*]
        //! ----------------------------------------------------------------------------
        if(!theAISShape->isVisible()) continue;

        const occHandle(AIS_ColoredShape) &underlyingShape = occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
        const occHandle(MeshVS_Mesh) &aMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.value(bodyIndex));

        Quantity_Color aColor;
        underlyingShape->Color(aColor);
        myAspect.SetColor(aColor);
        aMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes,meshNodesVisible);

        aMesh->SetDisplayMode(theMeshDisplayMode);
        occMeshContext->RecomputePrsOnly(aMesh,false);

        if(aMesh==occHandle(MeshVS_Mesh)()) exit(8887);
        if(aMesh.IsNull())
        {
            cout<<"occPreGLWidget::displayAllMeshes()->____the mesh object is NULL____"<<endl;
            exit(8888);
        }

        cout<<"occPreGLWidget::displayAllMeshes()->____showing a mesh with: "<<aMesh->GetDataSource()->GetAllElements().Extent()<<" elements____"<<endl;
        cout<<"occPreGLWidget::displayAllMeshes()->____"<<aMesh->GetDataSource()->GetAllNodes().Extent()<<" nodes____"<<endl;

        cout<<"occPreGLWidget::displayAllMeshes()->____trying to display the mesh"<<endl;
        occMeshContext->Display(aMesh,false);
        cout<<"occPreGLWidget::displayAllMeshes()->____mesh displayed"<<endl;
    }

    /*
    //! ---------------------------------------
    //! handle the display 3D prismatic meshes
    //! ---------------------------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractivePrismaticMeshes.begin(); it!=myMapOfInteractivePrismaticMeshes.end(); ++it)
    {
        const occHandle(AIS_InteractiveObject) &curInteractiveMesh = it.value();
        if(!curInteractiveMesh.IsNull()) continue;
        const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(myMapOfInteractiveShapes.value(it.key()));
        int bodyIndex = theAISShape->index();
        if(!theAISShape->isVisible()) continue;

        const occHandle(AIS_ColoredShape) &underlyingShape = occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
        const occHandle(MeshVS_Mesh) &aMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractivePrismaticMeshes.value(bodyIndex));

        Quantity_Color aColor;
        underlyingShape->Color(aColor);
        myAspect.SetColor(aColor);
        aMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes,meshNodesVisible);
        aMesh->SetDisplayMode(theMeshDisplayMode);

        occMeshContext->RecomputePrsOnly(aMesh,false);
        occMeshContext->Display(aMesh,false);
    }
*/

    //! -------------------
    //! update the viewers
    //! -------------------
    if(updateViewers)
    {
        occContext->UpdateCurrentViewer();
        occMeshContext->UpdateCurrentViewer();
    }
}

//! -------------------------------
//! function: setWorkingMode_Model
//! details:  slot
//! -------------------------------
void occPreGLWidget::setWorkingMode_Model()
{
    if(myCurWorkingMode!=curWorkingMode_onModel)
    {
        cout<<"occPreGLWidget::setWorkingMode_Model()->____ON MODEL____"<<endl;

        //! --------------------------
        //! handle auto-triangulation
        //! --------------------------
        if(!occContext->DefaultDrawer()->IsAutoTriangulation())
            occContext->DefaultDrawer()->SetAutoTriangulation(Standard_True);

        //! ----------------
        //! hide the meshes
        //! ----------------
        if(myCurWorkingMode == curWorkingMode_onMesh) this->hideAllMeshes();

        //! ------------------------------------
        //! remove transparency from the bodies
        //! ------------------------------------
        bool isEnabled = false;
        bool updateViewer = true;
        this->setTransparency(isEnabled,updateViewer);

        //! ------------------------
        //! handle the display mode
        //! ------------------------
        switch(myCurDisplayMode)
        {
        case CurDisplayMode_ShadedExterior:
            occGLWidget::setShadedExteriorView();
            break;
        case CurDisplayMode_ShadedExteriorAndEdges:
            occGLWidget::setShadedExteriorAndEdgesView();
            break;
        case CurDisplayMode_Wireframe:
            occGLWidget::setWireframeView();
            break;
        }

        //! ------------------------
        //! Set the selection style
        //! ------------------------
        //this->setSelectionStyle(Quantity_NOC_GREEN,float(0.1));

        //! -----------------------------------------------------
        //! Change the internal status - update the working mode
        //! -----------------------------------------------------
        myCurWorkingMode = curWorkingMode_onModel;

        //! ------------------------------
        //! enable the selection buttons
        //! ------------------------------
        emit requestEnableSelectionButtons();

        //! --------------------------
        //! hide the results tool bar
        //! --------------------------
        ResultsToolBar *theResultsToolBar = static_cast<ResultsToolBar*>(tools::getWidgetByName("resultsToolBar"));
        theResultsToolBar->setVisible(false);
    }
}

//! -------------------------------------------------
//! function: setWorkingMode_Mesh
//! details:  slot - set working mode "mesh"
//!           Don't use display/hide function inside
//! -------------------------------------------------
void occPreGLWidget::setWorkingMode_Mesh()
{
    if(myCurWorkingMode!=curWorkingMode_onMesh)
    {
        cout<<"occPreGLWidget::setWorkingMode_Mesh()->____ON MESH____"<<endl;

        //! ----------------------------------------
        //! handle auto-triangulation
        //! replace the triangulation with the mesh
        //! ----------------------------------------
        if(occContext->DefaultDrawer()->IsAutoTriangulation())
        {
            cout<<"occPreGLWidget::setWorkingMode_Mesh()->____auto triangulation OFF____"<<endl;
            occContext->DefaultDrawer()->SetAutoTriangulation(false);
        }
        //this->replaceTriangulation();

        //! ------------------------------------
        //! remove transparency from the bodies
        //! ------------------------------------
        bool isEnabled = false;
        bool updateViewer = true;
        this->setTransparency(isEnabled,updateViewer);

        //! ---------------------------------------------------------------------------------
        //! display also the geometry in order to allow selection, but in wireframe mode
        //! note: displayAllMeshes() internally reads the status variable "myCurDisplayMode"
        //! ---------------------------------------------------------------------------------
        switch(myCurDisplayMode)
        {
        case CurDisplayMode_ShadedExteriorAndEdges:
            occGLWidget::setShadedExteriorAndEdgesView();
            break;
        case CurDisplayMode_ShadedExterior:
            occGLWidget::setShadedExteriorView();
            break;
        case CurDisplayMode_Wireframe:
            occGLWidget::setWireframeView();
            break;
        }
        occContext->UpdateCurrentViewer();

        //! ------------------------
        //! set the selection style
        //! ------------------------
        //setSelectionStyle(Quantity_NOC_GREEN,float(0.1));

        //! ------------------------------------------------------------
        //! change the internal status and enable the selection buttons
        //! ------------------------------------------------------------
        myCurWorkingMode = curWorkingMode_onMesh;
        emit requestEnableSelectionButtons();

        //! --------------------------
        //! hide the results tool bar
        //! --------------------------
        ResultsToolBar *theResultsToolBar = static_cast<ResultsToolBar*>(tools::getWidgetByName("resultsToolBar"));
        theResultsToolBar->setVisible(false);
    }
}

//! ------------------------------------------
//! function: setWireframwView
//! details:  both for shapes both for meshes
//! ------------------------------------------
void occPreGLWidget::setWireframeView()
{
    cout<<"occPreGLWidget::setWireframeView->____function called____"<<endl;

    //! --------------------------------------------------------------
    //! save the previous display mode and change the internal status
    //! --------------------------------------------------------------
    myCurDisplayMode_old = myCurDisplayMode;
    myCurDisplayMode = CurDisplayMode_Wireframe;

    //! ---------------------
    //! wireframe for shapes
    //! ---------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
        //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    {
        const occHandle(AIS_ColoredShape) &theShape = occHandle(AIS_ColoredShape)::DownCast(it.value());
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //const occHandle(AIS_ColoredShape) &theShape = occHandle(AIS_ColoredShape)::DownCast(apair.second);
        if(!theShape.IsNull()) occContext->SetColor(theShape,Quantity_NOC_BLACK,Standard_False);
    }
    occContext->SetDisplayMode(AIS_WireFrame,Standard_False);

    if(myCurWorkingMode==curWorkingMode_onMesh)
    {
        //! ---------------------
        //! wireframe for meshes
        //! ---------------------
        this->displayAllMeshes();
    }
    occContext->UpdateCurrentViewer();
    occMeshContext->UpdateCurrentViewer();
}

//! --------------------------------
//! function: setShadedExteriorView
//! details:
//! --------------------------------
void occPreGLWidget::setShadedExteriorView()
{
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    {
        int k = it.key();
        const occHandle(AIS_ColoredShape) &theShape = occHandle(AIS_ColoredShape)::DownCast(it.value());
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //int k = apair.first;
        //const occHandle(AIS_ColoredShape) &theShape = occHandle(AIS_ColoredShape)::DownCast(apair.second);
        if(!theShape.IsNull()) occContext->SetColor(theShape,myShapeColor.getColor(k),Standard_False);
    }
    occGLWidget::setShadedExteriorView();
}

//! ----------------------------------------
//! function: setShadedExteriorAndEdgesView
//! details:
//! ----------------------------------------
void occPreGLWidget::setShadedExteriorAndEdgesView()
{
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    {
        const occHandle(AIS_ColoredShape) &theShape = occHandle(AIS_ColoredShape)::DownCast(it.value());
        int k = it.key();
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //int k = apair.first;
        //const occHandle(AIS_ColoredShape) &theShape = occHandle(AIS_ColoredShape)::DownCast(apair.second);
        if(!theShape.IsNull()) occContext->SetColor(theShape,myShapeColor.getColor(k),Standard_False);
    }
    occGLWidget::setShadedExteriorAndEdgesView();
}

//! ----------------------------------------
//! function: setWorkingMode_NamedSelection
//! details:  slot
//! -----------------------------------------
void occPreGLWidget::setWorkingMode_NamedSelection()
{
    //! hide the mesh, when needed
    if(myCurWorkingMode == curWorkingMode_onMesh) this->hideAllMeshes();

    //! set the selection style
    //setSelectionStyle(Quantity_NOC_GREEN,float(0.1));

    //! store the previous view mode
    myCurDisplayMode_old = myCurDisplayMode;

    //! Set the current working mode - update the working mode
    myCurWorkingMode = curWorkingMode_onNamedSelections;

    //! enable the selection buttons
    emit requestEnableSelectionButtons();

    //! --------------------------
    //! hide the results tool bar
    //! --------------------------
    ResultsToolBar *theResultsToolBar = static_cast<ResultsToolBar*>(tools::getWidgetByName("resultsToolBar"));
    theResultsToolBar->setVisible(false);
}

//! ------------------------------
//! function: reset custom colors
//! details:
//! ------------------------------
void occPreGLWidget::resetCustomColors(bool updateViewer)
{
    AIS_ListOfInteractive aList;
    occContext->ObjectsInside(aList,AIS_KOI_Shape,-1);
    if(!aList.IsEmpty())
    {
        for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
        {
            const occHandle(AIS_ColoredShape) &curAISShape = occHandle(AIS_ColoredShape)::DownCast(it.value());
            if(!curAISShape.IsNull())
            {
                curAISShape->ClearCustomAspects();
                occContext->RecomputePrsOnly(curAISShape,Standard_False,Standard_False);
            }
        }
    }
    //! ------------------
    //! update the viewer
    //! ------------------
    if(updateViewer) occContext->UpdateCurrentViewer();
}

//! ------------------------------------
//! function: set working mode contacts
//! details:
//! ------------------------------------
void occPreGLWidget::setWorkingMode_Contacts()
{
    cout<<"occPreGLWidget::setWorkingMode_Contacts()->____ON CONTACT____"<<endl;

    //! -----------------------------
    //! Set the current working mode
    //! -----------------------------
    myCurWorkingMode = curWorkingMode_onContact;

    if(myCurDisplayMode_old != CurDisplayMode_Wireframe)
    {
        myCurDisplayMode_old = myCurDisplayMode;
        switch(myCurDisplayMode_old)
        {
        case CurDisplayMode_ShadedExterior:
            myCurDisplayMode = CurDisplayMode_ShadedExterior;
            this->setShadedExteriorView();
            break;
        case CurDisplayMode_ShadedExteriorAndEdges:
            myCurDisplayMode = CurDisplayMode_ShadedExteriorAndEdges;
            this->setShadedExteriorAndEdgesView();
            break;
        }
    }

    //! --------------------------
    //! handle auto-triangulation
    //! --------------------------
    if(!occContext->DefaultDrawer()->IsAutoTriangulation())
        occContext->DefaultDrawer()->SetAutoTriangulation(Standard_True);

    //! ------------------------
    //! set the selection style
    //! ------------------------
    //setSelectionStyle(Quantity_NOC_GREEN,float(0.1));

    //! --------------------
    //! hide all the meshes
    //! --------------------
    this->hideAllMeshes();

    //! ---------------------------------
    //! remove highlight from the bodies
    //! ---------------------------------
    this->unhighlightBody(false);

    //! ------------------------
    //! reset the custom colors
    //! ------------------------
    this->resetCustomColors();

    //! ------------------------------------------
    //! show all the AIS_Shapes with transparency
    //! ------------------------------------------
    this->setTransparency(true,true);

    //! -----------------------------
    //! enable the selection buttons
    //! -----------------------------
    emit requestEnableSelectionButtons();

    //! --------------------------
    //! hide the results tool bar
    //! --------------------------
    ResultsToolBar *theResultsToolBar = static_cast<ResultsToolBar*>(tools::getWidgetByName("resultsToolBar"));
    theResultsToolBar->setVisible(false);
}

//! ------------------------
//! function: hideAllMeshes
//! details:  slot
//! ------------------------
void occPreGLWidget::hideAllMeshes(bool updateViewer)
{
    occMeshContext->RemoveAll(updateViewer);
}

//!------------------------------
//! function: hideSelectedBodies
//! details:
//!-------------------------------
void occPreGLWidget::hideSelectedBodies()
{
    //cout<<"occPreGLWidget::hideSelectedBodies()->____function called____"<<endl;
    for(occContext->InitSelected();occContext->MoreSelected();occContext->InitSelected())
    {
        //! -----------------------------------------------------------------
        //! change the shape visibility. Set the flag BEFORE calling Erase()
        //! -----------------------------------------------------------------
        const occHandle(AIS_ExtendedShape) &theShape = (occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive()));
        theShape->setShapeVisibility(Standard_False);

        //! ----------------------------------------------------------------------------
        //! erase a MeshVS_Mesh interactive object, if present (volume or surface mesh)
        //! ----------------------------------------------------------------------------
        int index = theShape->index();
        if(!myDS2->ArrayOfMeshDS.value(index).IsNull() || !myDS2->ArrayOfMeshDS2D.value(index).IsNull())
        {
            occMeshContext->Erase(myMapOfInteractiveMeshes.value(index),false);
            //occMeshContext->Erase(myMapOfInteractiveMeshes.at(index),false);
        }
        occContext->Erase(occContext->SelectedInteractive(),Standard_False);
    }
    occContext->CloseLocalContext(occContext->IndexOfCurrentLocal());
    myLocalCtxNumber=occContext->OpenLocalContext();

    //! -----------------------------------------------------------------
    //! reactivate the current "Standard selection mode", since when
    //! the local context is closed, the selection modes are deactivated
    //! -----------------------------------------------------------------
    this->reactivateCurrentStandardSelectionMode();
}

//! --------------------------------------------------------------------
//! function: buildMeshIOs
//! details:  build the mesh interactive objects. The mesh data sources
//!           are read from the meshDataBase
//! --------------------------------------------------------------------
void occPreGLWidget::buildMeshIOs()
{
    cout<<"occPreGLWidget::buildMeshIOs()->____function called____"<<endl;
    if(myDS2->ArrayOfMeshDS.size()==0) return;

    //! -----------------------------------------
    //! use material "GOLD" and change the color
    //! -----------------------------------------
    Graphic3d_MaterialAspect anAspect(Graphic3d_NOM_GOLD);

    for(int k=1;k<=myDS2->bodyMap.size();k++)
    {
        const TopoDS_Shape &shape = myDS2->bodyMap.value(k);
        int bodyIndex = myDS2->bodyMap.key(shape);

        occHandle(MeshVS_DataSource) aMeshDS;
        if(!myDS2->ArrayOfMeshDS.value(bodyIndex).IsNull())
        {
            //! both surface and volume mesh are available for the current body
            if(isMeshViewVolume==false)
            {
                //! use the surface mesh for generating the AI mesh object
                aMeshDS = myDS2->ArrayOfMeshDS2D.value(bodyIndex);
            }
            else
            {
                //! use the volume mesh datasource for generating the AI mesh object
                aMeshDS = myDS2->ArrayOfMeshDS.value(bodyIndex);
                //cout<<"____number of nodes: "<<aMeshDS->GetAllNodes().Extent()<<"____"<<endl;
                //cout<<"____number of elements: "<<aMeshDS->GetAllElements().Extent()<<"____"<<endl;
            }
        }
        else
        {
            if(!myDS2->ArrayOfMeshDS2D.value(bodyIndex).IsNull())
            {
                //! only the surface mesh is available for the current body
                aMeshDS = myDS2->ArrayOfMeshDS2D.value(bodyIndex);
            }
        }

        if(!aMeshDS.IsNull())
        {
            const occHandle(MeshVS_Mesh) &aMesh = new MeshVS_Mesh();
            aMesh->SetDataSource(aMeshDS);

            //! ---------------------
            //! presentation builder
            //! ---------------------
            occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(aMesh);
            if(aBuilder.IsNull()) exit(8888);
            aMesh->AddBuilder(aBuilder,Standard_False);

            //! ---------------------------------------------------------------------
            //! retrieve the color of the underlying shape and assign it to the mesh
            //! ---------------------------------------------------------------------
            const occHandle(AIS_ColoredShape) &underlyingShape = occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.value(k));

            Quantity_Color aColor;
            underlyingShape->Color(aColor);
            anAspect.SetColor(aColor);

            //! --------------------------------------------------------------------------
            //! Display nodes - use it for checking if the mesh is actually 2-nd order...
            //! --------------------------------------------------------------------------
            aMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, Standard_False);
            aMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, Standard_True);
            aMesh->GetDrawer()->SetInteger(Aspect_TOM_BALL,1);
            aMesh->GetDrawer()->SetColor(MeshVS_BP_NodalColor,NODE_MESH_COLOR);

            //! --------------------------------------------------------------------
            //! Draw the edges of the volume and surface mesh with different colors
            //! --------------------------------------------------------------------
            if(!myDS2->ArrayOfMeshDS.value(bodyIndex).IsNull())
            {
                aMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,VOLUME_MESH_EDGECOLOR);
            }
            else
            {
                aMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,SURFACE_MESH_ONLY_EDGECOLOR);
            }

            aMesh->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial, anAspect);
            aMesh->SetDisplayMode(MeshVS_DMF_Shading);
            //aMesh->SetDisplayMode(MeshVS_DMF_WireFrame);
            //aMesh->SetDisplayMode(MeshVS_DMF_Shrink);

            //! ----------------------------------------------------------------
            //! Insert the mesh object into the array of Mesh IO for displaying
            //! ----------------------------------------------------------------
            if(myMapOfInteractiveMeshes.value(bodyIndex).IsNull()) myMapOfInteractiveMeshes.insert(bodyIndex,aMesh);
        }
        else
        {
            cout<<"occPreGLWidget::buildMeshIOs()->____inserting NULL meshIO for shape nr: "<<k<<"____"<<endl;
            myMapOfInteractiveMeshes.insert(bodyIndex,occHandle(MeshVS_Mesh)());
        }
    }
    cout<<"occPreGLWidget::buildMeshIOs()->____mesh objects created____"<<endl;
}

//! -----------------------------------------------
//! function: buildPrismaticMeshIO
//! details:  build the prismatic mesh for preview
//! -----------------------------------------------
void occPreGLWidget::buildPrismaticMeshIO()
{
    cout<<"occPreGLWidget::buildPrismaticMeshIO()->____function called____"<<endl;

    //! --------------------
    //! use material "GOLD"
    //! --------------------
    Graphic3d_MaterialAspect aspect(Graphic3d_NOM_GOLD);

    for(QMap<int,occHandle(MeshVS_DataSource)>::iterator it= myDS2->ArrayOfMeshDS.begin(); it!=myDS2->ArrayOfMeshDS.end(); ++it)
    {
        int bodyIndex = it.key();

        //! ---------------------------------
        //! retrieve the 3D mesh data source
        //! ---------------------------------
        occHandle(Ng_MeshVS_DataSource3D) meshDS = occHandle(Ng_MeshVS_DataSource3D)::DownCast(it.value());

        //! --------------------------------------------------
        //! check if inflation is defined on the current body
        //! --------------------------------------------------
        bool hasBodyInflation = myDS2->HasPrismaticFaces.value(bodyIndex);
        if(hasBodyInflation)
        {
            if(!meshDS.IsNull())
            {
                cout<<"occPreGLWidget::buildPrismaticMeshIO()->____building the prismatic mesh IO for body: "<<bodyIndex<<"____"<<endl;

                //! ----------------------------
                //! the interactive mesh object
                //! ----------------------------
                const occHandle(MeshVS_Mesh) &aMesh = new MeshVS_Mesh();
                aMesh->SetDataSource(meshDS);

                occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(aMesh);

                //! ---------------------
                //! presentation builder
                //! ---------------------
                aMesh->AddBuilder(aBuilder,Standard_False);

                //! ---------------------------------------------------------------------
                //! retrieve the color of the underlying shape and assign it to the mesh
                //! ---------------------------------------------------------------------
                const occHandle(AIS_ColoredShape) &underlyingShape = occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));

                Quantity_Color aColor;
                underlyingShape->Color(aColor);
                aspect.SetColor(aColor);
                aMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, Standard_False);
                aMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, Standard_True);
                aMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,VOLUME_MESH_EDGECOLOR);
                aMesh->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial, aspect);
                aMesh->SetDisplayMode(MeshVS_DMF_Shading);

                //! ----------------------------------------------------------------
                //! Insert the mesh object into the array of Mesh IO for displaying
                //! ----------------------------------------------------------------
                myMapOfInteractivePrismaticMeshes.insert(bodyIndex,aMesh);
            }
        }
    }
}

//!--------------------------------------------------
//! function: show all the bodies
//! details:  overrides occGLWidget::showAllBodies()
//!--------------------------------------------------
void occPreGLWidget::showAllBodies()
{
    cout<<"occPreGLWidget::showAllBodies()->____function called____"<<endl;

    //! ---------------------------------------------------------------------------------
    //! The DisplayAll() works in the neutral point, so the local context must be closed
    //! 25/08/2017: check ... old note
    //! ---------------------------------------------------------------------------------
    occContext->CloseLocalContext(occContext->IndexOfCurrentLocal());

    //! ------------------------------------
    //! show only the "active" bodies
    //! the bodies are shown wireframe mode
    //! ------------------------------------
    if(myCurWorkingMode == curWorkingMode_onSolution) occContext->SetDisplayMode(AIS_WireFrame,Standard_True);

    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    {
        const occHandle(AIS_ExtendedShape) &AIS = occHandle(AIS_ExtendedShape)::DownCast(it.value());
        int bodyIndex = myDS2->bodyMap.key(AIS->Shape());
        if(myDS2->MapOfIsActive.value(bodyIndex)==true)
        {
            occContext->Display(AIS,false);
            AIS->setShapeVisibility(Standard_True);
        }
    }
    myLocalCtxNumber = occContext->OpenLocalContext();
    this->reactivateCurrentStandardSelectionMode();

    AIS_ListOfInteractive thelistOfDIsplayed;
    occContext->DisplayedObjects(AIS_KOI_Shape,0,thelistOfDIsplayed, Standard_False);
    AIS_ListIteratorOfListOfInteractive it;
    for(it.Initialize(thelistOfDIsplayed);it.More();it.Next())
    {
        const occHandle(AIS_ExtendedShape) &curAISShape = occHandle(AIS_ExtendedShape)::DownCast(it.Value());
        curAISShape->setShapeVisibility(Standard_True);
        switch(myCurWorkingMode)
        {
        case curWorkingMode_onContact:
            occContext->SetTransparency(curAISShape,TRANSPARENCY_IN_WORKING_MODE_CONTACT,Standard_False);
            break;
        case curWorkingMode_onModel:
            occContext->SetTransparency(curAISShape,TRANSPARENCY_IN_WORKING_MODE_MODEL,Standard_False);
            break;
        case curWorkingMode_onMesh:
            occContext->SetTransparency(curAISShape,TRANSPARENCY_IN_WORKING_MODE_MESH,Standard_False);
            break;
        case curWorkingMode_onSolution:
            occContext->SetTransparency(curAISShape,TRANSPARENCY_IN_WORKING_MODE_MODEL,Standard_False);
            break;
        default:
            occContext->SetTransparency(curAISShape,TRANSPARENCY_IN_WORKING_MODE_MODEL,Standard_False);
            break;
        }
    }
    occView->ZFitAll();
    occContext->UpdateCurrentViewer();
    //occMeshContext->UpdateCurrentViewer();

    if(curAction3D()!=CurAction3D_Nothing) occContext->SetAutomaticHilight(Standard_False);

    //! ----------------------------------------------------------------
    //! this is a temporary solution to use hide/show context services
    //! on the mesh, showing them together with the underlying geometry
    //! ----------------------------------------------------------------
    Graphic3d_NameOfMaterial theMat = Graphic3d_NOM_GOLD;
    switch(myCurWorkingMode)
    {
    case curWorkingMode_onMesh: this->displayAllMeshes(false,theMat); break;
    case curWorkingMode_onSolution: this->hideAllMeshes(); break;
    }
    cout<<"occPreGLWidget::showAllBodies()->____function exiting____"<<endl;
}

//! ---------------------------------------------------------
//! function: refreshMeshView
//! details:  recompute the presentation of the mesh objects
//!           and immediately refresh the presentation
//! ---------------------------------------------------------
void occPreGLWidget::refreshMeshView(bool onlyExterior)
{
    //! ---------------------------
    //! update the status variable
    //! ---------------------------
    isMeshViewVolume = (onlyExterior==true? false:true);

    //! -----------------------------------
    //! reset the mesh interactive objects
    //! QMap performs insert and replace
    //! -----------------------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveMeshes.begin(); it!=myMapOfInteractiveMeshes.end(); ++it)
    {
        myMapOfInteractiveMeshes.insert(it.key(),occHandle(MeshVS_Mesh)());
    }
    this->buildMeshIOs();
    this->showAllBodies();
}

//! --------------------------------
//! function: hideAllTheOtherBodies
//! details:
//! --------------------------------
void occPreGLWidget::hideAllTheOtherBodies()
{
    AIS_ListOfInteractive list;
    AIS_ListIteratorOfListOfInteractive listIt;

    //! -----------------------------------------------------------------------------
    //! the list of the selected AIS objects: it can contains duplicates
    //! because the interactive object can be retrieved also from subshape selection
    //! -----------------------------------------------------------------------------
    for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
    {
        occHandle(AIS_ExtendedShape) aShape = occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());
        list.Append(aShape);
    }
    listIt.Initialize(list);

    //! Build The list of the AIS shapes without duplicates
    AIS_ListOfInteractive finalList;

    //! Put the first shape into the list
    const occHandle(AIS_ExtendedShape) &theFirst = occHandle(AIS_ExtendedShape)::DownCast(list.First());
    finalList.Append(theFirst);

    AIS_ListIteratorOfListOfInteractive itFinalList;

    //! Elimination of the duplicated AIS_Shapes in the list
    for(;listIt.More();listIt.Next())
    {
        itFinalList.Initialize(finalList);
        const occHandle(AIS_ExtendedShape) &outer = occHandle(AIS_ExtendedShape)::DownCast(listIt.Value());
        for(int k=0;itFinalList.More();itFinalList.Next())
        {
            const occHandle(AIS_ExtendedShape) &inner = occHandle(AIS_ExtendedShape)::DownCast(itFinalList.Value());
            if(outer!= inner)
            {
                k++;
            }
            //! This is equivalent to a multiple logical AND
            if(k==finalList.Extent())finalList.Append(outer);
        }
    }

    //! ----------------------------------------------------
    //! list of displayed AIS_Shapes with signature "Shape"
    //! ----------------------------------------------------
    AIS_ListOfInteractive theListOfDisplayedShapes;
    occContext->DisplayedObjects(AIS_KOI_Shape, 0, theListOfDisplayedShapes, Standard_False);

    AIS_ListIteratorOfListOfInteractive it;
    AIS_ListOfInteractive theListToBeHidden;
    AIS_ListIteratorOfListOfInteractive itListToBeHidden;

    for(it.Initialize(theListOfDisplayedShapes);it.More();it.Next())
    {
        int k = 0;
        for(itFinalList.Initialize(finalList);itFinalList.More();itFinalList.Next())
        {
            if(it.Value()!=itFinalList.Value())
            {
                k++;
            }
            //! This is equivalent to a multiple logical AND
            if(k==finalList.Extent())
            {
                theListToBeHidden.Append(it.Value());
            }
        }
    }

    //! Finally hide the shapes
    for(itListToBeHidden.Initialize(theListToBeHidden);itListToBeHidden.More();itListToBeHidden.Next())
    {
        const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(itListToBeHidden.Value());

        //! ---------------------------------------------------------------------
        //! the main shape of the AIS
        //! difference with respect to the base class method: hide also the mesh
        //! ---------------------------------------------------------------------
        const TopoDS_Shape &theShape = theAISShape->Shape();
        int bodyIndex = myDS2->bodyMap.key(theShape);
        if(!myDS2->ArrayOfMeshDS.value(bodyIndex).IsNull())
        {
            const occHandle(MeshVS_Mesh) &theMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.value(bodyIndex));
            //const occHandle(MeshVS_Mesh) &theMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.at(bodyIndex));
            occMeshContext->Erase(theMesh,Standard_False);
        }

        //! update the flag
        theAISShape->setShapeVisibility(Standard_False);
        occContext->Erase(itListToBeHidden.Value(),Standard_False);
    }
    occContext->CloseLocalContext(occContext->IndexOfCurrentLocal());
    myLocalCtxNumber=occContext->OpenLocalContext();
    this->reactivateCurrentStandardSelectionMode();
}

//! -----------------------
//! function: showFaceMesh
//! details:
//! -----------------------
#ifndef USE_FACE_DS_BUILDER
void occPreGLWidget::showFaceMesh()
{
    //cout<<"occPreGLWidget::showFaceMesh()->____function called____"<<endl;
    this->hideAllMeshes();

    std::vector<std::pair<int,int>> vecPairs;
    occContext->InitSelected();
    if(occContext->MoreSelected())
    {
        if(occContext->SelectedShape().ShapeType()==TopAbs_FACE)
        {
            std::pair<int,int> pair;
            for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
            {
                const occHandle(AIS_ExtendedShape) &AIS_SHape = occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());
                const TopoDS_Shape &theParentShape =AIS_SHape->Shape();
                const TopoDS_Shape &theFace = occContext->SelectedShape();
                this->getTopologyNumber(theParentShape,theFace,pair.first,pair.second);
                vecPairs.push_back(pair);
            }
        }
    }

    occHandle(Ng_MeshVS_DataSourceFace) summedDS;
    QList<occHandle(Ng_MeshVS_DataSourceFace)> listOfMeshDS;
    for(int k=0; k<vecPairs.size(); k++)
    {
        int bodyIndex = vecPairs.at(k).first;
        int faceIndex = vecPairs.at(k).second;

        const occHandle(Ng_MeshVS_DataSourceFace) &curFaceDS = occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDS2->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
        if(curFaceDS.IsNull()) continue;
        listOfMeshDS.append(curFaceDS);

        cout<<"occPreGLWidget::showFaceMesh()->____number of elements: "<<curFaceDS->GetAllElements().Extent()<<"____"<<endl;
        cout<<"occPreGLWidget::showFaceMesh()->____number of nodes: "<<curFaceDS->GetAllNodes().Extent()<<"____"<<endl;
    }

    occHandle(Ng_MeshVS_DataSourceFace) aFaceDS = new Ng_MeshVS_DataSourceFace(listOfMeshDS);
    if(aFaceDS.IsNull())
    {
        QMessageBox::warning(this,APPNAME,"Error in merging the selected face meshes",QMessageBox::Ok);
        return;
    }

    aFaceDS->writeMesh("D:/testMesh.txt",2);

    //! the mesh interactive object
    occHandle(MeshVS_Mesh) aFaceMesh = new MeshVS_Mesh();
    aFaceMesh->SetDataSource(aFaceDS);

    //! create and add the presentation builder
    occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(aFaceMesh);
    aFaceMesh->AddBuilder(aBuilder,false);

    //! cosmetic for displaying the face mesh
    Graphic3d_MaterialAspect myAspect(Graphic3d_NOM_GOLD);
    myAspect.SetColor(Quantity_NOC_RED);

    aFaceMesh->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial,myAspect);
    aFaceMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes,true);
    aFaceMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges,true);
    aFaceMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
    aFaceMesh->SetDisplayMode(MeshVS_DMF_Shading);

    occMeshContext->Display(aFaceMesh,false);
    occContext->UpdateCurrentViewer();
    emptyTheSelection();
}
#endif

#ifdef USE_FACE_DS_BUILDER
#include <facedatasourcebuilder.h>
#include <indexedmapofmeshdatasources.h>
void occPreGLWidget::showFaceMesh()
{
    cout<<"occPreGLWidget::showFaceMesh()->____function called____"<<endl;
    QList<TopoDS_Face> listOfFaces;
    for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
    {
        const TopoDS_Shape &selectedShape = occContext->SelectedShape();
        if(selectedShape.ShapeType()==TopAbs_FACE)
        {
            //cout<<"occPreGLWidget::showFaceMesh()->____adding a face____"<<endl;
            listOfFaces<<TopoDS::Face(occContext->SelectedShape());
        }
    }
    if(listOfFaces.isEmpty()) return;

    cout<<"occPreGLWidget::showFaceMesh()->____number of selected faces: "<<listOfFaces.length()<<"____"<<endl;

    faceDataSourceBuilder aFaceDSBuilder;
    aFaceDSBuilder.setFaces(listOfFaces);
    aFaceDSBuilder.setDataBase(myDS2);

    IndexedMapOfMeshDataSources dataSources;
    aFaceDSBuilder.perform(dataSources);
    cout<<"occPreGLWidget::showFaceMesh()->____the data sources have been built____"<<endl;

    //! -----------------------------
    //! display the mesh datasources
    //! -----------------------------
    for(IndexedMapOfMeshDataSources::iterator it=dataSources.begin(); it!=dataSources.end(); ++it)
    {
        //int bodyIndex = it.key();
        const occHandle(MeshVS_DataSource) &faceDS = it.value();

        //! ----------------------------
        //! the mesh interactive object
        //! ----------------------------
        occHandle(MeshVS_Mesh) meshIO = new MeshVS_Mesh();
        meshIO->SetDataSource(faceDS);

        //! ----------------------------------------
        //! create and add the presentation builder
        //! ----------------------------------------
        occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(meshIO);
        meshIO->AddBuilder(aBuilder,false);

        //! --------------------------------------
        //! cosmetic for displaying the face mesh
        //! --------------------------------------
        Graphic3d_MaterialAspect aspect(Graphic3d_NOM_GOLD);
        aspect.SetColor(Quantity_NOC_RED);

        meshIO->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial,aspect);
        meshIO->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes,false);
        meshIO->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges,true);
        meshIO->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
        meshIO->SetDisplayMode(MeshVS_DMF_Shading);

        occMeshContext->Display(meshIO,false);
    }
    occMeshContext->UpdateCurrentViewer();
    this->emptyTheSelection();
}
#endif

//! -------------------------------
//! function: showEdgeMesh
//! details:  for testing purposes
//! -------------------------------
void occPreGLWidget::showEdgeMesh()
{
    cout<<"occPreGLWidget::showEdgeMesh()->____function called____"<<endl;
    this->hideAllMeshes();

    std::vector<std::pair<int,int>> vecPairs;
    occContext->InitSelected();
    if(occContext->MoreSelected())
    {
        std::pair<int,int> pair;
        if(occContext->SelectedShape().ShapeType()==TopAbs_EDGE)
        {
            for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
            {
                const occHandle(AIS_ExtendedShape) &AIS_SHape = occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());
                const TopoDS_Shape &theParentShape =AIS_SHape->Shape();
                const TopoDS_Shape &theChildShape = occContext->SelectedShape();
                int parentShapeIndex, subShapeIndex;
                this->getTopologyNumber(theParentShape,theChildShape,parentShapeIndex,subShapeIndex);
                pair.first=parentShapeIndex;
                pair.second=subShapeIndex;
                vecPairs.push_back(pair);
            }
        }
    }

    for(int k=0; k<vecPairs.size(); k++)
    {
        std::pair<int,int> curPair = vecPairs.at(k);
        int parentShapeIndex = curPair.first;
        int childShapeIndex = curPair.second;
        const occHandle(MeshVS_DataSource) &anEdgeDS = myDS2->ArrayOfMeshDSOnEdges.getValue(parentShapeIndex,childShapeIndex);
        occHandle(MeshVS_Mesh) anEdgeMesh = new MeshVS_Mesh();
        anEdgeMesh->SetDataSource(anEdgeDS);
        occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(anEdgeMesh);
        anEdgeMesh->AddBuilder(aBuilder,false);
        Graphic3d_MaterialAspect myAspect(Graphic3d_NOM_GOLD);
        myAspect.SetColor(Quantity_NOC_RED);
        anEdgeMesh->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, true);
        anEdgeMesh->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, true);
        anEdgeMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_RED);
        if(anEdgeMesh.IsNull())
        {
            cout<<"NULL edge mesh"<<endl;
        }
        else occMeshContext->Display(anEdgeMesh,false);
    }
    occMeshContext->UpdateCurrentViewer();
    emptyTheSelection();
}


//! ----------------------
//! function: getMeshSize
//! details:
//! ----------------------
void occPreGLWidget::getMeshSize()
{
    int NN =0;
    int NE =0;
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!= myMapOfInteractiveShapes.end(); ++it)
    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!= myMapOfInteractiveShapes.end(); ++it)
    {
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(apair.second);
        const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(it.value());
        int bodyIndex = theAISShape->index();
        //! ask for the number of nodes and elements only if the mesh data structure is not null
        if(!myDS2->ArrayOfMeshDS.value(bodyIndex).IsNull())
        {
            NN = NN + myDS2->ArrayOfMeshDS.value(bodyIndex)->GetAllNodes().Extent();
            NE = NE + myDS2->ArrayOfMeshDS.value(bodyIndex)->GetAllElements().Extent();
        }
    }
    QString msg = QString("Number of nodes: %1 \nNumber of elements: %2").arg(NN).arg(NE);
    QMessageBox::information(this,APPNAME,msg,QMessageBox::Ok);
}

//! ---------------------------
//! function: ShowContextMenu1
//! details:  slot
//! ---------------------------
void occPreGLWidget::ShowContextMenu1(const QPoint& pos)
{
    //! this avoids crashes when nothing has been loaded into the context
    if(myCurWorkingMode!=curWorkingMode_none)
    {
        //! clear the context menu
        myContextMenu->clear();

        //! position
        QPoint globalPos = this->mapToGlobal(pos);

        SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
        bool isEnabled = sm->isSomethingRunning()? false:true;
        QModelIndex curIndex = sm->myTreeView->currentIndex();

        if(curIndex.isValid())
        {
            //! -------------------------------
            //! retrieve the current node type
            //! -------------------------------
            SimulationNodeClass *node = curIndex.data(Qt::UserRole).value<SimulationNodeClass*>();
            SimulationNodeClass::nodeType type = node->getType();
            SimulationNodeClass::nodeType family = node->getFamily();

            switch(family)
            {
            case SimulationNodeClass::nodeType_structuralAnalysis:
                contextMenuBuilder::buildStructuralAnalysisContextMenu(myContextMenu,false,isEnabled);
                contextMenuBuilder::buildThermalAnalysisContextMenu(myContextMenu,true,isEnabled);
                break;

            case SimulationNodeClass::nodeType_meshControl:
            {
                contextMenuBuilder::buildMeshContextMenu(myContextMenu,false,isEnabled);

                //! ---------------
                //! add mesh tools
                //! ---------------
                if(isEnabled) this->buildMeshToolsContextMenu(myContextMenu);
            }
                break;

            case SimulationNodeClass::nodeType_geometry:
            {
                switch(type)
                {
                case SimulationNodeClass::nodeType_geometry:
                    contextMenuBuilder::buildModelRootContextMenu(myContextMenu,false,isEnabled);
                    break;
                case SimulationNodeClass::nodeType_geometryBody:
                    contextMenuBuilder::buildGeometryContextMenu(myContextMenu,false,isEnabled);
                    break;
                }
            }
                break;

            case SimulationNodeClass::nodeType_connection:
                contextMenuBuilder::buildContactContextMenu(myContextMenu,false,isEnabled);
                break;
            }
        }

        //! -------------------------------
        //! build the suppression menu
        //! -------------------------------
        if(isEnabled) this->buildSuppressionContextMenu();

        //! -----------
        //! basic menu
        //! -----------
        this->buildMinimalContexetMenu();

        //! -----------------------------------------------
        //! handling the selected action
        //!
        //! remote displacement                         28
        //! remote force                                29
        //! fixed support                               31
        //! cylindrical support                         32
        //! frictionless support                        33
        //! force                                       34
        //! moment                                      35
        //! pressure                                    36
        //! displacement                                37
        //! thermal condition                           38
        //! standard earth gravity                      39
        //! acceleration                                10
        //! rotational velocity                         11
        //!
        //! suppress                                    57
        //! unsuppress                                  58
        //! invert suppression set                      59
        //! suppress all other bodies                   60
        //! unsuppress all bodies                       61
        //! insert repair tool                          62
        //! hide body                                   63
        //! hide all other bodies                       64
        //!
        //! insert remote point root                    45
        //! insert remote point                         46
        //! insert coordinate system                    47
        //!
        //! flip contact                                48
        //! insert contact                              49
        //! insert contact group                        44
        //! generate automatic connections              43
        //!
        //! insert mesh method                          50
        //! insert body sizing                          51
        //! insert new mesh method                      52
        //! insert edge sizing                          59
        //! generate volume mesh                        54
        //! generate surface mesh                       55
        //! clear all meshes                            56
        //!
        //! insert named selection                      70
        //! insert face sizing                          72
        //! insert vertex sizing                        69
        //!
        //! insert imported body temperature dist       53
        //! insert initial stress                       80
        //! insert initial strain                       81
        //! insert external data                        71
        //! insert openfoam scalar field translator    100
        //! rename                                     101
        //! duplicate                                   98
        //! delete                                      99
        //!
        //! update post object                         200
        //! insert total displacement                  201
        //! insert directional displacement            206
        //! insert stress equivalent                   202
        //! insert maximum principal                   207
        //! insert middle principal                    208
        //! insert minimum principal                   209
        //!
        //! insert equivalent strain                   203
        //! insert strain intensity                    214
        //! insert max principal strain                211
        //! insert middle principal strain             212
        //! insert minimum principal strain            213
        //! evaluate result                            204
        //! clear generated data                       205
        //!
        //! show face mesh                             300
        //! show edge mesh                             301
        //! export STL                                 302
        //! export cloud                               303
        //! show healing elements                      304
        //! display curvature                          305
        //! discretize the face                        306
        //! show tetrahedra                            307
        //! show pyramids                              308
        //! show prism                                 309
        //! -----------------------------------------------
        QAction* selectedItem = myContextMenu->exec(globalPos);

        //! check if the selection is valid
        if(selectedItem)
        {
            cout<<"____action nr. "<<selectedItem->data().toInt()<<" called____"<<endl;
            switch(selectedItem->data().toInt())
            {
            //! For cases {0, 1, 2, 3, 4} which handle the change of the selection mode
            //! a signal is emitted: it is received by the MainWindow class, which in turn
            //! calls the slots toggle<..>SelectionMode(), which in turn calls the slot
            //! slot setSelectionMode(). Passing through the MainWindow is for
            //! handling the status (checked/unchecked) of the buttons in the toolbar
            case 0:
                //! selection mode vertex
                emit selectionModeVertex(true);
                break;
            case 1:
                //! selection mode edge
                emit selectionModeEdge(true);
                break;
            case 2:
                //! selection mode face
                emit selectionModeFace(true);
                break;
            case 3:
                //! selection mode solid
                emit selectionModeSolid(true);
                break;
            case 4:
                //! pick coordinates mode (the TopABS_FACE selection mode is activated)
                emit selectionModePickPointCoordinates(true);
                break;
            case 5:
                occView->SetProj(V3d_Xpos);
                this->FitAll();
                break;
            case 6:
                occView->SetProj(V3d_Xneg);
                this->FitAll();
                break;
            case 7:
                occView->SetProj(V3d_Ypos);
                this->FitAll();
                break;
            case 8:
                occView->SetProj(V3d_Yneg);
                this->FitAll();
                break;
            case 9:
                occView->SetProj(V3d_Zpos);
                this->FitAll();
                break;
            case 10:
                occView->SetProj(V3d_Zneg);
                this->FitAll();
                break;
            case 11:
                this->isometricView();
                break;
            case 12:
                //cout<<"occGLWidget::ShowContextMenu->____context menu request____ZOOM TO FIT (FIT ALL)____"<<endl;
                //! Please, do not use this->FitAll() because this internally changes the "myCurrentActio3D" value,
                //! In that case, when the fit has been performed, the current action 3D button (if anyone is pressed)
                //! is pressed, but it would not correspond to an activated panning, rotating, or window zooming mode
                occView->FitAll();

                switch(myCurAction3D)
                {
                case(CurAction3D_Rotation):
                    setAction3D_Rotation();
                    break;
                case(CurAction3D_Panning):
                    setAction3D_Pan();
                    break;
                case(CurAction3D_WindowZooming):
                    setAction3D_WindowZooming();
                    break;
                }
                break;
            case 13:
                //! important: emit signal AFTER
                this->hideSelectedBodies();
                emit requestSynchVisibility();
                break;
            case 14:
                //! important: emit signal AFTER
                hideAllTheOtherBodies();
                emit requestSynchVisibility();
                break;
            case 15:
                this->showAllBodies();
                //! important: emit signal AFTER
                emit requestSynchVisibility();
                break;
            case 38:
                //! thermal condition
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisThermalCondition);
                break;
            case 40:    //OK
                //! create an acceleration
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Acceleration);
                break;
            case 41:    //OK
                //! create an acceleration
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RotationalVelocity);
                break;
            case 44:    //OK
                //! insert a connection group
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_connectionGroup);
                break;
            case 49:    //OK
                //! insert a contact pair
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_connectionPair);
                break;
            case 54:    //OK
                //! generate the 3D mesh
                emit requestGenerateMesh(true);
                break;
            case 55:    //OK
                //! generate the surface mesh
                emit requestGenerateMesh(false);
                break;
            case 56:   //OK
                //! clear mesh from viewer
                this->clearMeshFromViewer();
                break;
            case 205:
                //! clear generated data
                emit requestClearGeneratedData();
                break;
            case 50:    //OK
                //! insert mesh method
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_meshBodyMeshMethod);
                break;
            case 51:    //OK
                //! insert body mesh sizing
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_meshBodyMeshControl);
                break;
            case 74:    //OK
                //!inset new mesh method
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_meshMethod);
                break;
            case 71:    //OK
                //! insert a mapper root
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_mapper);
                break;
            case 72:    //OK
                //! insert face mesh sizing;
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_meshFaceSize);
                break;
            case 76:
                //! export step file
                cout<<"____sending signal 76____"<<endl;
                emit requestExportStepFile();
                break;
            case 68:    //OK
                //! insert edge mesh sizing;
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_meshEdgeSize);
                break;
            case 69:    //OK
                //! insert vertex sizing;
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_meshVertexSize);
                break;
            case 57:    //OK
                //! suppress body
                emit requestChangeSuppressionStatus(Property::SuppressionStatus_Suppressed);
                break;
            case 23:
                //! suppress all other bodies
                suppressAllOtherBodies();
                break;
            case 24:
                //! unsuppress bodies
                emit requestChangeSuppressionStatus(Property::SuppressionStatus_Active);
                break;
            case 206:   //OK
                //! insert boundary condition displacement
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement);
                break;
            case 28:    //
                //! insert boundary condition remote displacement
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteDisplacement);
                break;
            case 29:
                //! insert boundary condition remote force
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_RemoteForce);
                break;
            case 31:
                //! insert boundary condition fixed support
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryContidion_FixedSupport);
                break;
            case 32:    //OK
                //! insert boundary condition cylindrical support
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_CylindricalSupport);
                break;
            case 33:    //OK
                //! insert boundary coondition frictionless support
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_FrictionlessSupport);
                break;
            case 34:    //OK
                //! insert boundary condition force
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Force);
                break;
            case 35:    //OK
                //! insert boundary condition moment
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Moment);
                break;
            case 36:    //OK
                //! insert boundary condition pressure
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Pressure);
                break;
            case 37:    //OK
                //! insert boundary condition displacement
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_structuralAnalysisBoundaryCondition_Displacement);
                break;
            case 90:    //OK
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisTemperature);
                break;
            case 91:    //OK
                emit requestCreateSimulationNode(SimulationNodeClass::nodeType_thermalAnalysisConvection);
                break;
            case 300:   //OK
                //! show face mesh
                this->showFaceMesh();
                break;
            case 301:   //OK
                //! show edge mesh
                this->showEdgeMesh();
                break;
            case 302:   //OK
                //! export STL
                this->exportSTL();
                break;
            case 303:   //OK
                //! export cloud
                this->exportCloud();
                break;
            case 304:   //OK
                //! show healing elements
                emit requestShowHealingElements();
                break;
            case 305:   //OK
                this->displayCurvatureMap();
                break;
            case 306:
            {
                TopoDS_Face theSelectedFace;
                for(occContext->InitSelected(); occContext->MoreSelected(); occContext->NextSelected())
                {
                    // take only the last face
                    theSelectedFace = TopoDS::Face(occContext->SelectedShape());
                }
                emit requestBuildFaceMesh(theSelectedFace);
            }
                break;

            case 307:
                //! ------------------
                //! display only tets
                //! ------------------
                for(occContext->InitSelected(); occContext->MoreSelected(); occContext->NextSelected())
                {
                    const TopoDS_Shape &aSolid = occContext->SelectedShape();
                    if(aSolid.ShapeType()==TopAbs_SOLID)
                    {
                        int bodyIndex = myDS2->bodyMap.key(aSolid);
                        const occHandle(Ng_MeshVS_DataSource3D) &mesh3D = occHandle(Ng_MeshVS_DataSource3D)::DownCast(
                                    myDS2->ArrayOfMeshDS.value(bodyIndex));
                        occHandle(Ng_MeshVS_DataSource3D) outputMesh;
                        MeshTools::filterVolumeElementsByType(mesh3D,TET,outputMesh);
                        if(!outputMesh.IsNull())
                        {
                            cout<<"____number of tets: "<<outputMesh->GetAllElements().Extent()<<"____"<<endl;
                            cout<<"____number of nodes: "<<outputMesh->GetAllNodes().Extent()<<"____"<<endl;
                        }
                        occPreGLWidget::displayMesh(occMeshContext,Quantity_NOC_RED,outputMesh);
                    }
                }
                this->clearSelected();
                break;

            case 308:
                //! ----------------------
                //! display only pyramids
                //! ----------------------
                for(occContext->InitSelected(); occContext->MoreSelected(); occContext->NextSelected())
                {
                    const TopoDS_Shape &aSolid = occContext->SelectedShape();
                    if(aSolid.ShapeType()==TopAbs_SOLID)
                    {
                        int bodyIndex = myDS2->bodyMap.key(aSolid);
                        const occHandle(Ng_MeshVS_DataSource3D) &mesh3D = occHandle(Ng_MeshVS_DataSource3D)::DownCast(
                                    myDS2->ArrayOfMeshDS.value(bodyIndex));
                        occHandle(Ng_MeshVS_DataSource3D) outputMesh;
                        MeshTools::filterVolumeElementsByType(mesh3D,PYRAM,outputMesh);
                        if(!outputMesh.IsNull())
                        {
                            cout<<"____number of pyramids: "<<outputMesh->GetAllElements().Extent()<<"____"<<endl;
                            cout<<"____number of nodes: "<<outputMesh->GetAllNodes().Extent()<<"____"<<endl;
                            occPreGLWidget::displayMesh(occMeshContext,Quantity_NOC_YELLOW,outputMesh);
                        }
                    }
                }
                this->clearSelected();
                break;

            case 309:
                //! --------------------
                //! display only prisms
                //! --------------------
                for(occContext->InitSelected(); occContext->MoreSelected(); occContext->NextSelected())
                {
                    const TopoDS_Shape &aSolid = occContext->SelectedShape();
                    if(aSolid.ShapeType()==TopAbs_SOLID)
                    {
                        int bodyIndex = myDS2->bodyMap.key(aSolid);
                        const occHandle(Ng_MeshVS_DataSource3D) &mesh3D = occHandle(Ng_MeshVS_DataSource3D)::DownCast(
                                    myDS2->ArrayOfMeshDS.value(bodyIndex));
                        occHandle(Ng_MeshVS_DataSource3D) outputMesh;
                        MeshTools::filterVolumeElementsByType(mesh3D,PRISM,outputMesh);
                        if(!outputMesh.IsNull())
                        {
                            cout<<"____number of prisms: "<<outputMesh->GetAllElements().Extent()<<"____"<<endl;
                            cout<<"____number of nodes: "<<outputMesh->GetAllNodes().Extent()<<"____"<<endl;
                            occPreGLWidget::displayMesh(occMeshContext,Quantity_NOC_VIOLET,outputMesh);
                        }
                    }
                }
                this->clearSelected();
                break;

            case 310:
                //! ------------------
                //! display only hexa
                //! ------------------
                for(occContext->InitSelected(); occContext->MoreSelected(); occContext->NextSelected())
                {
                    const TopoDS_Shape &aSolid = occContext->SelectedShape();
                    if(aSolid.ShapeType()==TopAbs_SOLID)
                    {
                        int bodyIndex = myDS2->bodyMap.key(aSolid);
                        const occHandle(Ng_MeshVS_DataSource3D) &mesh3D = occHandle(Ng_MeshVS_DataSource3D)::DownCast(
                                    myDS2->ArrayOfMeshDS.value(bodyIndex));
                        occHandle(Ng_MeshVS_DataSource3D) outputMesh;
                        MeshTools::filterVolumeElementsByType(mesh3D,HEXA,outputMesh);
                        if(!outputMesh.IsNull())
                        {
                            cout<<"____number of prisms: "<<outputMesh->GetAllElements().Extent()<<"____"<<endl;
                            cout<<"____number of nodes: "<<outputMesh->GetAllNodes().Extent()<<"____"<<endl;
                            occPreGLWidget::displayMesh(occMeshContext,Quantity_NOC_GREEN,outputMesh);
                        }
                    }
                }
            }
        }
    }
}

//! --------------------------------------
//! function: buildSuppressionContextMenu
//! details:
//! --------------------------------------
void occPreGLWidget::buildSuppressionContextMenu()
{   
    //! ------------------------------------------
    //! check if some bodies have been suppressed
    //! ------------------------------------------
    QWidget *smw = tools::getWidgetByName("simmanager");
    SimulationManager *sm = static_cast<SimulationManager*>(smw);
    QStandardItem *itemGeometryRoot = sm->getTreeItem(SimulationNodeClass::nodeType_geometry);
    int NbBodies = itemGeometryRoot->rowCount();
    bool somethingIsSuppressed = false;
    for(int i=0; i<NbBodies; i++)
    {
        QExtendedStandardItem* itemBody = static_cast<QExtendedStandardItem*>(itemGeometryRoot->child(i,0));
        SimulationNodeClass* curNode = itemBody->data(Qt::UserRole).value<SimulationNodeClass*>();
        Property::SuppressionStatus ss = curNode->getPropertyItem("Suppressed")->data(Qt::UserRole).value<Property>().getData().value<Property::SuppressionStatus>();

        if(ss==Property::SuppressionStatus_Suppressed)
        {
            somethingIsSuppressed = true;
            break;
        }
    }

    //! -------------------------------------
    //! check if something has been selected
    //! -------------------------------------
    occContext->InitSelected();
    if(occContext->MoreSelected())
    {
        contextMenuBuilder::addActionSuppress(myContextMenu);
        contextMenuBuilder::addActionSuppressAllOther(myContextMenu);
        if(somethingIsSuppressed) contextMenuBuilder::addActionUnsuppressAllBodies(myContextMenu);
    }
    myContextMenu->addSeparator();
}

//! ------------------------------
//! function: clearMeshFromViewer
//! details:
//! ------------------------------
void occPreGLWidget::clearMeshFromViewer()
{
    cerr<<"occPreGLWidget::clearMeshFromViewer()->____function called____"<<endl;
    //! if something has been selected delete the mesh on the selected bodies
    occContext->InitSelected();
    if(occContext->MoreSelected())
    {
        //! ----------------------------------------------------
        //! remove the mesh data source from the mesh data base
        //! ----------------------------------------------------
        for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
        {
            //! --------------------------
            //! the current selected body
            //! --------------------------
            const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());
            int bodyIndex = myDS2->bodyMap.key(theAISShape->Shape());

            this->invalidateMesh(bodyIndex);

            /*
            //! ----------------------------
            //! clear the mesh data sources
            //! ----------------------------
            myDS2->ArrayOfMeshDS.insert(bodyIndex, occHandle(MeshVS_DataSource)());
            myDS2->ArrayOfMeshDS2D.insert(bodyIndex, occHandle(MeshVS_DataSource)());

            //! ---------------------------------
            //! delete the face mesh datasources
            //! ---------------------------------
            for(int faceNr=1; faceNr<=myDS2->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent(); faceNr++)
                myDS2->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr, occHandle(MeshVS_DataSource)());

            //! ---------------------------------
            //! delete the edge mesh datasources
            //! ---------------------------------
            for(int edgeNr=1; edgeNr<=myDS2->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent(); edgeNr++)
                myDS2->ArrayOfMeshDSOnEdges.setValue(bodyIndex,edgeNr, occHandle(MeshVS_DataSource)());
            */

            //! ---------------------------
            //! set the flag to be updated
            //! ---------------------------
            myDS2->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,true);
        }

        //! ------------------------------------------
        //! remove from the viewer and delete the IOs
        //! ------------------------------------------
        for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
        {
            const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());
            int bodyIndex = theAISShape->index();

            //! clear the mesh MeshVS_Mesh from the viewer
            const occHandle(MeshVS_Mesh) &aMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.value(bodyIndex));
            //const occHandle(MeshVS_Mesh) &aMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.at(bodyIndex));
            occMeshContext->Remove(aMesh,Standard_False);

            //myMapOfInteractiveMeshes[bodyIndex] = occHandle(MeshVS_Mesh)();
            myMapOfInteractiveMeshes.insert(bodyIndex,occHandle(MeshVS_Mesh)());
        }
        emptyTheSelection();
    }
    else //! remove all the meshes on all the bodies
    {
        //! ----------------------------------------
        //! remove all the meshes from the database
        //! ----------------------------------------
        //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it= myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
        for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it= myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
        {
            //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
            //const occHandle(AIS_ExtendedShape) &theShape = occHandle(AIS_ExtendedShape)::DownCast(apair.second);
            const occHandle(AIS_ExtendedShape) &theShape = occHandle(AIS_ExtendedShape)::DownCast(it.value());
            int bodyIndex = theShape->index();

            this->invalidateMesh(bodyIndex);

            /*
            //! clear the mesh data sources
            myDS2->ArrayOfMeshDS.insert(bodyIndex, occHandle(MeshVS_DataSource)());
            myDS2->ArrayOfMeshDS2D.insert(bodyIndex, occHandle(MeshVS_DataSource)());

            for(int faceNr=1; faceNr<=myDS2->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent(); faceNr++)
            {
                myDS2->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr, occHandle(MeshVS_DataSource)());
            }
            for(int edgeNr=1; edgeNr<=myDS2->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent(); edgeNr++)
            {
                myDS2->ArrayOfMeshDSOnEdges.setValue(bodyIndex,edgeNr, occHandle(MeshVS_DataSource)());
            }
            */

            //! set the flag to be updated
            myDS2->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex,true);
        }

        //! -------------------------------------------------------------------------------
        //! remove the meshes from the viewer
        //! since the external mesh is incorporated into the view, the auto-triangulation
        //! is reactivated, and the presentation of the AIS_Shape is recomputed, using the
        //! auto-generated BRepMesh
        //! -------------------------------------------------------------------------------
        occContext->DefaultDrawer()->SetAutoTriangulation(Standard_True);

        for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it= myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
        //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it= myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
        {
            //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
            //int k = apair.first;
            //const occHandle(MeshVS_Mesh) &aMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.at(k));
            int k = it.key();
            const occHandle(MeshVS_Mesh) &aMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.value(k));
            occMeshContext->Remove(aMesh,Standard_False);

            //const occHandle(AIS_ExtendedShape) &theShape = occHandle(AIS_ExtendedShape)::DownCast(apair.second);
            const occHandle(AIS_ExtendedShape) &theShape = occHandle(AIS_ExtendedShape)::DownCast(it.value());

            int bodyIndex = theShape->index();

            //myMapOfInteractiveMeshes[bodyIndex] = occHandle(MeshVS_Mesh)();
            myMapOfInteractiveMeshes.insert(bodyIndex,occHandle(MeshVS_Mesh)());
            occContext->RecomputePrsOnly(theShape,Standard_False);
        }
    }

    //! update
    occMeshContext->UpdateCurrentViewer();
    occContext->UpdateCurrentViewer();

    //! diagnostic
    printSummary();
}

//! -------------------------
//! function: invalidateMesh
//! details:
//! -------------------------
void occPreGLWidget::invalidateMesh(int bodyIndex, Standard_Boolean updateViewer)
{
    //! ------------------------------------
    //! display the obsolete mesh in yellow
    //! ------------------------------------
    const occHandle(MeshVS_Mesh) &theCurIOMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.value(bodyIndex));
    //const occHandle(MeshVS_Mesh) &theCurIOMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.at(bodyIndex));
    if(!theCurIOMesh.IsNull())
    {
        theCurIOMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_RED);
        occMeshContext->RecomputePrsOnly(theCurIOMesh,Standard_False,Standard_False);
    }

    //! ---------------------------------------
    //! reset the mesh surface and volume mesh
    //! ---------------------------------------
    myDS2->ArrayOfMeshDS2D.insert(bodyIndex,occHandle(MeshVS_DataSource)());
    myDS2->ArrayOfMeshDS.insert(bodyIndex,occHandle(MeshVS_DataSource)());

    myDS2->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, Standard_True);
    int Nfaces = myDS2->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
    for(int faceNr=1; faceNr<= Nfaces; faceNr++)
    {
        myDS2->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(MeshVS_DataSource)());
    }

    int Nedges = myDS2->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent();
    for(int edgeNr=1; edgeNr<=Nedges; edgeNr++)
    {
        myDS2->ArrayOfMeshDSOnEdges.setValue(bodyIndex,edgeNr,occHandle(MeshVS_DataSource)());
    }
    if(updateViewer) occContext->UpdateCurrentViewer();
}

//! -------------------------
//! function: invalidateMesh
//! details:  overload
//! -------------------------
void occPreGLWidget::invalidateMeshes(const std::vector<int> &indexes)
{
    for(std::vector<int>::const_iterator anIter = indexes.cbegin(); anIter!=indexes.cend(); anIter++)
    {
        int bodyIndex = *anIter;
        const occHandle(MeshVS_Mesh) &theCurIOMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.value(bodyIndex));
        //const occHandle(MeshVS_Mesh) &theCurIOMesh = occHandle(MeshVS_Mesh)::DownCast(myMapOfInteractiveMeshes.at(bodyIndex));
        if(theCurIOMesh.IsNull()) continue;

        //! ------------------------------------
        //! display the obsolete mesh in yellow
        //! ------------------------------------
        theCurIOMesh->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_YELLOW);
        occMeshContext->RecomputePrsOnly(theCurIOMesh,Standard_False,Standard_False);

        //! ------------------------------
        //! reset surface and volume mesh
        //! ------------------------------
        myDS2->ArrayOfMeshDS.insert(bodyIndex,occHandle(Ng_MeshVS_DataSource3D)());
        myDS2->ArrayOfMeshDS2D.insert(bodyIndex,occHandle(Ng_MeshVS_DataSource2D)());

        //! ----------------------------------------------
        //! set the mesh of the body as "must be updated"
        //! ----------------------------------------------
        myDS2->ArrayOfMeshIsToBeUdpdated.insert(bodyIndex, Standard_True);

        //! --------------------------------------------------------------------------
        //! MeshVS_DataSource for the faces and edges: null MeshVS_DataSource handle
        //! -------------------------------------------------------------------------
        int Nfaces = myDS2->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
        for(int faceNr=0; faceNr<= Nfaces; faceNr++)
        {
            myDS2->ArrayOfMeshDSOnFaces.setValue(bodyIndex,faceNr,occHandle(Ng_MeshVS_DataSource2D)());
        }
        int Nedges = myDS2->MapOfBodyTopologyMap.value(bodyIndex).edgeMap.Extent();
        for(int edgeNr=1; edgeNr<=Nedges; edgeNr++)
        {
            myDS2->ArrayOfMeshDSOnEdges.setValue(bodyIndex,edgeNr,occHandle(Ng_MeshVS_DataSource1D)());
        }
    }
    occContext->UpdateCurrentViewer();
}

//! ------------------------------
//! function: invalidateAllMeshes
//! details:
//! ------------------------------
void occPreGLWidget::invalidateAllMeshes()
{
    for(int k=1;k<=myDS2->bodyMap.size();k++) this->invalidateMesh(k, Standard_False);
    occMeshContext->UpdateCurrentViewer();
}

//! ------------------------------------
//! function: printSummary
//! details:  summary after CAD loading
//! ------------------------------------
void occPreGLWidget::printSummary()
{
    cout<<endl<<endl;
    cout<<"------------------------------------------------------------";
    cout<<"------------------------------------------------------------------------------"<<endl;
    cout<<"BI\t"<<"BIFDS\t"<<"ShapeType\t"<<"DS\t"<<"IO\t"
       <<"Type\t"<<"2DEng\t"<<"Type\t"<<"3DEng\t"<<"Order\t"<<"NN2D"<<"\t"<<"NE2D"<<"\t"<<"NN3D"<<"\t"<<"NE3D"<<endl;
    cout<<endl;

    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!= myMapOfInteractiveShapes.end(); ++it)
    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!= myMapOfInteractiveShapes.end(); ++it)
    {
        //! interactive AIS shape
        const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(it.value());
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //const occHandle(AIS_ExtendedShape) &theAISShape = occHandle(AIS_ExtendedShape)::DownCast(apair.second);

        //! topology and geometry content
        const TopoDS_Shape &theShape = theAISShape->Shape();

        //! body index
        int bodyIndex = theAISShape->index();

        //! body index from the geometry data structure - must be the same of the previous
        int bodyIndex2 = myDS2->bodyMap.key(theShape);

        //! shape type
        Standard_CString shapeType;
        if(theShape.ShapeType()==TopAbs_SOLID)shapeType="3D body";

        //! mesh DS is valid?
        Standard_CString meshDSStatus;
        if(!myDS2->ArrayOfMeshDS.value(bodyIndex).IsNull())meshDSStatus="VALID";
        else meshDSStatus="NULL";

        //! mesh VS is valid?
        Standard_CString meshVSStatus;
        if(!myMapOfInteractiveMeshes.value(bodyIndex).IsNull()) meshVSStatus="VALID";
        //if(!myMapOfInteractiveMeshes.at(bodyIndex).IsNull()) meshVSStatus="VALID";
        else meshVSStatus="EMPTY";

        //! 2D & 3D mesh engine
        Standard_CString me2D = (myDS2->ArrayOfMesh2DEngine.value(bodyIndex) == Property::meshEngine2D_Netgen)? "Ng2D" : "OCC";
        Standard_CString me3D = (myDS2->ArrayOfMesh3DEngine.value(bodyIndex) == Property::meshEngine3D_Netgen)? "Ng3D" : "N.A.";

        //! mesh order: first/second
        Standard_CString mesh_Order = (myDS2->ArrayOfMeshOrder.value(bodyIndex) == Property::meshOrder_First)? "1-st" : "2-nd";

        //! surface mesh type
        //Standard_CString smt = (myDS2->ArrayOfSurfaceMeshType.value(bodyIndex) == Property::meshType_Surface_AllTrig)? "trig" : "quad";

        //! volume mesh type
        //Standard_CString vmt = (myDS2->ArrayOfVolumeMeshType.value(bodyIndex) == Property::meshType_Volume_AllTet)? "tet" : "hexa";

        //! number of surface nodes and elements
        int NP2D = 0;
        int NE2D = 0;
        int NP3D = 0;
        int NE3D = 0;

        //! number of nodes and elements
        if(!myDS2->ArrayOfMeshDS.value(bodyIndex).IsNull())
        {
            NP3D = NP3D + myDS2->ArrayOfMeshDS.value(bodyIndex)->GetAllNodes().Extent();
            NE3D = NE3D + myDS2->ArrayOfMeshDS.value(bodyIndex)->GetAllElements().Extent();
            myDS2->ArrayOfMeshDS.value(bodyIndex)->IsKind("Ng_MeshVS_DataSource2D");
        }
        //! dump to the standard output
        //cout<<bodyIndex<<"\t"<<bodyIndex2<<"\t"<<shapeType<<"\t\t"<<meshDSStatus<<"\t"<<meshVSStatus
        //   <<"\t"<<smt<<"\t"<<me2D<<"\t"<<vmt<<"\t"<<me3D<<"\t"<<mesh_Order<<"\t"<<NP2D<<"\t"<<NE2D<<"\t"<<NP3D<<"\t"<<NE3D<<endl;
        cout<<bodyIndex<<"\t"<<bodyIndex2<<"\t"<<shapeType<<"\t\t"<<meshDSStatus<<"\t"<<meshVSStatus
           <<"\t"<<me2D<<"\t"<<me3D<<"\t"<<mesh_Order<<"\t"<<NP2D<<"\t"<<NE2D<<"\t"<<NP3D<<"\t"<<NE3D<<endl;
    }
}

//! ---------------------------------
//! function: set the selection mode
//! details:
//! ---------------------------------
#include <Aspect_TypeOfHighlightMethod.hxx>
//void occPreGLWidget::setSelectionStyle(Quantity_NameOfColor color, float transparency)
//{
//    occHandle(Prs3d_Drawer) hs = new Prs3d_Drawer();
//    hs->SetColor(color);
//    hs->SetTransparency(transparency);
//    occContext->SetSelectionStyle(hs);
//}

//! -----------------------------------------------------------
//! function: removeObsoleteMeshes
//! details:  slot: remove the obsolete meshes from the viewer
//! -----------------------------------------------------------
void occPreGLWidget::removeObsoleteMeshes()
{
    //! ------------------------------------------
    //! remove from the display the invalid grids
    //! ------------------------------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    {
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        const occHandle(AIS_ExtendedShape) &curAIS = occHandle(AIS_ExtendedShape)::DownCast(it.value());
        //const occHandle(AIS_ExtendedShape) &curAIS = occHandle(AIS_ExtendedShape)::DownCast(apair.second);
        const TopoDS_Shape &theShape = curAIS->Shape();
        int bodyIndex = myDS2->bodyMap.key(theShape);

        Standard_Boolean isMeshToBeUpdated = myDS2->ArrayOfMeshIsToBeUdpdated.value(bodyIndex);
        if(isMeshToBeUpdated==true)
        {
            occMeshContext->Remove(myMapOfInteractiveMeshes.value(bodyIndex),false);
            myMapOfInteractiveMeshes.insert(bodyIndex,occHandle(MeshVS_Mesh)());

            //occMeshContext->Remove(myMapOfInteractiveMeshes.at(bodyIndex),false);
            //myMapOfInteractiveMeshes[bodyIndex] = occHandle(MeshVS_Mesh)();
        }
    }
}

//! ------------------------------------
//! function: suppress all other bodies
//! details:
//! ------------------------------------
void occPreGLWidget::suppressAllOtherBodies()
{
    std::vector<int> indexOfSelected;
    for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
    {
        const occHandle(AIS_ExtendedShape) &theAIS = occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());
        const TopoDS_Shape &theShape = theAIS->Shape();
        int bodyIndex = myDS2->bodyMap.key(theShape);
        indexOfSelected.push_back(bodyIndex);
    }
    for(int i=1; i<=myDS2->bodyMap.size(); i++)
    {
        int n=0;
        for(int k=0; k<indexOfSelected.size(); k++)
        {
            if(i!=indexOfSelected[k])n++;
        }
        if(n==indexOfSelected.size() && myDS2->MapOfIsActive.value(i)==false)
        {
            myDS2->MapOfIsActive.insert(i,false);
            occContext->Remove(myMapOfInteractiveShapes.value(i),Standard_False);
            //occContext->Remove(myMapOfInteractiveShapes.at(i),Standard_False);
        }
    }
}

//! ------------------------------
//! function: clear custom colors
//! details:
//! -------------------------------
void occPreGLWidget::clearAllCustomColors()
{
    cout<<"occPreGLWidget::clearAllCustomColors()->____function called____"<<endl;

    //! clean all the custom colors
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    {
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //const occHandle(AIS_ColoredShape) &AISShape = occHandle(AIS_ColoredShape)::DownCast(apair.second);
        const occHandle(AIS_ColoredShape) &AISShape = occHandle(AIS_ColoredShape)::DownCast(it.value());
        AISShape->ClearCustomAspects();
    }
    //! recompute the presentation
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    {
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //const occHandle(AIS_ColoredShape) &AISShape = occHandle(AIS_ColoredShape)::DownCast(apair.second);
        const occHandle(AIS_ColoredShape) &AISShape = occHandle(AIS_ColoredShape)::DownCast(it.value());
        occContext->RecomputePrsOnly(AISShape,Standard_False);
    }
    //! update the viewer
    occContext->UpdateCurrentViewer();
}

//! ----------------------------------
//! function: displayColoredSubshapes
//! details:
//! ----------------------------------
void occPreGLWidget::displayColoredSubshapes(const TopTools_ListOfShape &listOfShape, Quantity_NameOfColor aColor, bool cleanPrevious, bool updateViewer)
{
    cout<<"occPreGLWidget::displayColoredSubshapes()->____function called____"<<endl;

    this->resetCustomColors();
    if(listOfShape.Extent()!=0)
    {
        //! -----------------
        //! handle the scope
        //! -----------------
        if(listOfShape.First().ShapeType()==TopAbs_SOLID)
        {
            TopTools_ListIteratorOfListOfShape anIter;
            for(anIter.Initialize(listOfShape);anIter.More();anIter.Next())
            {
                int bodyIndex = myDS2->bodyMap.key(anIter.Value());
                //const occHandle(AIS_ColoredShape) &AISShape =  occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.at(bodyIndex));
                const occHandle(AIS_ColoredShape) &AISShape =  occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
                AISShape->SetCustomColor(anIter.Value(),aColor);
            }
        }
        else
        {
            TopTools_ListIteratorOfListOfShape anIter;
            int bodyIndex, subShapeIndex;
            TopAbs_ShapeEnum type;
            for(anIter.Initialize(listOfShape);anIter.More();anIter.Next())
            {
                myDS2->getSubShapeNr(anIter.Value(),bodyIndex,subShapeIndex,type);
                //const occHandle(AIS_ColoredShape) &AISShape =  occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.at(bodyIndex));
                const occHandle(AIS_ColoredShape) &AISShape =  occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
                AISShape->SetCustomColor(anIter.Value(),aColor);
            }
        }
    }
    else
    {
        //! ----------------------------
        //! clean all the custom colors
        //! ----------------------------
        //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
        for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it=myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
        {
            //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
            //const occHandle(AIS_ColoredShape) &AISShape =  occHandle(AIS_ColoredShape)::DownCast(apair.second);
            const occHandle(AIS_ColoredShape) &AISShape =  occHandle(AIS_ColoredShape)::DownCast(it.value());
            AISShape->ClearCustomAspects();
        }
    }
    //! ---------------------------
    //! recompute the presentation
    //! ---------------------------
    //for(std::map<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); ++it)
    {
        //std::pair<int,occHandle(AIS_InteractiveObject)> apair = *it;
        //const occHandle(AIS_ColoredShape) &AISShape = occHandle(AIS_ColoredShape)::DownCast(apair.second);
        const occHandle(AIS_ColoredShape) &AISShape = occHandle(AIS_ColoredShape)::DownCast(it.value());
        occContext->RecomputePrsOnly(AISShape,Standard_False);
    }
    //! ------------------
    //! update the viewer
    //! ------------------
    if(updateViewer) occContext->UpdateCurrentViewer();
}

//! --------------------------------
//! function: hideBody
//! details:  hide a list of bodies
//! --------------------------------
void occPreGLWidget::hideBody(const TColStd_ListOfInteger &listOfBodyNumbers)
{
    cout<<"occPreGLWidget::hideBody()->____function called____"<<endl;
    cout<<"occPreGLWidget::hideBody()->____number of bodies to hide: "<<listOfBodyNumbers.Extent()<<"____"<<endl;
    occContext->CloseLocalContext(occContext->IndexOfCurrentLocal());
    TColStd_ListIteratorOfListOfInteger anIter;
    for(anIter.Initialize(listOfBodyNumbers); anIter.More(); anIter.Next())
    {
        int bodyIndex = anIter.Value();
        const occHandle(AIS_ExtendedShape) &IO = occHandle(AIS_ExtendedShape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
        //const occHandle(AIS_ExtendedShape) &IO = occHandle(AIS_ExtendedShape)::DownCast(myMapOfInteractiveShapes.at(bodyIndex));
        IO->setShapeVisibility(Standard_False);
        //occContext->Remove(IO,false);
        occContext->Erase(IO, false);
    }
    myLocalCtxNumber=occContext->OpenLocalContext();

    //! -------------------------------------------------------------------------
    //! now the selection modes must be reactivated, because when the
    //! context is closed, the selection modes (and the selection list) are lost
    //! -------------------------------------------------------------------------
    this->reactivateCurrentStandardSelectionMode();
    occContext->UpdateCurrentViewer();
}

//! ----------------------------------------------
//! function - slot: show body
//! details: show a list of bodies only if active
//! ----------------------------------------------
void occPreGLWidget::showBody(const TColStd_ListOfInteger &listOfBodies)
{
    cout<<"occPreGLWidget::showBody()->____function called____"<<endl;
    occContext->CloseLocalContext(occContext->IndexOfCurrentLocal());
    TColStd_ListIteratorOfListOfInteger anIter;
    for(anIter.Initialize(listOfBodies); anIter.More(); anIter.Next())
    {
        int bodyIndex = anIter.Value();
        const occHandle(AIS_ExtendedShape) &IO = occHandle(AIS_ExtendedShape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
        //const occHandle(AIS_ExtendedShape) &IO = occHandle(AIS_ExtendedShape)::DownCast(myMapOfInteractiveShapes.at(bodyIndex));
        occContext->Display(IO,false);
        IO->setShapeVisibility(Standard_True);
    }
    myLocalCtxNumber=occContext->OpenLocalContext();

    //! -------------------------------------------------------------------------
    //! now the selection modes must be reactivated, because when the
    //! context is closed, the selection modes (and the selection list) are lost
    //! -------------------------------------------------------------------------
    this->reactivateCurrentStandardSelectionMode();
    occContext->UpdateCurrentViewer();
}

//! -------------------------------------
//! function: highlight body
//! details:  highlight a list of bodies
//! -------------------------------------
void occPreGLWidget::highlightBody(const QList<int> &listOfBodyNumbers)
{
    cout<<"occPreGLWidget::highlightBody()->____function called____"<<endl;
    this->unhighlightBody(false);

    //! -----------------------------
    //! highlight the list of bodies
    //! -----------------------------
    occHandle(Prs3d_Drawer) highlightDrawer = new Prs3d_Drawer();
    highlightDrawer->SetMethod(Aspect_TOHM_COLOR);
    highlightDrawer->SetColor(static_cast<Quantity_Color>(Quantity_NOC_GREEN));
    highlightDrawer->SetTransparency(0.0);

    occContext->SetHighlightStyle(Prs3d_TypeOfHighlight_LocalSelected,highlightDrawer);
    for(QList<int>::const_iterator it = listOfBodyNumbers.cbegin(); it!=listOfBodyNumbers.cend(); ++it)
    {
        int bodyIndex = *it;
        //const occHandle(AIS_Shape) &anAISShape = occHandle(AIS_Shape)::DownCast(myMapOfInteractiveShapes.at(bodyIndex));
        const occHandle(AIS_Shape) &anAISShape = occHandle(AIS_Shape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
        occContext->HilightWithColor(anAISShape,highlightDrawer,false);
    }
    occContext->UpdateCurrentViewer();
}

//! ------------------------------------------------------
//! function: unhighlightBody
//! details:  slot - remove the highlight from the bodies
//! ------------------------------------------------------
void occPreGLWidget::unhighlightBody(bool updateViewer)
{
    //! remove the highlight from all the bodies
    AIS_ListOfInteractive AISList;
    AIS_ListIteratorOfListOfInteractive iter;
    occContext->ObjectsInside(AISList,AIS_KOI_Shape,0);
    for(iter.Initialize(AISList);iter.More();iter.Next())
    {
        occContext->Unhilight(iter.Value(),false);
    }
    if(updateViewer==true) occContext->UpdateCurrentViewer();
}

//! -----------------------------------------------------------------
//! function: displayShapeCopy
//! details:  aim: showing a boundary condition on a subtopology
//!           using different colors. Showing a mesh control on an
//!           edge using solids segments. Showing an edge control on
//!           a vertex using a pinball (region of influence)
//! -----------------------------------------------------------------
void occPreGLWidget::displayShapeCopy(const TopTools_ListOfShape &list1,
                                      const TopTools_ListOfShape &list2,
                                      Quantity_NameOfColor color1,
                                      Quantity_NameOfColor color2,
                                      QVariant options)
{
    //cout<<"occPreGLWidget::displayShapeCopy()->____function called____"<<endl;

    const double tol = 1e-2;

    //! ------------------------------
    //! access the simulation manager
    //! ------------------------------
    SimulationManager *sm = static_cast<SimulationManager*>(tools::getWidgetByName("simmanager"));
    SimulationNodeClass *theCurNode = sm->getCurrentNode();

    TopoDS_Compound myCompound1, myCompound2;
    static occHandle(AIS_Shape) anOldShape1,anOldShape2;

    //! ------------------------------------------
    //! do not update the viewer (done at the end)
    //! ------------------------------------------
    if(!anOldShape1.IsNull()) occContext->Remove(anOldShape1,false);
    if(!anOldShape2.IsNull()) occContext->Remove(anOldShape2,false);

    if(!list1.IsEmpty())
    {
        //! --------------------------------------------
        //! get the type of shape contained in the list
        //! --------------------------------------------
        TopAbs_ShapeEnum typeList1 = list1.First().ShapeType();
        TopoDS_Builder myBuilder1;
        switch(typeList1)
        {
        case TopAbs_VERTEX:
        {
            double radius = options.toDouble();
            myBuilder1.MakeCompound(myCompound1);
            TopTools_ListIteratorOfListOfShape anIter;
            for(anIter.Initialize(list1);anIter.More();anIter.Next())
            {
                const TopoDS_Vertex &theVertex = TopoDS::Vertex(anIter.Value());
                const TopoDS_Shape &s = BRepPrimAPI_MakeSphere::BRepPrimAPI_MakeSphere(BRep_Tool::Pnt(theVertex),radius);
                myBuilder1.Add(myCompound1,s);
            }
            anOldShape1 = new AIS_ColoredShape(myCompound1);
            //anOldShape1->SetZLayer(Graphic3d_ZLayerId_TopOSD);

            //! -----------------------
            //! needed for shaded view
            //! -----------------------
            if(!occContext->DefaultDrawer()->IsAutoTriangulation()) occContext->DefaultDrawer()->SetAutoTriangulation(Standard_True);
            anOldShape1->SetTransparency(0.7);

            //! ----------------------
            //! hide sphere wireframe
            //! ----------------------
            anOldShape1->Attributes()->SetFaceBoundaryDraw(false);
            occContext->SetColor(anOldShape1,Quantity_Color(Quantity_NOC_RED),false);
            occContext->Display(anOldShape1,1,-1,false);
        }
        break;

        case TopAbs_EDGE:
        {
            cout<<"____handling edge____"<<endl;
            myBuilder1.MakeCompound(myCompound1);

            //! -----------------------------------------------------------------------
            //! in this case the QVariant "options" contains "true" if the calling
            //! functions deals with a mesh control item: in this case a certain
            //! number of cylinders, showing the mesh esge size is shown; otherwise
            //! the edge is copied and shown with the color of the boundary condititon
            //! -----------------------------------------------------------------------
            if(options.toBool())
            {
                cout<<"____handling edge: a mesh control___"<<endl;
                TopTools_ListIteratorOfListOfShape itEdges;
                for(itEdges.Initialize(list1);itEdges.More();itEdges.Next())
                {
                    cout<<"____accessing the edge BRep____"<<endl;
                    //! -------------------------------------------
                    //! access the BRep representation of the edge
                    //! -------------------------------------------
                    const TopoDS_Edge &edge = TopoDS::Edge(itEdges.Value());
                    if(BRep_Tool::Degenerated(edge)) continue;
                    BRepAdaptor_Curve BRepAdaptor(edge);
                    GeomAdaptor_Curve curve = BRepAdaptor.Curve();
                    //! -------------------------------------------
                    //! calculate the actual length
                    //! -------------------------------------------
                    CPnts_AbscissaPoint CP;
                    CP.Init(curve);

                    double L = CP.Length(curve,tol);
                    double newDeltaL;
                    int NbDivisions = 1;
                    //cout<<"____curve length: " <<L<<"____"<<endl;
                    //! ----------------------------------------------
                    //! retrieve the element size/number of divisions
                    //! from the GUI
                    //! ----------------------------------------------
                    double elementSize_GUI;
                    int typeOfSizing = theCurNode->getPropertyValue<int>("Sizing type");
                    switch(typeOfSizing)
                    {
                    case 0:
                    {
                        //cout<<"____handling case element size____"<<endl;
                        //! ------------------------------
                        //! the element size is specified
                        //! ------------------------------
                        QExtendedStandardItem *item = theCurNode->getPropertyItem("Element size");
                        elementSize_GUI = item->data(Qt::UserRole).value<Property>().getData().toDouble();
                        NbDivisions = int(L/elementSize_GUI);
                        if(NbDivisions ==0) NbDivisions++;
                        newDeltaL = L/NbDivisions;
                    }
                        break;

                    case 1:
                    {
                        //cout<<"____handling case number of divisions____"<<endl;
                        //! -------------------------------------
                        //! the number of divisions is specified
                        //! -------------------------------------
                        QExtendedStandardItem *item = theCurNode->getPropertyItem("Number of divisions");
                        int NbDivisions_GUI = item->data(Qt::UserRole).value<Property>().getData().toInt();
                        //cout<<"____number of divisions: "<<NbDivisions_GUI<<"____"<<endl;
                        if(NbDivisions_GUI==0) NbDivisions_GUI++;
                        newDeltaL = L/NbDivisions_GUI;
                        //cout<<"____new deltaL: "<<newDeltaL<<"____"<<endl;
                        NbDivisions = NbDivisions_GUI;
                    }
                        break;
                    }
                    //! --------------------
                    //! build the cylinders
                    //! --------------------
                    double s_old = BRepAdaptor.FirstParameter();
                    for(int i=0; i<NbDivisions; i++)
                    {
                        CP.Perform(newDeltaL,s_old,tol);
                        double s = CP.Parameter();
                        gp_Pnt P0 = BRepAdaptor.Value(s_old);
                        gp_Pnt P1 = BRepAdaptor.Value(s);

                        //! axis of the cylinder
                        gp_Vec vec(P0,P1);
                        gp_Dir dir(vec);
                        gp_Ax2 ax(P0,dir);
                        TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(ax,1,newDeltaL-0.3*newDeltaL).Shape();
                        myBuilder1.Add(myCompound1,aCyl);
                        s_old = s;
                        //cout<<i<<"____"<<"(x,y,z) = "<<"("<<P0.X()<<", "<<P0.Y()<<", "<<P0.Z()<<")"<<endl;
                    }
                }

                //! ----------------------------------
                //! display the sequence of cylinders
                //! use model z-layer
                //! ----------------------------------
                anOldShape1 = new AIS_MeshSegmentMarker(myCompound1);
                anOldShape1->SetColor(Quantity_Color(Quantity_NOC_YELLOW));

                if(!occContext->DefaultDrawer()->IsAutoTriangulation()) occContext->DefaultDrawer()->SetAutoTriangulation(Standard_True);
                anOldShape1->SetTransparency(0.0);

                //! ----------------------------------------
                //! hide wirefreame (and also the seam edge)
                //! ----------------------------------------
                anOldShape1->Attributes()->SetFaceBoundaryDraw(false);
            }
            else
            {
                //! ----------------------------------------------
                //! build a thick cylinder to overlap to the edge
                //! not working... to do
                //! ----------------------------------------------
                cout<<"____handling edge: not a mesh control___"<<endl;
                TopTools_ListIteratorOfListOfShape itEdges;
                for(itEdges.Initialize(list1); itEdges.More(); itEdges.Next())
                {
                    myBuilder1.Add(myCompound1,itEdges.Value());
                }
            }

            occContext->Display(anOldShape1,1,-1,false);
        }
            break;

            //! ---------------------------------------------------------------------------------------
            //! The default case is for TopAbs_SOLID (and COMPSOLID, and COMPOUND) and for TopAbs_FACE
            //! ---------------------------------------------------------------------------------------
        default:
        {
            TopTools_ListIteratorOfListOfShape anIter;
            myBuilder1.MakeCompound(myCompound1);
            for(anIter.Initialize(list1);anIter.More();anIter.Next())
            {
                myBuilder1.Add(myCompound1,anIter.Value());
            }
            anOldShape1 = new AIS_Shape(myCompound1);

            //! ----------------------
            //! use the model z-layer
            //! ----------------------
            //anOldShape1->SetZLayer(Graphic3d_ZLayerId_TopOSD);
            anOldShape1->SetColor(color1);
            anOldShape1->Attributes()->ShadingAspect()->SetTransparency(0.0,Aspect_TOFM_FRONT_SIDE);
            anOldShape1->Attributes()->ShadingAspect()->SetTransparency(1.0,Aspect_TOFM_BACK_SIDE);
            //! ----------------------------------------------------------------------------
            //! display: second argument: Shaded mode: third argument: no selection allowed
            //! does not update the viewer (done at the end)
            //! ----------------------------------------------------------------------------
            occContext->Display(anOldShape1,1,-1,false);
        }
            break;
        }
        //occContext->UpdateCurrentViewer();
    }

    if(!list2.IsEmpty())
    {
        //! --------------------------------------------
        //! get the type of shape contained in the list
        //! --------------------------------------------
        TopAbs_ShapeEnum typeList2 = list2.First().ShapeType();
        if(typeList2!=TopAbs_VERTEX)
        {
            TopoDS_Builder myBuilder2;
            TopTools_ListIteratorOfListOfShape anIter;
            myBuilder2.MakeCompound(myCompound2);
            for(anIter.Initialize(list2);anIter.More();anIter.Next())
            {
                myBuilder2.Add(myCompound2,anIter.Value());
            }
            anOldShape2 = new AIS_ColoredShape(myCompound2);
            anOldShape2->SetColor(color2);
            anOldShape2->Attributes()->ShadingAspect()->SetTransparency(0.0,Aspect_TOFM_FRONT_SIDE);
            anOldShape2->Attributes()->ShadingAspect()->SetTransparency(1.0,Aspect_TOFM_BACK_SIDE);
            //! ----------------------------------------------------------------------------
            //! display: second argument: Shaded mode: third argument: no selection allowed
            //! does not update the viewer (done at the end)
            //! ----------------------------------------------------------------------------
            occContext->Display(anOldShape2,1,-1,false);
        }
        else
        {
            ;
        }
    }
    occContext->UpdateCurrentViewer();
}


void occPreGLWidget::displayShapeCopy1(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color)
{
    cout<<"occPreGLWidget::displayShapeCopy1()->____function called____"<<endl;

    TopoDS_Compound myCompound1;
    static occHandle(AIS_Shape) anOldShape1;

    //! ------------------------------------------
    //! do not update the viewer (done at the end)
    //! ------------------------------------------
    if(!anOldShape1.IsNull()) occContext->Remove(anOldShape1,false);

    if(!listShapes.IsEmpty())
    {
        //! --------------------------------------------
        //! get the type of shape contained in the list
        //! --------------------------------------------
        TopAbs_ShapeEnum typeList1 = listShapes.First().ShapeType();
        TopoDS_Builder myBuilder1;
        switch(typeList1)
        {
        case TopAbs_VERTEX:
        {
            ;
        }
        break;

        default:
        {
            TopTools_ListIteratorOfListOfShape anIter;
            myBuilder1.MakeCompound(myCompound1);
            for(anIter.Initialize(listShapes);anIter.More();anIter.Next())
            {
                myBuilder1.Add(myCompound1,anIter.Value());
            }
            anOldShape1 = new AIS_Shape(myCompound1);

            //! ----------------------
            //! use the model z-layer
            //! ----------------------
            //anOldShape1->SetZLayer(Graphic3d_ZLayerId_TopOSD);
            anOldShape1->SetColor(color);
            anOldShape1->Attributes()->ShadingAspect()->SetTransparency(0.0,Aspect_TOFM_FRONT_SIDE);
            anOldShape1->Attributes()->ShadingAspect()->SetTransparency(1.0,Aspect_TOFM_BACK_SIDE);
            //! ----------------------------------------------------------------------------
            //! display: second argument: Shaded mode: third argument: no selection allowed
            //! does not update the viewer (done at the end)
            //! ----------------------------------------------------------------------------
            occContext->Display(anOldShape1,1,-1,false);
        }
            break;
        }
        occContext->UpdateCurrentViewer();
    }
}


void occPreGLWidget::displayShapeCopy2(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color)
{
    //cout<<"occPreGLWidget::displayShapeCopy()->____function called____"<<endl;

    TopoDS_Compound myCompound1;
    static occHandle(AIS_Shape) anOldShape1;

    //! ------------------------------------------
    //! do not update the viewer (done at the end)
    //! ------------------------------------------
    if(!anOldShape1.IsNull()) occContext->Remove(anOldShape1,false);

    if(!listShapes.IsEmpty())
    {
        //! --------------------------------------------
        //! get the type of shape contained in the list
        //! --------------------------------------------
        TopAbs_ShapeEnum typeList1 = listShapes.First().ShapeType();
        TopoDS_Builder myBuilder1;
        switch(typeList1)
        {
        case TopAbs_VERTEX:
        {
            ;
        }
        break;

        default:
        {
            TopTools_ListIteratorOfListOfShape anIter;
            myBuilder1.MakeCompound(myCompound1);
            for(anIter.Initialize(listShapes);anIter.More();anIter.Next())
            {
                myBuilder1.Add(myCompound1,anIter.Value());
            }
            anOldShape1 = new AIS_Shape(myCompound1);

            //! ----------------------
            //! use the model z-layer
            //! ----------------------
            //anOldShape1->SetZLayer(Graphic3d_ZLayerId_TopOSD);
            anOldShape1->SetColor(color);
            anOldShape1->Attributes()->ShadingAspect()->SetTransparency(0.0,Aspect_TOFM_FRONT_SIDE);
            anOldShape1->Attributes()->ShadingAspect()->SetTransparency(1.0,Aspect_TOFM_BACK_SIDE);
            //! ----------------------------------------------------------------------------
            //! display: second argument: Shaded mode: third argument: no selection allowed
            //! does not update the viewer (done at the end)
            //! ----------------------------------------------------------------------------
            occContext->Display(anOldShape1,1,-1,false);
        }
            break;
        }
        occContext->UpdateCurrentViewer();
    }
}

//! ----------------------
//! function: key pressed
//! details:
//! ----------------------
void occPreGLWidget::keyPressEvent(QKeyEvent *theKey)
{
    switch(theKey->key())
    {
    case Qt::Key_F3: extendSelectionToAjacent(); break;
    case Qt::Key_F4:
    {
        cout<<"____F4 pressed____"<<endl;
    }
        break;
    case Qt::Key_F7: occView->FitAll(); break;
    case Qt::Key_F10: emit highlightmeshface(); break;
    case Qt::Key_F11: emit highlightmeshedge(); break;

        //switch(myCurAction3D)
        //{
        //case(CurAction3D_Rotation): setAction3D_Rotation(); break;
        //case(CurAction3D_Panning): setAction3D_Pan(); break;
        //case(CurAction3D_WindowZooming): setAction3D_WindowZooming(); break;
        //}
        //break;

    case Qt::Key_F8:
    {
        cout<<"____F8 pressed: displaying triangulation____"<<endl;
        for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
        {
            const occHandle(AIS_Shape) &anAIS_Shape = occHandle(AIS_Shape)::DownCast(occContext->SelectedInteractive());
            const TopoDS_Shape &curShape = anAIS_Shape->Shape();
            this->displayTriangulation(curShape);
        }
        occContext->UpdateCurrentViewer();
        this->clearSelected();
    }
        break;
    case Qt::Key_F9:
    {
        this->hideSelectedBodies();
        emit requestSynchVisibility();
    }
        break;
    case (Qt::Key_B): if(theKey->modifiers()==Qt::ControlModifier) emit selectionModeSolid(true); break;
    case (Qt::Key_F): if(theKey->modifiers()==Qt::ControlModifier) emit selectionModeFace(true); break;
    case (Qt::Key_E): if(theKey->modifiers()==Qt::ControlModifier) emit selectionModeEdge(true); break;
    case (Qt::Key_P): if(theKey->modifiers()==Qt::ControlModifier) emit selectionModeVertex(true); break;
    case(Qt::Key_C): if(theKey->modifiers()==Qt::ControlModifier) break;
    case(Qt::Key_N): if(theKey->modifiers()==Qt::ControlModifier) break;
    case(Qt::Key_M): if(theKey->modifiers()==Qt::ControlModifier) break;
    case (Qt::Key_A):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            const Graphic3d_RenderingParams &renPar = occView->RenderingParams();
            if(renPar.Method==Graphic3d_RM_RASTERIZATION)
            {
                occView->ChangeRenderingParams().NbMsaaSamples=128;
                occView->ChangeRenderingParams().Method = Graphic3d_RM_RAYTRACING;
                occView->ChangeRenderingParams().IsShadowEnabled = Standard_True;
                QMessageBox::information(this,APPNAME,"Raytracing, Antialiasing, Shadows ON",QMessageBox::Ok);
            }
            else
            {
                occView->ChangeRenderingParams().NbMsaaSamples=0;
                occView->ChangeRenderingParams().Method = Graphic3d_RM_RASTERIZATION;
                occView->ChangeRenderingParams().IsShadowEnabled = Standard_False;
                QMessageBox::information(this,"APPNAME","Raytracing, Antialiasing, Shadows OFF",QMessageBox::Ok);
            }
            occView->Redraw();
        }
        break;

    case(Qt::Key_S):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            if(occView->ShadingModel()==V3d_GOURAUD)
            {
                occView->SetShadingModel(V3d_PHONG);
                QMessageBox::information(this,APPNAME,"Shading mode Phong",QMessageBox::Ok);
            }
            else
            {
                occView->SetShadingModel(V3d_GOURAUD);
                QMessageBox::information(this,APPNAME,"Shading mode Gouraud",QMessageBox::Ok);

            }
            occView->Update();
        }
    }
}

//! ----------------------------------------------------------------------
//! function: replaceTriangulation
//! details:  replace the triangulation with an externally generated mesh
//! ----------------------------------------------------------------------
void occPreGLWidget::replaceTriangulation()
{
    cout<<"occPreGLWidget::replaceTriangulation()->____function called____"<<endl;

    //! -----------------------------------------------------------------------
    //! retrieve the AIS_Shape with signature "0" ("Shape") within the context
    //! -----------------------------------------------------------------------
    AIS_ListOfInteractive listOfAIS;
    AIS_ListIteratorOfListOfInteractive anIter;
    occContext->ObjectsInside(listOfAIS, AIS_KOI_Shape, 0);
    for(anIter.Initialize(listOfAIS); anIter.More(); anIter.Next())
    {
        //! -----------------------------------------------------------
        //! check if the body has an externally generated surface mesh
        //! -----------------------------------------------------------
        const occHandle(AIS_Shape) &theCurAIS = occHandle(AIS_Shape)::DownCast(anIter.Value());
        const TopoDS_Shape &curShape = theCurAIS->Shape();
        int bodyIndex = myDS2->bodyMap.key(curShape);

        if(!myDS2->ArrayOfMeshDS2D.value(bodyIndex).IsNull())
        {
            cout<<"occPreGLWidget::replaceTriangulation()->____rebuild using an externally generated mesh____"<<endl;

            //! -----------------------------
            //! clean the BRep triangulation
            //! -----------------------------
            BRepTools::Clean(curShape);
            int NbFaces = myDS2->MapOfBodyTopologyMap.value(bodyIndex).faceMap.Extent();
            for(int faceNr=1; faceNr<=NbFaces; faceNr++)
            {
                const occHandle(Ng_MeshVS_DataSourceFace) &faceMesh =
                        occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDS2->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceNr));

                if(!faceMesh.IsNull())
                {
                    occHandle(Poly_Triangulation) theTriangulation;
                    MeshTools::toPolyTriangulation(faceMesh, theTriangulation);
                    TopoDS_Face curFace = TopoDS::Face(myDS2->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindKey(faceNr));
                    //! ----------------------------
                    //! update the mesh of the face
                    //! ----------------------------
                    BRepMesh_ShapeTool::AddInFace(curFace, theTriangulation);
                }
            }
        }
        else
        {
            //! -------------------------------------------------------------
            //! no FEM mesh for replacing the tessellation: check if
            //! the shape has a triangulation (this could have been deleted)
            //! -------------------------------------------------------------
            if(!BRepTools::Triangulation(curShape,1e10))
            {
                //! the BRep mesh was previously cleaned: rebuild the visualization mesh
                double theLinDeflection = 0.1;
                bool isRelative = true;
                double theAngDeflection = 0.3;  //! 17.1 deg
                bool isInParallel = true;
                BRepMesh_IncrementalMesh(curShape, theLinDeflection, isRelative, theAngDeflection, isInParallel);
            }
        }
    }
}

//! -------------------------------
//! function: displayTriangulation
//! details:
//! -------------------------------
void occPreGLWidget::displayTriangulation(const TopoDS_Shape &aShape)
{
    BRep_Builder builder;
    TopoDS_Compound Comp;
    builder.MakeCompound(Comp);
    for (TopExp_Explorer ex(aShape,TopAbs_FACE) ; ex.More(); ex.Next())
    {
        TopoDS_Face F = TopoDS::Face(ex.Current());
        //cout<<"retrieving triangulation from face Nr: "<<myDS2->MapOfBodyTopologyMap.value(bodyIndex).faceMap.FindIndex(ex.Current())<<endl;
        TopLoc_Location L;
        Handle (Poly_Triangulation) facing = BRep_Tool::Triangulation(F,L);
        TColgp_Array1OfPnt tab(1,(facing->NbNodes()));
        tab = facing->Nodes();
        Poly_Array1OfTriangle tri(1,facing->NbTriangles());
        tri = facing->Triangles();
        for (Standard_Integer i=1;i<=(facing->NbTriangles());i++)
        {
            Poly_Triangle trian = tri.Value(i);
            Standard_Integer index1,index2,index3,M,N;
            trian.Get(index1,index2,index3);
            for (Standard_Integer j=1;j<=3;j++)
            {
                switch (j)
                {
                case 1: M = index1; N = index2; break;
                case 2: N = index3; break;
                case 3: M = index2; break;
                }
                BRepBuilderAPI_MakeEdge ME(tab.Value(M),tab.Value(N));
                if (ME.IsDone())
                {
                    builder.Add(Comp,ME.Edge());
                }
            }
        }
    }
    occHandle(AIS_Shape) c = new AIS_Shape(Comp);
    occContext->Display(c,Standard_True);
    occContext->UpdateCurrentViewer();
}

//! --------------------
//! function: esportSTL
//! detail:
//! --------------------
void occPreGLWidget::exportSTL()
{
    //cout<<"occPreGLWidget::exportSTL()->____function called____"<<endl;
    exportingTools::exportSTL(occContext, myDS2);
}

//! ---------------------------------------------------------------
//! function: exportCloud
//! details:  write the nodes coordinates of the mesh of the
//!           selected entities - now working with faces and edges
//! ---------------------------------------------------------------
void occPreGLWidget::exportCloud()
{
    //cout<<"occPreGLWidget::exportCloud()->____function called____"<<endl;
    exportingTools::exportCloud(occContext, myDS2);
}

//! ----------------------
//! function: exportCloud
//! details:
//! ----------------------
bool occPreGLWidget::event(QEvent *event)
{
    //cout<<"occPreGLWidget::event()->____event: "<<event->type()<<"____"<<endl;
    if(event->type()==QBackgroundEvent::type())
    {
        cout<<"occPreGLWidget::event()->____backgound event change received____"<<endl;
        QBackgroundEvent *bkgEvent = static_cast<QBackgroundEvent*>(event);
        int R1 = bkgEvent->R1();
        int G1 = bkgEvent->G1();
        int B1 = bkgEvent->B1();
        int R2 = bkgEvent->R2();
        int G2 = bkgEvent->G2();
        int B2 = bkgEvent->B2();
        int gradient = bkgEvent->getGradient();
        this->setBackgroundColor(R1/255.0,G1/255.0,B1/255.0,R2/255.0,G2/255.0,B2/255.0,gradient);
        return true;
    }
    else
    {
        //! Make sure the rest of events are handled
        return occGLWidget::event(event);
    }
}

//! -----------------------
//! function: focusInEvent
//! details:
//! -----------------------
void occPreGLWidget::focusInEvent(QFocusEvent *focusEvent)
{
    Q_UNUSED(focusEvent)

    cerr<<"___focus on widget: "<<this->objectName().toStdString()<<"____"<<endl;
    cerr<<"___"<<focusEvent->MouseButtonPress<<"____"<<endl;
    emit requestChangeDelegateContext(this->getContext());
}

//! -------------------------------
//! function: displayCurvatureMap
//! details:  for testing purposes
//! -------------------------------
void occPreGLWidget::displayCurvatureMap()
{
    cout<<"occPreGLWidget::displayCurvatureMap()->____function called____"<<endl;
    this->hideAllMeshes();

    occContext->InitSelected();
    if(occContext->MoreSelected())
    {
        //! --------------------------------
        //! retrieve the selected face tags
        //! --------------------------------
        std::vector<std::pair<int,int>> vecPairs;
        if(occContext->SelectedShape().ShapeType()==TopAbs_FACE)
        {
            std::pair<int,int> pair;
            for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
            {
                const occHandle(AIS_ExtendedShape) &AIS_SHape = occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive());
                const TopoDS_Shape &theParentShape =AIS_SHape->Shape();
                const TopoDS_Shape &theFace = occContext->SelectedShape();
                this->getTopologyNumber(theParentShape,theFace,pair.first,pair.second);
                vecPairs.push_back(pair);
            }
        }

        this->emptyTheSelection();

        //! ------------------------------
        //! pile up the mesh data sources
        //! ------------------------------
        QList<occHandle(Ng_MeshVS_DataSourceFace)> faceMeshDSlist;
        for(int k=0; k<vecPairs.size(); k++)
        {
            int bodyIndex = vecPairs.at(k).first;
            int faceIndex = vecPairs.at(k).second;

            const occHandle(Ng_MeshVS_DataSourceFace) &curFaceMesh =
                    occHandle(Ng_MeshVS_DataSourceFace)::DownCast(myDS2->ArrayOfMeshDSOnFaces.getValue(bodyIndex,faceIndex));
            faceMeshDSlist<<curFaceMesh;
        }

        //! -------------------------------
        //! sum the face mesh data sources
        //! -------------------------------
        occHandle(Ng_MeshVS_DataSourceFace) summedDS = new Ng_MeshVS_DataSourceFace(faceMeshDSlist);

        occHandle(MeshVS_Mesh) coloredMesh;
        int mode = 1;   //! Gaussian curvature
        summedDS->computeDiscreteCurvature(mode);
        MeshTools::buildColoredMesh(summedDS,summedDS->myCurvature,coloredMesh,0,6.28,10,true);
        //MeshTools::buildColoredMesh(summedDS,summedDS->myCurvatureGradient,coloredMesh,0,6.28,10,true);

        occMeshContext->Display(coloredMesh,false);
        occMeshContext->UpdateCurrentViewer();
    }
}

//! ----------------------
//! function: displayMesh
//! details:
//! ----------------------
void occPreGLWidget::displayMesh(const occHandle(AIS_InteractiveContext) &aCTX,
                                 Quantity_NameOfColor aColor,
                                 const occHandle(MeshVS_DataSource) &aMeshDS)
{
    cout<<"occPreGLWidget::displayMesh()->____function called____"<<endl;

    //! ----------------------------
    //! the mesh interactive object
    //! ----------------------------
    occHandle(MeshVS_Mesh) aMeshIO = new MeshVS_Mesh();
    aMeshIO->SetDataSource(aMeshDS);

    //! ----------------------------------------
    //! create and add the presentation builder
    //! ----------------------------------------
    occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(aMeshIO);
    aMeshIO->AddBuilder(aBuilder,Standard_False);

    //! --------------------------------------
    //! cosmetic for displaying the face mesh
    //! --------------------------------------
    Graphic3d_MaterialAspect myAspect(Graphic3d_NOM_GOLD);
    myAspect.SetColor(static_cast<Quantity_Color>(aColor));

    aMeshIO->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial, myAspect);
    aMeshIO->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes, Standard_False);
    aMeshIO->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges, Standard_True);
    aMeshIO->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
    aMeshIO->SetDisplayMode(MeshVS_DMF_Shading);

    aCTX->Display(aMeshIO,Standard_True);
    cout<<"occPreGLWidget::displayMesh()->____mesh displayed____"<<endl;
}

//! ------------------------------------------
//! function: updateViewerAfterDataBasechange
//! details:
//! ------------------------------------------
void occPreGLWidget::updateViewerAfterDataBaseChange()
{
    cout<<"occPreGLWidget::updateViewerAfterDataBaseChange()->____function called____"<<endl;
    if(myDS2!=NULL)
    {
        //! ------------------------------------------------
        //! clear the maps of interactive shapes and meshes
        //! ------------------------------------------------
        myMapOfInteractiveShapes.clear();
        myMapOfInteractiveShapes.clear();
        //! -----------------------------------------------
        //! remove the interactive objects from the viewer
        //! -----------------------------------------------
        occContext->RemoveAll(true);
        occMeshContext->RemoveAll(true);

        //! ------------------------------------------
        //! create the interactive shapes and display
        //! ------------------------------------------
        this->createInteractiveShapes();
        this->displayCAD();
    }
}

//! ------------------------------
//! function: showMeshDataSources
//! details:  for developers(?)
//! ------------------------------
void occPreGLWidget::showMeshDataSources(const IndexedMapOfMeshDataSources &indexedMapOfDS)
{
    cout<<"occPreGLWidget::showMeshDataSources()->____function called____"<<endl;
    this->emptyTheSelection();
    for(IndexedMapOfMeshDataSources::const_iterator it = indexedMapOfDS.cbegin(); it!=indexedMapOfDS.cend(); ++it)
    {
        const occHandle(MeshVS_DataSource) &aDS = it.value();
        occPreGLWidget::displayMesh(occContext,Quantity_NOC_GREEN,aDS);
    }
}

//! ----------------------
//! function: displayMesh
//! details:  utility
//! ----------------------
#include <Geom_Plane.hxx>
void occPreGLWidget::displayMesh(const occHandle(MeshVS_DataSource) &aMeshDS,
                                 Quantity_NameOfColor aColorName,
                                 bool showMeshEdges)
{
    cout<<"occPreGLWidget::displayMesh()->____function called____"<<endl;

    occHandle(MeshVS_Mesh) meshIO = new MeshVS_Mesh();
    meshIO->SetDataSource(aMeshDS);
    occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(meshIO);
    meshIO->AddBuilder(aBuilder,false);
    Graphic3d_MaterialAspect aspect(Graphic3d_NOM_GOLD);
    aspect.SetColor(aColorName);

    meshIO->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial,aspect);
    meshIO->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes,false);
    meshIO->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges,showMeshEdges);
    meshIO->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
    meshIO->SetDisplayMode(MeshVS_DMF_Shading);
    occContext->Display(meshIO,true);
}

//! ----------------------------------------------------------------------------------
//! function: addClipPlane
//! details:  overrides the base class function. Here the clipped objects
//!           are retrieved not using the context functions (such ::ObjectsInside
//!           because there is not an AIS_KOI for MeshVS_Mesh), but using the
//!           private map members myMapOfInteractiveShapes & myMapOfInteractiveMeshes
//! ----------------------------------------------------------------------------------
void occPreGLWidget::addClipPlane(double A, double B, double C, double D, int ID, bool isOn)
{
    cout<<"occPreGLWidget::addClipPlane()->____function called____"<<endl;

    //! ------------------------------------------------------------------
    //! add the clip plane only if it is not already contained in the map
    //! ------------------------------------------------------------------
    occHandle(Graphic3d_ClipPlane) aClipPlane;
    if(myMapOfClipPlanes.value(ID,aClipPlane).IsNull())
    {
        aClipPlane = new Graphic3d_ClipPlane();
        myMapOfClipPlanes.insert(ID,aClipPlane);
    }
    else
    {
        aClipPlane = myMapOfClipPlanes.value(ID);
    }
    aClipPlane->SetEquation (gp_Pln (A,B,C,D));

    //! ------------
    //! set capping
    //! ------------
    bool cappingOn = false;
    if(myCurWorkingMode == curWorkingMode_onModel || myCurWorkingMode == curWorkingMode_onContact) cappingOn = true;
    aClipPlane->SetCapping(cappingOn);

    //! ---------------------------------------
    //! add the clip plane to the AIS_Shape(s)
    //! ---------------------------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); it++)
    {
        const occHandle(AIS_InteractiveObject) &curShapeObject = it.value();
        if(curShapeObject.IsNull()) continue;
        if(!curShapeObject->ClipPlanes().IsNull()) curShapeObject->ClipPlanes()->Clear();
        for(QMap<int,occHandle(Graphic3d_ClipPlane)>::iterator it = myMapOfClipPlanes.begin(); it!=myMapOfClipPlanes.end(); it++)
        {
            const occHandle(Graphic3d_ClipPlane) &curPlane = it.value();
            curShapeObject->AddClipPlane(curPlane);
        }
    }

    //! ---------------------------------------------------
    //! add the clip plane to the interactive mesh objects
    //! ---------------------------------------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveMeshes.begin(); it!=myMapOfInteractiveMeshes.end(); it++)
    {
        occHandle(AIS_InteractiveObject) curMeshObject = it.value();
        if(curMeshObject.IsNull()) continue;
        if(!curMeshObject->ClipPlanes().IsNull()) curMeshObject->ClipPlanes()->Clear();
        for(QMap<int,occHandle(Graphic3d_ClipPlane)>::iterator it = myMapOfClipPlanes.begin(); it!=myMapOfClipPlanes.end(); it++)
        {
            const occHandle(Graphic3d_ClipPlane) &curPlane = it.value();
            curMeshObject->AddClipPlane(curPlane);
        }
    }

    //! -------------------------------------
    //! add the clip plane to the whole view
    //! left here for documentation
    //! -------------------------------------
    //occView->AddClipPlane(aClipPlane);

    //! ------------------------------------------------
    //! activate the clipping plane and update the view
    //! ------------------------------------------------
    aClipPlane->SetOn(isOn);
    occView->Redraw();

    //! --------------------------------------
    //! add/replace the clip plane to the map
    //! --------------------------------------
    myMapOfClipPlanes.insert(ID,aClipPlane);
    cout<<"____final number of clip planes: "<<myMapOfClipPlanes.size()<<"____"<<endl;

    occHandle(AIS_Plane) anAISPlane;
    if(!myMapOfHandlePlanes.contains(ID))
    {
        occHandle(Geom_Plane) aGeomPlane = new Geom_Plane(A,B,C,D);
        anAISPlane = new AIS_Plane(aGeomPlane);
        myMapOfHandlePlanes.insert(ID,anAISPlane);
        occContext->Display(anAISPlane,AIS_Shaded);
    }
    else
    {
        occHandle(Geom_Plane) aGeomPlane = new Geom_Plane(A,B,C,D);
    }
}

//! ----------------------------------------------------------
//! function: updateClipPlanes
//! details:  remove all the clip planes from all the objects
//!           and activate the "activeClipPlanes"
//! ----------------------------------------------------------
void occPreGLWidget::updateClipPlanes(const std::vector<int> &activeClipPlanes)
{
    if(activeClipPlanes.size()==0) return;

    //! ---------------------------------------------------
    //! remove the clip plane from each of the AIS_Shape
    //! remove the clip plane from the MeshVS_Mesh objects
    //! ---------------------------------------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); it++)
    {
        const occHandle(AIS_InteractiveObject) &curShapeObject = it.value();
        if(curShapeObject.IsNull()) continue;
        const occHandle(Graphic3d_SequenceOfHClipPlane) &shapeClipPlanes = curShapeObject->ClipPlanes();
        for(int i = shapeClipPlanes->Lower(); i<=shapeClipPlanes->Upper(); i++)
        {
            const occHandle(Graphic3d_ClipPlane) &shapeClipPlane = shapeClipPlanes->Value(i);
            int clipPlaneNb = myMapOfClipPlanes.key(shapeClipPlane);
            if(std::find(activeClipPlanes.begin(),activeClipPlanes.end(),clipPlaneNb)!=activeClipPlanes.end()) shapeClipPlane->SetOn(true);
            else shapeClipPlane->SetOn(false);
        }
    }
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveMeshes.begin(); it!=myMapOfInteractiveMeshes.end(); it++)
    {
        const occHandle(AIS_InteractiveObject) &curMeshObject = it.value();
        if(curMeshObject.IsNull()) continue;
        const occHandle(Graphic3d_SequenceOfHClipPlane) &meshObjectClipPlanes = curMeshObject->ClipPlanes();

        if(meshObjectClipPlanes.IsNull())
        {
            //! ------------------------------------------------------
            //! the mesh interactive object has no clip plane defined
            //! ------------------------------------------------------
            //cout<<"____the mesh objects have not planes: adding____"<<endl;
            for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator itMesh=myMapOfInteractiveMeshes.begin(); itMesh!=myMapOfInteractiveMeshes.end(); itMesh++)
            {
                const occHandle(AIS_InteractiveObject) &aMeshObject = it.value();
                for(QMap<int,occHandle(Graphic3d_ClipPlane)>::iterator itCP = myMapOfClipPlanes.begin(); itCP != myMapOfClipPlanes.end(); itCP++)
                {
                    const occHandle(Graphic3d_ClipPlane) &curClipPlane = itCP.value();
                    int clipPlaneNr = itCP.key();
                    if(std::find(activeClipPlanes.cbegin(),activeClipPlanes.cend(),clipPlaneNr)==activeClipPlanes.cend()) curClipPlane->SetOn(true);
                    else curClipPlane->SetOn(false);
                    aMeshObject->AddClipPlane(curClipPlane);
                }
            }
            return;
        }
        for(int i = meshObjectClipPlanes->Lower(); i<=meshObjectClipPlanes->Upper(); i++)
        {
            const occHandle(Graphic3d_ClipPlane) &shapeClipPlane = meshObjectClipPlanes->Value(i);
            int clipPlaneNb = myMapOfClipPlanes.key(shapeClipPlane);
            if(std::find(activeClipPlanes.begin(),activeClipPlanes.end(),clipPlaneNb)!=activeClipPlanes.end()) shapeClipPlane->SetOn(true);
            else shapeClipPlane->SetOn(false);
        }
    }
}

//! ------------------------------------------
//! function: removeClipPlane
//! details:  overrides the base class method
//! ------------------------------------------
void occPreGLWidget::removeClipPlane(int ID)
{
    cout<<"occPreGLWidget::removeClipPlane()->____removing clip plane ID: "<<ID<<"____"<<endl;

    const occHandle(Graphic3d_ClipPlane) &clipPlane = myMapOfClipPlanes.value(ID);

    //! -------------------------------------------------
    //! remove the clip plane from each of the AIS_Shape
    //! -------------------------------------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveShapes.begin(); it!=myMapOfInteractiveShapes.end(); it++)
    {
        const occHandle(AIS_InteractiveObject) &curShapeObject = it.value();
        if(curShapeObject.IsNull()) continue;
        curShapeObject->RemoveClipPlane(clipPlane);
    }
    //! --------------------------------------------------
    //! remove the clip plane frm the MeshVS_Mesh objects
    //! --------------------------------------------------
    for(QMap<int,occHandle(AIS_InteractiveObject)>::iterator it = myMapOfInteractiveMeshes.begin(); it!=myMapOfInteractiveMeshes.end(); it++)
    {
        const occHandle(AIS_InteractiveObject) &curMeshObject = it.value();
        if(curMeshObject.IsNull()) continue;
        curMeshObject->RemoveClipPlane(clipPlane);
    }

    //! ---------------------------------------------
    //! remove the clip plane ID from the whole view
    //! left here for documentation
    //! ---------------------------------------------
    //occView->RemoveClipPlane(clipPlane);

    //! -------
    //! redraw
    //! -------
    occView->Redraw();

    //! -----------------------------------
    //! remove the clip plane from the map
    //! -----------------------------------
    myMapOfClipPlanes.remove(ID);
    cout<<"____final number of clip planes: "<<myMapOfClipPlanes.size()<<"____"<<endl;
}

//! ------------------------------------------
//! function: updateClipPlaneTranslation
//! details:  overrides the base class method
//! ------------------------------------------
void occPreGLWidget::updateClipPlaneTranslation(int ID, int zVal, const QVector<double> &coeffs)
{
    //cout<<"occGLWidget::updateClipPlaneTranslation()->____function called. ID: "<<ID<<" zVal: "<<zVal<<"____"<<endl;

    double a = coeffs.at(0);
    double b = coeffs.at(1);
    double c = coeffs.at(2);
    double d = coeffs.at(3);

    gp_Pln aPlane(a,b,c,d);
    gp_Ax1 planeAxis = aPlane.Axis();
    gp_Dir translationDirection = planeAxis.Direction();
    gp_Vec translationVector(translationDirection);
    double deltaZ = double(zVal)*1000.0/1000.0;

    translationVector.Scale(deltaZ);

    const occHandle(Graphic3d_ClipPlane) &curClipPlane = myMapOfClipPlanes.value(ID);
    curClipPlane->SetEquation(aPlane.Translated(translationVector));
    occView->Redraw();
}

//! ----------------------------
//! function: buildSlicedMeshIO
//! details:
//! ----------------------------
void occPreGLWidget::buildSlicedMeshIO(const QMap<int,occHandle(MeshVS_DataSource)> &slicedMeshDS)
{
    cout<<"occPreGLWidget::buildSlicedMeshIO()->____function called____"<<endl;
    for(QMap<int,occHandle(MeshVS_DataSource)>::const_iterator it= slicedMeshDS.cbegin(); it!=slicedMeshDS.cend(); it++)
    {
        occHandle(Ng_MeshVS_DataSource3D) aSlicedMeshDS = occHandle(Ng_MeshVS_DataSource3D)::DownCast(it.value());
        if(aSlicedMeshDS.IsNull()) continue;

        //! ----------------------------------------
        //! retrieve the color of the current shape
        //! ----------------------------------------
        int bodyIndex = it.key();
        Quantity_NameOfColor colorName = this->shapeColor(bodyIndex).Name();

        occHandle(MeshVS_Mesh) meshIO = new MeshVS_Mesh();
        meshIO->SetDataSource(aSlicedMeshDS);
        occHandle(MeshVS_MeshPrsBuilder) aBuilder = new MeshVS_MeshPrsBuilder(meshIO);
        meshIO->AddBuilder(aBuilder,false);

        Graphic3d_MaterialAspect aspect(Graphic3d_NOM_GOLD);
        aspect.SetColor(colorName);

        meshIO->GetDrawer()->SetMaterial(MeshVS_DA_FrontMaterial,aspect);
        meshIO->GetDrawer()->SetBoolean(MeshVS_DA_DisplayNodes,false);
        meshIO->GetDrawer()->SetBoolean(MeshVS_DA_ShowEdges,true);
        meshIO->GetDrawer()->SetColor(MeshVS_DA_EdgeColor,Quantity_NOC_BLACK);
        meshIO->SetDisplayMode(MeshVS_DMF_Shading);

        mySlicedMeshesIO.insert(bodyIndex,meshIO);
        occContext->Display(meshIO,true);
    }
    cout<<"occPreGLWidget::buildSlicedMeshIO()->____exiting function____"<<endl;
}

//! ----------------------------
//! function: eraseSlicedMeshes
//! details:
//! ----------------------------
void occPreGLWidget::eraseSlicedMeshes()
{
    bool updateViewer = false;
    for(QMap<int,occHandle(MeshVS_Mesh)>::iterator it=mySlicedMeshesIO.begin(); it!=mySlicedMeshesIO.end(); ++it)
    {
        occContext->Erase(it.value(),updateViewer);
    }
}

//! ----------------------------
//! function: applyCustomColors
//! details:
//! ----------------------------
#include <StdSelect_BRepOwner.hxx>
void occPreGLWidget::applyCustomColors(const QMap<GeometryTag,TopoDS_Shape> &subShapeMaps, Quantity_NameOfColor aColor, bool updateViewer)
{
    if(subShapeMaps.isEmpty()) return;
    //if(AIS_DisplayMode(occContext->DisplayMode())==AIS_Shaded)
    //{
        for(QMap<GeometryTag,TopoDS_Shape>::const_iterator it = subShapeMaps.cbegin(); it!=subShapeMaps.cend(); it++)
        {
            const GeometryTag &aTag = it.key();
            const TopoDS_Shape &curSubShape = it.value();
            int index = aTag.parentShapeNr;
            const occHandle(AIS_ColoredShape) &aColored = occHandle(AIS_ColoredShape)::DownCast(myMapOfInteractiveShapes.value(index));
            aColored->SetCustomColor(curSubShape,aColor);
            occContext->RecomputePrsOnly(aColored,true,true);
        }
    //}

    /*
    static occHandle(AIS_Shape) tempShape;
    if(AIS_DisplayMode(occContext->DisplayMode())==AIS_WireFrame)
    {
        TopoDS_Compound aComp;
        TopoDS_Builder aBuilder;
        aBuilder.MakeCompound(aComp);
        for(QMap<GeometryTag,TopoDS_Shape>::const_iterator it = subShapeMaps.cbegin(); it!=subShapeMaps.cend(); it++)
        {
            const TopoDS_Shape &curSubShape = it.value();
            aBuilder.Add(aComp,curSubShape);
        }
        tempShape = new AIS_Shape(aComp);
        occContext->SetColor(tempShape,aColor,false);
        occContext->Display(tempShape,AIS_Shaded,-1,false,false);
    }
    AIS_ListOfInteractive l;
    occContext->ObjectsInside(l,AIS_KOI_Shape,0);
    cout<<"@ ----------------------"<<endl;
    cout<<"@ objects inside: "<<l.Size()<<endl;
    cout<<"@ ----------------------"<<endl;
    */
    if(updateViewer==true) occContext->UpdateCurrentViewer();
}
