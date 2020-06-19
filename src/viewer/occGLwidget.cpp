#define ValZWMin 1
#define TRANSPARENCY 0.0
#define limitAngle 30.1     //! the limit angle for adjacent selection

//! ----------------
//! custom includes
//! ----------------
#include "occGLwidget.h"
#include "ais_extendedshape.h"
#include "geomtoolsclass.h"
#include "tools.h"
#include "mydefines.h"
#include "ais_cormarker.h"
#include "ais_spheremarker.h"
#include "markers.h"

//! ---
//! Qt
//! ---
#include <QFileDialog>
#include <QMouseEvent>
#include <QRubberBand>
#include <QStyleFactory>
#include <QMenu>
#include <QActionGroup>
#include <QDebug>
#include <QCursor>
#include <QDateTime>
#include <QCursor>
#include <QBitmap>
#include <QPixmap>
#include <QMessageBox>

//! ----
//! OCC
//! ----
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopOpeBRepBuild_Tools.hxx>
#include <TopExp.hxx>
#include <TColStd_ListOfInteger.hxx>
#include <TColStd_ListIteratorOfListOfInteger.hxx>
#include <TopoDS_ListOfShape.hxx>
#include <TopoDS_ListIteratorOfListOfShape.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <AIS_RubberBand.hxx>
#include <Graphic3d_RenderingParams.hxx>
#include <Graphic3d_RenderingMode.hxx>
#include <Graphic3d_AspectFillArea3d.hxx>
#include <Graphic3d_Vec2.hxx>
#include <V3d_TypeOfShadingModel.hxx>
#include <OpenGl_GraphicDriver.hxx>
#include <OpenGl_Caps.hxx>
#include <Aspect_TypeOfHighlightMethod.hxx>
#include <AIS_StatusOfPick.hxx>
#include <AIS_Shape.hxx>
#include <AIS_ListOfInteractive.hxx>
#include <AIS_ListIteratorOfListOfInteractive.hxx>
#include <Aspect_TypeOfTriedronPosition.hxx>
#include <Prs3d_ShadingAspect.hxx>
#include <GProp_GProps.hxx>
#include <BRep_Tool.hxx>
#include <BRepGProp.hxx>
#include <gp_Pnt.hxx>
#include <Prs3d_LineAspect.hxx>
#include <ProjLib.hxx>
#include <GC_MakeLine.hxx>
#include <GeomAPI_IntCS.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <ElSLib.hxx>
#include <Graphic3d_Camera.hxx>
#include <gp_Dir.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <Font_FontAspect.hxx>
#include <Graphic3d_TransformPers.hxx>
#include <Geom_Surface.hxx>
#include <Quantity_Color.hxx>
#include <Quantity_TypeOfColor.hxx>
#include <Aspect_GradientFillMethod.hxx>
#include <Prs3d_IsoAspect.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <Geom_Axis2Placement.hxx>
#include <Prs3d_ArrowAspect.hxx>
#include <Prs3d_DatumAspect.hxx>
#include <AIS_Trihedron.hxx>
#include <Aspect_WidthOfLine.hxx>
#include <AIS_KindOfInteractive.hxx>
#include <Graphic3d_MaterialAspect.hxx>
#include <Graphic3d_NameOfMaterial.hxx>

static double rx = 0.0;
static double ry = 0.0;
static int sx = 0.0;
static int sy = 0.0;
static bool zRotation = Standard_False;

//! ------------------------------------
//! function: GetGraphicDriver
//! detalis:  return the graphic driver
//! ------------------------------------
static occHandle(Graphic3d_GraphicDriver)& GetGraphicDriver()
{
    static occHandle(OpenGl_GraphicDriver) aGraphicDriver;
    //aGraphicDriver->ChangeOptions().ffpEnable=Standard_True;
    return aGraphicDriver;
}

//! --------------------
//! function: timeStamp
//! details:
//! --------------------
static QString timeStamp()
{
    QDateTime dateTime;
    QString dateFormat = "dd/MM/yyyy";
    QString dateString = dateTime.currentDateTime().toString(dateFormat);
    QString timeFormat = "hh:mm";
    QString timeString = dateTime.currentDateTime().toString(timeFormat);
    QString timeStamp;
    timeStamp.append("Model\n").append("Date: ").append(dateString).append("\n").append("Time: ").append(timeString);
    return timeStamp;
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
occGLWidget::occGLWidget(QWidget* parent):QGLWidget(parent),
    myLocalCtxNumber(0),
    myCurDisplayMode(CurDisplayMode_Wireframe),
    myCurDisplayMode_old(CurDisplayMode_Unset),
    myCurAction3D(CurAction3D_Nothing),
    myCurAction3D_old(CurAction3D_Nothing), //! experimental
    myAllowSinglePick(Standard_False),
    myCurGlobalSelectionMode(CurGlobalSelectionMode_Single),
    myCurSelectionMode(CurSelection_Nothing),
    myXmax(0),
    myXmin(0),
    myYmax(0),
    myYmin(0),
    myRectBand(NULL),
    myDisplayQuality(displayQuality_Medium),
    myRotationPointType(RotationPointType_gravity),
    myCurClipPlaneID(-1) //! experimental
{
    cout<<"occGLWidget::occGLWidget()->____CONSTRUCTOR CALLED____"<<endl;

    //! ------------------------
    //! activate mouse tracking
    //! ------------------------
    this->setMouseTracking(true);

    myContextMenu = new QMenu(this);

    myCursorModeMenu = new QMenu("Cursor mode",this);
    myViewModeMenu = new QMenu("View", this);

    //! deactivate the context menu request
    //this->setContextMenuPolicy(Qt::CustomContextMenu);
    //connect(this,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(ShowContextMenu(const QPoint&)));

    //! -------------------
    //! create the actions
    //! -------------------
    this->createTheActions();

    //! ----------------------------
    //! Set the initial cursor type
    //! ----------------------------
    this->setCursor(Qt::ArrowCursor);

    //! ------------------------------------------------
    //! the camera target - (0, 0, 0) at the beginning,
    //! when no object is diplayed
    //! ------------------------------------------------
    myCOR.SetX(0);
    myCOR.SetY(0);
    myCOR.SetZ(0);

    //! -----------------------------------------------
    //! COV (center of view if the camera) has changed
    //! -----------------------------------------------
    connect(this,SIGNAL(CORchanged(gp_Pnt,bool)),this,SLOT(updateCOR(gp_Pnt,bool)));

    //! ------------------------
    //! the graphic window text
    //! ------------------------
    myTextLabel = new AIS_TextLabel;
    occHandle(Graphic3d_TransformPers) trs = new Graphic3d_TransformPers(Graphic3d_TMF_2d);
    trs->SetCorner2d(Aspect_TOTP_LEFT_UPPER);
    Graphic3d_Vec2i offset(20,50);
    trs->SetOffset2d(offset);

    myTextLabel->SetZLayer(Graphic3d_ZLayerId_TopOSD);
    myTextLabel->SetTransformPersistence(trs);

    Standard_CString theContent = timeStamp().toStdString().c_str();
    myTextLabel->SetText(theContent);
    myTextLabel->SetColor(Quantity_NOC_BLACK);
    myTextLabel->SetFontAspect(Font_FA_Regular);
    myTextLabel->SetFont("Arial");
    myTextLabel->SetHeight(14);

    myFloatingLabel = new AIS_TextLabel;
    myFloatingLabel->SetZLayer(Graphic3d_ZLayerId_TopOSD);
    myFloatingLabel->SetColor(Quantity_NOC_BLACK);
    myFloatingLabel->SetFontAspect(Font_FA_Regular);
    myFloatingLabel->SetFont("Arial");
    myFloatingLabel->SetHeight(14);
    myFloatingLabel->SetText("");

    /*
    static double rx = 0.;
    static double ry = 0.;
    static int sx = 0;
    static int sy = 0;
    static Standard_Boolean zRotation = Standard_False;
    */
}

//! ---------------------
//! function: destructor
//! ---------------------
occGLWidget::~occGLWidget()
{    
    cout<<"occGLWidget::~occGLWidget()->____DESTRUCTOR CALLED____"<<endl;
}

//! ----------------------------------------------------------------
//! function: createTheActions
//! details:  to be inserted in the context menu of the widget.
//!           no signal-slot connection is implemented and action
//!           Events - triggered() are handled using a switch-case
//!           structure (see the occGLWidget:showContextMenu(), and
//!           shortcuts through keyboard events. It's working fine
//! ----------------------------------------------------------------
void occGLWidget::createTheActions()
{
    //! create the actions
    selectionModeActionsGroup = new QActionGroup(this);
    selectVertex = new QAction(QIcon(":/icons/icon_select vertex.png"),"Vertex (Ctrl+P)",selectionModeActionsGroup);
    selectVertex->setCheckable(true);
    selectVertex->setData(0);

    selectEdge = new QAction(QIcon(":/icons/icon_select edge.png"),"Edge (Ctrl+E)",selectionModeActionsGroup);
    selectEdge->setCheckable(true);
    selectEdge->setData(1);

    selectFace = new QAction(QIcon(":/icons/icon_select face.png"),"Face (Ctrl+F)",selectionModeActionsGroup);
    selectFace->setCheckable(true);
    selectFace->setData(2);

    selectSolid = new QAction(QIcon(":/icons/icon_select solid.png"),"Solid  (Ctrl+B)",selectionModeActionsGroup);
    selectSolid->setCheckable(true);
    selectSolid->setData(3);

    selectPickPointCoordinates = new QAction(QIcon(":/icons/icon_pick point coordinates.png"),"Hit point coordinates",selectionModeActionsGroup);
    selectPickPointCoordinates->setCheckable(true);
    selectPickPointCoordinates->setData(4);

    //! adds the actions to the actions group
    selectionModeActionsGroup->addAction(selectVertex);
    selectionModeActionsGroup->addAction(selectEdge);
    selectionModeActionsGroup->addAction(selectFace);
    selectionModeActionsGroup->addAction(selectSolid);
    selectionModeActionsGroup->addAction(selectPickPointCoordinates);

    //! sets the exclusive mode of selection
    selectionModeActionsGroup->setExclusive(true);

    //! view actions
    viewFront= new QAction(QIcon(":/icons/icon_front view.png"),"Front",this);
    viewFront->setData(5);
    viewBack = new QAction(QIcon(":/icons/icon_back view.png"),"Back",this);
    viewBack->setData(6);
    viewRight = new QAction(QIcon(":/icons/icon_right view.png"),"Right",this);
    viewRight->setData(7);
    viewLeft = new QAction(QIcon(":/icons/icon_left view.png"),"Left",this);
    viewLeft->setData(8);
    viewTop = new QAction(QIcon(":/icons/icon_top view.png"),"Top",this);
    viewTop->setData(9);
    viewBottom = new QAction(QIcon(":/icons/icon_bottom view.png"),"Bottom",this);
    viewBottom->setData(10);
    viewIsometric = new QAction(QIcon(":/icons/icon_isometric view.png"),"Isometric view",this);
    viewIsometric->setData(11);
    actionZoomToFit = new QAction(QIcon(":/icons/icon_fit all.png"),"Zoom to Fit (F7)",this);
    actionZoomToFit->setData(12);

    //! show - hide bodies
    actionHideBody = new QAction(QIcon(":/icons/icon_lamp OFF.png"),"Hide body (F9)",this);
    actionHideBody->setData(13);

    //! action hide all other bodies
    actionHideAllOtherBodies = new QAction(QIcon(":/icons/icon_lamp OFF.png"),"Hide all other bodies",this);
    actionHideAllOtherBodies->setData(14);

    //! action show all bodies
    actionShowAllBodies = new QAction(QIcon(":/icons/icon_lamp ON.png"),"Show all bodies",this);
    actionShowAllBodies->setData(15);

    //! action select all
    actionSelectAll = new QAction(QIcon(":/icons/icon_select all.png"), QString("Select all"),this);
    actionSelectAll->setData(16);
}

//! ---------------
//! function: init
//! details:
//! ---------------
#include "wbtrihedron.h"
void occGLWidget::init()
{    
    //! --------------------------------
    //! Create Aspect_DisplayConnection
    //! --------------------------------
    occHandle(Aspect_DisplayConnection) aDisplayConnection = new Aspect_DisplayConnection();

    //! -----------------------
    //! Get the graphic driver
    //! -----------------------
    if (GetGraphicDriver().IsNull())
    {
        GetGraphicDriver() = new OpenGl_GraphicDriver(aDisplayConnection);
        GetGraphicDriver()->EnableVBO(true);
    }

    //! ----------------------
    //! Get the window handle
    //! ----------------------
    WId window_handle = (WId) winId();

    //! ------------------
    //! Create the window
    //! ------------------
    occHandle(WNT_Window) wind = new WNT_Window((Aspect_Handle) window_handle);

    //! ------------------------------
    //! Create V3dViewer and V3d_View
    //! ------------------------------
    occViewer = new V3d_Viewer(GetGraphicDriver()); // OCC 7.1.0
    occView = occViewer->CreateView();
    occView->SetWindow(wind);
    occView->ChangeRenderingParams().IsAntialiasingEnabled=true;
    if(!wind->IsMapped()) wind->Map();

    //! -----------------------------------------------------------
    //! Set background color: read it from the "settings.txt" file
    //! -----------------------------------------------------------
    QString fileSettingName = QString(SYSTEM_PROGRAM_DATA)+"/WB/settings.txt";
    cout<<"occGLWidget::init()->____"<<fileSettingName.toStdString()<<"____"<<endl;

    ifstream fileSettings;
    fileSettings.open(fileSettingName.toStdString());
    if(fileSettings.is_open())
    {
        //! ----------------------------------------------
        //! the settings file has been succesfully opened
        //! ----------------------------------------------
        int gradient, r1,g1,b1,r2,g2,b2;
        char tmp[64];        
        std::string val;

        int N=0;

        //! jump over line 1
        std::getline(fileSettings,val);

        //! read the type of gradient
        std::getline(fileSettings,val);
        N = sscanf(val.c_str(),"%s%d",tmp,&gradient);
        cout<<tmp<<" "<<gradient<<endl;

        //! read the first colot
        std::getline(fileSettings,val);
        sscanf(val.c_str(),"%s%d%d%d",tmp,&r1,&b1,&g1);
        cout<<tmp<<" "<<r1<<" "<<g1<<" "<<b1<<endl;

        //! read the second color
        std::getline(fileSettings,val);
        sscanf(val.c_str(),"%s%d%d%d",tmp,&r2,&b2,&g2);
        cout<<tmp<<" "<<r2<<" "<<g2<<" "<<b2<<endl;

        this->setBackgroundColor(r1/255.0,g1/255.0,b1/255.0,r2/255.0,g2/255.0,b2/255.0,gradient);
    }
    else
    {
        //! settings.txt has not been opened: use a default background
        occView->SetBackgroundColor(Quantity_NOC_WHITE);
    }

    occView->MustBeResized();

    //! ------------------------------
    //! create AIS_InteractiveContext
    //! ------------------------------
    occContext = new AIS_InteractiveContext(occViewer);

    occHandle(Prs3d_LineAspect) aHiddenLineAspect = new Prs3d_LineAspect(Quantity_NOC_GRAY,Aspect_TOL_DASH,1.0);
    occContext->DefaultDrawer()->SetHiddenLineAspect(aHiddenLineAspect);

    //! --------------------
    //! the selection style
    //! --------------------
    //const occHandle(Prs3d_Drawer) selectionStyle = new Prs3d_Drawer();
    //occHandle(Prs3d_ShadingAspect) shadingAspect = new Prs3d_ShadingAspect();
    //selectionStyle->SetColor(static_cast<Quantity_NameOfColor>(Quantity_NOC_RED));
    //selectionStyle->SetShadingAspect(shadingAspect);
    //occContext->SetSelectionStyle(selectionStyle);

    //! -----------------------------------------------
    //! create an additional context for the mesh view
    //! -----------------------------------------------
    occMeshContext = new AIS_InteractiveContext(occViewer);

    //! --------------------------------------------------
    //! create an additional context for the results view
    //! --------------------------------------------------
    occPostContext = new AIS_InteractiveContext(occViewer);

    //! -------------------------------
    //! Set up lights - default lights
    //! -------------------------------
    occViewer->SetDefaultLights();
    occViewer->SetLightOn();

    //! --------------
    //! Set the triad
    //! --------------
    occView->TriedronDisplay(Aspect_TOTP_RIGHT_LOWER,Quantity_NOC_BLACK,0.1,V3d_ZBUFFER);

    //! -------------------------------------------------------------------------
    //! Create an orthographic View in this Viewer (actually this is the default)
    //! -------------------------------------------------------------------------
    occView->Camera()->SetProjectionType(Graphic3d_Camera::Projection_Orthographic);

    //! Update the visualization in this View
    occView->Redraw();

    //! The text label and the (NULL) floating (x, y, z) label
    occContext->Display(myFloatingLabel,true);

    //! a z-layer
    occViewer->AddZLayer(my_zLayer);
}

//! ----------------------
//! function: paint event
//! details:
//! ----------------------
void occGLWidget::paintEvent(QPaintEvent* e)
{
    Q_UNUSED(e);

    if (occContext.IsNull())
    {
        this->init();
    }
    Standard_CString theContent = timeStamp().toStdString().c_str();
    myTextLabel->SetText(theContent);
    occView->Redraw();
}

//! -------------------------------
//! function: write the text label
//! details:
//! -------------------------------
void occGLWidget::writeLabel(const QString& theTextLabel, int xPos, int yPos, occHandle(Graphic3d_TransformPers) &trs)
{
    Standard_CString content = theTextLabel.toStdString().c_str();
    myFloatingLabel->SetText(content);

    Graphic3d_Vec2i offset(xPos+10,yPos-10);
    trs->SetOffset2d(offset);
}

//! -----------------------
//! function: resize event
//! details:
//! -----------------------
void occGLWidget::resizeEvent(QResizeEvent* e)
{
    Q_UNUSED(e);
    if(!occView.IsNull())
    {
        occView->MustBeResized();
        occContext->RecomputePrsOnly(myTextLabel,Standard_False,Standard_False);
    }
}

//! -----------------------------------
//! function: set a gradient bakground
//! details:
//! -----------------------------------
void occGLWidget::setBackgroundColor(double R1, double G1, double B1, double R2, double G2, double B2, int tof)
{
    cout<<"occGLWidget::setBackgroundColor()->____function called____"<<endl;
    cout<<"occGLWidget::setBackgroundColor()->____1st("<<R1<<","<<G1<<","<<B1<<")_2nd("<<R2<<","<<G2<<","<<B2<<")_gradient: "<<tof<<"____"<<endl;
    Quantity_Color c1(R1,G1,B1,Quantity_TOC_RGB);
    Quantity_Color c2(R2,G2,B2,Quantity_TOC_RGB);
    if(tof==0)
    {
        cout<<"occGLWidget::setBackgroundColor()->____setting uniform background color____"<<endl;
        occView->SetBackgroundColor(c1);
    }
    else
    {
        cout<<"occGLWidget::setBackgroundColor()->____setting gradient background color____"<<endl;
        Aspect_GradientFillMethod typeOfGradient = static_cast<Aspect_GradientFillMethod>(tof);
        occView->SetBgGradientColors(c1,c2,typeOfGradient);
    }
}

//! ---------------------------
//! function: mouse well event
//! details:
//! ---------------------------
void occGLWidget::onMouseWheel(const int theFlags, const int theDelta, const QPoint thePoint)
{
    Q_UNUSED(theFlags);
    Standard_Integer aFactor=16;
    Standard_Integer aX=thePoint.x();
    Standard_Integer aY=thePoint.y();

    if (theDelta>0)
    {
        aX += aFactor;
        aY += aFactor;
    }
    else
    {
        aX -= aFactor;
        aY -= aFactor;
    }
    occView->Zoom(thePoint.x(),thePoint.y(),aX,aY);
}

//! ----------------------------
//! function: mouse press event
//! details:
//! ----------------------------
void occGLWidget::mousePressEvent(QMouseEvent *e)
{
    //cout<<"occGLWidget::mousePressEvent()->____mouse button pressed____"<<endl;
    //! if the left button has been pressed
    if (e->button() == Qt::LeftButton)
    {
        if(myCurAction3D==CurAction3D_Panning) setCursor(Qt::SizeAllCursor);
        onLButtonDown((e->buttons() | e->modifiers()), e->pos());
    }
    //! if the mid button has been pressed
    else if (e->button() == Qt::MidButton)
    {
        onMButtonDown((e->buttons() | e->modifiers()), e->pos());
    }
    //! if the right button has been pressed
    else if (e->button() == Qt::RightButton)
    {
        onRButtonDown((e->buttons() | e->modifiers()), e->pos());
    }
}

//! ----------------------------
//! function: mouse wheel event
//! details:
//! ----------------------------
void occGLWidget::wheelEvent(QWheelEvent *e)
{
    onMouseWheel(e->buttons(),e->delta(),e->pos());
}

//! ------------------------------
//! function: mouse release event
//! details:
//! ------------------------------
void occGLWidget::mouseReleaseEvent(QMouseEvent *e)
{
    if (e->button() == Qt::LeftButton)
    {
        //! ----------------------------------------------------------
        //! back from a pan operation - returns to the default cursor
        //! ----------------------------------------------------------
        if(myCurAction3D==CurAction3D_Panning ||
                myCurAction3D==CurAction3D_Rotation ||
                myCurAction3D==CurAction3D_WindowZooming)
        {
            setCursor(Qt::ArrowCursor);
        }

        //! experimental
        if(myCurAction3D==CurAction3D_PlaneDrag)
        {
            myCurAction3D = myCurAction3D_old;
        }

        onLButtonUp(e->buttons() | e->modifiers(), e->pos());
    }
    else if (e->button() == Qt::MidButton)
    {
        onMButtonUp(e->buttons() | e->modifiers(), e->pos());
    }
    else if (e->button() == Qt::RightButton)
    {
        onRButtonUp(e->buttons() | e->modifiers(), e->pos());
    }
}

//! ------------------------
//! function: onLButtonDown
//! details:
//! ------------------------
void occGLWidget::onLButtonDown(const int theFlags,const QPoint thePoint)
{
    Q_UNUSED(theFlags);

    //! --------------------------------
    //! Save the current mouse position
    //! --------------------------------
    myXmin = thePoint.x();
    myYmin = thePoint.y();
    myXmax = thePoint.x();
    myYmax = thePoint.y();

    if(myCurAction3D == CurAction3D_Rotation)
    {
        //! -------------------------------------------
        //! RotationPointType_gravity at the beginning
        //! -------------------------------------------
        this->startRotation(thePoint.x(),thePoint.y(),myRotationPointType,myCOR);
    }

    //! ------------------------
    //! Change the cursor type
    //! ------------------------
    QBitmap mask(":/cursors/mask.bmp");

    if(myCurAction3D==CurAction3D_Rotation)
    {
        QPixmap cursor(":/cursors/cursor_rotate.png",".png",Qt::AutoColor);
        cursor.setMask(mask);
        setCursor(cursor);
    }
    else if(myCurAction3D==CurAction3D_WindowZooming)
    {
        QPixmap cursor(":/cursors/icon_zoom.png",".png",Qt::AutoColor);
        cursor.setMask(mask);
        setCursor(cursor);
    }
    //! experimental
    else if(myCurAction3D==CurAction3D_PlaneDrag)
    {
        occContext->MoveTo(thePoint.x(),thePoint.y(),occView,true);
        occContext->Select(true);
        cout<<"____"<<occContext->SelectedInteractive()->get_type_name()<<"____"<<endl;
    }
}

//! -----------------------------
//! function: onMButtonDown
//! details:  mid button pressed
//! -----------------------------
void occGLWidget::onMButtonDown(const int theFlags,const QPoint thePoint)
{
    Q_UNUSED(theFlags);
    Q_UNUSED(thePoint);
}

//! -------------------------------
//! function: right button pressed
//! -------------------------------
void occGLWidget::onRButtonDown(const int theFlags,const QPoint thePoint)
{
    Q_UNUSED(theFlags);
    Q_UNUSED(thePoint);

    //! Save the current mouse position
    myXmin = thePoint.x();
    myYmin = thePoint.y();
    myXmax = thePoint.x();
    myYmax = thePoint.y();
}

//! ----------------------
//! function: onLButtonUp
//! details:
//! ----------------------
void occGLWidget::onLButtonUp(const int theFlags,const QPoint thePoint)
{
    Q_UNUSED(theFlags);

    //static occHandle(AIS_Shape) detectedShape;

    //! -----------------------------------------------------------
    //! Handling the 3D actions window zooming, rotation, panning
    //! The automatic highlighting has been previously deactivated
    //! -----------------------------------------------------------
    AIS_StatusOfPick PS;
    switch(myCurAction3D)
    {
    case(CurAction3D_WindowZooming):
        myXmax = thePoint.x();
        myYmax = thePoint.y();
        if((abs(myXmin - myXmax)>ValZWMin) || (abs(myYmin - myYmax)>ValZWMin))
        {
            occView->WindowFitAll(myXmin, myYmin, myXmax, myYmax);
        }
        break;

    case(CurAction3D_Rotation): case (CurAction3D_Panning):

        myXmax = thePoint.x();
        myYmax = thePoint.y();
        if(myLocalCtxNumber>0)
        {
            if(thePoint.x()==myXmax && thePoint.y()==myYmin)
            {
                //! One single click when the rotation mode is active.
                //! Effect: the center ("At" point) of the camera is changed
                //! If the "air" is clicked the new At point is the origin
                //! (0, 0, 0) (the default one?)
                //! If the a model point P is clicked this is the new "At point"
                occContext->MoveTo(myXmax, myYmin, occView,true);
                if(!occContext->DetectedShape().IsNull())
                {
                    //! ------------------------------------
                    //! click on a (visible face of a) body
                    //! ------------------------------------
                    gp_Pnt newCOR = hitPoint(thePoint.x(), thePoint.y(), occContext->DetectedShape());
                    emit CORchanged(newCOR, true);
                }
                else
                {
                    //! -------------------------
                    //! click on the empty space
                    //! -------------------------
                    double X_gravity, Y_gravity, Z_gravity;                    
                    X_gravity = Y_gravity = Z_gravity = 0.0;

                    gp_Pnt newCOR(X_gravity, Y_gravity, Z_gravity);
                    emit CORchanged(newCOR, false);
                }
            }
        }
        break;

    default:
        break;
    }

    //! ---------------------------------------------
    //! The global selection mode is "Single select"
    //! ---------------------------------------------
    if(myCurGlobalSelectionMode==CurGlobalSelectionMode_Single)
    {
        if(myCurSelectionMode!=CurSelection_Nothing &&
                myCurSelectionMode!=CurSelection_PointCoordinatesPicking &&
                myAllowSinglePick==Standard_False)
        {
            if(thePoint.x()==myXmin && thePoint.y()==myYmin)
            {
                if (theFlags & Qt::ControlModifier)
                {
                    occContext->MoveTo(thePoint.x(),thePoint.y(),occView,false);
                    PS=occContext->ShiftSelect(true);
                }
                else
                {
                    occContext->MoveTo(thePoint.x(),thePoint.y(),occView,false);
                    PS=occContext->Select(true);
                }
                //cout<<"occGLWidget::onLButtonUp->____STATUS OF PICK "<<PS<<"____"  <<endl;
            }
            //! Emits the selectionChanged()
            emit selectionChanged();
        }
        //! if nothing has been picked reset the permanent message
        if(PS==1) emit statusBarMessage("");
    }
    else if(myCurGlobalSelectionMode==CurGlobalSelectionMode_Multiple)
    {
        //! --------------------------------------------------
        //! The global selection mode is "multiple selection"
        //! a non zero selection area exists
        //! --------------------------------------------------
        if(abs(myXmin-myXmax)>ValZWMin || abs(myYmin-myYmax)>ValZWMin)
        {
            //! actually selects the shapes within the box
            PS=occContext->Select(myXmin, myYmin, myXmax, myYmax, occView, true);
            //! Emits selectionChanged()
            emit selectionChanged();
        }

        //! in case of one single click deselect all
        if(thePoint.x()==myXmin && thePoint.y()==myYmin)
        {
            //if(this->hasFocus())
                occContext->ClearSelected(true);
            emit statusBarMessage("");
        }
    }

    //! hide the rubber band, if present
    if (myRectBand)
    {
        myRectBand->hide();
        occView->RedrawImmediate();
    }
}

//! ----------------------
//! function: onMButtonUp
//! details:
//! ----------------------
void occGLWidget::onMButtonUp(const int theFlags,const QPoint thePoint)
{
    Q_UNUSED(thePoint);
    Q_UNUSED(theFlags);
}

//! ----------------------
//! function: onRButtonUp
//! details:
//! ----------------------
void occGLWidget::onRButtonUp(const int theFlags, const QPoint thePoint)
{
    Q_UNUSED(thePoint);
    Q_UNUSED(theFlags);
}

//! ----------------------
//! function: onMouseMove
//! details:
//! ----------------------
void occGLWidget::onMouseMove(const int theFlags, QPoint thePoint)
{
    double x_cur = thePoint.x();
    double y_cur = thePoint.y();

    occContext->MoveTo(x_cur,y_cur,occView, true);

    //! -------------
    //! experimental
    //! -------------
    //double Xp,Yp,Zp;
    //occView->Convert(x_cur,y_cur,Xp,Yp,Zp);
    //cout<<"____("<<Xp<<", "<<Yp<<", "<<Zp<<")____"<<endl;
    //! -----------------
    //! end experimental
    //! -----------------

    static TopoDS_Shape theOldDetectedShape;

    //! The 3D actions: rotation, panning, window zooming
    if (theFlags & Qt::LeftButton)
    {
        switch (myCurAction3D)
        {
        case CurAction3D_Rotation:
        {
            this->rotate(x_cur,y_cur,myRotationPointType,myCOR);
        }
            break;

        case CurAction3D_Panning:
            occView->Pan(thePoint.x()-myXmax, myYmax-y_cur);
            myXmax = x_cur;
            myYmax = y_cur;
            break;

        case CurAction3D_WindowZooming:
            myXmax = x_cur;
            myYmax = y_cur;
            drawRubberBand(myXmin, myYmin, x_cur, y_cur);
            break;

            //! -------------
            //! experimental
            //! -------------
        case CurAction3D_PlaneDrag:
        {
            int curClipPlaneID = this->getCurrentClipPlaneID();
            int delta = int(sqrt(pow(x_cur-myXmin,2)+pow(y_cur-myYmin,2)));
            cout<<"____dragging clip plane ID: "<<curClipPlaneID<<" delta: "<<delta<<"____"<<endl;
            myXmax = x_cur;
            myYmax = y_cur;

            //cesere
            //! --------------------------------
            //! retrieve the current clip plane
            //! --------------------------------
            const occHandle(Graphic3d_ClipPlane) &curClipPlane = myMapOfClipPlanes.value(curClipPlaneID);
            Graphic3d_ClipPlane::Equation planeEquation = curClipPlane->GetEquation();
            double a = *planeEquation.GetData();
            double b = *(planeEquation.GetData()+1);
            double c = *(planeEquation.GetData()+2);
            double d = *(planeEquation.GetData()+3);

            gp_Pln aPlane(a,b,c,d);
            gp_Ax1 planeAxis = aPlane.Axis();
            gp_Dir translationDirection = planeAxis.Direction();
            gp_Vec translationVector(translationDirection);
            double deltaZ = double(delta)*1000.0/1000.0;
            translationVector.Scale(deltaZ);
            aPlane.Translate(translationVector);

            aPlane.Coefficients(a,b,c,d);
            occHandle(Geom_Plane) geomPlane = new Geom_Plane(a,b,c,d);
            myMapOfHandlePlanes.value(curClipPlaneID)->SetComponent(geomPlane);
            occContext->Redisplay(myMapOfHandlePlanes.value(curClipPlaneID),true,false);
            //curClipPlane->SetEquation(aPlane);

            occView->Redraw();
        }
            break;
        }

        //! Handling current global selection modes multiple and single
        //! CurGlobalSelectionMode_Multiple -> draw a rectangle
        //! CurGlobalSelectionMode_Single -> handles "drag selection"
        switch(myCurGlobalSelectionMode)
        {
        case CurGlobalSelectionMode_Multiple:
            if(myCurAction3D == CurAction3D_Nothing)
            {
                //cout<<"doing multiple selection"<<endl;
                myXmax = x_cur;
                myYmax = y_cur;
                drawRubberBand(myXmin, myYmin, x_cur, y_cur);
            }
            break;

        case CurGlobalSelectionMode_Single:
            //!------------------------------------------------------------------------------
            //! "if" statement explanation:
            //! 1)  A local context must be open, otherwise the "occContext->DetectedShape()"
            //!     would make the code crash
            //! 2)  No action 3D must be active (no panning, no window zooming, no rotating)
            //! 3)  The current selection mode cannot be point picking, otherwise the entity
            //!     would become "selected" when the mouse is moving over them
            /// -----------------------------------------------------------------------------
            if(occContext->IndexOfCurrentLocal()>0 && myCurAction3D==CurAction3D_Nothing
                    && myCurSelectionMode!=CurSelection_PointCoordinatesPicking)
            {
                if(occContext->DetectedShape()!=theOldDetectedShape)
                {
                    //cout<<"new shape detected "<<endl;

                    //! fills a list of the previously selected shapes
                    TopoDS_ListOfShape listOfTopoDS_Shapes;
                    for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
                    {
                        listOfTopoDS_Shapes.Append(occContext->SelectedShape());
                    }
                    //cout<<"list of previously selected entities extent "<<listOfTopoDS_Shapes.Extent()<<endl;

                    TopoDS_ListIteratorOfListOfShape theIt;
                    theIt.Initialize(listOfTopoDS_Shapes);
                    for(int k=0;theIt.More();theIt.Next())
                    {
                        if(theIt.Value()!=occContext->DetectedShape())
                        {
                            k++;
                        }
                        if(k==listOfTopoDS_Shapes.Extent())  // equivalent to a multiple logical "And"
                        {
                            occContext->ShiftSelect(true);
                            emit selectionChanged();
                        }
                    }
                    theOldDetectedShape=occContext->DetectedShape();
                }
                //! This handles the case in which the selection is initially empty
                //! selecting what has been detected under the mouse pointer
                else if(occContext->NbSelected()==0)
                {
                    occContext->Select(true);
                    emit selectionChanged();
                }
            }
            break;

        default:
            break;
        }
    }
    //! -------------------------------------------------------------------
    //! If the mouse is moving and the pick point coordinate selection mode
    //! is active, write the coordinates of the picked point on the screen.
    //! Don't do it if the action3D is panning or rotating
    //! -------------------------------------------------------------------
    if(myCurSelectionMode == CurSelection_PointCoordinatesPicking &&
            myCurAction3D!=CurAction3D_Rotation && myCurAction3D!=CurAction3D_Panning
            && myCurGlobalSelectionMode!=CurGlobalSelectionMode_Multiple)
    {
        //! moreover a local context is open
        if(myLocalCtxNumber!=0)
        {
            if(!occContext->DetectedShape().IsNull())
            {
                gp_Pnt pickedPoint = hitPoint(x_cur, y_cur, occContext->DetectedShape());

                //!-------------------------------------
                //! write the coordinates on the screen
                //!-------------------------------------
                const QString &theTextLabel = QString("x = %1\ny = %2\nz = %3").
                        arg(pickedPoint.X()).arg(pickedPoint.Y()).arg(pickedPoint.Z());

                Standard_CString content = theTextLabel.toStdString().c_str();
                myFloatingLabel->SetText(content);

                //! set the persistance type
                occHandle(Graphic3d_TransformPers) trs1 = new Graphic3d_TransformPers(Graphic3d_TMF_2d);

                //! the handle point
                trs1->SetCorner2d(Aspect_TOTP_LEFT_UPPER);

                //! apply the transformation
                Graphic3d_Vec2i offset(x_cur+10,y_cur-10);
                trs1->SetOffset2d(offset);
                myFloatingLabel->SetTransformPersistence(trs1);

                //! removed since unnecessary [?]: the recompute presentation method slows down
                //! the update of the label content
                //! occContext->RecomputePrsOnly(myFloatingLabel, Standard_True, Standard_False);
            }
            else //! no face is detected -> erase the label
            {
                //! no point picked: write on the screen only the label with time and date
                myFloatingLabel->SetText("");
            }
            occContext->RecomputePrsOnly(myFloatingLabel,Standard_True, Standard_False);
        }
    }
}

//! -------------------------
//! function: mouseMoveEvent
//! details:
//! -------------------------
void occGLWidget::mouseMoveEvent(QMouseEvent *e)
{
    onMouseMove(e->buttons(),e->pos());
}

//! --------------------------------
//! function: mouseDoubleClickEvent
//! details:
//! --------------------------------
void occGLWidget::mouseDoubleClickEvent(QMouseEvent *e)
{
    Q_UNUSED(e)
    ;
}


//! ------------------------
//! function: keyPressEvent
//! details:
//! ------------------------
void occGLWidget::keyPressEvent(QKeyEvent *theKey)
{
    switch(theKey->key())
    {
    case Qt::Key_F3:
        cout<<"occGLWidget::keyPressEvent->____F3 pressed____"<<endl;
        extendSelectionToAjacent();
        break;

    case Qt::Key_F4:
        cout<<"occGLWidget::keyPressEvent->____F4 pressed____"<<endl;
        break;

    case Qt::Key_F7:
        //cout<<"occGLWidget::keyPressEvent->____F7 pressed____"<<endl;
        occView->FitAll();
        break;

    case Qt::Key_F10:
        cout<<"Key F2 pressed - highlighting a mesh face"<<endl;
        emit highlightmeshface();
        //emit highlightmeshedge();
        break;

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

         default:  // for the compiler
            break;
        }

        break;

    case Qt::Key_F8:
        //cout<<"occGLWidget::keyPressEvent->____F8 pressed____"<<endl;
        break;

    case Qt::Key_F9:
        //cout<<"occGLWidget::keyPressEvent->____F9 pressed____"<<endl;
        this->hideSelectedBodies();
        break;

    case (Qt::Key_B):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            //cout<<"occGLWidget::keyPressEvent->____Ctrl+B pressed____"<<endl;
            emit selectionModeSolid(true);
        }
        break;

    case (Qt::Key_F):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            //cout<<"occGLWidget::keyPressEvent->____Ctrl+F pressed____"<<endl;
            emit selectionModeFace(true);
        }
        break;

    case (Qt::Key_E):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            //cout<<"occGLWidget::keyPressEvent->____Ctrl+E pressed____"<<endl;
            emit selectionModeEdge(true);
        }
        break;

    case (Qt::Key_P):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            //cout<<"occGLWidget::keyPressEvent->____Ctrl+P pressed____"<<endl;
            emit selectionModeVertex(true);
        }
        break;
    case(Qt::Key_C):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            //cout<<"occGLWidget::keyPressEvent->____Ctrl+C pressed____"<<endl;
        }
        break;
    case(Qt::Key_N):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            //cout<<"occGLWidget::keyPressEvent->____Ctrl+N pressed____"<<endl;
        }
        break;
    case(Qt::Key_M):
        if(theKey->modifiers()==Qt::ControlModifier)
        {
            //cout<<"occGLWidget::keyPressEvent->____Ctrl+M pressed____"<<endl;
        }
        break;
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
            occView->Redraw();
        }
    default:
        break;
    }
}

//! -------------------------
//! function: drawRubberBand
//! details:
//! -------------------------
void occGLWidget::drawRubberBand(const int minX,const int minY,const int maxX,const int maxY)
{

    QRect aRect;

    // sets the rectangle
    (minX<maxX)?(aRect.setX(minX)):(aRect.setX(maxX));
    (minY<maxY)?(aRect.setY(minY)):(aRect.setY(maxY));
    aRect.setWidth(abs(maxX-minX));
    aRect.setHeight(abs(maxY-minY));

    if (!myRectBand)
    {
        myRectBand = new QRubberBand(QRubberBand::Rectangle, this);
        // frame in MS Windows style
        myRectBand->setStyle(QStyleFactory::create("windows"));
    }
    myRectBand->setGeometry(aRect);
    myRectBand->show();
    /*
    occHandle(AIS_RubberBand) myRB = new AIS_RubberBand(Quantity_NOC_AZURE,Aspect_TOL_SOLID,
                                                     Quantity_NOC_AZURE,0.75,1.0);

    myRB->SetRectangle(minX,minY,maxX,maxY);
    occContext->Display(myRB);
    occView->RedrawImmediate();
    */
}

//! --------------------------------
//! function: setAction3D_PlaneDrag
//! details:
//! --------------------------------
void occGLWidget::setAction3D_PlaneDrag()
{
    myCurAction3D_old = myCurAction3D;
    myCurAction3D =CurAction3D_PlaneDrag;
    switch(myCurSelectionMode)
    {
    case(CurSelection_Vertex): occContext->DeactivateStandardMode(TopAbs_VERTEX); break;
    case(CurSelection_Edge): occContext->DeactivateStandardMode(TopAbs_EDGE); break;
    case(CurSelection_Face): occContext->DeactivateStandardMode(TopAbs_FACE); break;
    case(CurSelection_Solid): occContext->DeactivateStandardMode(TopAbs_SOLID); break;
    case(CurSelection_PointCoordinatesPicking): break;
    }
}

//! --------------------------------------------------
//! function: setAction3D_Rotation
//! details:  set "rotation" as current the 3D action
//! --------------------------------------------------
void occGLWidget::setAction3D_Rotation()
{
    myCurAction3D=CurAction3D_Rotation;
    myAllowSinglePick = Standard_True;

    switch(myCurSelectionMode)
    {
    case(CurSelection_Vertex):
        occContext->DeactivateStandardMode(TopAbs_VERTEX);
        break;
    case(CurSelection_Edge):
        occContext->DeactivateStandardMode(TopAbs_EDGE);
        break;
    case(CurSelection_Face):
        // do nothing
        break;
    case(CurSelection_Solid):
        occContext->DeactivateStandardMode(TopAbs_SOLID);
        break;
    case(CurSelection_PointCoordinatesPicking):
        // do nothing
        break;
    default:
        break;
    }
    occContext->ActivateStandardMode(TopAbs_FACE);
    occContext->SetAutomaticHilight(Standard_False);
}

//! ---------------------------------------------
//! function: setAction3D_Pan
//! details:  set "pan" as current the 3D action
//! ---------------------------------------------
void occGLWidget::setAction3D_Pan()
{
    myCurAction3D=CurAction3D_Panning;
    myAllowSinglePick = Standard_True;

    switch(myCurSelectionMode)
    {
    case(CurSelection_Vertex):
        occContext->DeactivateStandardMode(TopAbs_VERTEX);
        break;
    case(CurSelection_Edge):
        occContext->DeactivateStandardMode(TopAbs_EDGE);
        break;
    case(CurSelection_Face):
        // do nothing
        break;
    case(CurSelection_Solid):
        occContext->DeactivateStandardMode(TopAbs_SOLID);
        break;
    case(CurSelection_PointCoordinatesPicking):
        // do nothing
        break;
    default:
        break;
    }
    occContext->ActivateStandardMode(TopAbs_FACE);
    occContext->SetAutomaticHilight(Standard_False);
}

//! --------------------------------------------------------
//! function: setAction3D_WindowZooming
//! details:  set "window zooming" as current the 3D action
//! --------------------------------------------------------
void occGLWidget::setAction3D_WindowZooming()
{
    myCurAction3D=CurAction3D_WindowZooming;
    myAllowSinglePick = Standard_False;

    switch(myCurSelectionMode)
    {
    case(CurSelection_Vertex):
        occContext->DeactivateStandardMode(TopAbs_VERTEX);
        break;
    case(CurSelection_Face):
        occContext->DeactivateStandardMode(TopAbs_FACE);
        break;
    case(CurSelection_Edge):
        occContext->DeactivateStandardMode(TopAbs_EDGE);
        break;
    case(CurSelection_Solid):
        occContext->DeactivateStandardMode(TopAbs_SOLID);
        break;
    case(CurSelection_PointCoordinatesPicking):
        occContext->DeactivateStandardMode(TopAbs_FACE);
        break;
    default:
        break;
    }
    occContext->SetAutomaticHilight(Standard_False);
    myCurSelectionMode=CurSelection_Nothing;
}

//! -----------------
//! function: FitAll
//! details:
//! -----------------
void occGLWidget::FitAll()
{
    myCurAction3D=CurAction3D_FitAll;
    // this is equivalent to myCurrentAction3D=curAction3D_Nothing
    myAllowSinglePick = Standard_False;
    occView->ZFitAll();
    occView->FitAll();
}

//! ------------------------
//! function: isometricView
//! details:
//! ------------------------
void occGLWidget::isometricView()
{
    occView->SetProj(V3d_XposYposZpos);
}

//! ------------------------------------
//! function: setSelectionMode
//! details:  activate a selection mode
//! ------------------------------------
#include <Aspect_TypeOfFacingModel.hxx>
#include <Graphic3d_MaterialAspect.hxx>
void occGLWidget::setSelectionMode(CurSelectionMode selectionMode)
{
    //! EXPERIMENTAL
    //! if(occContext->IndexOfCurrentLocal()<1)myLocalCtxNumber==occContext->OpenLocalContext();

    //! during the selection the 3D operations are not allowed
    myCurAction3D=CurAction3D_Nothing;

    //! After pressing a 3D view operation button the current selection
    //! mode is not changed, so the previous selection is kept active
    //! If, returning from a 3D view operation, the selection
    //! mode is changed, the previous selection is cleared
    if(selectionMode!=myCurSelectionMode)
    {
        //cout<<"occGLWidget::setSelectionMode->____CURRENT SELECTION MODE CHANGED____"<<endl;
        emptyTheSelection();
    }

    //! ---------------------------------------------
    //! activate one of the standard selection modes
    //! try and error lead to the following setting
    //! for Prs3D_TypeOgHighlight
    //! ---------------------------------------------
    Prs3d_TypeOfHighlight toh = Prs3d_TypeOfHighlight_LocalSelected;
    Aspect_TypeOfHighlightMethod atoh = Aspect_TOHM_COLOR;
    switch(selectionMode)
    {
    case CurSelection_Solid:
    {
        myCurSelectionMode = CurSelection_Solid;
        myAllowSinglePick = Standard_False;
        occContext->ActivateStandardMode(TopAbs_SOLID);
        occContext->DeactivateStandardMode(TopAbs_FACE);
        occContext->DeactivateStandardMode(TopAbs_EDGE);
        occContext->DeactivateStandardMode(TopAbs_VERTEX);
        if(!occContext->AutomaticHilight()) occContext->SetAutomaticHilight(Standard_True);

        occHandle(Prs3d_Drawer) selectionDrawer = new Prs3d_Drawer();
        selectionDrawer->SetDisplayMode(AIS_Shaded);
        selectionDrawer->SetMethod(atoh);
        selectionDrawer->SetMethod(Aspect_TOHM_COLOR);
        selectionDrawer->SetColor(static_cast<Quantity_NameOfColor>(Quantity_NOC_GREEN));
        occContext->SetSelectionStyle(selectionDrawer);
    }
        break;

    case CurSelection_Face:
    {
        myCurSelectionMode = CurSelection_Face;
        myAllowSinglePick = Standard_False;
        occContext->DeactivateStandardMode(TopAbs_SOLID);
        occContext->ActivateStandardMode(TopAbs_FACE);
        occContext->DeactivateStandardMode(TopAbs_EDGE);
        occContext->DeactivateStandardMode(TopAbs_VERTEX);
        if(!occContext->AutomaticHilight()) occContext->SetAutomaticHilight(Standard_True);

        occHandle(Prs3d_Drawer) selectionDrawer = new Prs3d_Drawer();
        selectionDrawer->SetDisplayMode(AIS_Shaded);
        selectionDrawer->SetMethod(Aspect_TOHM_COLOR);
        selectionDrawer->SetColor(static_cast<Quantity_NameOfColor>(Quantity_NOC_GREEN));

        occContext->SetHighlightStyle(Prs3d_TypeOfHighlight_LocalSelected,selectionDrawer);
        occContext->SetSelectionStyle(selectionDrawer);
    }
        break;

    case CurSelection_Edge:
    {
        myCurSelectionMode = CurSelection_Edge;
        myAllowSinglePick = Standard_False;
        occContext->DeactivateStandardMode(TopAbs_SOLID);
        occContext->DeactivateStandardMode(TopAbs_FACE);
        occContext->ActivateStandardMode(TopAbs_EDGE);
        occContext->DeactivateStandardMode(TopAbs_VERTEX);
        if(!occContext->AutomaticHilight())  occContext->SetAutomaticHilight(Standard_True);

        occHandle(Prs3d_Drawer) selectionDrawer = new Prs3d_Drawer();
        selectionDrawer->SetDisplayMode(AIS_Shaded);
        selectionDrawer->SetColor(static_cast<Quantity_NameOfColor>(Quantity_NOC_GREEN));
        occContext->SetHighlightStyle(Prs3d_TypeOfHighlight_LocalSelected,selectionDrawer);
        occContext->SetSelectionStyle(selectionDrawer);
    }
        break;

    case CurSelection_Vertex:
    {
        myCurSelectionMode = CurSelection_Vertex;
        myAllowSinglePick = Standard_False;
        occContext->DeactivateStandardMode(TopAbs_SOLID);
        occContext->DeactivateStandardMode(TopAbs_FACE);
        occContext->DeactivateStandardMode(TopAbs_EDGE);
        occContext->ActivateStandardMode(TopAbs_VERTEX);
        if(!occContext->AutomaticHilight()) occContext->SetAutomaticHilight(Standard_True);

        occHandle(Prs3d_Drawer) selectionDrawer = new Prs3d_Drawer();
        selectionDrawer->SetDisplayMode(AIS_Shaded);
        selectionDrawer->SetColor(static_cast<Quantity_NameOfColor>(Quantity_NOC_GREEN));
        occContext->SetSelectionStyle(selectionDrawer);
    }
        break;

    case CurSelection_PointCoordinatesPicking:
    {
        myCurSelectionMode = CurSelection_PointCoordinatesPicking;
        myAllowSinglePick = Standard_False;
        occContext->DeactivateStandardMode(TopAbs_SOLID);
        occContext->ActivateStandardMode(TopAbs_FACE);
        occContext->DeactivateStandardMode(TopAbs_EDGE);
        occContext->DeactivateStandardMode(TopAbs_VERTEX);
        occContext->SetAutomaticHilight(Standard_False);
    }
        break;

    default: break;
    }
}

//! ---------------------------------------------------
//! function: emptyTheSelection
//! details:  the function does not change the context
//! ---------------------------------------------------
void occGLWidget::emptyTheSelection()
{
    if(occContext->IndexOfCurrentLocal()>0)
    {
        cout<<"occGLWidget::emptyTheSelection()->____function called____"<<endl;
        // deselect the subshapes
        for(occContext->InitSelected();occContext->MoreSelected();occContext->InitSelected())
        {
            occContext->ClearSelected(Standard_False);
        }
    }
    else
    {
        cout<<"occGLWidget::emptyTheSelection()->____function called____"<<endl;
        for(occContext->InitCurrent();occContext->MoreCurrent();occContext->InitCurrent())
        {
            occContext->ClearCurrents(Standard_False);
        }
    }
    occContext->UpdateCurrentViewer();
}

//! ----------------------------------------
//! function: curAction3D
//! details:  returns the current 3D action
//! ----------------------------------------
CurAction3D occGLWidget::curAction3D()
{
    return myCurAction3D;
}

//! ---------------------------------------------
//! function: curSelectionMode
//! details:  returns the current selection mode
//! ---------------------------------------------
TopAbs_ShapeEnum occGLWidget::curSelectionMode()
{
    TopAbs_ShapeEnum value;
    switch (myCurSelectionMode)
    {
    case CurSelection_Solid: value = TopAbs_SOLID; break;
    case CurSelection_Face: value = TopAbs_FACE; break;
    case CurSelection_Edge: value = TopAbs_EDGE; break;
    case CurSelection_Vertex: value = TopAbs_VERTEX; break;
    default: break;
    }
    return value;
}

//! -----------------------------------------------------------------
//! function: sets sthe global selection mode (single or box select)
//! details: change in the internal status of the widget
//! -----------------------------------------------------------------
void occGLWidget::setGlobalCurSelectionMode(int index)
{
    switch (index)
    {
    case 0:
        myCurGlobalSelectionMode = CurGlobalSelectionMode_Single;
        emit statusBarMessage("Single select");
        break;

    case 1:
        myCurGlobalSelectionMode = CurGlobalSelectionMode_Multiple;
        emit statusBarMessage("Multiple select");
        break;

    default:
        break;
    }
}

//! ---------------------------
//! function: setWireFrameView
//! details:
//! ---------------------------
void occGLWidget::setWireframeView()
{
    //! save the previous display mode
    myCurDisplayMode_old = myCurDisplayMode;

    //! change the internal status
    myCurDisplayMode = CurDisplayMode_Wireframe;

    occHandle(Prs3d_LineAspect) la = new Prs3d_LineAspect(Quantity_NOC_BLACK,Aspect_TOL_DASH,2.0);
    occContext->DefaultDrawer()->SetLineAspect(la);
    occContext->SetDisplayMode(AIS_WireFrame,Standard_False);
    occContext->UpdateCurrentViewer();
}

//! ----------------------------------------
//! function: setShadedExteriorAndEdgesView
//! details:
//! ----------------------------------------
void occGLWidget::setShadedExteriorAndEdgesView()
{
    //if(myCurDisplayMode!=myCurDisplayMode_old)
    //{
        //! --------------------------------------------------------------
        //! save the previous display mode and change the internal status
        //! --------------------------------------------------------------
        myCurDisplayMode_old = myCurDisplayMode;
        myCurDisplayMode = CurDisplayMode_ShadedExteriorAndEdges;

        //! ------------------------------
        //! the face boundaries are shown
        //! ------------------------------
        occContext->DefaultDrawer()->SetFaceBoundaryDraw(Standard_True);
        occContext->DefaultDrawer()->FaceBoundaryAspect()->SetColor(Quantity_NOC_BLACK);
        occContext->DefaultDrawer()->FaceBoundaryAspect()->SetWidth(1.0);
        occContext->SetDisplayMode(AIS_Shaded,Standard_False);

        //! ------------------------------------------------------------
        //! list of visible AIS_Shapes: signature "0" means that only
        //! the AIS_Shape(s) having signature "Shape" are accounted for
        //! ------------------------------------------------------------
        AIS_ListOfInteractive theListOfDisplayedShapes;
        occContext->DisplayedObjects(AIS_KOI_Shape, 0, theListOfDisplayedShapes, Standard_False);

        AIS_ListIteratorOfListOfInteractive it;
        for(it.Initialize(theListOfDisplayedShapes);it.More();it.Next())
        {
            const occHandle(AIS_Shape) &AIS = occHandle(AIS_Shape)::DownCast(it.Value());
            occContext->RecomputePrsOnly(AIS,Standard_False,Standard_False);
        }
        occContext->UpdateCurrentViewer();
    //}
}

//! ------------------------------------
//! function: setShadedExteriorView
//! details:  needs updateCurrentViewer
//! ------------------------------------
void occGLWidget::setShadedExteriorView()
{
    //if(myCurDisplayMode!=myCurDisplayMode_old)
    //{
        //! --------------------------------------------------------------
        //! save the previous display mode and change the internal status
        //! --------------------------------------------------------------
        myCurDisplayMode_old = myCurDisplayMode;
        myCurDisplayMode = CurDisplayMode_ShadedExterior;

        //! the edges are not shown
        occContext->DefaultDrawer()->SetFaceBoundaryDraw(Standard_False);
        occContext->SetDisplayMode(AIS_Shaded,Standard_False);

        //! -------------------------------------------------------------------
        //! list of visible AIS_Shapes: signature "0" means that only the only
        //! the AIS_Shape(s) having signature "Shape" are accounted for
        //! -------------------------------------------------------------------
        AIS_ListOfInteractive theListOfDisplayedShapes;
        occContext->DisplayedObjects(AIS_KOI_Shape, 0, theListOfDisplayedShapes, Standard_False);

        AIS_ListIteratorOfListOfInteractive it;
        for(it.Initialize(theListOfDisplayedShapes);it.More();it.Next())
        {
            occContext->RecomputePrsOnly(it.Value(),Standard_False,Standard_False);
        }
        occContext->UpdateCurrentViewer();
    //}
}

//! ---------------------------------------------------------------- //
//! function: properties of the selection                            //
//! details: this function is called when the selection is achanged  //
//!          At the end of the computations (volume, area, ...) the  //
//!          statusBarMessage(QString) carrying the info is emitted  //
//! ---------------------------------------------------------------- //
void occGLWidget::selectionProperties()
{
    QString message;
    GProp_GProps props;

    TColStd_ListOfInteger listOfActiveModes;
    TColStd_ListIteratorOfListOfInteger listIt;

    //! loam - "l"ist "o"f "a"ctive "m"odes is a copy
    //! of "listOfActiveModes" introduced for using
    //! the method "contains()" of the class "QList"
    //! "TColStd_ListOfInteger" does not contain that method)
    QList<int> loam;

    //! Searches for the active modes
    //! First of all checks if a local context is open, otherwise
    //! "occContext->ActivateddStandardModes()" causes a crash
    if(occContext->IndexOfCurrentLocal()>0)
    {
        //! builds the list of the active modes
        listOfActiveModes = occContext->ActivatedStandardModes();
        for(listIt.Initialize(listOfActiveModes);listIt.More();listIt.Next())
        {
            // cout<<"____active mode "<<listIt.Value()<<endl;
            loam.append(listIt.Value());
        }
    }

    //! check if the selection is not empty
    occContext->InitSelected();
    if(occContext->MoreSelected())
    {
        //! 1 - vertices - TopAbs_VERTEX
        if(loam.contains(1))
        {
            //cout<<"vertex"<<endl;
            if(occContext->NbSelected()==1)   //! one single point selected - retrieve the coordinates
            {
                //cout<<"occPreGLWidget::selectionProperties->____A SINGLE POINT HAS BEEN SELECTED____"<<endl;
                occContext->InitSelected();
                const TopoDS_Shape &selectedShape = occContext->SelectedShape();
                const gp_Pnt &pnt = BRep_Tool::Pnt(TopoDS::Vertex(selectedShape));
                Standard_Real x,y,z;
                pnt.Coord(x,y,z);
                message = QString("(x, y, z) = (%1, %2, %3").arg(z).arg(y).arg(z);
                //cout<<"occPreGLWidget::selectionProperties->____X = "<<x<<" Y = "<<y<<" Z = "<<z<<"____"<<endl;
            }
            else if(occContext->NbSelected()==2) // two points selected - calculate the distance
            {
                occContext->InitSelected();
                //! retrieves the first point
                const TopoDS_Shape &selectedShape1=occContext->SelectedShape();
                const gp_Pnt &pnt1 = BRep_Tool::Pnt(TopoDS::Vertex(selectedShape1));
                //! retrieves the second point
                occContext->NextSelected();
                const TopoDS_Shape &selectedShape2=occContext->SelectedShape();
                const gp_Pnt &pnt2 = BRep_Tool::Pnt(TopoDS::Vertex(selectedShape2));
                Standard_Real distance = pnt1.Distance(pnt2);
                message = QString("Distance = %1").arg(distance);
                //cout<<"occPreGLWidget::selectionProperties->____distance: "<<distance<<"____"<<endl;
            }
            else
            {
                message = QString("%1 vertex selected").arg(occContext->NbSelected());
            }
        }
        //! 2 - edges - TopAbs_EDGE or 3 - wires - TopAbs_WIRE
        else if (loam.contains(2) || loam.contains(3))
        {
            //cout<<"edge or wire"<<endl;
            Standard_Real len;
            for(len =0.0, occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
            {
                BRepGProp::LinearProperties(occContext->SelectedShape(),props);
                len=len+props.Mass();
            }
            message = QString("%1 edges selected - total length = %2")
                    .arg(occContext->NbSelected()).arg(len);
            //cout<<"occPreGLWidget::selectionProperties->____length: "<<len<<"____"<<endl;
        }
        //! 4 - face - TopAbs_FACE
        else if (loam.contains(4))
        {
            //cout<<"face"<<endl;
            Standard_Real area;
            for(area = 0.0, occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
            {
                BRepGProp::SurfaceProperties(occContext->SelectedShape(),props);
                area=area+props.Mass();
            }
            message = QString("%1 surfaces selected - total area = %2")
                    .arg(occContext->NbSelected()).arg(area);
            //cout<<"occPreGLWidget::selectionProperties->____area: "<<area<<"____"<<endl;
        }
        //! 6 - solid - TopAbs_SOLID or 5 - shell - TopAbs_SHELL
        else if (loam.contains(5)  || loam.contains(6))
        {
            //cout<<"shell or volume"<<endl;
            if(occHandle(AIS_Shape)::DownCast(occContext->SelectedInteractive())
                ->Shape().ShapeType()==TopAbs_SOLID)
            {
                //cout<<"solid selected"<<endl;
                Standard_Real volume;
                for(volume = 0.0,occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
                {
                    BRepGProp::VolumeProperties(occContext->SelectedShape(),props);
                    volume=volume+props.Mass();
                }
                message = QString("%1 volumes selected - total volume = %3")
                        .arg(occContext->NbSelected()).arg(volume);
                //cout<<"occPreGLWidget::selectionProperties->____volume: "<<volume<<"____"<<endl;
            }
            else if(occHandle(AIS_Shape)::DownCast(occContext->SelectedInteractive())
                ->Shape().ShapeType()==TopAbs_SHELL)
            {
                //cout<<"shell"<<endl;
                Standard_Real area;
                for(area = 0.0, occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
                {
                    BRepGProp::SurfaceProperties(occContext->SelectedShape(),props);
                    area=area+props.Mass();
                }
                message = QString("%1 shells selected - total area = %2")
                        .arg(occContext->NbSelected()).arg(area);
                //cout<<"occPreGLWidget::selectionProperties->____shell area: "<<area<<"____"<<endl;
            }
        }
    }
    //! emits a system message with the geometry info of the selection
    emit statusBarMessage(message);
}

//! -----------------------------
//! function: hideSelectedBodies
//! details:
//! -----------------------------
void occGLWidget::hideSelectedBodies()
{
    //! the for cycle erases also the temporary object used for highlithing the selection
    for(occContext->InitSelected();occContext->MoreSelected();occContext->InitSelected())
    {
        //! the shape visibility flag change must be placed BEFORE the Erase() method
        (occHandle(AIS_ExtendedShape)::DownCast(occContext->SelectedInteractive()))->setShapeVisibility(Standard_False);
        occContext->Erase(occContext->SelectedInteractive(),Standard_False);
    }
    occContext->CloseLocalContext(occContext->IndexOfCurrentLocal());
    myLocalCtxNumber=occContext->OpenLocalContext();

    //! ---------------------------------------------------------------
    //! Reactivate the current selection mode (when the context
    //! is closed the selection modes and the selection list are lost)
    //! ---------------------------------------------------------------
    this->reactivateCurrentStandardSelectionMode();
}

//! ------------------------
//! function: showAllBodies
//! details:
//! ------------------------
void occGLWidget::showAllBodies()
{
    occContext->CloseLocalContext(occContext->IndexOfCurrentLocal());
    //! The DisplayAll() works in the neutral point, so the local context must be closed

    occContext->DisplayAll(false);
    myLocalCtxNumber = occContext->OpenLocalContext(true);

    //! ---------------------------------------------------------------
    //! Reactivate the current selection mode (when the context
    //! is closed the selection modes and the selection list are lost)
    //! ---------------------------------------------------------------
    this->reactivateCurrentStandardSelectionMode();

    //! ---------------------------------------------
    //! Mark all the AIS_ExtendedShape(s) as visible
    //! ---------------------------------------------
    AIS_ListOfInteractive thelistOfDIsplayed;
    Standard_Integer objectSignature = 0;
    occContext->DisplayedObjects(AIS_KOI_Shape, objectSignature, thelistOfDIsplayed, Standard_False);

    for(AIS_ListIteratorOfListOfInteractive it(thelistOfDIsplayed);it.More();it.Next())
    {
        const occHandle(AIS_ExtendedShape) &curAISShape = occHandle(AIS_ExtendedShape)::DownCast(it.Value());
        curAISShape->setShapeVisibility(Standard_True);
        curAISShape->SetTransparency(TRANSPARENCY);
    }
    occView->ZFitAll();
    occContext->UpdateCurrentViewer();

    //! If a showAllBodies() is called when an action3D is active
    //! the automatic highlight should be deactivated
    if(myCurAction3D!=CurAction3D_Nothing)occContext->SetAutomaticHilight(Standard_False);
}

//! -------------------------------------------------------------
//! function: hide all the other bodies
//! details:  among the visible shapes, hide only the unselected
//!--------------------------------------------------------------
void occGLWidget::hideAllTheOtherBodies()
{
    AIS_ListOfInteractive finalList;
    for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
    {
        const occHandle(AIS_Shape) &curAISShape = occHandle(AIS_Shape)::DownCast(occContext->SelectedInteractive());
        if(!finalList.Contains(curAISShape)) finalList.Append(curAISShape);
    }

    //! -----------------------------
    //! list of the displayed shapes
    //! -----------------------------
    AIS_ListOfInteractive theListOfDisplayedShapes;
    Standard_Integer objectSignature = 0;
    occContext->DisplayedObjects(AIS_KOI_Shape, objectSignature, theListOfDisplayedShapes, Standard_False);

    AIS_ListOfInteractive theListToBeHidden;

    for(AIS_ListIteratorOfListOfInteractive it(theListOfDisplayedShapes);it.More();it.Next())
    {
        int k = 0;
        for(AIS_ListIteratorOfListOfInteractive itFinalList(finalList);itFinalList.More();itFinalList.Next())
        {
            if(it.Value()!=itFinalList.Value()) k++;

            //! equivalent to a multiple logical AND
            if(k==finalList.Extent()) theListToBeHidden.Append(it.Value());
        }
    }

    //! ------------------------
    //! finally hide the shapes
    //! ------------------------
    for(AIS_ListIteratorOfListOfInteractive itListToBeHidden(theListToBeHidden);itListToBeHidden.More();itListToBeHidden.Next())
    {
        occHandle(AIS_ExtendedShape)::DownCast(itListToBeHidden.Value())->setShapeVisibility(Standard_False);
        occContext->Erase(itListToBeHidden.Value(),Standard_False);
    }

    occContext->CloseLocalContext(occContext->IndexOfCurrentLocal());
    myLocalCtxNumber=occContext->OpenLocalContext();
    //!occContext->UpdateCurrentViewer();

    //! Now the selection modes must be reactivated, because when the
    //! context is closed, the selection modes are lost (and the list
    //! of the selected objects is empty)
    this->reactivateCurrentStandardSelectionMode();
}

//! -------------------------------------------------
//! function: reactivateCurrentStandardSelectionMode
//! details:
//! -------------------------------------------------
void occGLWidget::reactivateCurrentStandardSelectionMode()
{
    //cout<<"occGLWidget::reactivateCurrentStandardSelectionMode()->____function called____"<<endl;
    switch(myCurSelectionMode)
    {
    case CurSelection_Solid:
        occContext->ActivateStandardMode(TopAbs_SOLID);
        occContext->SetAutomaticHilight(true);
        //setTypeOfHighlight(Standard_True);
        break;
    case CurSelection_Face:
        occContext->ActivateStandardMode(TopAbs_FACE);
        occContext->SetAutomaticHilight(true);
        //setTypeOfHighlight(Standard_False);
        break;
    case CurSelection_Edge:
        occContext->ActivateStandardMode(TopAbs_EDGE);
        occContext->SetAutomaticHilight(true);
        //setTypeOfHighlight(Standard_False);
        break;
    case CurSelection_Vertex:
        occContext->ActivateStandardMode(TopAbs_VERTEX);
        occContext->SetAutomaticHilight(true);
        //setTypeOfHighlight(Standard_False);
        break;
    default:
        break;
    }
}

//! -----------------------------------
//! function: buildMinimalContexetMenu
//! details:
//! -----------------------------------
void occGLWidget::buildMinimalContexetMenu()
{
    //! adds the "isometric view" action
    myContextMenu->addAction(viewIsometric);
    //! adds the "Zoom to Fit" action
    myContextMenu->addAction(actionZoomToFit);
    //! adds a separator
    myContextMenu->addSeparator();

    //! ---------------------------------------------------------------- //
    //! Something is selected in the viewer -> enriches the context menu //
    //! The number of displayed and hidden shapes must be calculated     //
    //! before building the submenu                                      //
    //! ---------------------------------------------------------------- //
    AIS_ListOfInteractive theListOfDisplayedShapes, theListOfHiddenShapes;

    Standard_Integer objectSignature = -1;
    occContext->DisplayedObjects(AIS_KOI_Shape, objectSignature, theListOfDisplayedShapes, Standard_False);
    occContext->ErasedObjects(AIS_KOI_Shape, objectSignature, theListOfHiddenShapes);

    if(theListOfDisplayedShapes.Extent()==0)  // all the bodies are hidden
    {
        myContextMenu->addAction(actionShowAllBodies);
        myContextMenu->addSeparator();
    }
    else if(theListOfHiddenShapes.Extent()>=1  && theListOfDisplayedShapes.Extent()>=1)  // at least one body is shown
    {
        myContextMenu->addAction(actionShowAllBodies);
        //! if something is selected gives the possibility to hide it, or to invert the visibility
        if(occContext->NbSelected()!=0)
        {
            myContextMenu->addAction(actionHideBody);
            myContextMenu->addAction(actionHideAllOtherBodies);
        }
        myContextMenu->addSeparator();
    }
    else if(theListOfHiddenShapes.Extent()==0)  // all the bodies are shown
    {
        //! if something is selected, give the possibility to hide it, or to invert the selection
        if(occContext->NbSelected()!=0)
        {
            myContextMenu->addAction(actionHideBody);
            myContextMenu->addAction(actionHideAllOtherBodies);
        }
    }

    //! --------------------------------
    //! setup of the "Cursor Mode" menu
    //! --------------------------------
    QList <QAction *> listOfSelectionModes;
    listOfSelectionModes.append(selectVertex);
    listOfSelectionModes.append(selectEdge);
    listOfSelectionModes.append(selectFace);
    listOfSelectionModes.append(selectSolid);

    //! If something is displayed the "pick point coordinates" row
    //! is added to the context sub-menu

    AIS_ListOfInteractive listOfDisplayed;

    //! the signature criterion is switched off
    occContext->DisplayedObjects(AIS_KOI_Shape, -1, listOfDisplayed, Standard_False);

    if(listOfDisplayed.Extent()>=1)
    {
        listOfSelectionModes.append(selectPickPointCoordinates);
    }
    myCursorModeMenu->addActions(listOfSelectionModes);

    //! If a selection mode has been activated through the MainWindow
    //! toolbar, here, in the context (sub)menu, the corresponding button
    //! should be in the "checked" status
    switch(myCurSelectionMode)
    {
    case CurSelection_Vertex:
        selectVertex->setChecked(true);
        break;
    case CurSelection_Edge:
        selectEdge->setChecked(true);
        break;
    case CurSelection_Face:
        selectFace->setChecked(true);
        break;
    case CurSelection_Solid:
        selectSolid->setChecked(true);
        break;
    case CurSelection_PointCoordinatesPicking:
        selectPickPointCoordinates->setChecked(true);
        break;
    default:
        break;
    }

    //! However, if a 3D action is active,
    //! the previous buttons must be unchecked
    if(myCurAction3D!=CurAction3D_Nothing)
    {
        selectVertex->setChecked(false);
        selectEdge->setChecked(false);
        selectFace->setChecked(false);
        selectSolid->setChecked(false);
        selectPickPointCoordinates->setChecked(false);
    }
    myContextMenu->addMenu(myCursorModeMenu);

    //! -----------------------------
    //! setup of the "View" sub-menu
    //! -----------------------------
    myViewModeMenu->setIcon(QIcon(":/icons/icon_view.png"));

    QList<QAction *> listOfActions;
    listOfActions.append(viewFront);
    listOfActions.append(viewBack);
    listOfActions.append(viewRight);
    listOfActions.append(viewLeft);
    listOfActions.append(viewTop);
    listOfActions.append(viewBottom);
    myViewModeMenu->addActions(listOfActions);

    //! add the "View" sub-menu
    myContextMenu->addMenu(myViewModeMenu);

    //! add a separator
    myContextMenu->addSeparator();

    //! add "Select all" action
    myContextMenu->addAction(actionSelectAll);
}

//! --------------------------
//! function: ShowContextMenu
//! details:
//! --------------------------
void occGLWidget::ShowContextMenu(const QPoint& pos)
{
    QPoint globalPos = this->mapToGlobal(pos);
    //! since the context menu is defined outside this function
    //! so it must be re-initialized

    //! adds the "isometric view" action
    myContextMenu->addAction(viewIsometric);
    //! adds the "Zoom to Fit" action
    myContextMenu->addAction(actionZoomToFit);
    //! adds a separator
    myContextMenu->addSeparator();

    //! ---------------------------------------------------------------- //
    //! Something is selected in the viewer -> enriches the context menu //
    //! The number of displayed and hidden shapes must be calculated     //
    //! before building the submenu                                      //
    //! ---------------------------------------------------------------- //
    AIS_ListOfInteractive theListOfDisplayedShapes, theListOfHiddenShapes;

    Standard_Integer objectSignature = -1;
    occContext->DisplayedObjects(AIS_KOI_Shape, objectSignature, theListOfDisplayedShapes, Standard_False);
    occContext->ErasedObjects(AIS_KOI_Shape, objectSignature, theListOfHiddenShapes);

    if(theListOfDisplayedShapes.Extent()==0)  //! all the bodies are hidden
    {
        myContextMenu->addAction(actionShowAllBodies);
        myContextMenu->addSeparator();
    }
    else if(theListOfHiddenShapes.Extent()>=1  && theListOfDisplayedShapes.Extent()>=1)  // at least one body is shown
    {
        myContextMenu->addAction(actionShowAllBodies);
        //! if something is selected gives the possibility to hide it, or to invert the visibility
        if(occContext->NbSelected()!=0)
        {
            myContextMenu->addAction(actionHideBody);
            myContextMenu->addAction(actionHideAllOtherBodies);
        }
        myContextMenu->addSeparator();
    }
    else if(theListOfHiddenShapes.Extent()==0)  //! all the bodies are shown
    {
        //! if something is selected gives the possibility to hide it, or to invert the selection
        if(occContext->NbSelected()!=0)
        {
            myContextMenu->addAction(actionHideBody);
            myContextMenu->addAction(actionHideAllOtherBodies);
        }
    }

    //! --------------------------------
    //! setup of the "Cursor Mode" menu
    //! --------------------------------
    QList <QAction *> listOfSelectionModes;
    listOfSelectionModes.append(selectVertex);
    listOfSelectionModes.append(selectEdge);
    listOfSelectionModes.append(selectFace);
    listOfSelectionModes.append(selectSolid);

    //! If something is displayed the "pick point coordinates" row
    //! is added to the context sub-menu

    AIS_ListOfInteractive listOfDisplayed;

    //! the signature criterion is switched off
    Standard_Integer objectSignature1 = -1;
    occContext->DisplayedObjects(AIS_KOI_Shape, objectSignature1, listOfDisplayed, Standard_False);

    if(listOfDisplayed.Extent()>=1)
    {
        listOfSelectionModes.append(selectPickPointCoordinates);
    }

    myCursorModeMenu->addActions(listOfSelectionModes);

    //! If a selection mode has been activated through the MainWindow
    //! toolbar, here, in the context (sub)menu, the corresponding button
    //! should be in the "checked" status
    switch(myCurSelectionMode)
    {
    case CurSelection_Vertex:
        selectVertex->setChecked(true);
        break;
    case CurSelection_Edge:
        selectEdge->setChecked(true);
        break;
    case CurSelection_Face:
        selectFace->setChecked(true);
        break;
    case CurSelection_Solid:
        selectSolid->setChecked(true);
        break;
    case CurSelection_PointCoordinatesPicking:
        selectPickPointCoordinates->setChecked(true);
        break;
    default:  // for the compiler
        break;
    }

    //! However, if a 3D action is active,
    //! the previous buttons must be unchecked
    if(myCurAction3D!=CurAction3D_Nothing)
    {
        selectVertex->setChecked(false);
        selectEdge->setChecked(false);
        selectFace->setChecked(false);
        selectSolid->setChecked(false);
        selectPickPointCoordinates->setChecked(false);
    }

    myContextMenu->addMenu(myCursorModeMenu);

    //! -----------------------------
    //! setup of the "View" sub-menu
    //! -----------------------------
    myViewModeMenu->setIcon(QIcon(":/icons/icon_view.png"));

    QList<QAction *> listOfActions;
    listOfActions.append(viewFront);
    listOfActions.append(viewBack);
    listOfActions.append(viewRight);
    listOfActions.append(viewLeft);
    listOfActions.append(viewTop);
    listOfActions.append(viewBottom);
    myViewModeMenu->addActions(listOfActions);

    //! add the "View" sub-menu
    myContextMenu->addMenu(myViewModeMenu);

    //! -----------------------------
    //! handling the selected action
    //! -----------------------------
    QAction* selectedItem = myContextMenu->exec(globalPos);

    //! check if the selection is valid
    if(selectedItem)
    {
        switch(selectedItem->data().toInt())
        {
        //! For cases {0, 1, 2, 3, 4} which handle the change of the selection mode
        //! a signal is emitted: it is received by the MainWindow class, which in turn
        //! calls the slots toggle<..>SelectionMode(), which in turn calls the slot
        //! slot setSelectionMode(). Passing through the MainWindow is for
        //! handling the status (checked/unchecked) of the buttons in the toolbar
        case 0:
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____SELECT VERTEX____"<<endl;
            emit selectionModeVertex(true);
            break;
        case 1:
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____SELECT EDGES____"<<endl;
            emit selectionModeEdge(true);
            break;
        case 2:
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____SELECT FACES____"<<endl;
            emit selectionModeFace(true);
            break;
        case 3:
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____SELECT SOLIDS____"<<endl;
            emit selectionModeSolid(true);
            break;
        case 4:
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____PICK POINT COORDINATES____"<<endl;
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
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____ISOMETRIC VIEW____"<<endl;
            this->isometricView();
            break;
        case 12:
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____ZOOM TO FIT (FIT ALL)____"<<endl;
            //! please, do not use this->FitAll() because this internally changes the "myCurrentActio3D" value,
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
             default:  // for the compiler
                break;
            }
            break;
        case 13:
            cout<<"occGLWidget::ShowContextMenu->____context menu request____Hide body____"<<endl;
            emit requestSynchVisibility(false);
            hideSelectedBodies();
            break;
        case 14:
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____Hide all other bodies___"<<endl;
            hideAllTheOtherBodies();
            break;
        case 15:
            //cout<<"occGLWidget::ShowContextMenu->____context menu request____Show all bodies____"<<endl;
            showAllBodies();
            break;
        case 16:
            cout<<"occGLWidget::ShowContextMenu->____context menu request____Select all____"<<endl;
            this->selectAll();
            break;
        default:
            break;
        }
    }
}

//! ----------------------------------------------------------------
//! function: hitPoint
//! details:  return the position of a point picked on the 3D model
//! ----------------------------------------------------------------
gp_Pnt occGLWidget::hitPoint(long x, long y, TopoDS_Shape shape)
{
    gp_Pnt resultPoint;

    double xEye, yEye, zEye, xAt, yAt, zAt;
    occView->Eye(xEye, yEye, zEye);
    occView->At(xAt, yAt, zAt);

    gp_Pnt EyePoint(xEye, yEye, zEye);
    gp_Pnt AtPoint(xAt, yAt, zAt);

    gp_Vec EyeVector(EyePoint, AtPoint);
    gp_Dir EyeDir(EyeVector);

    gp_Pln PlaneOfView = gp_Pln(AtPoint, EyeDir);

    Standard_Real theX, theY, theZ;
    occView->Convert(x, y, theX, theY, theZ);
    gp_Pnt ConvertedPoint (theX, theY, theZ);

    gp_Pnt2d ConvertedPointOnPlane = ProjLib::Project(PlaneOfView, ConvertedPoint);
    gp_Pnt shapePoint = ElSLib::Value(ConvertedPointOnPlane.X(), ConvertedPointOnPlane.Y(), PlaneOfView);
    resultPoint = shapePoint;   // initializes with a very far point from the camera

    GC_MakeLine line(shapePoint, EyeDir);

    TopExp_Explorer exp;
    TopAbs_State aState;

    for (exp.Init(shape,TopAbs_FACE);exp.More();exp.Next())
    {
        TopoDS_Face face = TopoDS::Face(exp.Current());
        BRepAdaptor_Surface surface(face);

        const GeomAdaptor_Surface& geomAdapSurf = surface.Surface();
        const occHandle(Geom_Surface)& geomSurf = geomAdapSurf.Surface();

        GeomAPI_IntCS inCS;

        const occHandle(Geom_Curve) &aLineCurve = line.Value();
        inCS.Perform(aLineCurve, geomSurf);
        if (inCS.IsDone())
        {
            if (inCS.NbPoints()!=0)
            {
                shapePoint = gp_Pnt(inCS.Point(1).XYZ());
                ShapeAnalysis_Surface shapeAnalysis(geomSurf);
                gp_Pnt2d shapePoint2D = shapeAnalysis.ValueOfUV(shapePoint, Precision::Confusion());
                BRepClass_FaceClassifier aClassifier(face, shapePoint2D, Precision::Confusion());
                aState = aClassifier.State();
                if((aState==TopAbs_ON) || (aState==TopAbs_IN))
                {
                    if(resultPoint.Distance(EyePoint)>shapePoint.Distance(shapePoint))
                    {
                        resultPoint = shapePoint;
                    }
                }
            }
        }
    }
    //cout<<"(x, y, z) = ("<<resultPoint.X()<<", "<<resultPoint.Y()<<" ,"<<resultPoint.Z()<<")"<<endl;
    return resultPoint;
}

//!---------------------------------------------------------------------//
//! function: updateCOR                                                 //
//! details:  update the definition of the center of rotation           //
//!           update also the view: if a point has been picked up on    //
//!           a model face, set display the model at the center of      //
//!           the viewer. If a point on the air has been picked, set    //
//!           both the origin and both the rotation center to a default //
//!           value (for the moment O(0,0,0))                           //
//! --------------------------------------------------------------------//
void occGLWidget::updateCOR(gp_Pnt newCOR, bool isOnFace)
{
    //cout<<"occGLWidget::updateCOR()->____function called____"<<endl;
    static bool isTheLastClickOnAir = false;

    //! ------------------------------
    //! update the center of rotation
    //! ------------------------------
    myCOR = newCOR;

    if(isOnFace)
    {
        //cout<<"occGLWidget::updateCOR()->____new COR("<<newCOR.X()<<", "<<newCOR.Y()<<", "<<newCOR.Z()<<")____"<<endl;
        gp_Dir cameraDirection = occView->Camera()->Direction();

        //! --------------------------
        //! updates the camera target
        //! --------------------------
        occView->Camera()->SetCenter(newCOR);
        occView->Camera()->SetDirection(cameraDirection);
        occView->Redraw();

        isTheLastClickOnAir = false;
        this->displayRotationCenter(myCOR, true);

        myRotationPointType = RotationPointType_selected;
    }
    else
    {
        //cout<<"occGLWidget::updateCOR()->____no point hit: the COV is on air____"<<endl;
        //cout<<"occGLWidget::updateCOR()->____new COR("<<newCOR.X()<<", "<<newCOR.Y()<<", "<<newCOR.Z()<<")____"<<endl;

        //! ------------------------------------------------------
        //! If the also the previous click was on air, do nothing
        //! Otherwise set the default camera "At" point
        //! ------------------------------------------------------
        gp_Dir cDirection;
        cDirection = occView->Camera()->Direction();

        //! ---------------------------------------------
        //! Return back to the previous rotation center
        //! ---------------------------------------------
        occView->Camera()->SetCenter(myCOR);
        occView->Camera()->SetDirection(cDirection);
        occView->Redraw();
        this->displayRotationCenter(myCOR, false);

        myRotationPointType = RotationPointType_gravity;
    }
}

//! --------------------------------------------------------
//! function: unsetSelectionModes
//! details:  needed when initializing (or re-initializing)
//! --------------------------------------------------------
void occGLWidget::unsetSelectionModes()
{
    myCurGlobalSelectionMode = CurGlobalSelectionMode_Single;
    myCurSelectionMode = CurSelection_Nothing;

    occContext->DeactivateStandardMode(TopAbs_SOLID);
    occContext->DeactivateStandardMode(TopAbs_FACE);
    occContext->DeactivateStandardMode(TopAbs_EDGE);
    occContext->DeactivateStandardMode(TopAbs_VERTEX);
}

//! --------------------------------------------------------
//! function: linearDeviation
//! details:  needed when initializing (or re-initializing)
//! --------------------------------------------------------
void occGLWidget::unsetViewOperations()
{
    myCurAction3D=CurAction3D_Nothing;
}

//! ----------------------------------
//! function: set display quality
//! details:  re-tessellate the shape
//! ----------------------------------
void occGLWidget::setDisplayQuality(displayQuality DQ)
{
    //! ----------------------------
    //! changes the internal status
    //! ----------------------------
    myDisplayQuality = DQ;

    AIS_ListOfInteractive listOfAIS;
    AIS_ListIteratorOfListOfInteractive anIter;

    //! ------------------------------------
    //! signature "0" is only for AIS_Shape
    //! ------------------------------------
    occContext->ObjectsInside(listOfAIS,AIS_KOI_Shape,0);

    //! ----------------------------------------
    //! retrieve the current meshing parameters
    //! ----------------------------------------
    double deviationAngle = occContext->DefaultDrawer()->DeviationAngle();
    double linearDeviation = occContext->DefaultDrawer()->DeviationCoefficient();

    double ad,ld;

    //! ------------------------
    //! change the tassellation
    //! ------------------------
    switch(myDisplayQuality)
    {
    case displayQuality_Low: { ad = 0.5; ld = 0.5; } break;
    case displayQuality_Medium: { ad = deviationAngle; ld = linearDeviation; } break;
    case displayQuality_High: { ad = 0.1; ld = 0.1; } break;
    }

    for(anIter.Initialize(listOfAIS); anIter.More(); anIter.Next())
    {
        occHandle(AIS_Shape) AISshape = occHandle(AIS_Shape)::DownCast(anIter.Value());
        TopoDS_Shape aShape = AISshape->Shape();
        BRepTools::Clean(aShape);
        BRepMesh_IncrementalMesh aMesher(aShape,ld,true,ad,true,false);
        occContext->RecomputePrsOnly(AISshape,false,false);
        Q_UNUSED(aMesher)
    }
    occContext->UpdateCurrentViewer();
}

//! -----------------------------------------------
//! function: extend to adjacent - one body - slot
//! details:
//! -----------------------------------------------

//! --------------------------------------------------------
//! function: getNormal
//! details:  normal at the center of the (U, V) projection
//! --------------------------------------------------------
static void getNormal(const TopoDS_Face& aFace, gp_Dir& aDNS)
{
    gp_Pnt aPoint;
    gp_Vec aD1U, aD1V;
    const occHandle(Geom_Surface) &aSurface = BRep_Tool::Surface(aFace);
    Standard_Real U1, U2, V1, V2;
    aSurface->Bounds(U1, U2, V1, V2);
    Standard_Real U = 0.5*(U1+U2);
    Standard_Real V = 0.5*(V1+V2);
    //! computes the first derivative (gradient)
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

//! -----------------------------------
//! function: extendSelectionToAjacent
//! details:  private slot
//! -----------------------------------
void occGLWidget::extendSelectionToAjacent()
{
    cout<<"occGLWidget::extendSelectionToAjacent()->____function called____"<<endl;

    TopTools_ListOfShape selection;

    //! ------------------------------------------------
    //! put the selected interactive shapes into a list
    //! ------------------------------------------------
    AIS_ListOfInteractive list;
    for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
    {
        list.Append(occContext->SelectedInteractive());
    }

    AIS_ListIteratorOfListOfInteractive it;
    for(it.Initialize(list);it.More();it.Next())
    {
        const occHandle(AIS_Shape) &theAISShape = occHandle(AIS_Shape)::DownCast(it.Value());
        const TopoDS_Shape &theShape = theAISShape->Shape();

        TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
        edgeFaceMap.Clear();
        TopExp::MapShapesAndAncestors(theShape,TopAbs_EDGE,TopAbs_FACE,edgeFaceMap);

        if(occContext->NbSelected()==1)
        {
            occContext->InitSelected();
            if(occContext->MoreSelected())
            {
                //! the current selected face
                TopoDS_Face theSelectedFace = TopoDS::Face(occContext->SelectedShape());
                selection.Append(theSelectedFace);

                //! scans the edges of the current face
                TopExp_Explorer edgeExp(theSelectedFace,TopAbs_EDGE);
                for(int k=1;edgeExp.More();edgeExp.Next(),k++)
                {
                    //! the current edge
                    TopoDS_Edge theCurEdge= TopoDS::Edge(edgeExp.Current());

                    //! the normal to the selected face calculated on the current edge
                    gp_Vec NSF;
                    TopOpeBRepBuild_Tools::GetNormalToFaceOnEdge(theSelectedFace,theCurEdge, NSF);

                    //! find the adjacent face
                    TopoDS_Shape theAdjacentFace;

                    //! get the adjacent face
                    bool faceFound = TopOpeBRepBuild_Tools::GetAdjacentFace(theSelectedFace,theCurEdge,edgeFaceMap,theAdjacentFace);

                    if(faceFound)
                    {
                        gp_Vec NAF;
                        TopOpeBRepBuild_Tools::GetNormalToFaceOnEdge(TopoDS::Face(theAdjacentFace),theCurEdge, NAF);
                        Standard_Real angle = NAF.Angle(NSF)*180/M_PI;

                        //! check normals
                        bool inRange = (angle<=limitAngle || angle>=(180-limitAngle));
                        if(inRange)
                        {
                            //!cout<<"Deviation angle: "<<theDeviationAngle<<" -> Face added"<<endl;
                            occContext->AddOrRemoveSelected(theAdjacentFace,true);
                            selection.Append(theAdjacentFace);
                        }
                    }
                }
            }
            cout<<"____number of selected entities: "<<selection.Extent()<<"____"<<endl;
        }
        // --- l'errore è qui --- non si possono più usare le funzioni del contesto

        else if (occContext->NbSelected()>1)
        {
            BRep_Builder theBuilder;
            TopoDS_Compound theFaceCompound;
            theBuilder.MakeCompound(theFaceCompound);

            for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
            {
                theBuilder.Add(theFaceCompound,occContext->SelectedShape());
            }

            TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMapOfFaceCompound;
            edgeFaceMapOfFaceCompound.Clear();

            //! This map is used for obtaining the free edges of a multiple face selection
            //! The mentioned faces can be connected (shell) or not (this is why I use a
            //! TopoDS_Compound for joining the selection of the faces
            TopExp::MapShapesAndAncestors(theFaceCompound,TopAbs_EDGE,TopAbs_FACE,edgeFaceMapOfFaceCompound);

            //cout<<"______error here______"<<endl;

            //! free boundaries of the face compound
            TopExp_Explorer edgeExp(theFaceCompound,TopAbs_EDGE);
            TopTools_ListOfShape freeEdgeList;
            int k1=0;
            for(;edgeExp.More();edgeExp.Next())
            {
                if(edgeFaceMapOfFaceCompound.FindFromKey(edgeExp.Current()).Extent()==1)
                {
                    freeEdgeList.Append(edgeExp.Current());
                    k1++;
                }
            }

            //cout<<"Number of free edges of the faces compound: "<<k1<<endl;

            int k=0;
            TopTools_ListOfShape list;
            for(occContext->InitSelected();occContext->MoreSelected();occContext->NextSelected())
            {
                //! the current selected face
                TopoDS_Face theCurSelectedFace = TopoDS::Face(occContext->SelectedShape());

                TopExp_Explorer edgeExp(occContext->SelectedShape(),TopAbs_EDGE);
                for(;edgeExp.More();edgeExp.Next())
                {
                    const TopoDS_Edge &theCurEdge = TopoDS::Edge(edgeExp.Current());
                    if(freeEdgeList.Contains(edgeExp.Current()))
                    {
                        k++;
                        //! find the adjacent face
                        TopoDS_Shape theAdjacentFace;
                        bool faceFound = TopOpeBRepBuild_Tools::GetAdjacentFace(occContext->SelectedShape(),edgeExp.Current(),edgeFaceMap,theAdjacentFace);

                        if(faceFound)
                        {
                            //cout<<"adjacent face found"<<endl;
                            //! the normal to the selected face calculated at the common edge
                            gp_Vec NSF;
                            TopOpeBRepBuild_Tools::GetNormalToFaceOnEdge(theCurSelectedFace,theCurEdge, NSF);
                            //! the normal to the adjacent face calculated at the common edge
                            gp_Vec NAF;
                            TopOpeBRepBuild_Tools::GetNormalToFaceOnEdge(TopoDS::Face(theAdjacentFace),theCurEdge, NAF);
                            Standard_Real theDeviationAngle = NAF.Angle(NSF)*180/M_PI;
                            //! check normals
                            bool inRange = (theDeviationAngle<=limitAngle || theDeviationAngle>=(180-limitAngle))? true: false;
                            if(inRange)
                            {
                                //cout<<"Deviation angle: "<<theDeviationAngle2<<" -> Face added"<<endl;
                                list.Append(theAdjacentFace);
                            }
                            else
                            {
                                //cout<<"Deviation angle: "<<theDeviationAngle2<<" -> Face rejected"<<endl;
                            }
                        }
                        else
                        {
                            //cout<<"adjacent face not found"<<endl;
                        }
                    }
                }
            }
            //!------------------------------------------------------------------

            if(list.Extent()>0)
            {
                //cout<<"list extent: "<<list.Extent()<<endl;

                //! elimination of the duplicated adjacent faces in the list
                TopTools_ListIteratorOfListOfShape listIt;
                TopTools_ListOfShape finalList;
                TopTools_ListIteratorOfListOfShape itFinalList;
                finalList.Append(list.First());
                listIt.Initialize(list);
                for(;listIt.More();listIt.Next())
                {
                    const TopoDS_Shape &outer = listIt.Value();
                    int k = 0;
                    itFinalList.Initialize(finalList);
                    for(; itFinalList.More(); itFinalList.Next())
                    {
                        const TopoDS_Shape &inner = itFinalList.Value();
                        if(outer!= inner)
                        {
                            k++;
                        }
                        //! equivalent to a multiple logical AND
                        if(k==finalList.Extent())finalList.Append(outer);
                    }
                }
                //cout<<"final list extent: "<<finalList.Extent()<<endl;
                //! finally add the faces to the selection
                for(itFinalList.Init(finalList);itFinalList.More();itFinalList.Next())
                {
                    occContext->AddOrRemoveSelected(itFinalList.Value(),true);
                }
            }
            //!------------------------------------------------------------------
        }

    }
    cout<<"Number of selected entities: "<<occContext->NbSelected()<<endl;
    emit selectionChanged();
}

//! ------------------------------------
//! function: get the current action 3D
//! details:
//! ------------------------------------
CurAction3D occGLWidget::currentAction3D()
{
    return myCurAction3D;
}

//! --------------------------
//! function: get the context
//! details:  access function
//! --------------------------
const occHandle(AIS_InteractiveContext)& occGLWidget::getContext() const
{
    return occContext;
}

//! ------------------
//! function: getView
//! details:
//! ------------------
const occHandle(V3d_View)& occGLWidget::getView() const
{
    return occView;
}

//! --------------------------------
//! function - slot: clear selected
//! details:
//! --------------------------------
void occGLWidget::clearSelected()
{
    occContext->ClearSelected(Standard_True);
}

//! ---------------------------
//! function: displayTrihedron
//! details:  slot
//! ---------------------------
void occGLWidget::displayTrihedron(QVector<double> Origin, QVector<QVector<double>> directionalData, int axisLength)
{
    cout<<"occGLWidget::displayTrihedron()->____function called____"<<endl;
    this->removeTrihedron();

    gp_Pnt P(Origin.at(0),Origin.at(1),Origin.at(2));
    QVector<double> xAxisData = directionalData.at(0);
    QVector<double> zAxisData = directionalData.at(2);

    gp_Dir N(zAxisData.at(0),zAxisData.at(1),zAxisData.at(2));
    gp_Dir V(xAxisData.at(0),xAxisData.at(1),xAxisData.at(2));
    gp_Ax2 theCoordinateSystem(P,N,V);

    occHandle(Geom_Axis2Placement) thePlacement =  new Geom_Axis2Placement(theCoordinateSystem);
    occHandle(AIS_Trihedron) theTrihedron = new AIS_Trihedron(thePlacement);
    theTrihedron->SetTextColor(Quantity_NOC_BLACK);

    if(axisLength==0)
    {
        double X,Y;
        occView->Size(X,Y);
        axisLength = int(sqrt(X*Y)/TRIHEDRON_AXIS_LENGTH_DIAGONAL_FACTOR)+1;
    }

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

    //! --------
    //! z-layer
    //! --------
    occContext->SetZLayer(theTrihedron,Graphic3d_ZLayerId_TopOSD);

    //! --------
    //! display
    //! --------
    occContext->Display(theTrihedron,0,-1,true);
}

//! --------------------------
//! function: removeTrihedron
//! details:
//! --------------------------
void occGLWidget::removeTrihedron()
{
    occContext->CloseLocalContext();

    //! remove the previous, if present
    AIS_ListOfInteractive aList;
    occContext->ObjectsInside(aList,AIS_KOI_Datum,-1);
    AIS_ListIteratorOfListOfInteractive anIter;
    for(anIter.Initialize(aList);anIter.More();anIter.Next()) occContext->Remove(anIter.Value(),false);
    occContext->UpdateCurrentViewer();

    myLocalCtxNumber = occContext->OpenLocalContext();
    this->reactivateCurrentStandardSelectionMode();
}

//! --------------------------
//! function: saveImageToDisk
//! details:  slot
//! --------------------------
void occGLWidget::saveImageToDisk()
{
    //cout<<"occGLWidget()::saveImageToDisk()->____function called____"<<endl;
    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(this,"Save picture as ",tools::getWorkingDir(),
                                                    "*.png;;*.jpg;;*.bmp;;*.gif",&selectedFilter,0);
    if(!fileName.isEmpty())
    {
        //cout<<"occGLWidget()::saveImageToDisk()->____"<<fileName.toStdString()<<"____"<<endl;
        occView->Dump(fileName.toStdString().c_str());
    }
}

//! -------------------------------------------------------------
//! function: selectAll
//! details:  select all according to the current selection mode
//!           it works only on visible shapes
//! -------------------------------------------------------------
void occGLWidget::selectAll()
{
    cout<<"occGLWidget::selectAll()->____function called____"<<endl;
    AIS_ListOfInteractive theListOfDisplayedAIS_Objects;
    AIS_ListIteratorOfListOfInteractive it;

    //! retrieve the displayed AIS_Shape(s)
    //! the signature criterion is not accounted for (-1)
    occContext->DisplayedObjects(AIS_KOI_Shape, -1, theListOfDisplayedAIS_Objects, Standard_False);

    switch(myCurSelectionMode)
    {
    case CurSelection_Solid:
    {
        for(it.Initialize(theListOfDisplayedAIS_Objects); it.More(); it.Next())
        {
            const occHandle(AIS_Shape) &theCurAIS = occHandle(AIS_Shape)::DownCast(it.Value());
            const TopoDS_Shape &theShape = theCurAIS->Shape();
            occContext->AddOrRemoveSelected(theShape,Standard_False);
        }
        occContext->UpdateCurrentViewer();
    }
        break;

    case CurSelection_Face:
    {
        for(it.Initialize(theListOfDisplayedAIS_Objects); it.More(); it.Next())
        {
            const occHandle(AIS_Shape) &theCurAIS = occHandle(AIS_Shape)::DownCast(it.Value());
            const TopoDS_Shape &theShape = theCurAIS->Shape();
            TopExp_Explorer anExp;
            //! add the faces of the current shape
            for(anExp.Init(theShape,TopAbs_FACE);anExp.More();anExp.Next())
            {
                occContext->AddOrRemoveSelected(anExp.Current(),Standard_False);
            }
        }
        occContext->UpdateCurrentViewer();
    }
        break;

    case CurSelection_Edge:
    {
        for(it.Initialize(theListOfDisplayedAIS_Objects); it.More(); it.Next())
        {
            const occHandle(AIS_Shape) &theCurAIS = occHandle(AIS_Shape)::DownCast(it.Value());
            const TopoDS_Shape &theShape = theCurAIS->Shape();
            TopExp_Explorer anExp;
            //! add the faces of the current shape
            for(anExp.Init(theShape,TopAbs_EDGE);anExp.More();anExp.Next())
            {
                occContext->AddOrRemoveSelected(anExp.Current(),Standard_False);
            }
        }
        occContext->UpdateCurrentViewer();
    }
        break;

    case CurSelection_Vertex:
    {
        for(it.Initialize(theListOfDisplayedAIS_Objects); it.More(); it.Next())
        {
            const occHandle(AIS_Shape) &theCurAIS = occHandle(AIS_Shape)::DownCast(it.Value());
            const TopoDS_Shape &theShape = theCurAIS->Shape();
            TopExp_Explorer anExp;
            //! add the faces of the current shape
            for(anExp.Init(theShape,TopAbs_VERTEX);anExp.More();anExp.Next())
            {
                occContext->AddOrRemoveSelected(anExp.Current(),Standard_False);
            }
        }
        occContext->UpdateCurrentViewer();
    }
        break;

    default:
        break;
    }
}

//! ------------------------
//! function: startRotation
//! details:
//! ------------------------
void occGLWidget::startRotation(int x,int y,RotationPointType theRotationPointType, const gp_Pnt& theSelectedPoint)
{
    if (!occView.IsNull())
    {
        switch (theRotationPointType)
        {
        case RotationPointType_gravity:
            occView->StartRotation(x,y,0.45);
            break;
        case RotationPointType_selected:
            //cout<<"____start rotation about point ("<<theSelectedPoint.X()<<", "<<theSelectedPoint.Y()<<", "<<theSelectedPoint.Z()<<")____"<<endl;
            sx = x; sy = y;

            double X,Y;
            occView->Size(X,Y);
            rx = Standard_Real(occView->Convert(X));
            ry = Standard_Real(occView->Convert(Y));

            occView->Rotate(0.,0.,0.,theSelectedPoint.X(),theSelectedPoint.Y(), theSelectedPoint.Z(),Standard_True);

            double zRotationThreshold;
            zRotation = Standard_False;
            zRotationThreshold = 0.45;
            if( zRotationThreshold > 0. )
            {
                Standard_Real dx = Abs(sx - rx/2.0);
                Standard_Real dy = Abs(sy - ry/2.0);
                Standard_Real dd = zRotationThreshold * (rx + ry)/2.0;
                if(dx>dd || dy> dd) zRotation = Standard_True;
            }
            break;
        default:
            break;
        }
        // VSR: 10.06.2015: next line commented out - causes ugly blinking on starting rotation with Perspective projection mode
        //activeView()->DepthFitAll();
    }
}


//! -----------------
//! function: rotate
//! details:
//! -----------------
void occGLWidget::rotate(int x, int y, RotationPointType theRotationPointType, const gp_Pnt& theSelectedPoint)
{
    if (!occView.IsNull())
    {
        switch (theRotationPointType)
        {
        case RotationPointType_gravity:
            occView->Rotation(x,y);
            break;
        case RotationPointType_selected:
            double dx, dy, dz;
            if(zRotation)
            {
                dz = atan2(Standard_Real(x)-rx/2., ry/2.-Standard_Real(y)) - atan2(sx-rx/2.,ry/2.-sy);
                dx = dy = 0.;
            }
            else
            {
                dx = (Standard_Real(x) - sx) * M_PI/rx;
                dy = (sy - Standard_Real(y)) * M_PI/ry;
                dz = 0.0;
            }

            occView->Rotate(dx, dy, dz,
                            theSelectedPoint.X(),theSelectedPoint.Y(), theSelectedPoint.Z(),
                            Standard_False);
            break;

        default:
            break;
        }
        //emit vpTransformed( this );
    }
    //  setZSize( getZSize() );
}

//! --------------------------------
//! function: displayRotationCenter
//! details:
//! --------------------------------
void occGLWidget::displayRotationCenter(gp_Pnt COR, bool isOnFace)
{
    cerr<<"occGLWidget::displayRotationCenter()->____function called____"<<endl;

    static occHandle(AIS_CORMarker) CORMarker;
    if(!CORMarker.IsNull())
    {
        cerr<<"index of the current context: "<<occContext->IndexOfCurrentLocal()<<endl;
        if(occContext->IndexOfCurrentLocal()!=0)
        {
            cerr<<"occGLWidget::displayRotationCenter()->____calling Erase()____"<<endl;
            occContext->Erase(CORMarker,true);
        }
        else
        {
            cerr<<"occGLWidget::displayRotationCenter()->____calling Remove()____"<<endl;
            occContext->Remove(CORMarker,true);
        }
    }

    if(isOnFace)
    {
        /*
        //! --------------------------------------------------------
        //! In this position:
        //! when the click does not intercept a face, the previous
        //! rotation center is used.
        //! Out of "if(isOnFace){ ... }":
        //! when the click does not intercept a face, the (0,0,0)
        //! rotation center is used.
        //! --------------------------------------------------------
        if(!CORMarker.IsNull())
        {
            cerr<<"index of the current context: "<<occContext->IndexOfCurrentLocal()<<endl;
            if( occContext->IndexOfCurrentLocal()!=0)
            {
                cerr<<"____calling Erase()____"<<endl;
                occContext->Erase(CORMarker,false);
            }
            else
            {
                cerr<<"____calling Remove()____"<<endl;
                occContext->Remove(CORMarker,false);
            }
        }
        */

        cout<<"occGLWidget::displayRotationCenter()->____newCOR("<<COR.X()<<", "<<COR.Y()<<", "<<COR.Z()<<")____"<<endl;

        //! -------------------------------
        //! display the center of rotation
        //! -------------------------------
        double X,Y;
        occView->Size(X,Y);
        double L = sqrt(X*Y);
        double radius = L/100.0;
        const TopoDS_Shape &sphere = BRepPrimAPI_MakeSphere(COR,radius);
        CORMarker = new AIS_CORMarker(sphere);
        CORMarker->Attributes()->SetFaceBoundaryDraw(false);

        //! ---------------------
        //! do not draw isolines
        //! ---------------------
        occHandle(Prs3d_IsoAspect) isoLineAspect = new Prs3d_IsoAspect(Quantity_NOC_GRAY,Aspect_TOL_SOLID,1.0,1);
        isoLineAspect->SetNumber(0);
        CORMarker->Attributes()->SetUIsoAspect(isoLineAspect);
        CORMarker->Attributes()->SetVIsoAspect(isoLineAspect);

        //! ------------
        //! shaded view
        //! ------------
        if(!occContext->DefaultDrawer()->IsAutoTriangulation()) occContext->DefaultDrawer()->SetAutoTriangulation(true);
        CORMarker->SetDisplayMode(AIS_Shaded);
        CORMarker->SetColor(Quantity_NOC_RED);
        occContext->Display(CORMarker,1,-1,false);
    }
    else
    {
        //! -------------------------------------------------------------
        //! if not on face erase the rotation center: same block of code
        //! -------------------------------------------------------------
        if(!CORMarker.IsNull())
        {
            cerr<<"not on face"<<endl;
            cerr<<"index of the current context: "<<occContext->IndexOfCurrentLocal()<<endl;
            if(occContext->IndexOfCurrentLocal()!=0)
            {
                cerr<<"occGLWidget::displayRotationCenter()->____calling Erase()____"<<endl;
                occContext->Erase(CORMarker,true);
            }
            else
            {
                cerr<<"occGLWidget::displayRotationCenter()->____calling Remove()____"<<endl;
                occContext->Remove(CORMarker,true);
            }
        }
    }
    //! ------------------
    //! update the viewer
    //! ------------------
    occContext->UpdateCurrentViewer();
}

//! -------------------------
//! function: displayShape()
//! details:
//! -------------------------
void occGLWidget::displaySphericalMarker(gp_Pnt P)
{
    static occHandle(AIS_Shape) anAIS_Shape;
    if(!anAIS_Shape.IsNull()) occContext->Remove(anAIS_Shape,true);

    if(!occContext->DefaultDrawer()->IsAutoTriangulation()) occContext->DefaultDrawer()->SetAutoTriangulation(true);

    //! --------------
    //! sphere marker
    //! --------------
    double X,Y;
    occView->Size(X,Y);
    double L = sqrt(X*Y);
    double radius = L/125.0;
    occHandle(AIS_Shape) aShape = markers::buildSphereMarker(P,radius);

    //! ---------------------
    //! do not show isolines
    //! ---------------------
    occHandle(Prs3d_IsoAspect) isoLineAspect = new Prs3d_IsoAspect(Quantity_NOC_GRAY,Aspect_TOL_SOLID,1.0,1);
    isoLineAspect->SetNumber(0);
    aShape->Attributes()->SetUIsoAspect(isoLineAspect);
    aShape->Attributes()->SetVIsoAspect(isoLineAspect);

    //! ----------------
    //! face boundaries
    //! ----------------
    aShape->Attributes()->SetFaceBoundaryDraw(false);

    occContext->Display(aShape,1,-1,false);
    occContext->UpdateCurrentViewer();

    anAIS_Shape = aShape;
}

//! -------------------------------------------------------
//! function: setTransparency
//! details:  enable/disable custom transparency on shapes
//!           AIS_KOI_Shape,signature=0
//! -------------------------------------------------------
void occGLWidget::setTransparency(bool isActive, bool updateViewer, double level)
{
    if(!isActive) level = 0.0;

    //! -------------------------------------------------------------
    //! retrieve the AIS_Shapes with signature "0" (Shape signature)
    //! -------------------------------------------------------------
    AIS_ListOfInteractive AISList;
    AIS_ListIteratorOfListOfInteractive iter;    
    occContext->ObjectsInside(AISList,AIS_KOI_Shape,0);
    for(iter.Initialize(AISList);iter.More();iter.Next())
    {
        const occHandle(AIS_Shape) &shape = occHandle(AIS_Shape)::DownCast(iter.Value());
        occContext->SetTransparency(shape,level,false);
    }
    if(updateViewer) occContext->UpdateCurrentViewer();
}


//! ---------------------------------------------
//! function: hideAllMarkers
//! details:  hide all markers, including triads
//! ----------------------------------------------
void occGLWidget::hideAllMarkers(bool updateViewer)
{
    cout<<"occGLWidget::hideAllMarkers()->____function called____"<<endl;

    occContext->CloseLocalContext(-1,false);
    AIS_ListOfInteractive theListOfIO;
    AIS_ListIteratorOfListOfInteractive it;
    QList<int> listOfSignatures;
    listOfSignatures<<CUSTOM_SIGNATURE_COR_MARKER<<
                      CUSTOM_SIGNATURE_MESH_SEGMENT_MARKER<<
                      CUSTOM_SIGNATURE_SPHERE_MARKER<<
                      CUSTOM_SIGNATURE_ARROW_MARKER<<
                      CUSTOM_SIGNATURE_DOUBLE_ARROW_MARKER<<
                      CUSTOM_SIGNATURE_CUSTOM_TRIHEDRON_MARKER<<
                      CUSTOM_SIGNATURE_CURVED_ARROW;

    //! -------------
    //! hide markers
    //! -------------
    for(QList<int>::iterator listIt = listOfSignatures.begin(); listIt!=listOfSignatures.end(); ++listIt)
    {
        int signature = *listIt;
        occContext->ObjectsByDisplayStatus(AIS_KOI_Shape,signature,AIS_DS_Displayed,theListOfIO);
        for(it.Initialize(theListOfIO); it.More(); it.Next())
        {
            occContext->Erase(it.Value(),false);
        }
    }

    //! ------------------------------------------
    //! hide also the triads (type AIS_KOI_Datum)
    //! ------------------------------------------
    //cout<<"____removing system of references____"<<endl;
    occContext->ObjectsByDisplayStatus(AIS_KOI_Datum,-1,AIS_DS_Displayed,theListOfIO);
    for(it.Initialize(theListOfIO); it.More(); it.Next())
    {
        occContext->Remove(it.Value(),false);
    }
    if(updateViewer) occContext->UpdateCurrentViewer();

    myLocalCtxNumber = occContext->OpenLocalContext();
    this->reactivateCurrentStandardSelectionMode();
}

//! ----------------------------------------------------------------
//! function: addClipPlane
//! details:  add a clip plane to the viewer, using the assigned ID
//!           and the coefficients {A,B,C,D} if the clip plane ID
//!           already exists, update it
//! ----------------------------------------------------------------
void occGLWidget::addClipPlane(double A, double B, double C, double D, int ID, bool isOn)
{
    //! ------------------------------------------------------------------
    //! add the clip plane only if it is not already contained in the map
    //! ------------------------------------------------------------------
    occHandle(Graphic3d_ClipPlane) aClipPlane;
    if(myMapOfClipPlanes.value(ID,aClipPlane).IsNull())
    {
        aClipPlane = new Graphic3d_ClipPlane();
    }
    else
    {
        aClipPlane = myMapOfClipPlanes.value(ID);
    }
    cout<<"occGLWidget::addClipPlane()->____initial number of clip planes: "<<myMapOfClipPlanes.size()<<"____"<<endl;

    //! --------------------------------------
    //! change equation of the clipping plane
    //! --------------------------------------
    aClipPlane->SetEquation (gp_Pln (A,B,C,D));

    //! ------------
    //! set capping
    //! ------------
    aClipPlane->SetCapping(false);

    //! ------------------------------------------------
    //! add the clip plane only to the AIS_Shape(s)
    //! a new clip plane is created if not existing;
    //! if already existing the clip plane is retrieved
    //! and set into the state "On"
    //! ------------------------------------------------
    AIS_ListOfInteractive listOfAISShapes;
    occContext->ObjectsInside(listOfAISShapes,AIS_KOI_Shape,0);
    for(AIS_ListIteratorOfListOfInteractive it(listOfAISShapes); it.More(); it.Next())
    {
        const occHandle(AIS_Shape) &curShapeObject = occHandle(AIS_Shape)::DownCast(it.Value());
        if(curShapeObject.IsNull()) continue;
        if(curShapeObject->ClipPlanes().IsNull())
        {
            curShapeObject->AddClipPlane(aClipPlane);
        }
        else
        {
            const occHandle(Graphic3d_ClipPlane) &curClipPlane = curShapeObject->ClipPlanes()->Value(ID);
            if(curClipPlane.IsNull()) curShapeObject->AddClipPlane(aClipPlane);
            else curShapeObject->ClipPlanes()->Value(ID)->SetOn(true);
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
    //aClipPlane->SetOn(isOn);
    occView->Redraw();

    //! --------------------------------------
    //! add/replace the clip plane to the map
    //! --------------------------------------
    myMapOfClipPlanes.insert(ID,aClipPlane);
    cout<<"occGLWidget::addClipPlane()->____final number of clip planes: "<<myMapOfClipPlanes.size()<<"____"<<endl;
}

//! --------------------------
//! function: removeClipPlane
//! details:
//! --------------------------
void occGLWidget::removeClipPlane(int ID)
{
    cout<<"occGLWidget::removeClipPlane()->____removing clip plane ID: "<<ID<<"____"<<endl;

    //! ----------------------------------------------
    //! remove thee clip plane ID from all the shapes
    //! ----------------------------------------------
    AIS_ListOfInteractive listOfAISShapes;
    occContext->ObjectsInside(listOfAISShapes,AIS_KOI_Shape,0);
    for(AIS_ListIteratorOfListOfInteractive it(listOfAISShapes); it.More(); it.Next())
    {
        const occHandle(AIS_Shape) &curShape = occHandle(AIS_Shape)::DownCast(it.Value());
        curShape->RemoveClipPlane(myMapOfClipPlanes.value(ID));
    }

    //! ---------------------------------------------
    //! remove the clip plane ID from the whole view
    //! ---------------------------------------------
    //occView->RemoveClipPlane(myMapOfClipPlanes.value(ID));

    //! -------
    //! redraw
    //! -------
    occView->Redraw();

    //! --------------------
    //! remove from the map
    //! --------------------
    myMapOfClipPlanes.remove(ID);
}

//! -------------------------------
//! function: getClipPlanesNbLimit
//! details:
//! -------------------------------
Standard_Integer occGLWidget::getClipPlanesNbLimit()
{
    return occViewer->Driver()->InquirePlaneLimit();
}

//! -------------------------
//! function: setClipPlaneOn
//! details:
//! -------------------------
void occGLWidget::setClipPlaneOn(int ID, bool isOn)
{
    const occHandle(Graphic3d_ClipPlane)& aClipPlane = myMapOfClipPlanes.value(ID);
    //const occHandle(Graphic3d_ClipPlane)& aClipPlane = myMapOfClipPlanes.at(ID);
    aClipPlane->SetOn(isOn);
    occView->Redraw();
}

//! -------------------------------------
//! function: updateClipPlaneTranslation
//! details:
//! -------------------------------------
void occGLWidget::updateClipPlaneTranslation(int ID, int zVal, const QVector<double> &coeffs)
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

//! ------------------------
//! function: hideAllBodies
//! details:  slot
//! ------------------------
void occGLWidget::hideAllBodies()
{
    occContext->CloseAllContexts(false);
    AIS_ListOfInteractive AISList;
    occContext->ObjectsByDisplayStatus(AIS_KOI_Shape,-1,AIS_DS_Displayed,AISList);
    AIS_ListIteratorOfListOfInteractive it;
    for(it.Initialize(AISList);it.More();it.Next())
    {
        const occHandle(AIS_ExtendedShape) &curAISShape = occHandle(AIS_ExtendedShape)::DownCast(it.Value());
        curAISShape->setShapeVisibility(false);
        occContext->Erase(curAISShape,false);
    }
    myLocalCtxNumber = occContext->OpenLocalContext(true);
    this->reactivateCurrentStandardSelectionMode();
}

//! -------------------------
//! function: ConvertToPlane
//! details:
//! -------------------------
Standard_Boolean occGLWidget::ConvertToPlane(const Standard_Integer Xs,
                                             const Standard_Integer Ys,
                                             Standard_Real& X,
                                             Standard_Real& Y,
                                             Standard_Real& Z,
                                             Standard_Boolean usePrecision)
{
    Standard_Real Xp = Xs, Yp = Ys;
    Standard_Real Xv,Yv,Zv;
    Standard_Real Vx,Vy,Vz;
    gp_Pln aPlane(occView->Viewer()->PrivilegedPlane());
    occView->Convert(Xp,Yp,Xv,Yv,Zv);
    occView->Proj(Vx,Vy,Vz);
    gp_Lin aLine(gp_Pnt(Xv,Yv,Zv), gp_Dir(Vx,Vy,Vz));
    IntAna_IntConicQuad theIntersection (aLine, aPlane, Precision::Angular());
    if (theIntersection.IsDone())
    {
        if (!theIntersection.IsParallel())
        {
            if (theIntersection.NbPoints() > 0)
            {
                gp_Pnt theSolution(theIntersection.Point(1));
                X = theSolution.X();
                Y = theSolution.Y();
                Z = theSolution.Z();
                if (usePrecision)
                {
                    //const double myPrecision = 0.01;
                    //X = (X 0. ? 1 : 0.)) * floor((abs(X)) / myPrecision) * myPrecision;
                    //Y = (Y 0. ? 1 : 0.)) * floor((abs(Y)) / myPrecision) * myPrecision;
                    //Z = (Z 0. ? 1 : 0.)) * floor((abs(Z)) / myPrecision) * myPrecision;
                }
                return Standard_True;
            }
        }
    }
    return Standard_False;
}
