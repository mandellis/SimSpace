#ifndef OCCGLWIDGET_H
#define OCCGLWIDGET_H

#define ValZWMin 1

//! ---
//! Qt
//! ---
#include <QGLWidget>
#include <QMap>

//! ----
//! OCC
//! ----
#include <OpenGl_GraphicDriver.hxx>
#include <Aspect_Handle.hxx>
#include <AIS_InteractiveContext.hxx>
#include <Aspect_DisplayConnection.hxx>
#include <WNT_Window.hxx>
#include <V3d_View.hxx>
#include <AIS_Shape.hxx>
#include <AIS_TextLabel.hxx>
#include <TopTools_ListOfShape.hxx>
#include <Graphic3d_ClipPlane.hxx>
#include <Graphic3d_SequenceOfHClipPlane.hxx>
#include <AIS_Plane.hxx>
#include <IntAna_IntConicQuad.hxx>

//! ----------------
//! custom includes
//! ----------------
#include "mydefines.h"
#include "actions3d.h"
#include "selectionmodes.h"
#include "displaymode.h"
#include "displayquality.h"
#include "ais_customsignatures.h"

#define HIGHLIGHT_COLOR Quantity_NOC_RED
#define TRANSPARENCY_HIGHLIGHT 0.50
#define TRIHEDRON_AXIS_LENGTH_DIAGONAL_FACTOR 20.0

class QRubberBand;
class QMenu;
class QPoint;
class QActionGroup;
class TopoDS_Face;

class occGLWidget: public QGLWidget
{
    Q_OBJECT

protected:

    //! rotation point type
    enum RotationPointType
    {
        RotationPointType_gravity,
        RotationPointType_selected
    };

    RotationPointType myRotationPointType;

    //! the 3D actions
    CurAction3D myCurAction3D;
    CurAction3D myCurAction3D_old;

    //! the selection mode
    CurSelectionMode myCurSelectionMode;

    //! the selection type
    SelectionType myCurSelectionType;

    //! the global selection mode - can be "single selection" or "multiple selection"
    CurGlobalSelectionMode myCurGlobalSelectionMode;    

    //! z-layer for model view
    Standard_Integer my_zLayerID_model;

    //! zlayer
    Standard_Integer my_zLayer;

    //! set transparency
    void setTransparency(bool isActive, bool updateViewer=false, double level = 0.5);

    //! current clip plane ID - experimental
    int myCurClipPlaneID;

public:

    //! constructor
    occGLWidget(QWidget* parent=0);

    //! destructor
    virtual ~occGLWidget();

    //! get the context
    const occHandle(AIS_InteractiveContext)& getContext() const;

    //! get the view
    const occHandle(V3d_View)& getView() const;

    //! view mode
    CurDisplayMode myCurDisplayMode, myCurDisplayMode_old;

    //! An option of curAction3D_Rotation
    Standard_Boolean myAllowSinglePick;

    //! maximum and minimum mouse pointer coordinates
    int myXmin, myYmin, myXmax, myYmax;

    //! rectangle for the mouse selection
    QRubberBand *myRectBand;

    //! return the current state of 3D action
    CurAction3D curAction3D();

    //! return the current selection mode
    TopAbs_ShapeEnum curSelectionMode();

    //! set the selection mode
    virtual void setSelectionMode(CurSelectionMode selectionMode);

    //! get the current selection mode
    CurSelectionMode getCurrentSelectionMode(){ return myCurSelectionMode; }

    //! create the actions
    void createTheActions();

    //! return the current 3D action
    CurAction3D currentAction3D();

    //! add a clip plane
    virtual void addClipPlane(double A, double B, double C, double D, int ID, bool isOn = true);

    //! remove a clip plane
    virtual void removeClipPlane(int ID);

    //! switch off a clip plane
    void setClipPlaneOn(int ID, bool isOn);

    //! update clip plane definition
    virtual void occGLWidget::updateClipPlaneCoefficients(int ID, const QVector<double> &coeffs);

    //! set the current clip plane
    virtual void setCurrentClipPlane(int curClipPlaneID)
    {
        myCurClipPlaneID = curClipPlaneID;
    }

    //! get the current clip plane ID - experimental
    virtual int getCurrentClipPlaneID() {return myCurClipPlaneID; }

    //! get clip planes Nb limit
    Standard_Integer getClipPlanesNbLimit();

    //! map of clipping planes
    QMap<int,occHandle(Graphic3d_ClipPlane)> myMapOfClipPlanes;

protected:

    //! active opened context - "0" is neutral, ">0" is for selection
    Standard_Integer myLocalCtxNumber;

    //! the occ viewer
    occHandle(V3d_Viewer) occViewer;

    //! the occ view
    occHandle(V3d_View) occView;

    //! the occ context
    occHandle(AIS_InteractiveContext) occContext;

    //! init
    virtual void init();

    //! clear the "current" or "selected" shapes
    void clearGeometrySelection();

    //! the left upper label content
    occHandle(AIS_TextLabel) myTextLabel;

    //! the floating label
    occHandle(AIS_TextLabel) myFloatingLabel;

    //! my camera target
    gp_Pnt myCOR;

    //! start rotation
    void startRotation(int x, int y, RotationPointType theRotationPointType, const gp_Pnt& theSelectedPoint);

    //! rotate
    void rotate(int x, int y, RotationPointType theRotationPointType, const gp_Pnt& theSelectedPoint);

    //! my display quality
    displayQuality myDisplayQuality;

    //! Paint events
    virtual void paintEvent(QPaintEvent *e);
    virtual void resizeEvent(QResizeEvent *e);

    //! Mouse events
    virtual void mousePressEvent(QMouseEvent *e);
    virtual void mouseReleaseEvent(QMouseEvent *e);
    virtual void mouseMoveEvent(QMouseEvent *e);
    virtual void mouseDoubleClickEvent(QMouseEvent *e);
    virtual void wheelEvent(QWheelEvent *e);
    virtual void onMouseMove(const int theFlags, QPoint thePoint);

    //! Keyboard events
    virtual void keyPressEvent(QKeyEvent *theKey);

    //! Button events
    virtual void onLButtonDown(const int theFlags,const QPoint thePoint);
    virtual void onMButtonDown(const int theFlags,const QPoint thePoint);
    virtual void onRButtonDown(const int theFlags,const QPoint thePoint);
    virtual void onLButtonUp(const int theFlags,const QPoint thePoint);
    virtual void onMButtonUp(const int theFlags,const QPoint thePoint);
    virtual void onRButtonUp(const int theFlags,const QPoint thePoint);
    virtual void onMouseWheel(const int theFlags,const int theDelta,const QPoint thePoint);
    //virtual void onMouseDoubleClick(const int theFlags, const QPoint thePoint);

    //! Rubber band
    void drawRubberBand(const int minX,const int minY,const int maxX,const int maxY);

    //! pick point
    gp_Pnt hitPoint(long x, long y, TopoDS_Shape shape);

    //! convert to plane
    Standard_Boolean ConvertToPlane(const Standard_Integer Xs,
                                    const Standard_Integer Ys,
                                    Standard_Real& X,
                                    Standard_Real& Y,
                                    Standard_Real& Z,
                                    Standard_Boolean usePrecision);

    //! Context menu - it is defined as protected member:
    //! The derived class can enrich this menu by adding
    //! additional actions
    QMenu *myContextMenu;
    QMenu *myCursorModeMenu;
    QMenu *myViewModeMenu;

    //! Actions - selection modes. Actually the "pick point coordinates" mode
    //! is not a "selection" action, but when this function is activated
    //! the other selection mode are deactivates
    QActionGroup *selectionModeActionsGroup;
    QAction *selectVertex;
    QAction *selectEdge;
    QAction *selectFace;
    QAction *selectSolid;
    QAction *selectPickPointCoordinates;

    //! ---------------------
    //! Actions - view modes
    //! ---------------------
    QAction *viewFront;
    QAction *viewBack;
    QAction *viewRight;
    QAction *viewLeft;
    QAction *viewTop;
    QAction *viewBottom;
    QAction *viewIsometric;
    QAction *actionZoomToFit;

    //! -------------------
    //! Show - hide bodies
    //! -------------------
    QAction *actionHideBody;
    QAction *actionHideAllOtherBodies;
    QAction *actionShowAllBodies;
    QAction *actionSelectAll;

public slots:

    //! sets the 3D action
    //! change the internal state of the widget
    void setAction3D_Rotation();
    void setAction3D_Pan();
    void setAction3D_WindowZooming();

    //! fit all
    void FitAll();

    //! select all: select all according to the current selection mode
    void selectAll();

    //! sets the view mode wireframe
    virtual void setWireframeView();

    //! sets the view mode shaded
    virtual void setShadedExteriorView();

    //! sets the view mode shaded and wireframe
    virtual void setShadedExteriorAndEdgesView();

    //! set selection type
    virtual void setSelectionType(SelectionType &aSelectionType) { myCurSelectionType = aSelectionType; }

    //! sets the global selection mode
    virtual void setGlobalCurSelectionMode(int);

    //! reactivate the current standard selection mode
    virtual void reactivateSelectionMode();

    //! set a gradient bakground
    void setBackgroundColor(double R1, double G1, double B1, double R2, double G2, double B2, int tof);

    //! sets the display quality
    void setDisplayQuality(displayQuality);

    //! display a trihedron
    void displayTrihedron(QVector<double> Origin, QVector<QVector<double> > directionalData, int axisLength=-1);

    //! remove the trihedra from the viewer
    void removeTrihedron();

    //! extends the selection to the adjacent - apply on a single body for the moment
    void extendSelectionToAjacent();

    //! save image to disk
    void saveImageToDisk();

    //! display shape
    void displaySphericalMarker(gp_Pnt C);

    //! hide all markers, including triads
    void hideAllMarkers(bool updateViewer = true);

    //! reset
    virtual void reset();

protected slots:

    //! geometric properties of a selection
    virtual void computeSelectionProperties();

    //! build a basic context menu structure
    void buildMinimalContexetMenu();

    //! hide all bodies
    virtual void hideAllBodies();

    //! hide selected bodies
    virtual void hideSelectedBodies();

    //! show all bodies
    virtual void showAllBodies();

    //! hide all the other bodies
    virtual void hideAllTheOtherBodies();

    //! set isometric view
    void isometricView();    

    //! unset the selection modes
    void unsetSelectionModes();

    //! unset the view operations
    void unsetViewOperations();

    //! write a timestamp
    void writeLabel(const QString& theTextLabel, int xPos, int yPos, opencascade::handle<Graphic3d_TransformPers> &);

    //! context menu
    virtual void ShowContextMenu(const QPoint&);

    //! highlight a list of shapes
    //! void highlightShapes(const TopTools_ListOfShape &aShapeToHighlight);

    //! clear selection
    void clearSelected();

signals:

    //! the selection is changed - emits a signal
    void selectionChanged();

    //! emits a statusBarMessage - to be written (e.g.) into the main window status bar
    void statusBarMessage(QString);

    //! The selection mode is changed within the widget by the context menu.
    //! When a selection mode is activated in this way it is also necessary
    //! to toggle the corresponding button within the mainwindow toolbar,
    //! and untoggle the other. In order to do that a signal is emitted, which
    //! is connnected to the MainWindow::toggleSelectionMode(curSelectionMode)
    void selectionModeVertex(bool);
    void selectionModeEdge(bool);
    void selectionModeFace(bool);
    void selectionModeSolid(bool);
    void selectionModePickPointCoordinates(bool);

    //! the center of view has been changed
    void CORchanged(gp_Pnt newCOR, bool isOnFace);

    ///! under testing - contruction
    void highlightmeshface();
    void highlightmeshedge();

    //! signal: the view mode has been changed
    void viewModeChanged(CurDisplayMode);

protected slots:

    //! update the camera target
    void updateCOR(gp_Pnt newCOR, bool isOnFace);

    //! display center of rotation
    void displayRotationCenter(gp_Pnt newOrigin, bool isOnFace);    

public:

    //! get clip planes
    QMap<int,occHandle(Graphic3d_ClipPlane)> getClipPlanes() { return myMapOfClipPlanes; }

    //! get scene bounding box
    void getSceneBoundingBox(double &lx, double &ly, double &lz);

    //! get selection type
    SelectionType getSelectionType() const { return myCurSelectionType; }

signals:

    //! request synch visibility
    void requestSynchVisibility(bool newIsVisible);
};

#endif // OCCGLWIDGET_H
