#ifndef OCCPREGLWIDGET_H
#define OCCPREGLWIDGET_H

//! ----------------
//! custom includes
//! ----------------
#include "occGLwidget.h"
#include "ais_extendedshape.h"
#include "ArrayOfColors.h"
#include "workingmode.h"
#include "simulationnodeclass.h"
#include "myenumvariables.h"
#include "ais_colorscaleextended.h"
#include <indexedmapofmeshdatasources.h>
#include "geometrytag.h"

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>
#include <AIS_SequenceOfInteractive.hxx>
#include <NCollection_Array1.hxx>
#include <TColStd_ListOfAsciiString.hxx>
#include <MeshVS_Mesh.hxx>
#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <MeshVS_Mesh.hxx>

//! ---
//! Qt
//! ---
#include <QFocusEvent>
#include <QMap>

//! ----
//! C++
//! ----
#include <set>
#include <map>

//! ---------------------
//! for testing purposes
//! ---------------------
#include <NCollection_Array1.hxx>
#include <ng_meshvs_datasourceface.h>

class QActionGroup;
class writeLabelClass;
class TopoDS_Shape;
class simulationDataBase;

class occPreGLWidget: public occGLWidget
{
    Q_OBJECT

protected:

    virtual void mousePressEvent(QMouseEvent *e) override;
    virtual void mouseReleaseEvent(QMouseEvent *e) override;
    virtual void mouseMoveEvent(QMouseEvent *e) override;
    virtual void onLButtonDown(const int theFlags, const QPoint thePoint) override;
    virtual void onRButtonDown(const int theFlags, const QPoint thePoint) override;
    virtual void onMButtonDown(const int theFlags, const QPoint thePoint) override;
    virtual void onLButtonUp(const int theFlags, const QPoint thePoint) override;
    virtual void onRButtonUp(const int theFlags, const QPoint thePoint) override;
    virtual void onMButtonUp(const int theFlags, const QPoint thePoint) override;
    virtual void onMouseMove(const int theFlags, QPoint thePoint) override;

protected:

    //! the occ context for the mesh view
    occHandle(AIS_InteractiveContext) occMeshContext;

    //! init
    virtual void init() override;

    //! paint event
    virtual void paintEvent(QPaintEvent *e) override;

    //! ... experimental ...
    virtual void focusInEvent(QFocusEvent *focusEvent);

    //! The mesh view mode
    bool isMeshViewVolume;

    //! The working mode
    curWorkingMode myCurWorkingMode;

    //! Map of interactive shapes
    QMap<int,occHandle(AIS_InteractiveObject)> myMapOfInteractiveShapes;

    //! Map of interactive meshes
    QMap<int,occHandle(AIS_InteractiveObject)> myMapOfInteractiveMeshes;

    //! Map of prismatic meshes
    QMap<int,occHandle(AIS_InteractiveObject)> myMapOfInteractivePrismaticMeshes;

    //! Shape colors: an array of colors for the shapes
    ArrayOfColors myShapeColor;

    //! Get the size of the mesh: total number of nodes and elements
    //! at the end of the meshing process
    void getMeshSize();

    //! The target of the camera
    Standard_Real myCOVX, myCOVY, myCOVZ;

    //! build the suppression menu
    void buildSuppressionContextMenu();

    //! The actions for the context menu
    QAction *actionGenerateMeshOnSelectedBodies;
    QAction *actionPreviewSurfaceOnSelectedBodies;
    QAction *actionClearGeneratedDataOnSelectedBodies;

    //! print a summary
    void printSummary();

public:

    //! constructor
    occPreGLWidget(QWidget *parent=0);

    //! destructor
    virtual ~occPreGLWidget();

    //! set the selection mode
    virtual void setSelectionMode(CurSelectionMode selectionMode) override;

    //! Display the CAD model
    void displayCAD(bool onlyLoad=false);

    //! set the simulation database
    void setSimulationdataBase(simulationDataBase* theDB) { myDS2 = theDB; }

    //! create the interactive shapes
    void createInteractiveShapes();

    //! return the mesh context
    occHandle(AIS_InteractiveContext) getMeshContext() {return occMeshContext; }

    //! get the interactive objects (AIS_Shape)
    QMap<int,occHandle(AIS_InteractiveObject)> getInteractiveObjects() { return myMapOfInteractiveShapes; }

    //! get the mesh objects (MeshVS_Mesh)
    QMap<int,occHandle(AIS_InteractiveObject)> getMeshObjects() { return myMapOfInteractiveMeshes; }

    //! get current working mode
    inline curWorkingMode getWorkingMode() { return myCurWorkingMode; }

    //! reset the mesh and shape maps
    void resetShapeMaps()
    {
        myMapOfInteractiveMeshes.clear();
        myMapOfInteractiveShapes.clear();

        //! clean the viewer
        getContext()->RemoveAll(true);
    }

    //! create and add a clip plane to the viewer
    virtual void addClipPlane(double A, double B, double C, double D, int ID, bool isOn = true) override;

    //! remove a clip plane
    virtual void removeClipPlane(int ID) override;

    //! update clip plane translation
    virtual void updateClipPlaneCoefficients(int ID, const QVector<double> &coeffs) override;

    //! update clip planes
    void updateClipPlanes(const std::vector<int> &activeClipPlanes);

    //! shape color
    Quantity_Color shapeColor(int bodyIndex)
    {
        const occHandle(AIS_Shape) &aShape = occHandle(AIS_Shape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
        Quantity_Color acolor;
        occContext->Color(aShape,acolor);
        return acolor;
    }

    //! set all elements visible
    void setAllElementsVisible();

    //! clip mesh
    void clipMesh();

public:

    //! The simulation database
    simulationDataBase *myDS2;

protected:

    //! topology numbers
    void getTopologyNumber(const TopoDS_Shape &theParentShape, const TopoDS_Shape &, int &parentShapeIndex, int &subShapeIndex);

    //! Set selection mode green solid
    //void setSelectionStyle(Quantity_NameOfColor color, float transparency);

    //! event
    virtual bool event(QEvent *event);

public slots:

    //! sets the vireframe view
    virtual void setWireframeView() override;

    //! sets the view mode shaded
    virtual void setShadedExteriorView() override;

    //! sets the view mode shaded and wireframe
    virtual void setShadedExteriorAndEdgesView() override;

    //! build the mesh interactive object
    void buildMeshIOs();

    //! build the sliced mesh interactive object
    void buildSlicedMeshIO(const QMap<int, opencascade::handle<MeshVS_DataSource> > &slicedMeshDS);

    //! set hidden elements
    void setHiddenElements(const std::map<int,occHandle(TColStd_HPackedMapOfInteger)> &hiddenElements);

    //! build the prismatic mesh interactive object
    void buildPrismaticMeshIO();

    //! reset the custom colors
    void resetCustomColors(bool updateViewer=false);

    //! display all the meshes
    void displayAllMeshes(bool meshNodesVisible=false, Graphic3d_NameOfMaterial theMaterial=Graphic3d_NOM_GOLD, bool updateViewers=true);

    //! hide all the meshes
    void hideAllMeshes(bool updateViewer=true);

    //! refresh mesh view
    virtual void refreshMeshView(bool onlyExterior);

    //! Reset all
    virtual void reset() override;

    //! Activates the operating mode "model"
    void setWorkingMode_Model();

    //! Activates the operating mode "mesh"
    void setWorkingMode_Mesh();

    //! Activate the operating mode "named selection"
    void setWorkingMode_NamedSelection();

    //! Activates the operating mode "mesh"
    void setWorkingMode_Contacts();

    //! Activate the operating mode "solution"
    //void setWorkingMode_Solution();

    //! set selection mode mesh
    void setMeshSelectionMode();

    //! set selection mode geometry
    void setGeometrySelectionMode();

    //! invalidate mesh
    void invalidateMesh(int bodyIndex, Standard_Boolean updateViewer=Standard_False);
    void invalidateMeshes(const std::vector<int> &indexes);
    void invalidateAllMeshes();

    //! display a contact pair
    void displayShapeCopy(const TopTools_ListOfShape &list1, const TopTools_ListOfShape &list2,
                          Quantity_NameOfColor color1, Quantity_NameOfColor color2, QVariant options=QVariant());

    void displayShapeCopy1(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color);
    void displayShapeCopy2(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color);

    //! display a list of shapes with a custom color
    void displayColoredSubshapes(const TopTools_ListOfShape &listOfShape = TopTools_ListOfShape(),
                                 Quantity_NameOfColor aColor=Quantity_NOC_YELLOW,
                                 bool cleanPrevious=true,
                                 bool updateViewer=true);

    //! clear all the custom colors (all the custom colors from all the interactive objects
    void clearAllCustomColors();

    //! hide body - used this method only to handle the suppression status
    virtual void  hideBody(const TColStd_ListOfInteger &listOfBodyNumbers);

    //! hide body - used this method only to handle the suppression status
    //virtual void  showBody(const std::vector<int> &bodiesVector);
    virtual void  showBody(const TColStd_ListOfInteger &listOfBodies);

    //! remove obsolete meshes
    void removeObsoleteMeshes();

    //! highlight bodies
    void highlightBody(const QList<int> &listOfBodyNumbers);

    //! unhighlight visible bodies
    void unhighlightBody(bool updateViewer=false);

    //! build triangulation
    void replaceTriangulation();

    //! display surface triangulation (for testing purposes)
    void displayTriangulation(const TopoDS_Shape &aShape);

    //! ----------------------------------
    //! change the transparency of bodies
    //! ----------------------------------
    void changeBodyTransparency(const std::vector<int> indexes, float transparency)
    {
        for(std::vector<int>::const_iterator it = indexes.cbegin(); it!= indexes.cend(); ++it)
        {
            int bodyIndex = *it;
            const occHandle(AIS_Shape) &theAIS_Shape = occHandle(AIS_Shape)::DownCast(myMapOfInteractiveShapes.value(bodyIndex));
            //const occHandle(AIS_Shape) &theAIS_Shape = occHandle(AIS_Shape)::DownCast(myMapOfInteractiveShapes.at(bodyIndex));
            occContext->SetTransparency(theAIS_Shape,transparency,false);
        }
        occContext->UpdateCurrentViewer();
    }

    //! ------------------------
    //! update the mesh context
    //! ------------------------
    void updateMeshContext()
    {
        cout<<"occPreGLWidget::updateMeshContext()->____function called: updating mesh view____"<<endl;
        AIS_ListOfInteractive meshList;
        //! -----------------------------------------------------------------------------------------------------------
        //! From OCC documentation: AIS_InteractiveContext::ObjectsInside()
        //! fills <aListOfIO> with objects of a particular Type and Signature with no consideration of display status.
        //! By Default, <WhichSignature> = -1 means control only on <WhichKind>.
        //! If <WhichKind> = AIS_KOI_None and <WhichSignature> = -1, all the objects are put into the list.
        //! -----------------------------------------------------------------------------------------------------------
        occMeshContext->ObjectsInside(meshList);
        AIS_ListIteratorOfListOfInteractive it;

        //! -------------------------------
        //! this is slow ... try to change
        //! -------------------------------
        for(it.Initialize(meshList); it.More(); it.Next())
        {
            const occHandle(MeshVS_Mesh) &aMeshVS_Mesh = occHandle(MeshVS_Mesh)::DownCast(it.Value());
            occMeshContext->Redisplay(aMeshVS_Mesh,true,false);
        }
        occMeshContext->UpdateCurrentViewer();
    }

    //! for testing purposes
    void displayCurvatureMap();


    //! apply custom colors
    void applyCustomColors(const QMap<GeometryTag,TopoDS_Shape> &subShapeMaps, Quantity_NameOfColor aColor, bool updateViewer = true);

private slots:

    //! show subshape mesh
    void showFaceMesh();
    void showEdgeMesh();

    void exportSTL();
    void exportCloud();

protected slots:

    void ShowContextMenu1(const QPoint&);

    //! Click on a topology -> returns the number
    void printTopologyNumber();

    //! Click on a topology -> returns the number
    void identifyTheSelection(int &vertexIndex, int &edgeIndex, int &faceIndex, int &solidIndex);

    //! Clears data on the selected bodies (meshes)
    void clearMeshFromViewer();

    //! Hide body
    void hideSelectedBodies() override;

    //! Show all bodies
    void showAllBodies() override;

    //! Hide all the other bodies
    void hideAllTheOtherBodies() override;

    //! suppress all other bodies
    void suppressAllOtherBodies();


public slots:

    //! -------------------------------------
    //! update viewer after data base change
    //! -------------------------------------
    void updateViewerAfterDataBaseChange();

    //! -----------------------
    //! show mesh data sources
    //! -----------------------
    void showMeshDataSources(const IndexedMapOfMeshDataSources &indexedMapOfDS);

    //! -----------------------------------------
    //! set transparency on working mode contact
    //! -----------------------------------------
    void setTransparency_onWorkingModeContact(double aLevel)
    {
        if(myCurWorkingMode == curWorkingMode_onContact) this->setTransparency(true,true,aLevel);
    }

private:

    //! Keyboard handling
    virtual void keyPressEvent(QKeyEvent *theKey) override;

    //! build mesh toos context menu
    void buildMeshToolsContextMenu(QMenu *aContextMenu);

signals:

    //! system messages update
    void newPanelMessage(QString);

    //! emit working mode changed
    void workingModeChanged(enum curWorkingMode);

    //! emit boundary condition changed
    void boundaryConditionsChanged(bool);

    //! new mesh control inserted
    void meshControlsChanged(bool);

    void requestSynchVisibility();
    void requestGenerateMesh(bool isVolume);
    void requestShowHealingElements();
    void requestCreateSimulationNode(SimulationNodeClass::nodeType aNodeType);
    void requestChangeSuppressionStatus(Property::SuppressionStatus newSuppressionStatus);
    void requestDisableSelectionButtons();
    void requestEnableSelectionButtons();
    void requestClearGeneratedData();
    void requestChangeDelegateContext(const occHandle(AIS_InteractiveContext) &aCTX);

    void requestBuildFaceMesh(const TopoDS_Face &aFace);
    void requestExportStepFile();

public:

    //! display Mesh
    void displayMesh(const occHandle(MeshVS_DataSource) &aMeshDS,
                     Quantity_NameOfColor aColorName = Quantity_NOC_ALICEBLUE,
                     bool showMeshEdges = true);

    static void displayMesh(const occHandle(AIS_InteractiveContext) &aCTX, Quantity_NameOfColor aColor, const occHandle(MeshVS_DataSource) &aMeshDS);
};

#endif // OCCPREGLWIDGET_H
