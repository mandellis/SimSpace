#ifndef SIMULATIONMANAGER_H
#define SIMULATIONMANAGER_H

//! ---
//! Qt
//! ---
#include <QWidget>
#include <QTreeView>
#include <QItemSelectionModel>
#include <QItemSelection>
#include <QList>

#include <AIS_InteractiveContext.hxx>
#include <AIS_Shape.hxx>
#include <TColStd_ListOfAsciiString.hxx>
#include <TopoDS_Compound.hxx>

//! -----------------------------------
//! for the experimental custom mesher
//! -----------------------------------
#include <TopoDS_Face.hxx>

//! ----
//! C++
//! ----
#include <vector>

//! ----------------
//! custom includes
//! ----------------
#include "simulationnodeclass.h"
#include <simulationdatabase.h>
#include "writesolverfileclass.h"
#include "writelabelclass.h"
#include "serializerclass.h"
#include "deserializerclass.h"
#include "listofmesh.h"
#include "postobject.h"
#include "frdreader.h"
#include "postengine.h"
#include "qoccprogressindicator.h"
#include "solutioninfo.h"
#include <indexedmapofmeshdatasources.h>
#include "detailviewer.h"
#include <userMessage.h>
#include <qhistogramdata.h>
#include "resultpresentation.h"

#ifdef COSTAMP_VERSION
#include <QProcess>
#endif

class inputFileGenerator;
class MeshingServer;
class QStandardItemModel;
class QStandardItem;
class QMenu;
class QTimer;

class CCXSolverManager;

class SimulationManager: public QWidget
{
    Q_OBJECT

public:

    QTreeView *myTreeView;

private:

    QStandardItemModel *myModel;
    QMenu *myContextMenu;

private:

    QStandardItem *myActiveAnalysisBranch;

    QStandardItem *RootItem;
    QStandardItem *Geometry_RootItem;
    QStandardItem *CoordinateSystems_RootItem;
    QStandardItem *RemotePoint_RootItem;
    QStandardItem *Mesh_RootItem;
    QStandardItem *Connections_RootItem;
    QStandardItem *NamedSelection_RootItem;
    QStandardItem *ThermalAnalysis_RootItem;
    QStandardItem *StaticAnalysis_RootItem;
    //QStandardItem *Results_RootItem;

    serializerClass *mySerializer;
    deserializerClass *myDeserializer;
    postEngine *myPostEngine;
    inputFileGenerator* myInputFileGenerator;

    //! experimental
    CCXSolverManager *theCCXSolverManager;

    //! experimental
    FrdReader *myFrdReader;


    //! internal timer
    QTimer *myTimer;

    //! meshing status
    bool myIsMeshingRunning;

    //! simulation status
    bool myIsCalculationRunning;

    //! detail viewer
    DetailViewer *myDetailViewer;

    //! current project directory
    //QString myCurrentProjectDir;

    //! root of the current running analysis
    QStandardItem *myCurrentRunningAnalysis;

private:

    //! the interactive context
    occHandle(AIS_InteractiveContext) myCTX;

    //! the text writer
    writeLabelClass *theTextWriter;

    //! get value of a property
    template<class T>
    inline T getPropertyValue(SimulationNodeClass *aNode, const QString &propertyName)
    {
        return aNode->getPropertyItem(propertyName)->data(Qt::UserRole).value<Property>().getData().value<T>();
    }

private:

    //! create the content
    void createContent();

    //! transfer mesh nodes
    void transferMeshNodes();

    //! get geometry selection
    ListOfShape getGeometrySelection() const;

    //! the simulation data base
    simulationDataBase *mySimulationDataBase;

    //! the selection model
    QItemSelectionModel *mySelectionModel;

    //! the meshing server
    MeshingServer *myMeshingServer;

private:

    void getTreeItemsRecursively(QStandardItemModel* model, QList<QExtendedStandardItem*> &items, QModelIndex parent = QModelIndex());
    //int calculateStartColumn() const;

    TopoDS_Shape fromTagToShape(const GeometryTag &aTag);
    TopTools_ListOfShape fromTagToShape(const QVector<GeometryTag> &vecLoc);

    //! clear generated data
    void clearGeneratedData();

    //! rename item based on definition
    void renameItemBasedOnDefinition();

    //! delete all children items
    void deleteAllChildrenItems();

    //! delete data sources from model items
    void deleteDataSourcesFromModel();

private slots:

    //! on mouse enter - UNUSED
    void onMouseEnter() { cout<<"++++"<<endl; myTreeView->resizeColumnToContents(0); }

    //! insert an item into the menu
    void handleItem(int type);

    //! mesh control item changed (suppressed, removed, modified):
    //! the corresponding mesh should be invalidated
    void handleMeshItemChange(QStandardItem *item);

    //! experimental
    void handleItemChange(QStandardItem*);

    //! show the context menu
    void showContextMenu(const QPoint& pos);

    //! build custom menu
    void buildCustomMenu(const QModelIndex &modelIndex);

    //! delete an item from the menu (not a root item)
    void deleteItem(QList<QModelIndex> indexesList = QList<QModelIndex>());

    //! highlighter
    void highlighter(QModelIndex modelIndex=QModelIndex());

    //! prepare for meshing
    TopTools_ListOfShape prepareForMeshing();

    //! build the mesh
    void buildMesh(bool isVolumeMesh);

    //! update the mesh statistics
    void updateMeshStatistics();

    //! suppress-unsuppress an item
    void changeNodeSuppressionStatus(Property::SuppressionStatus newSuppressionStatus);

    //! unsuppress all the geometry items (geometry nodes)
    void unsuppressAllBodies();

    //! invert suppression set
    void invertSuppressionSet();

    //! handle visibility change coming from detail viewer
    void handleVisibilityChange(bool newIsVisible);

    //! change the element control
    void ChangeElementControl();

    //! change the scoping method
    //void changeScopingMethod();

    //! change scope color
    //void changeColor();

    //! create local simulation node
    //void createSimulationNode(SimulationNodeClass::nodeType type, QVariant addOptions=QVariant());

    //! update the solver type
    //void updateSolverType();

    //! update node model
    //void updateNodeModel();

    //! rename item
    void renameItem();

    //! update a node name, using the name of the item
    void updateNodeName();

    //! experimental - to do
    //! create boundary condition text descriptor
    QString createItemDescriptor() const;

    void updateTimer();


public slots:

    //! write solver input file (CCX) - initialized with a default parameter
    void writeSolverInputFile();

    //! start analysis
    bool startAnalysis(const QString &projectDataDir);

    //! update load tables
    void resizeTabularData();

    //! handle tabular data
    void HandleTabularData();

    //! handle load definitions
    void handleFilmCoefficientLoadDefinitionChanged(const QString &textData);
    void handleReferenceTemperatureLoadDefinitionChanged(const QString &textData);
    void handleLoadMagnitudeDefinitionChanged(const QString &textData);
    void handleLoadXDefinitionChanged(const QString &textData);
    void handleLoadYDefinitionChanged(const QString &textData);
    void handleLoadZDefinitionChanged(const QString &textData);

    //! configure post engine
    void configurePostEngine();

    //! start post engine
    void startPostEngine();

    //! configure and start
    void configureAndStartPostEngine();

    //! compute and display mesh metric
    void computeAndDisplayMeshMetric();

public:

    //! ------------------------
    //! constructor - default -
    //! ------------------------
    SimulationManager(QWidget *parent=0);

    //! --------------
    //! constructor I
    //! --------------
    SimulationManager(const occHandle(AIS_InteractiveContext)& aCTX, QWidget *parent=0);

    //! -----------
    //! destructor
    //! -----------
    virtual ~SimulationManager()
    {
        cout<<"SimulationManager::~SimulationManager()->____DESTRUCTOR CALLED____"<<endl;
#ifdef COSTAMP_VERSION
        //! --------------------------------------------
        //! closing all the time step builder instances
        //! --------------------------------------------
        for(int i=0; i<timeStepBuilderPids.size(); i++)
        {
            cout<<"SimulationManager::~SimulationManager()->____KILLING TIME STEP BUILDER "<<i<<"____"<<endl;
            QProcess::execute(QString("taskkill /PID %1").arg(timeStepBuilderPids[i]));
        }
#endif
    }

    //! set the context
    void setContext(const occHandle(AIS_InteractiveContext) &aCTX);

    //! set the selection model
    void setSelectionModel();

    //! set the database
    void setDataBase(simulationDataBase *aDB);

    //! load the CAD model
    bool loadCADModel(const QString &fileName, TopoDS_Compound &shapeFromReader,
                      QList<QString> &listOfNames,
                      const occHandle(QOccProgressIndicator) &aProgressIndicator);

    //! create simulation data base
    void createSimulationDataBase(const TopoDS_Shape &shapeFromReader,
                                  const QString &fileName,
                                  const QList<QString> &listOfNames);

    //! create an empty simulation data base
    void createSimulationDataBaseEmpty();

    //! get the database
    simulationDataBase* getDataBase() { return mySimulationDataBase; }

    //! build database from disk
    void buildDataBaseFromDisk(const QString &fileName);

    //! get the current node
    SimulationNodeClass* getCurrentNode();

    //! get the model
    QStandardItemModel *getModel() const { return myModel; }

    //! get status
    bool isSomethingRunning()
    {
        if(myIsCalculationRunning == true || myIsMeshingRunning == true) return true;
        return false;
    }

    //! ----------------------
    //! set the detail viewer
    //! ----------------------
    void setDetailViewer(DetailViewer *theDetailViewer) { myDetailViewer = theDetailViewer; }

    //! ---------------------------------------
    //! set and get the active analysis branch
    //! ---------------------------------------
    void setTheActiveAnalysisBranch();
    QStandardItem* getActiveAnalysisBranch();

    //! setCurrentProjectDir
    //void setCurrentProjectDir(QString aDir) { myCurrentProjectDir = aDir;}

protected:

    //! event filter
    bool eventFilter(QObject *object, QEvent *event);

public slots:

    //! -----------------------
    //! create simulation node
    //! -----------------------
    void createSimulationNode(SimulationNodeClass::nodeType type, QVariant addOptions=QVariant());

    //! -----------------------
    //! handle meshing results
    //! -----------------------
    void handleMeshingResults(bool isMeshingSuccessfull);

    //! change scope color
    void changeColor();

    //! clear the tree
    void clearTree()
    {
        myModel->clear();
    }

    //! ----------------------------------------
    //! for the connection with the meshToolBar
    //! ----------------------------------------
    void buildVolumeMesh() { this->buildMesh(true); }
    void buildSurfaceMesh() { this->buildMesh(false); }

    //! -----------------
    //! export STEP file
    //! -----------------
    void exportSTEPFile();
    void exportBREPFile();

signals:

    //! working mode changed
    void requestSetWorkingMode(int);

    //! request mesh invalidate
    void requestMeshInvalidate(std::vector<int> meshIndexes);

    //! request for mesh generation
    void generateMesh(bool isVolume);

    //! clear mesh
    void requestClearMesh();

    //! clear the property viewer
    void clearPropertyViewer();

    //! ------------------------------------------------------------------------
    //! request custom color: use with default arguments to clean custom colors
    //! and update the viewer
    //! ------------------------------------------------------------------------
    void requestCustomColor(const TopTools_ListOfShape &shapes = TopTools_ListOfShape(),
                            Quantity_NameOfColor = Quantity_NameOfColor(),
                            bool clearPrevious = true,
                            bool updateViewer = true);

    //! request a selection mode - the request is sent to the mainWindow,
    //! by consequence the buttons status is automatically synchronized
    void request3DBodySelectionMode(bool);
    void request2DBodySelectionMode(bool);
    void request1DBodySelectionMode(bool);
    void request0DBodySelectionMode(bool);

    void requestRemoveObsoleteMeshes();
    void requestBuildMeshIOs();
    void requestHideMeshes();
    void requestHideSlicedMeshes();
    void requestShowMeshes(bool areMeshNodesVisible);
    void requestShowBody(TColStd_ListOfInteger ListOfBodyNumbers);
    void updatedetailViewer(QModelIndex);
    void requestHighlightBody(const QList<int> listOfBodies);
    void requestUnhighlightBodies(bool updateViewer);
    void requestHideBody(TColStd_ListOfInteger listOfBodies);
    void requestShowAllBodies();
    void requestReactivateCurrentStandardSelectionMode();
    void requestDisplayShapeCopy(const TopTools_ListOfShape &list1, const TopTools_ListOfShape &list2,
                                 Quantity_NameOfColor color1, Quantity_NameOfColor color2, QVariant options = QVariant());

    void requestDisplayShapeCopy1(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color);
    void requestDisplayShapeCopy2(const TopTools_ListOfShape &listShapes, Quantity_NameOfColor color);

    /*
    void requestDisplayMaster(const TopTools_ListOfShape &list1, const TopTools_ListOfShape &list2,
                                 Quantity_NameOfColor color1, Quantity_NameOfColor color2, QVariant options = QVariant());

    void requestDisplaySlave(const TopTools_ListOfShape &list1, const TopTools_ListOfShape &list2,
                                 Quantity_NameOfColor color1, Quantity_NameOfColor color2, QVariant options = QVariant());
                                 */

    void requestTabularData(QModelIndex index);

    //! this signal is emitted for activating the "Shape selector"
    void requestStartEditingScope();

    void requestStartEditingMagnitude();
    void requestStartEditingXcomponent();

    //! request the viewer to show a trihedron
    void requestDisplayTrihedron(QVector<double>,QVector<QVector<double>>,int size=0);

    //! -------------------
    //! tabular data stuff
    //! -------------------
    void requestShowColumns(QList<int> columnsToShow);
    void requestShowGraph(CustomTableModel *tabData, QList<int> columnsToShow);
    void requestClearGraph();
    void requestHideAllColumnsAndHeaders();     //! unused
    void requestHideFirstRow();
    void requestShowFirstRow();

    void requestCreateColorBox(double min, double max, int Nintervals);
    //void requestDisplayResult(const postObject& aPostObject);
    void requestDisplayResult(postObject& aPostObject);
    //void requestDisplayResult(postObject aPostObject);
    void requestHideAllResults();
    void requestHideSingleResult(const postObject& aPostObject);
    void requestSetActiveCentralTab(const QString& widgetName);
    void requestUpdateConvergenceViewer(const QList<solutionInfo> &solutionInfoList);

    //! ----------------------------------------------------
    //! display mesh quality data using a mainwindow widget
    //! ----------------------------------------------------
    void requestDisplayMeshMetricHystogram(const histogramData &hysData);

private:

    QList<QStandardItem *> ItemListFromListOfShape(TopTools_ListOfShape *listOfShapes);
    //QExtendedStandardItem* getTreeItem(nodeType theNodeType);
    int getInsertionRow() const;

    //SimulationNodeClass* getAnalysisSettingsNodeFromCurrentItem() const;
    //QExtendedStandardItem *getAnalysisSettingsItemFromCurrentItem() const;

    void duplicateItem(QExtendedStandardItem *item=Q_NULLPTR);
    void swapContact();

public:

    QExtendedStandardItem* getTreeItem(SimulationNodeClass::nodeType theNodeType);
    QList<QExtendedStandardItem *> getAllTreeItemOfType(SimulationNodeClass::nodeType theNodeType);
    SimulationNodeClass* getAnalysisSettingsNodeFromCurrentItem() const;
    QExtendedStandardItem *getAnalysisSettingsItemFromCurrentItem() const;
    //int calculateStartColumn() const;

public slots:

    //! synchronize the "Visible" control with the OCC viewer
    void synchVisibility();

    //! save to disk
    void saveSimulationDataBase(const QString &savingDir, const QString &fileName);

    //! node list builder
    QList<SimulationNodeClass*> nodeListBuilder(const QString &savedFilesDir);

    //! new definition of the "Scope" property
    //void updateTags();

    //! start interpolation
    void interpolate();

    //! interpolate one time
    //void interpolatePrivate(QString filePath = QString(), bool multiInterpolate=false);

    //! --------------------------
    //! interpolate several times
    //! --------------------------
    void interpolatePrivate(int mode);

    //! translate open foam scalar data
    bool translateOpenFoamScalarData();

    //! handle global mesh control change
    void handleGlobalMeshControlChange();

    //! handle solution component change: a solution component has been changed.
    //! Remove the data from the item and clear the colored mesh
    void handleSolutionComponentChanged();

    //! show healing elements
    void showHealingElements();

    //! update results presentation
    void updateResultsPresentation();

    //! read results file
    void readResultsFile(const QString &fileName, const QString &solutionDataDir);

    //! experimental
    void customMesherBuildFaceMesh(const TopoDS_Face &aFace);

private:

    //! retrieve the result contained into an item in the form of MeshVS_Mesh object
    bool retrieveCurrentItemResult(postObject &aPostObject);
    QList<postObject> retrieveAllResults();

    //! --
    void callPostEngineEvaluateResult();

    //! evaluate the result in post item
    //void callPostEngineEvaluateResult_private(SimulationNodeClass *curNode, bool immediatelyDisplay=true);
    void callPostEngineEvaluateResult_private(QStandardItem *curItem, bool immediatelyDisplay=true);

    //! evaulate all results
    void evaluateAllResults();

    //! update remote point absolute coordinates
    void updateRemotePointAbsCoordinates();

    //! create automatic connections
    void createAutomaticConnections();

private slots:

    //! retrieve solver info
    void retrieveSolverInfo();

#ifdef COSTAMP_VERSION
    //! start time step builder
    void COSTAMP_startTimeStepBuilder();

private slots:

    //! change the time history file location
    void COSTAMP_changeTimeHistoryFile(const QString &fileLocation);

    //! add process parameters
    bool COSTAMP_addProcessParameters() ;

private:

     //! instances of the time step builder
     std::vector<qint64> timeStepBuilderPids;

     //! times selected from time step builder
     std::vector<double> tSbList;
#endif

signals:

    //! experimental - request an update of the mesh context
    void requestUpdateMeshView();

    //! handle double view port for contacts
    void requestShowDoubleViewPort(bool isVisible);

    //! hide a list of bodies from the contact/master viewer
    void requestHideBodiesFromViewers(const TColStd_ListOfInteger &bodyListNbs);

    //! display a list of bodies on the master view port
    void requestShowBodyOnMasterViewPort(const TColStd_ListOfInteger &listOfActiveBodies, const TColStd_ListOfInteger &listOfNotActiveBodies);

    //! display a list of bodies on the slave view port
    void requestShowBodyOnSlaveViewPort(const TColStd_ListOfInteger &listOfActiveBodies, const TColStd_ListOfInteger &listOfNotActiveBodies);

    //! request display a spherical marker (i.e. for showing a remote point)
    void requestDisplaySphericalMarker(gp_Pnt P);

    //! request hide all markers
    void requestHideAllMarkers(bool updateViewer=true);

    //! request refresh viewer
    void requestRefreshViewer();

public:

    //! build and preview the prismatic layers
    bool previewPrismaticLayer();

    //! show a mesh - utility
    void displayFaceMesh(const occHandle(MeshVS_DataSource) &aMeshDS,
                                         Quantity_NameOfColor aColorName = Quantity_NOC_AQUAMARINE1,
                                         bool showMeshEdges = true);

private:

    //! show marker if any
    void displayMarker();

    //! experimental
    void buildMeshIO();

    //! experimental
    void resetAndUpdateModel();

    //! experimental
    void generateBoundaryConditionsMeshDS(bool computeDual);

    //! replicate bolts
    void replicateBolt();

    //! add parent time tag
    void addParentTimeTag(SimulationNodeClass *aNode);

signals:

    //! request update viewport
    void requestUpdateViewport();

    //! show mesh data sources
    void requestShowMeshDataSources(const IndexedMapOfMeshDataSources &meshDataSources);

    //! set the clipper data base
    void requestSetClipperDataBase(meshDataBase *theDB);

    //! display text on console
    void requestDisplayTextOnConsole(const QString &text);

    //! request reset custom colors
    void requestResetCustomColors(bool updateViewer = true);

    //! request apply custom colors
    void requestApplyCustomColor(const QMap<GeometryTag,TopoDS_Shape> &subShapesMap, Quantity_NameOfColor aColor, bool updateViewer = true);
};

#endif // SIMULATIONMANAGER_H
