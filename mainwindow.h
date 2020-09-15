#ifndef MAINWINDOW_H
#define MAINWINDOW_H

//! -------------------------------------------
//! redefinition of opencascade Handle() macro
//! -------------------------------------------
#include "occhandle.h"

//! ----
//! OCC
//! ----
#include <OCCLicense_Activate.hxx>

//! ----------------
//! custom includes
//! ----------------
#include "mydefines.h"
#include "actions3d.h"
#include "displayquality.h"
#include "displaymode.h"
#include "workingmode.h"
#include "myenumvariables.h"
#include "openfoamreader.h"
#include "optionsWidget/optionswidget.h"
#include "qprogressevent.h"
#include "global.h"
#include <qhistogram.h>

//! ----
//! OCC
//! ----
#include <AIS_InteractiveContext.hxx>

//! ---
//! Qt
//! ---
#include <QMainWindow>
#include <QAction>

class occGLWidget;
class occPreGLWidget;
class occPostWidget;

class QToolBar;
class QComboBox;
class QLabel;
class QTextEdit;
class QTabWidget;
class QPlainTextEdit;
class QDockWidget;
class QTableView;
class QPushButton;
class QPushButtonExtended;
class geometryManager;
class SimulationManager;
class MeshManager;
class SimulationManager;
class DetailViewer;
class TableWidget;

class TabularDataViewerClass1;

class QProgressBar;
class systemConsole;
class QTabWidgetExtended;

class ConvergenceDataChart;
class ConvergenceDataChart1;

class QProgressIndicator;
class ResultsToolBar;
class MeshToolBar;
class dockableViewPort;
class clipTool;
class memoryProfiler;
class userMessagesModel;
class QHistogram;

namespace Ui
{
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    //! constructor
    explicit MainWindow(QWidget *parent = 0);

    //! destructor
    ~MainWindow();

private:

    Ui::MainWindow *ui;

    //! --------------
    //! debug console
    //! --------------
    systemConsole *myDebugConsole;

    //! -------------------------------
    //! memory usage: this is dockable
    //! -------------------------------
    memoryProfiler *myMemoryProfiler;

    //! ------------------------------------------------------------
    //! convergence data viewer - it displays the convergence trend
    //! ------------------------------------------------------------
    ConvergenceDataChart *myConvergenceDataChart;
    ConvergenceDataChart1 *myConvergenceDataChart1;

    //! ---------------------
    //! the worksheet viewer
    //! ---------------------
    systemConsole *mySimulationMonitor;

    //! -------------------
    //! mesh metric viewer
    //! -------------------
    QHistogram *myHistogram;

    //! ----------------------------------------------------------
    //! user messages
    //! results of the simulation, of the meshing operations, ...
    //! ----------------------------------------------------------
    userMessagesModel* myUserMessages;

    //! -------------
    //! progress bar
    //! -------------
    QProgressBar *myProgressBar;

    //! -------------------------------------
    //! File, Geometry, Tools, Solution menu
    //! -------------------------------------
    QMenu *FileMenu;
    QMenu *GeometryMenu;

    QMenu *MeshMenu;
    QMenu *MeshMenuInsert;

    QMenu *ToolsMenu;
    QMenu *SolutionMenu;

    QMenu *ViewMenu;
    QMenu *viewSubMenuDisplayQuality;
    QMenu *viewSubMenuWindows;

    //! ------------------------
    //! The mainwindow toolbars
    //! ------------------------
    QToolBar* projectToolbar;
    QToolBar* viewAndSelectionToolbar;
    QToolBar* transformationsToolBar;
    QToolBar* solutionToolBar;
    MeshToolBar *meshToolBar;
    ResultsToolBar *resultsToolBar;

    //! -----------------------
    //! The permanent messages
    //! -----------------------
    QLabel *statusLabel;

    //! ----------------
    //! Graphic windows
    //! ----------------
    //occGLWidget *myMainOCCViewer;
    //occPreGLWidget *myMainOCCViewer;
    occPostWidget *myMainOCCViewer;
    dockableViewPort *myDockableMasterViewPort;
    dockableViewPort *myDockableSlaveViewPort;

    //! -----------------------------
    //! actions: "File" menu actions
    //! -----------------------------
    QAction *actionOpenProject;
    QAction *actionSaveProject;
    QAction *actionSaveProjectAs;
    QAction *actionImport;
    QAction *actionClose;
    QAction *actionExit;

    //! -----------------------------
    //! actions: "View" menu actions
    //! -----------------------------
    QActionGroup *groupOfActions_DisplayModes;
    QAction *actionShadedExteriorAndEdges;
    QAction *actionShadedExterior;
    QAction *actionWireframe;

    //! --------------------------------------------
    //! "View" subMenu  - "Display quality" actions
    //! --------------------------------------------
    QActionGroup *displayQualityActionGroup;
    QAction *actionHigh;
    QAction *actionMedium;
    QAction *actionCoarse;

    QAction *actionShowTabularData;
    QAction *actionShowSimulationManager;
    QAction *actionShowDetailViewer;
    QAction *actionShowGraphViewer;
    QAction *actionShowDebugWindow;
    QAction *actionShowSectionPlanes;
    QAction *actionShowMemoryProfiler;
    QAction *actionShowPhytonConsole;

    //! -----------------------------
    //! actions: "Mesh" menu actions
    //! -----------------------------
    QAction* actionPreviewMesh;
    QAction *actionGenerateMesh;
    QAction *actionLoadMesh;
    QAction *actionClearMesh;

    //! -----------------------------------------------
    //! "Mesh" subMenu - "Insert mesh control" actions
    //! -----------------------------------------------
    QAction *actionInsertMethod;
    QAction *actionInsertBodySizing;
    QAction *actionInsertFaceSizing;
    QAction *actionInsertEdgeSizing;
    QAction *actionInsertVertexSizing;
    QAction *actionInsertPrismaticLayer;

    //! -------------------
    //! "Solution" actions
    //! -------------------
    QAction *actionSolve;
    QAction *actionWriteInputFile;
    QAction *actionReadResults;

    //! -----
    QAction *actionImageToFile;
    QAction *actionColorBox;

    //! --------------------------
    //! "View operations" actions
    //! --------------------------
    QAction* actionToggle3Drotation;
    QAction* actionTogglePan;
    QAction* actionToggleWindowZoom;
    QAction* actionFitAll;

    //! --------------------------------
    //! Used by the "Context menu" only
    //! --------------------------------
    QActionGroup *groupOfViewActions;

    //! ---------------------------------
    //! action: "Selection mode" actions
    //! ---------------------------------
    QAction* actionToggleSolidSelect;
    QAction* actionToggleFaceSelect;
    QAction* actionToggleEdgeSelect;
    QAction* actionToggleVertexSelect;
    QAction* actionTogglePickPointCoordinates;

    //! -------------------------
    //! action type of selection
    //! -------------------------
    QAction *actionSingleSelect;
    QAction *actionMultipleSelect;
    QAction *actionSelectAll;

    //! ------------------------
    //! action extend selection
    //! ------------------------
    QAction *actionExtendToAdjacent;
    QAction *actionExtendToLimits;
    QAction *actionExtendToConnections;

    //! --------------------
    //! action translations
    //! --------------------
    QAction* actionTranslationX;
    QAction* actionTranslationY;
    QAction* actionTranslationZ;
    QAction* actionRotationX;
    QAction* actionRotationY;
    QAction* actionRotationZ;
    QAction* actionDeleteTransformation;

    //! --------
    //! options
    //! --------
    QAction *actionOptions;

    //! ---------------------------------------------------
    //! The selection combobox - single or multiple select
    //! ---------------------------------------------------
    QPushButtonExtended *myPushButtonTypeOfSelection;
    QPushButtonExtended *myImageButton;

    //! ----------------------
    //! Extend selection type
    //! ----------------------
    QPushButton* myExtendSelectionButton;

    //! --------------------------------
    //! Import status - status variable
    //! --------------------------------
    geometryImportStatus myGeometryImportStatus;

    //! -------------------------------------------
    //! the current project name - status variable
    //! -------------------------------------------
    QString myCurrentProjectName;

    //! ------------------------------------
    //! Working directory - status variable
    //! ------------------------------------
    QString myWorkingDir;

    //! -----------------------
    //! the simulation manager
    //! -----------------------
    SimulationManager *mySimulationManager;

    //! ------------------------------------------
    //! simulation manager non dockable container
    //! ------------------------------------------
    QTabWidget *mySimulationManagerContainer;

    //! ------------------
    //! the detail viewer
    //! ------------------
    DetailViewer *myDetailViewer;

    //! -------------------------------------
    //! detail viewer non dockable container
    //! -------------------------------------
    QTabWidget *myDetailViewerContainer;

    //! -----------------------------
    //! dock container for hystogram
    //! -----------------------------
    QDockWidget *myHistogramDockContainer;

    //! ------------------------------------------
    //! graph viewer - loads vs time using charts
    //! ------------------------------------------
    TabularDataViewerClass1 *myTabularDataGraphViewer1;

    //! -------------
    //! tabular data
    //! -------------
    TableWidget *myTabularDataView;

    //! -----------------------
    //! the progress indicator
    //! -----------------------
    QProgressIndicator *myProgressIndicator;

    //! ----------
    //! clip tool
    //! ----------
    clipTool *myClipTool;

    //! -----------------------
    //! the central tab widget
    //! -----------------------
    QTabWidgetExtended *myCentralTabWidget;

    //! -------------
    //! dock widgets
    //! -------------
    QDockWidget *myTabularDataDock;
    QDockWidget *myTabularDataViewerDock;
    QDockWidget *mySimulationManagerDock;
    QDockWidget *myDetailViewerDock;
    QDockWidget *myDebugConsoleDock;
    QDockWidget *myClipToolDock;

private:

    //! event
    virtual bool MainWindow::event(QEvent *event);

    //! create the menus
    void createMenu();

    //! create dock widgets
    void createDockWidgets();

    //! creates the actions
    void createActions();

    //! create the toolbars
    void createToolBars();

    //! create the connections
    void setUpConnections();

    //! this handle the logic of the mainwindos buttons (automatic release-enable/press-disable)
    void handleViewAndSelectionButtons(QAction *action);

    //! write solver input file
    //void writeSolverInputFile();

private slots:

    //! create item mesh method - wrapper
    void createItemMeshMethod();

    //! create item mesh body sizing - wrapper
    void createItemMeshBodySizing();

    //! create item mesh face sizing - wrapper
    void createItemMeshFaceSizing();

    //! create item mesh edge sizing - wrapper
    void createItemMeshEdgeSizing();

    //! create item mesh vertex sizing - wrapper
    void createItemMeshVertexSizing();

    //! create item mesh prismatic layer - wrapper
    void createItemMeshPrismaticLayer();

    //! Handle the "Select all" request from the toolbar
    void HandleSelectAll();

    //! handle color bx visibility
    void handleColorBoxVisibility(bool isVisible);

    //! close application
    void myClose();

    //! save project
    void saveProject();

    //! save project as - experimental
    void saveProjectAs();

    //! open project - experimental
    void openProject();

    //! import geometry file - dialog
    void importFile(QString &fileName = QString(""));

    //! write solver input file - dialog
    void writeSolverInputFile();

    //! closes the current session
    void closeTheModel();

    //! activates the rotation mode
    void toggle3Drotation(bool isActivated);

    //! activates the pan mode
    void togglePan(bool isActivated);

    //! activates the box zoom
    void toggleWindowZoom(bool isActivated);

    //! fit the model
    void fitAll(bool isActivated);

    //! activate the TopoDS_Solid standard selection mode
    void toggleSolidSelectionMode(bool isActivated);

    //! activates the TopoDS_Face standard selection mode
    void toggleFaceSelectionMode(bool isActivated);

    //! activates the TopoDS_Edge standard selection mode
    void toggleEdgeSelectionMode(bool isActivated);

    //! activates the TopoDS_Vertex standard selection mode
    void toggleVertexSelectionMode(bool isActivated);

    //! activates the point coordinates picking mode
    void togglePointCoordinatesPickingMode(bool isActivated);

    //! extend to adjacent
    void callExtendSelectionToAdjacent();
    void callExtendSelectionToLimits();
    void callExtendSelectionToConnections();

    //! display quality - method calling public slot of the occPreGLWidget
    void setHighQualityDisplay();

    //! display quality - method calling public slot of the occPreGLWidget
    void setMediumQualityDisplay();

    //! display quality - method calling public slot of the occPreGLWidget
    void setLowQualityDisplay();

    //! set the working mode
    void setWorkingMode(int workingModeNumber);

    //! add a coordinate system transformation
    void addCoordinateSystemTranslationX();
    void addCoordinateSystemTranslationY();
    void addCoordinateSystemTranslationZ();
    void addCoordinateSystemRotationX();
    void addCoordinateSystemRotationY();
    void addCoordinateSystemRotationZ();
    void removeTransformation();

    void setSelectionModeBox();
    void setSelectionModeSingle();
    void setSelectionModeGeometry();
    void setSelectionModeMesh();

    void showOptionsWidget();

    void readResultsFile();

public slots:

    //! start analysis
    void startAnalysis();

    //! handle working mode menu - a temporary solution
    //void checkWorkingModeItem(curWorkingMode curWorkMode);

    //! handle the view mode menu - a temporary solution
    void checkViewModeItem(CurDisplayMode theCurDisplayMode);

    //! reactivate all the selection buttons
    void enableSelectionButtons();

    //! disable all the selection buttons
    void disableSelectionButtons();

    //! disable/enable a button
    void disableButton(const QString &buttonName, bool disable=true);

    //! double click on shape selector
    void startEditingScope();

    //! start editing the magnitude
    void startEditingMagnitude();

    //! start editing Xcomponent
    void startEditingXcomponent();

    //! experimental
    void setDoubleViewPortVisible(bool isVisible);

    //! experimental
    void clearMasterSlaveViewPorts(const TColStd_ListOfInteger &bodyListNbs);

    //! experimental
    void showBodiesOnMasterViewPort(const TColStd_ListOfInteger &listOfActiveBodies, const TColStd_ListOfInteger &listOfNotActiveBodies);

    //! experimental
    void showBodiesOnSlaveViewPort(const TColStd_ListOfInteger &listOfActiveBodies, const TColStd_ListOfInteger &listOfNotActiveBodies);

    //! update viewport
    void updateViewport();

    //! display mesh metric
    void displayMeshMetric(const histogramData& hysData);

public:

    //! get the interactive context
    const occHandle(AIS_InteractiveContext)& getContext() const;

    //! get the simulation manager
    SimulationManager* getSimulationManager() { return mySimulationManager; }

    //! get the progress indicator
    QProgressIndicator* getProgressIndicator() { return myProgressIndicator; }

protected:

    //! close event
    void closeEvent(QCloseEvent *event);

private:

    //! create a working dir
    bool createWorkingDir();

private slots:

    //! update the settings file
    void updateSettings();

signals:

    //! extend selection
    void requestExtendSelectionToAdjacent();
    void requestExtendSelectioToLimits();
    void requestExtendSelectionToConnections();
};

#endif // MAINWINDOW_H
