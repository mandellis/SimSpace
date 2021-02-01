//! -------
//! quazip
//! -------
#include <ext/quazip/inc/quazip.h>
#include <ext/quazip/inc/JlCompress.h>

//! ----------------
//! custom includes
//! ----------------
#include "ui_mainwindow.h"
#include "mydefines.h"
#include "mainwindow.h"
#include "occPreGLwidget.h"
#include "occpostwidget.h"
#include "simulationdatabase.h"
#include "simulationmanager.h"
#include "simulationnodeclass.h"
#include "detailviewer.h"
#include "shapeselector.h"
#include "generaldelegate.h"

#include "src/gui/tabularData/tabulardataviewerclass1.h"

#include "src/gui/tabularData/tablewidget.h"
#include "qpushbuttonextended.h"
#include "src/utils/tools.h"
#include "optionsWidget/optionswidget.h"
#include "systemConsole/systemconsole.h"
#include "src/utils/ccout.h"
#include "qconsoleevent.h"
#include "qtabwidgetextended.h"

#include "convergencedatachart.h"

#include "systemConsole/SimulationMonitor.h"
#include "qprogressindicator.h"
#include "qoccprogressindicator.h"
#include "resultstoolbar.h"
#include "meshtoolbar.h"
#include "dockableviewport.h"
#include "src/utils/cliptool/clipTool.h"
#include "src/memory/memoryprofiler.h"
#include <qhistogram.h>
#include <src/utils/tools.h>
//! ----
//! C++
//! ----
#include <iostream>
#include <fstream>

//! ---
//! Qt
//! ---
#include <QTime>
#include <QBitmap>
#include <QFileDialog>
#include <QToolBar>
#include <QComboBox>
#include <QPushButton>
#include <QLabel>
#include <QMessageBox>
#include <QDockWidget>
#include <QPlainTextEdit>
#include <QTabWidget>
#include <QProgressBar>
#include <QDebug>

using namespace std;

//! -----------------------------
//! function: thirdPartyActivate
//! details:
//! -----------------------------
static void thirdPartyActivate()
{
    cout<<"thirdPartyActivate()->____Activating third party products____"<<endl;

    //! activate OCC products                                       //
    //Standard_CString emeshActivationKey ="        PUZZLE_DIE         0        0      0    0 3f85a38fd58cab86c9f3831e6025957ba6f63394b841eaf92f88da8e13ddd1f7 Co.Stamp s.r.l. / Master key";
    //Standard_CString omfActivationKey =  "          PUZZLE_DIE         0        0      0    0 8d19e677b16fb8e2e0c2301df07656dee5d2ba44cc03eebeeeca55d64ed3ed2c Co.Stamp s.r.l. / Master key";

    Standard_CString emeshKey ="        PUZZLE_DIE         0        0      0    0 5ba43f346b9dfbd4ece27b57443b309961676c9bb98db05e041047c22e93a907 Co.Stamp s.r.l. / Master key";
    Standard_CString omfKey =  "          PUZZLE_DIE         0        0      0    0 613699cc791fa76a08d0c5f8534c2db8eb89565b94cfcdcbc600ba4639221061 Co.Stamp s.r.l. / Master key";

    OCCLicense_Activate("EMESH-7.3",emeshKey);
    OCCLicense_Activate("OMF-7.3",omfKey);
    //char filePath[256];
    //sprintf(filePath,"D:/Work/Qt/pro_25.0_OCC7.3.0/License_09092019/License.txt");
    //OCCLicense_LoadKeyFile(filePath,true);

    //! -------------
    //! user message
    //! -------------
    userMessage m;
    m.isDone = true;
    m.message = "Activating OCC third party libraries";
    Global::status().myMessages->appendMessage(m);
}

//! ---------------------------------------------------------
//! function: createWorkingDir
//! details:  create a working directory and open a log file
//! ---------------------------------------------------------
bool MainWindow::createWorkingDir()
{
    //! Once the application has been opened it reads the application settings from
    //! the file "settings.txt", which for the moment is the directory "C:/ProgramData/WB"
    //! [1] If the directory "WB" does not exist, create it, then enter it.
    //! [2] Search for the file "settings.txt": if it does not exist create it and write
    //!
    //!     - Workdir   <path of the Workdir>
    //!     - ...
    //!     - ...
    //!
    //! [3] read the file "../../[]/settings.txt" and init the "myWorkingDir"
    //! the location where to store the settings.txt file: typically a system directory

    QDir::setCurrent(QString(SYSTEM_PROGRAM_DATA));
    QDir curDirS = QDir::current();
    curDirS.setPath(QDir::current().absolutePath());
    cout<<"MainWindow::createWorkingDir()->____application initial current directory: "<<curDirS.absolutePath().toStdString()<<"____"<<endl;

    if(curDirS.cd("WB")==false)
    {
        cout<<"MainWindow::createWorkingDir()->____WB directory not existing: creating____"<<endl;
        curDirS.mkdir("WB");
        curDirS.cd("WB");
    }
    cout<<"MainWindow::createWorkingDir()->____current directory: "<<curDirS.absolutePath().toStdString()<<"____"<<endl;

    QString settingsFilePath = curDirS.absolutePath()+"/"+"settings.txt";
    char settingsFilePathC[512];
    sprintf(settingsFilePathC,settingsFilePath.toStdString().c_str());

    cout<<"MainWindow::createWorkingDir()->____setting.txt file path: "<<settingsFilePathC<<"____"<<endl;

    FILE *settings = fopen(settingsFilePathC,"r");

    char wd[512],tmp[512];

    if((settings)==NULL)
    {
        cout<<"MainWindow::createWorkingDir()->____the settings file does not exist: creating____";
        settings = fopen(settingsFilePathC,"w");
        //! write a value into the file
        fprintf(settings,"workingdir\t%s\n",curDirS.absolutePath().append("/Workdir").toStdString().c_str());
        //! default bkg type: vertical gradient
        fprintf(settings,"defaultBkgGradient\t%d\n",int(2));
        //! default first color: (DEF_R1,DEF_G1,DEF_B1)
        fprintf(settings,"defaultFirstColor\t%d\t%d\t%d\n",DEF_R1,DEF_G1,DEF_B1);
        //! default second color: (DEF_R2,DEF_G2,DEF_B2)
        fprintf(settings,"defaultSecondColor\t%d\t%d\t%d\n",DEF_R2,DEF_G2,DEF_B2);
        fclose(settings);
    }
    else
    {
        cout<<"MainWindow::createWorkingDir()->____the settings file exists: reading settings____";
        fscanf(settings,"%s%s",tmp,wd);
        printf("working dir: %s\n",wd);
        fclose(settings);
    }

    settings = fopen(settingsFilePathC,"r");
    fscanf(settings,"%s%s",tmp,wd);
    fclose(settings);

    //! --------------------------------------
    //! set the application working directory
    //! --------------------------------------
    QDir::setCurrent(QString::fromLatin1(wd));
    myWorkingDir = QDir::current().absolutePath();

    //! -------------
    //! user message
    //! -------------
    userMessage m;
    m.isDone = true;
    m.message = QString("Working directory: \"")+myWorkingDir+"\" set";
    Global::status().myMessages->appendMessage(m);

    cout<<"MainWindow::createWorkingDir()->____Working dir: "<<myWorkingDir.toStdString()<<"____"<<endl;
    return true;
}

//! -----------------------------------
//! function: updateSettings
//! details:  update the settings file
//! -----------------------------------
void MainWindow::updateSettings()
{
    cout<<"MainWindow::updateSettings->____function called____"<<endl;
    QDir curDirS;
    curDirS.setCurrent(QString(SYSTEM_PROGRAM_DATA).append("/WB"));
    QString fs = curDirS.absolutePath().append("/settings.txt");
    FILE *settings = fopen(fs.toStdString().c_str(),"w");
    if(settings!=NULL)
    {
        int defaultBackGroundGradient = 2;
        int r1 = DEF_R1;
        int g1 = DEF_G1;
        int b1 = DEF_B1;
        int r2 = DEF_R2;
        int g2 = DEF_G2;
        int b2 = DEF_B2;

        fprintf(settings,"workingdir\t%s\n"
                         "defaultBkgGradient\t%d\n"
                         "defaultFirstColor\t%d\t%d\t%d\n"
                         "defaultSecondColor\t%d\t%d\t%d\n",
                myWorkingDir.toStdString().c_str(),
                defaultBackGroundGradient,
                r1,g1,b1,
                r2,g2,b2);
    }
    fclose(settings);
}

//! --------------------------
//! function: event
//! details:  handling events
//! --------------------------
bool MainWindow::event(QEvent *event)
{
    if (event->type() == QProgressEvent::type())
    {
        //cout<<"MainWindow::event()->____Progress event received____"<<endl;
        QProgressEvent *aProgressEvent = static_cast<QProgressEvent*>(event);
        switch(aProgressEvent->action())
        {
        case QProgressEvent_Init:
        {
            //cout<<"MainWindow::event()->____Initialize the progress bar____"<<endl;
            int min = aProgressEvent->min();
            int max = aProgressEvent->max();
            myProgressBar->setFormat(aProgressEvent->getMessage()+" %p");
            myProgressBar->setRange(min, max);
            myProgressBar->setValue(min);
            myProgressBar->show();
        }
            break;

        case QProgressEvent_Reset:
        {
            //cout<<"MainWindow::event()->____Reset the progress bar____"<<endl;
            myProgressBar->setFormat("Meshing finished");
            myProgressBar->setRange(0,100);
            myProgressBar->setValue(0);
            myProgressBar->hide();
        }
            break;

        case QProgressEvent_Update:
        {
            //cout<<"MainWindow::event()->____Update the progress bar____"<<endl;
            int val = aProgressEvent->val();
            myProgressBar->setFormat(aProgressEvent->getMessage()+" %p");
            myProgressBar->setValue(val);
        }
            break;
        }
        return true;
    }
    else
    {
        //! Make sure the rest of events are handled
        return QMainWindow::event(event);
    }
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------

MainWindow::MainWindow(QWidget *parent):QMainWindow(parent),ui(new Ui::MainWindow)
{
    cout<<"MainWindow::MainWindow()->____constructor called____"<<endl;

    //! -----------------------------------------------------
    //! set the ownership of the top and bottom left corners
    //! -----------------------------------------------------
    setCorner(Qt::TopLeftCorner, Qt::LeftDockWidgetArea);
    setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);

    //! --------------------------------------------------------
    //! user messages
    //! once created the pointer is passed to the Global object
    //! --------------------------------------------------------
    myUserMessages = new userMessagesModel(this);
    Global::status().myMessages = myUserMessages;

    //! ----------------------------
    //! create a progress indicator
    //! ----------------------------
    myProgressBar = new QProgressBar(this);
    myProgressBar->setFixedWidth(300);
    myProgressBar->hide();
    myProgressBar->setTextVisible(true);
    myProgressBar->setStyleSheet(QString::fromUtf8("text-align: center;"));

    //! -------------
    //! experimental
    //! -------------
    //optionsWidget *options = new optionsWidget(this);

    //! -----------------------
    //! create the working dir
    //! -----------------------
    this->createWorkingDir();

    //! -------------------------------
    //! activation of the OCC products
    //! -------------------------------
    thirdPartyActivate();

    //! -------------
    //! setup the Ui
    //! -------------
    ui->setupUi(this);

    //! -----------------
    //! set window title
    //! -----------------
    this->setWindowTitle(APPNAME);

    //! ----------------------------------
    //! geometry import - status variable
    //! ----------------------------------
    myGeometryImportStatus = geometryImportStatus_NotLoaded;

    //! ----------------------------------------
    //! the current save file - status variable
    //! ----------------------------------------
    myCurrentProjectName ="";

    //! ---------------------------------------------------------------------------
    //! the progress indicator
    //! horizontal style, two progress bars
    //! the widget has no parent, so it must be deleted manually in the destructor
    //! ---------------------------------------------------------------------------
    bool hasAdditionalBar = true;
    myProgressIndicator = new QProgressIndicator(hasAdditionalBar);

    //! ----------------------------------
    //! tab widget as central widget: add
    //! 1st) the main OCC viewer
    //! 2nd) the worksheet viewer
    //! 3rd) the convergence data viewer
    //! ----------------------------------
    myCentralTabWidget = new QTabWidgetExtended(this);
    myCentralTabWidget->setTabShape(QTabWidget::Triangular);
    myCentralTabWidget->setTabPosition(QTabWidget::South);
    this->setCentralWidget(myCentralTabWidget);

    //! -----------------------------------------------------
    //! 1st) tab of "theCentralTagWidget": "myMainOCCViewer"
    //! -----------------------------------------------------
    myMainOCCViewer = new occPostWidget(this);
    myMainOCCViewer->setObjectName("maingwindow");
    //myMainOCCViewer->setFocusPolicy(Qt::ClickFocus);
    myCentralTabWidget->addTab(myMainOCCViewer,"3D viewer");

    //! ---------------------------------------------------------
    //! 2nd) tab of "theCentralTagWidget": "mySimulationMonitor"
    //! ---------------------------------------------------------
    mySimulationMonitor= new SimulationMonitor(true,false,this);
    mySimulationMonitor->setObjectName("worksheetViewer");
    //mySimulationMonitor->setContinuoslyLogging(true,tools::getWorkingDir().append("/SolutionData/SolverOutput.txt"));
    myCentralTabWidget->addTab(mySimulationMonitor,"Worksheet");

    //! ------------------------------------------------------------
    //! 3rd) tab of "theCentralTagWidget": "myConvergenceDataChart"
    //! ------------------------------------------------------------
    //myConvergenceDataChart = new ConvergenceDataChart(this);
    //myConvergenceDataChart->setObjectName("ConvergenceDataChart");

    myConvergenceDataChart = new ConvergenceDataChart(this);
    myConvergenceDataChart->setObjectName("ConvergenceDataChart");

    myCentralTabWidget->addTab(myConvergenceDataChart,"Convergence");

    //! -----------------------------------------------------------
    //! the tabular data graph view - loads vs time through charts
    //! ------------------------------------------------------------
    myTabularDataGraphViewer1 = new TabularDataViewerClass1(this);

    //! -----------------------------
    //! graph viewer and dock window
    //! -----------------------------
    myTabularDataViewerDock = new QDockWidget("Graph viewer", this);
    //myTabularDataViewerDock->setWidget(myTabularDataGraphViewer);
    myTabularDataViewerDock->setWidget(myTabularDataGraphViewer1);
    this->addDockWidget(Qt::BottomDockWidgetArea, myTabularDataViewerDock);
    myTabularDataViewerDock->setVisible(true);

    //! --------------------------------------------------------
    //! the custom table model viewer - loads in a table view
    //! the view model "myUserMessages" has been created before
    //! --------------------------------------------------------
    myTabularDataView = new TableWidget(this);
    myTabularDataView->setMessagesModel(myUserMessages);

    //! ---------------------------------------
    //! the hystogram for the mesh metric view
    //! ---------------------------------------
    //myHistogram = new QHistogram(this);

    /*
    //! ---------------------
    //! test system messages
    //! ---------------------
    for(int n=0; n<10; n++)
    {
        userMessage testMessage(false,QString("test message %1").arg(n+1));
        Global::status().myMessages->appendMessage(testMessage);
    }
    //! ----------------
    //! end of the test
    //! ----------------
    */

    //! -------------------------------------------------------------------
    //! table widget and dock widget
    //! "TableWidget" is the "View" of the tabular data "customTableModel"
    //! -------------------------------------------------------------------
    myTabularDataDock = new QDockWidget("Tabular data viewer", this);
    myTabularDataDock->setWidget(myTabularDataView);
    this->addDockWidget(Qt::BottomDockWidgetArea,myTabularDataDock);
    myTabularDataDock->setVisible(true);

    //! -----------------------
    //! the simulation manager
    //! -----------------------
    mySimulationManager = new SimulationManager(this);
    mySimulationManager->setObjectName("simmanager");

    //! ------------------------------------
    //!  ... and its non-dockable container
    //! ------------------------------------
    mySimulationManagerContainer = new QTabWidget(this);
    mySimulationManagerContainer->setTabPosition(QTabWidget::South);
    mySimulationManagerContainer->setTabShape(QTabWidget::Triangular);
    mySimulationManagerContainer->addTab(mySimulationManager,"Tree view");

    //! -------------------------------
    //! simulation manager dock widget
    //! -------------------------------
    mySimulationManagerDock = new QDockWidget("Simulation manager", this);
    mySimulationManagerDock->setWidget(mySimulationManagerContainer);
    this->addDockWidget(Qt::LeftDockWidgetArea,mySimulationManagerDock);
    mySimulationManagerDock->setVisible(true);

    //! ------------------
    //! the detail viewer
    //! ------------------
    myDetailViewer = new DetailViewer(this);
    myDetailViewer->setObjectName("detailViewer");
    mySimulationManager->setDetailViewer(myDetailViewer);

    //! -------------------------------
    //! ... its non-dockable container
    //! -------------------------------
    myDetailViewerContainer = new QTabWidget(this);
    myDetailViewerContainer->setTabShape(QTabWidget::Triangular);
    myDetailViewerContainer->setTabPosition(QTabWidget::South);
    myDetailViewerContainer->setTabShape(QTabWidget::Triangular);
    myDetailViewerContainer->addTab(myDetailViewer,"Details");

    //! -------------------
    //! detail viewer dock
    //! -------------------
    myDetailViewerDock = new QDockWidget("Detail viewer", this);
    myDetailViewerDock->setWidget(myDetailViewerContainer);
    this->addDockWidget(Qt::LeftDockWidgetArea, myDetailViewerDock);
    myDetailViewerDock->setVisible(true);

    //! ----------------------------------
    //! clip tool and its dockable window
    //! ----------------------------------
    myClipTool = new clipTool(this);

    QWidget *container = new QWidget(this);
    QHBoxLayout *h = new QHBoxLayout;
    container->setLayout(h);
    h->addWidget(myClipTool);

    myClipToolDock = new QDockWidget("Clipping planes");
    myClipToolDock->setWidget(container);
    this->addDockWidget(Qt::LeftDockWidgetArea, myClipToolDock);
    myClipToolDock->setVisible(false);

    //! -----------
    //! status bar
    //! -----------
    statusLabel = new QLabel(this);
    statusLabel->setText(APPNAME);
    statusBar()->addPermanentWidget(statusLabel);

    //! ---------------------------------------------------------------------------
    //! uncomment if you want to put the progress bar into the bottom-right corner
    //! ---------------------------------------------------------------------------
    //statusBar()->addPermanentWidget(myProgressIndicator);
    //statusBar()->addPermanentWidget(myProgressBar);

    //! -------------------------------------------------
    //! creates the actions, menu, toolbar, dock windows
    //! -------------------------------------------------
    this->createActions();
    this->createMenu();
    this->createToolBars();
    this->createDockWidgets();
    this->setUpConnections();

    //! ----------------------
    //! set size and position
    //! ----------------------
    this->resize(1600,900);

    //! -------------------------------------
    //! create an empty simulation data base
    //! -------------------------------------
    mySimulationManager->createSimulationDataBaseEmpty();
}

//! ---------------------
//! function: destructor
//! function: details
//! ---------------------
MainWindow::~MainWindow()
{
    cout<<"MainWindow::~MainWindow()->____DESTRUCTOR CALLED____"<<endl;

    //! ------------------------------
    //! delete the progress indicator
    //! ------------------------------
    if(myProgressIndicator!=Q_NULLPTR) delete myProgressIndicator;

    //! -------------------------------
    //! update the "settings.txt" file
    //! -------------------------------
    this->updateSettings();
    delete ui;
}

//! ----------------------------
//! function: createDockWidgets
//! details:
//! ----------------------------
void MainWindow::createDockWidgets()
{
    cout<<"MainWindow::createDockWidgets()->____function called____"<<endl;

    //! --------------------------------------------------------------
    //! debug console dock: the debug console is created as "running"
    //!
    //! constructor of systemConsole:
    //! systemConsole {
    //! bool isCreatedAsRunning = true,
    //! bool isMenuBarVisible = true,
    //! bool isContinuosLog = false,
    //! const QString &myLogFile = QString(),
    //! QWidget *parent = 0};
    //! --------------------------------------------------------------
    myDebugConsoleDock = new QDockWidget("Debug console",this);

#ifndef DEBUG_VERSION
    myDebugConsole = new systemConsole(true,false,false,"",this);
#else
    QString logFileName = tools::getWorkingDir().append("/log.txt");
    myDebugConsole = new systemConsole(true,false,true,logFileName,this);
#endif

    myDebugConsole->setObjectName("debugConsole");
    myDebugConsoleDock->setWidget(myDebugConsole);
    this->addDockWidget(Qt::RightDockWidgetArea,myDebugConsoleDock);
    myDebugConsoleDock->setVisible(false);

    //! --------------------------------
    //! memory profiler dockable widget
    //! --------------------------------
    myMemoryProfiler = new memoryProfiler(2500,this);
    myMemoryProfiler->setWindowTitle("Memory profiler");
    this->addDockWidget(Qt::RightDockWidgetArea,myMemoryProfiler);
    myMemoryProfiler->setVisible(false);

    //! ----------------------------------------------------
    //! viewports for handling contacts
    //! note: "doubleViewPort" class inherits "QDockWidget"
    //! ----------------------------------------------------
    myDockableMasterViewPort= new dockableViewPort(this);
    myDockableMasterViewPort->setObjectName("topViewPort");
    myDockableMasterViewPort->setWindowTitle("Master viewer");
    this->addDockWidget(Qt::RightDockWidgetArea,myDockableMasterViewPort);
    myDockableMasterViewPort->setVisible(true);

    myDockableSlaveViewPort= new dockableViewPort(this);
    myDockableSlaveViewPort->setObjectName("bottomViewPort");
    myDockableSlaveViewPort->setWindowTitle("Slave viewer");
    this->addDockWidget(Qt::RightDockWidgetArea,myDockableSlaveViewPort);
    myDockableSlaveViewPort->setVisible(true);

    //! ------------------------------------
    //! Hystogram plot for mesh metric plot
    //! ------------------------------------
    //myHistogramDockContainer = new QDockWidget(this);
    //myHistogramDockContainer->setWindowTitle("Mesh metric");
    //myHistogramDockContainer->setWidget(myHistogram);
    //this->addDockWidget(Qt::RightDockWidgetArea,myHistogramDockContainer);
    //myHistogramDockContainer->setFloating(true);
    //myHistogramDockContainer->setVisible(false);
}

//! ------------------------------------
//! function: clearMasterSlaveViewPorts
//! details:
//! ------------------------------------
void MainWindow::clearMasterSlaveViewPorts(const TColStd_ListOfInteger &bodyListNbs)
{
    myDockableMasterViewPort->hideBody(bodyListNbs);
    myDockableSlaveViewPort->hideBody(bodyListNbs);
}

//! -------------------------------------
//! function: showBodiesOnMasterViewPort
//! details:
//! -------------------------------------
void MainWindow::showBodiesOnMasterViewPort(const TColStd_ListOfInteger &listOfActiveBodies, const TColStd_ListOfInteger &listOfNotActiveBodies)
{
    if(myMainOCCViewer->getWorkingMode()!=curWorkingMode_onContact) return;

    //! ----------------------
    //! selectable AIS_Shapes
    //! ----------------------
    TColStd_ListIteratorOfListOfInteger it;
    for(it.Initialize(listOfActiveBodies); it.More(); it.Next())
    {
        int bodyIndex = it.Value();
        const occHandle(AIS_ExtendedShape)& curAIS = occHandle(AIS_ExtendedShape)::DownCast(myDockableMasterViewPort->getViewPort()->getInteractiveObjects().value(bodyIndex));
        myDockableMasterViewPort->getContext()->SetDisplayMode(curAIS,AIS_Shaded,false);
        myDockableMasterViewPort->getContext()->Display(curAIS,false);
    }
    //! ------------------------
    //! unselectable AIS_Shapes
    //! ------------------------
    TColStd_ListIteratorOfListOfInteger it1;
    for(it1.Initialize(listOfNotActiveBodies); it1.More(); it1.Next())
    {
        int bodyIndex = it1.Value();
        const occHandle(AIS_ExtendedShape) &curAIS = occHandle(AIS_ExtendedShape)::DownCast(myDockableMasterViewPort->getViewPort()->getInteractiveObjects().value(bodyIndex));
        //const occHandle(AIS_ExtendedShape) &curAIS = occHandle(AIS_ExtendedShape)::DownCast(myDockableMasterViewPort->getViewPort()->getInteractiveObjects().at(bodyIndex));
        myDockableMasterViewPort->getContext()->SetDisplayMode(curAIS,AIS_WireFrame,false);
        myDockableMasterViewPort->getContext()->Display(curAIS,false);
        myDockableMasterViewPort->getContext()->Deactivate(curAIS);
    }
    myDockableMasterViewPort->FitAll();
    myDockableMasterViewPort->getContext()->UpdateCurrentViewer();
}

//! ------------------------------------
//! function: showBodiesOnSlaveViewPort
//! details:
//! ------------------------------------
void MainWindow::showBodiesOnSlaveViewPort(const TColStd_ListOfInteger &listOfActiveBodies, const TColStd_ListOfInteger &listOfNotActiveBodies)
{
    if(myMainOCCViewer->getWorkingMode()!=curWorkingMode_onContact) return;

    //! ----------------------
    //! selectable AIS_Shapes
    //! ----------------------
    TColStd_ListIteratorOfListOfInteger it;
    for(it.Initialize(listOfActiveBodies); it.More(); it.Next())
    {
        int bodyIndex = it.Value();
        const occHandle(AIS_ExtendedShape) &curAIS = occHandle(AIS_ExtendedShape)::DownCast(myDockableSlaveViewPort->getViewPort()->getInteractiveObjects().value(bodyIndex));
        //const occHandle(AIS_ExtendedShape) &curAIS = occHandle(AIS_ExtendedShape)::DownCast(myDockableSlaveViewPort->getViewPort()->getInteractiveObjects().at(bodyIndex));
        myDockableSlaveViewPort->getContext()->SetDisplayMode(curAIS,AIS_Shaded,false);
        myDockableSlaveViewPort->getContext()->Display(curAIS,false);
    }

    //! ------------------------
    //! unselectable AIS_Shapes
    //! ------------------------
    TColStd_ListIteratorOfListOfInteger it1;
    for(it1.Initialize(listOfNotActiveBodies); it1.More(); it1.Next())
    {
        int bodyIndex = it1.Value();
        const occHandle(AIS_ExtendedShape) &curAIS = occHandle(AIS_ExtendedShape)::DownCast(myDockableSlaveViewPort->getViewPort()->getInteractiveObjects().value(bodyIndex));
        //const occHandle(AIS_ExtendedShape) &curAIS = occHandle(AIS_ExtendedShape)::DownCast(myDockableSlaveViewPort->getViewPort()->getInteractiveObjects().at(bodyIndex));
        myDockableSlaveViewPort->getContext()->SetDisplayMode(curAIS,AIS_WireFrame,false);
        myDockableSlaveViewPort->getContext()->Display(curAIS,false);
        myDockableSlaveViewPort->getContext()->Deactivate(curAIS);
    }
    myDockableSlaveViewPort->FitAll();
    myDockableSlaveViewPort->getContext()->UpdateCurrentViewer();
}

//! ---------------------
//! function: createMenu
//! details:
//! ---------------------
void MainWindow::createMenu()
{
    //! --------------------
    //! "File" menu - setup
    //! --------------------
    FileMenu = menuBar()->addMenu("File");
    FileMenu->addAction(actionOpenProject);
    FileMenu->addAction(actionSaveProject);
    FileMenu->addAction(actionSaveProjectAs);
    FileMenu->addAction(actionImport);
    FileMenu->addAction(actionClose);
    FileMenu->addAction(actionExit);

    //! ------------------------
    //! "Geometry" menu - setup
    //! ------------------------
    GeometryMenu = menuBar()->addMenu("Geometry");

    //! --------------------
    //! "Mesh" menu - setup
    //! --------------------
    MeshMenu = menuBar()->addMenu("Mesh");

    //! ----------------------
    //! Mesh submenu "Insert"
    //! ----------------------
    MeshMenuInsert = MeshMenu->addMenu("Insert");
    MeshMenuInsert->setIcon(QIcon(":/icons/icon_insert.png"));
    MeshMenu->addSeparator();

    MeshMenuInsert->addAction(actionInsertMethod);
    MeshMenuInsert->addSeparator();
    MeshMenuInsert->addAction(actionInsertBodySizing);
    MeshMenuInsert->addAction(actionInsertFaceSizing);
    MeshMenuInsert->addAction(actionInsertEdgeSizing);
    MeshMenuInsert->addAction(actionInsertVertexSizing);
    MeshMenuInsert->addSeparator();
    MeshMenuInsert->addAction(actionInsertPrismaticLayer);

    MeshMenu->addAction(actionPreviewMesh);
    MeshMenu->addAction(actionGenerateMesh);
    MeshMenu->addSeparator();
    MeshMenu->addAction(actionLoadMesh);
    MeshMenu->addSeparator();
    MeshMenu->addAction(actionClearMesh);

    //! -------------
    //! "Tools" menu
    //! -------------
    ToolsMenu = menuBar()->addMenu("Tools");
    ToolsMenu->addAction(actionWriteInputFile);
    ToolsMenu->addAction(actionReadResults);
    ToolsMenu->addSeparator();
    ToolsMenu->addAction(actionOptions);

    //! --------------------
    //! "View" menu - setup
    //! --------------------
    ViewMenu = menuBar()->addMenu("View");
    ViewMenu->addAction(actionShadedExteriorAndEdges);
    ViewMenu->addAction(actionShadedExterior);
    ViewMenu->addAction(actionWireframe);

    //! ---------------------------------
    //! "View" subMenu "Display quality"
    //! ---------------------------------
    viewSubMenuDisplayQuality = new QMenu("Display quality",ViewMenu);
    ViewMenu->addMenu(viewSubMenuDisplayQuality);
    viewSubMenuDisplayQuality->addAction(actionHigh);
    viewSubMenuDisplayQuality->addAction(actionMedium);
    viewSubMenuDisplayQuality->addAction(actionCoarse);

    //! -------------------------
    //! "View" subMenu "Windows"
    //! -------------------------
    viewSubMenuWindows = new QMenu("Windows",ViewMenu);
    ViewMenu->addMenu(viewSubMenuWindows);
    viewSubMenuWindows->addAction(actionShowSimulationManager);
    viewSubMenuWindows->addAction(actionShowDetailViewer);
    viewSubMenuWindows->addAction(actionShowTabularData);
    viewSubMenuWindows->addAction(actionShowGraphViewer);    
    viewSubMenuWindows->addAction(actionShowDebugWindow);
    viewSubMenuWindows->addAction(actionShowSectionPlanes);
    viewSubMenuWindows->addAction(actionShowMemoryProfiler);
    viewSubMenuWindows->addAction(actionShowPhytonConsole);
}

//! -----------------------------
//! function: create the actions
//! details:
//! -----------------------------
void MainWindow::createActions()
{
    //!---------------------
    //! "File" menu actions
    //!---------------------
    actionOpenProject = new QAction("Open project",this);
    actionOpenProject->setIconVisibleInMenu(true);
    actionOpenProject->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_O));
    actionOpenProject->setIcon(QIcon(":/icons/icon_open project.png"));
    connect(actionOpenProject,SIGNAL(triggered(bool)),this,SLOT(openProject()));

    actionSaveProject = new QAction("Save project",this);
    actionSaveProject->setIconVisibleInMenu(true);
    actionSaveProject->setIcon(QIcon(":/icons/icon_save project.png"));
    actionSaveProject->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_S));
    connect(actionSaveProject,SIGNAL(triggered()),this,SLOT(saveProject()));

    actionSaveProjectAs = new QAction("Save project as",this);
    actionSaveProjectAs->setIconVisibleInMenu(true);
    actionSaveProjectAs->setIcon(QIcon(":/icons/icon_save project as.png"));
    actionSaveProjectAs->setShortcut(QKeySequence(Qt::SHIFT + Qt::Key_S));
    connect(actionSaveProjectAs,SIGNAL(triggered()),this,SLOT(saveProjectAs()));

    actionImport = new QAction("Import external geometry",this);
    actionImport->setIconVisibleInMenu(true);
    actionImport->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_I));
    actionImport->setIcon(QIcon(":/icons/icon_import geometry file.png"));
    connect(actionImport,SIGNAL(triggered()),this,SLOT(importFile()));

    actionExit = new QAction("Exit",this);
    actionExit->setIconVisibleInMenu(true);
    actionExit->setIcon(QIcon(":/icons/icon_exit.png"));
    actionExit->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_Q));
    connect(actionExit,SIGNAL(triggered()),this,SLOT(close()));

    actionClose = new QAction("Close the project",this);
    actionClose->setIconVisibleInMenu(true);
    actionClose->setIcon(QIcon(":/icons/icon_discard.png"));
    actionClose->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_X));
    connect(actionClose,SIGNAL(triggered(bool)),this,SLOT(closeTheModel()));

    //!-------------------------
    //! "Geometry" menu actions
    //!-------------------------

    //!---------------------
    //! "Mesh" menu actions
    //! --------------------
    actionPreviewMesh = new QAction("Preview the surface mesh", this);
    actionPreviewMesh->setIcon(QIcon(":/icons/icon_surface mesh.png"));
    connect(actionPreviewMesh,SIGNAL(triggered(bool)),mySimulationManager,SLOT(buildSurfaceMesh()));

    //! Generate volume mesh
    actionGenerateMesh = new QAction("Generate the mesh",this);
    actionGenerateMesh->setIcon(QIcon(":/icons/icon_volume mesh.png"));
    connect(actionGenerateMesh,SIGNAL(triggered(bool)),mySimulationManager,SLOT(buildVolumeMesh()));

    //! Clear mesh
    actionClearMesh = new QAction("Clear all the meshes",this);
    actionClearMesh->setIcon(QIcon(":/icons/icon_clear data.png"));
    connect(actionClearMesh,SIGNAL(triggered(bool)),myMainOCCViewer,SLOT(clearMeshFromViewer()));

    //! Meshing method
    actionInsertMethod = new QAction("Insert method",this);
    actionInsertMethod->setIcon(QIcon(":/icons/icon_mesh method.png"));
    connect(actionInsertMethod,SIGNAL(triggered(bool)),this,SLOT(createItemMeshMethod()));

    //! Body sizing
    actionInsertBodySizing  = new QAction("Body sizing",this);
    actionInsertBodySizing->setIcon(QIcon(":/icons/icon_volume mesh.png"));
    connect(actionInsertBodySizing,SIGNAL(triggered(bool)),this,SLOT(createItemMeshBodySizing()));

    //! Face sizing
    actionInsertFaceSizing = new QAction("Face sizing",this);
    actionInsertFaceSizing->setIcon(QIcon(":/icons/icon_mesh face sizing.png"));
    connect(actionInsertFaceSizing,SIGNAL(triggered(bool)),this,SLOT(createItemMeshFaceSizing()));

    //! Edge sizing
    actionInsertEdgeSizing = new QAction("Edge sizing",this);
    actionInsertEdgeSizing->setIcon(QIcon(":/icons/icon_mesh edge sizing.png"));
    connect(actionInsertEdgeSizing,SIGNAL(triggered(bool)),this,SLOT(createItemMeshEdgeSizing()));

    //! Vertex sizing
    actionInsertVertexSizing = new QAction("Vertex sizing",this);
    actionInsertVertexSizing->setIcon(QIcon(":/icons/icon_point.png"));
    connect(actionInsertVertexSizing,SIGNAL(triggered(bool)),this,SLOT(createItemMeshVertexSizing()));

    //! Prismatic layer
    actionInsertPrismaticLayer = new QAction("Insert prismatic layer",this);
    actionInsertPrismaticLayer->setIcon(QIcon(":/icons/icon_prismatic layer.png"));
    connect(actionInsertPrismaticLayer,SIGNAL(triggered(bool)),this,SLOT(createItemMeshPrismaticLayer()));

    //! load mesh from file
    actionLoadMesh = new QAction("Load mesh from file",this);
    actionLoadMesh->setIcon(QIcon(":/icons/icon_open from file.png"));

    //!-------------------------
    //! "Solution" menu actions
    //!-------------------------
    actionSolve = new QAction("Solve",this);
    actionSolve->setIcon(QIcon(":/icons/icon_solve.png"));
    connect(actionSolve,SIGNAL(triggered(bool)),this,SLOT(startAnalysis()));

    actionWriteInputFile = new QAction("Write solver input file",this);
    actionWriteInputFile->setIcon(QIcon(":/icons/icon_write solver input file.png"));
    connect(actionWriteInputFile,SIGNAL(triggered(bool)),this,SLOT(writeSolverInputFile()));

    actionReadResults = new QAction("Read results file");
    actionReadResults->setIcon(QIcon(":/icons/icon_reading.png"));
    connect(actionReadResults,SIGNAL(triggered(bool)),this,SLOT(readResultsFile()));

    //! image button
    actionImageToFile = new QAction("Save image",this);
    actionImageToFile->setIcon(QIcon(":/icons/icon_image to file.png"));
    connect(actionImageToFile,SIGNAL(triggered(bool)),myMainOCCViewer,SLOT(saveImageToDisk()));

    //! color box
    actionColorBox = new QAction("Show/hide color box",this);
    actionColorBox->setIcon(QIcon(":/icons/icon_colorbox.png"));
    actionColorBox->setCheckable(true);
    connect(actionColorBox,SIGNAL(toggled(bool)),this,SLOT(handleColorBoxVisibility(bool)));

    actionOptions = new QAction("Edit options", this);
    actionOptions->setIcon(QIcon(":/icons/icon_options.png"));
    connect(actionOptions,SIGNAL(triggered(bool)),this,SLOT(showOptionsWidget()));

    //!---------------------
    //! "View" menu actions
    //!---------------------
    groupOfActions_DisplayModes= new QActionGroup(this);

    //! first view mode - shaded exterior and wireframe
    actionShadedExteriorAndEdges = new QAction("Shaded Exterior and Edges", groupOfActions_DisplayModes);
    actionShadedExteriorAndEdges->setCheckable(true);

    //! set the shaded exterior and edge display mode as initially active
    actionShadedExteriorAndEdges->setChecked(true);
    connect(actionShadedExteriorAndEdges,SIGNAL(triggered()),myMainOCCViewer,SLOT(setShadedExteriorAndEdgesView()));

    //! second view mode - shaded exterior only
    actionShadedExterior = new QAction("Shaded Exterior", groupOfActions_DisplayModes);
    actionShadedExterior->setCheckable(true);
    connect(actionShadedExterior,SIGNAL(triggered()),myMainOCCViewer,SLOT(setShadedExteriorView()));

    //! third mode - wireframe view
    actionWireframe = new QAction("Wireframe", groupOfActions_DisplayModes);
    actionWireframe->setCheckable(true);
    connect(actionWireframe,SIGNAL(triggered()),myMainOCCViewer,SLOT(setWireframeView()));

    //! The following grouping sets the "exclusive mode"
    groupOfActions_DisplayModes->addAction(actionShadedExteriorAndEdges);
    groupOfActions_DisplayModes->addAction(actionShadedExterior);
    groupOfActions_DisplayModes->addAction(actionWireframe);
    groupOfActions_DisplayModes->setExclusive(true);

    //! "View" subMenu  - "Display quality" actions
    //! create a group of actions
    displayQualityActionGroup= new QActionGroup(this);
    actionHigh = new QAction("High", displayQualityActionGroup);
    actionHigh->setCheckable(true);
    connect(actionHigh,SIGNAL(triggered(bool)),this,SLOT(setHighQualityDisplay()));

    actionMedium = new QAction("Medium", displayQualityActionGroup);
    actionMedium->setCheckable(true);
    connect(actionMedium,SIGNAL(triggered()),this,SLOT(setMediumQualityDisplay()));

    actionCoarse = new QAction("Coarse", displayQualityActionGroup);
    actionCoarse->setCheckable(true);
    connect(actionCoarse,SIGNAL(triggered()),this,SLOT(setLowQualityDisplay()));

    //! fills the actions group
    displayQualityActionGroup->addAction(actionHigh);
    displayQualityActionGroup->addAction(actionMedium);
    displayQualityActionGroup->addAction(actionCoarse);
    displayQualityActionGroup->setExclusive(true);

    //! -----
    //! view
    //! -----
    actionShowGraphViewer = new QAction("Graph viewer",this);
    actionShowGraphViewer->setCheckable(true);
    actionShowGraphViewer->setChecked(true);

    actionShowTabularData = new QAction("Tabular data",this);
    actionShowTabularData->setCheckable(true);
    actionShowTabularData->setChecked(true);

    actionShowSimulationManager = new QAction("Simulation manager", this);
    actionShowSimulationManager->setCheckable(true);
    actionShowSimulationManager->setChecked(true);

    actionShowDetailViewer = new QAction("Detail viewer", this);
    actionShowDetailViewer->setCheckable(true);
    actionShowDetailViewer->setChecked(true);

    actionShowDebugWindow = new QAction("Debug window",this);
    actionShowDebugWindow->setCheckable(true);
    actionShowDebugWindow->setChecked(true);

    actionShowSectionPlanes = new QAction("Section planes",this);
    actionShowSectionPlanes->setCheckable(true);
    actionShowSectionPlanes->setChecked(false);

    actionShowMemoryProfiler = new QAction("Memory profiler",this);
    actionShowMemoryProfiler->setCheckable(true);
    actionShowMemoryProfiler->setChecked(false);

    actionShowPhytonConsole = new QAction("Phyton console",this);
    actionShowPhytonConsole->setCheckable(true);
    actionShowPhytonConsole->setChecked(false);

    //! "View" toolbar actions
    actionFitAll = new QAction("Toggle fit all",this);
    actionFitAll->setIcon(QIcon(":/icons/icon_fit all.png"));
    actionFitAll->setCheckable(true);
    connect(actionFitAll,SIGNAL(toggled(bool)),this,SLOT(fitAll(bool)));

    actionToggle3Drotation = new QAction("Toggle rotation",this);
    actionToggle3Drotation->setIcon(QIcon(":/icons/icon_rotate.png"));
    actionToggle3Drotation->setCheckable(true);
    connect(actionToggle3Drotation,SIGNAL(toggled(bool)),this,SLOT(toggle3Drotation(bool)));

    actionTogglePan = new QAction("Toggle pan",this);
    actionTogglePan->setIcon(QIcon(":/icons/icon_pan.png"));
    actionTogglePan->setCheckable(true);
    connect(actionTogglePan,SIGNAL(toggled(bool)),this,SLOT(togglePan(bool)));

    actionToggleWindowZoom = new QAction("Toggle zoom",this);
    actionToggleWindowZoom->setIcon(QIcon(":/icons/icon_zoom.png"));
    actionToggleWindowZoom->setCheckable(true);
    connect(actionToggleWindowZoom,SIGNAL(toggled(bool)),this,SLOT(toggleWindowZoom(bool)));

    //! "Selection mode" actions
    actionToggleSolidSelect = new QAction("Toggle solid selection",this);
    actionToggleSolidSelect->setIcon(QIcon(":/icons/icon_select solid.png"));
    actionToggleSolidSelect->setCheckable(true);
    connect(actionToggleSolidSelect,SIGNAL(toggled(bool)),this,SLOT(toggleSolidSelectionMode(bool)));

    actionToggleFaceSelect = new QAction("Toggle face selection",this);
    actionToggleFaceSelect->setIcon(QIcon(":/icons/icon_select face.png"));
    actionToggleFaceSelect->setCheckable(true);
    connect(actionToggleFaceSelect,SIGNAL(toggled(bool)),this,SLOT(toggleFaceSelectionMode(bool)));

    actionToggleEdgeSelect = new QAction("Toggle egde selection",this);
    actionToggleEdgeSelect->setIcon(QIcon(":/icons/icon_select edge.png"));
    actionToggleEdgeSelect->setCheckable(true);
    connect(actionToggleEdgeSelect,SIGNAL(toggled(bool)),this,SLOT(toggleEdgeSelectionMode(bool)));

    actionToggleVertexSelect = new QAction("Toggle vertex selection",this);
    actionToggleVertexSelect->setIcon(QIcon(":/icons/icon_select vertex.png"));
    actionToggleVertexSelect->setCheckable(true);
    connect(actionToggleVertexSelect,SIGNAL(toggled(bool)),this,SLOT(toggleVertexSelectionMode(bool)));

    actionTogglePickPointCoordinates = new QAction("Pick point coordinates",this);
    actionTogglePickPointCoordinates->setIcon(QIcon(":/icons/icon_pick point coordinates.png"));
    actionTogglePickPointCoordinates->setCheckable(true);
    connect(actionTogglePickPointCoordinates,SIGNAL(toggled(bool)),this,SLOT(togglePointCoordinatesPickingMode(bool)));

    //! action extend selection
    actionExtendToAdjacent = new QAction(QIcon(":/icons/icon_extend to adjacent.png"),"Extend to adjacent");
    connect(actionExtendToAdjacent,SIGNAL(triggered(bool)),this,SLOT(callExtendSelectionToAdjacent()));
    actionExtendToLimits = new QAction(QIcon(":/icons/icon_extend to limit.png"),"Extend to limits");
    connect(actionExtendToLimits,SIGNAL(triggered(bool)),this,SLOT(callExtendSelectionToLimits()));
    actionExtendToConnections = new QAction(QIcon(":/icons/icon_extend to limit.png"),"Extend to connections");
    connect(actionExtendToConnections,SIGNAL(triggered(bool)),this,SLOT(callExtendSelectionToConnections()));

    //! actions tranformations
    QBitmap mask(":/icons/icon_mask.bmp");
    actionTranslationX = new QAction("Translate X",this);
    actionTranslationY = new QAction("Translate Y",this);
    actionTranslationZ = new QAction("Translate Z",this);
    actionRotationX = new QAction("Rotation X",this);
    actionRotationY = new QAction("Rotation Y",this);
    actionRotationZ = new QAction("Rotation Z",this);

    actionDeleteTransformation = new QAction("Delete transformation",this);

    actionTranslationX->setData(1);
    actionTranslationY->setData(2);
    actionTranslationZ->setData(3);
    actionRotationX->setData(4);
    actionRotationY->setData(5);
    actionRotationZ->setData(6);
    actionDeleteTransformation->setData(7);

    QPixmap pixmapTranslationX(":/icons/icon_translation x.png",".png",Qt::AutoColor);
    QPixmap pixmapTranslationY(":/icons/icon_translation y.png",".png",Qt::AutoColor);
    QPixmap pixmapTranslationZ(":/icons/icon_translation z.png",".png",Qt::AutoColor);
    QPixmap pixmapRotationX(":/icons/icon_rotation x.png",".png",Qt::AutoColor);
    QPixmap pixmapRotationY(":/icons/icon_rotation y.png",".png",Qt::AutoColor);
    QPixmap pixmapRotationZ(":/icons/icon_rotation z.png",".png",Qt::AutoColor);

    pixmapTranslationX.setMask(mask);
    pixmapTranslationY.setMask(mask);
    pixmapTranslationZ.setMask(mask);
    pixmapRotationX.setMask(mask);
    pixmapRotationY.setMask(mask);
    pixmapRotationZ.setMask(mask);

    QIcon iconTranslationX(pixmapTranslationX);
    QIcon iconTranslationY(pixmapTranslationY);
    QIcon iconTranslationZ(pixmapTranslationZ);
    QIcon iconRotationX(pixmapRotationX);
    QIcon iconRotationY(pixmapRotationY);
    QIcon iconRotationZ(pixmapRotationZ);

    actionTranslationX->setIcon(iconTranslationX);
    actionTranslationY->setIcon(iconTranslationY);
    actionTranslationZ->setIcon(iconTranslationZ);
    actionRotationX->setIcon(iconRotationX);
    actionRotationY->setIcon(iconRotationY);
    actionRotationZ->setIcon(iconRotationZ);
    actionDeleteTransformation->setIcon(QIcon(":/icons/icon_delete.png"));

    connect(actionTranslationX,SIGNAL(triggered(bool)),this,SLOT(addCoordinateSystemTranslationX()));
    connect(actionTranslationY,SIGNAL(triggered(bool)),this,SLOT(addCoordinateSystemTranslationY()));
    connect(actionTranslationZ,SIGNAL(triggered(bool)),this,SLOT(addCoordinateSystemTranslationZ()));
    connect(actionRotationX,SIGNAL(triggered(bool)),this,SLOT(addCoordinateSystemRotationX()));
    connect(actionRotationY,SIGNAL(triggered(bool)),this,SLOT(addCoordinateSystemRotationY()));
    connect(actionRotationZ,SIGNAL(triggered(bool)),this,SLOT(addCoordinateSystemRotationZ()));
    connect(actionDeleteTransformation,SIGNAL(triggered(bool)),SLOT(removeTransformation()));
}

//! -------------------------
//! function: createToolbars
//! details:
//! -------------------------
void MainWindow::createToolBars()
{
    //! add the "project" toolbar
    projectToolbar = addToolBar("Project toolbar");
    projectToolbar->setAllowedAreas(Qt::TopToolBarArea);
    projectToolbar->setVisible(true);
    projectToolbar->setMovable(true);

    projectToolbar->addAction(actionImport);
    projectToolbar->addAction(actionOpenProject);
    projectToolbar->addAction(actionSaveProject);
    projectToolbar->addAction(actionSaveProjectAs);

    //! add the "view and selection" toolbar to the main window
    viewAndSelectionToolbar = addToolBar("View and selection toolbar");
    viewAndSelectionToolbar->setAllowedAreas(Qt::TopToolBarArea);
    viewAndSelectionToolbar->setVisible(true);
    viewAndSelectionToolbar->setMovable(true);

    //! add the select geometry/select mesh button
    QPushButtonExtended *buttonToggleSelectGeometryMesh = new QPushButtonExtended(this);
    buttonToggleSelectGeometryMesh->setIcon(QIcon(":/icons/icon_select geometry.png"));
    viewAndSelectionToolbar->addWidget(buttonToggleSelectGeometryMesh);

    QMenu *menuToggleSelectGeometryMesh = new QMenu(this);
    buttonToggleSelectGeometryMesh->setMenu(menuToggleSelectGeometryMesh);

    QAction *actionToggleSelectGeometry = new QAction(QIcon(":/icons/icon_select geometry.png"),"Select geometry",this);
    connect(actionToggleSelectGeometry,SIGNAL(triggered(bool)),this,SLOT(setSelectionModeGeometry()));

    QAction *actionToggleSelectMesh = new QAction(QIcon(":/icons/icon_select mesh.png"),"Select mesh",this);
    connect(actionToggleSelectMesh,SIGNAL(triggered(bool)),this,SLOT(setSelectionModeMesh()));

    menuToggleSelectGeometryMesh->addAction(actionToggleSelectGeometry);
    menuToggleSelectGeometryMesh->addAction(actionToggleSelectMesh);

    //! add the pick point coordinates
    viewAndSelectionToolbar->addAction(actionTogglePickPointCoordinates);

    //! add the "view" actions to the "view and selection" toolbar
    viewAndSelectionToolbar->addAction(actionToggle3Drotation);
    viewAndSelectionToolbar->addAction(actionTogglePan);
    viewAndSelectionToolbar->addAction(actionToggleWindowZoom);
    viewAndSelectionToolbar->addAction(actionFitAll);

    //! add a separator
    viewAndSelectionToolbar->addSeparator();

    //! add the "type of selection" menu button to the toolbar
    myPushButtonTypeOfSelection = new QPushButtonExtended(this);
    myPushButtonTypeOfSelection->setIcon(QIcon(":/icons/icon_single select.png"));
    viewAndSelectionToolbar->addWidget(myPushButtonTypeOfSelection);
    QMenu *menuButtonTypeOfSelection = new QMenu(this);
    myPushButtonTypeOfSelection->setMenu(menuButtonTypeOfSelection);

    actionSingleSelect = new QAction(QIcon(":/icons/icon_single select.png"),"Single select");
    connect(actionSingleSelect,SIGNAL(triggered(bool)),this,SLOT(setSelectionModeSingle()));

    actionMultipleSelect = new QAction(QIcon(":/icons/icon_multiple select.png"),"Box select");
    connect(actionMultipleSelect,SIGNAL(triggered(bool)),this,SLOT(setSelectionModeBox()));

    actionSelectAll = new QAction(QIcon(":/icons/icon_select all.png"),"Salect all");
    connect(actionSelectAll,SIGNAL(triggered(bool)),this,SLOT(HandleSelectAll()));

    menuButtonTypeOfSelection->addAction(actionSingleSelect);
    menuButtonTypeOfSelection->addAction(actionMultipleSelect);
    menuButtonTypeOfSelection->addAction(actionSelectAll);

    //! add the "selection mode" actions to the "view and selection" toolbar
    viewAndSelectionToolbar->addAction(actionToggleSolidSelect);
    viewAndSelectionToolbar->addAction(actionToggleFaceSelect);
    viewAndSelectionToolbar->addAction(actionToggleEdgeSelect);
    viewAndSelectionToolbar->addAction(actionToggleVertexSelect);

    //! add separator
    viewAndSelectionToolbar->addSeparator();

    //! extend selection
    myExtendSelectionButton = new QPushButton();
    myExtendSelectionButton->setIcon(QIcon(":/icons/icon_extend to adjacent.png"));
    QMenu *menuButton = new QMenu(this);
    menuButton->addAction(actionExtendToAdjacent);
    menuButton->addAction(actionExtendToLimits);
    menuButton->addAction(actionExtendToConnections);
    myExtendSelectionButton->setMenu(menuButton);
    viewAndSelectionToolbar->addWidget(myExtendSelectionButton);

    //! -------------------------
    //! add the solution toolbar
    //! -------------------------
    solutionToolBar = addToolBar("Solution toolbar");
    solutionToolBar->setAllowedAreas(Qt::TopToolBarArea);
    solutionToolBar->setVisible(true);
    solutionToolBar->setMovable(true);

    //! add the actions
    solutionToolBar->addAction(actionSolve);
    solutionToolBar->addAction(actionWriteInputFile);

    //! ------------------
    //! save image button
    //! ------------------
    myImageButton = new QPushButtonExtended();
    myImageButton->setIcon(QIcon(":/icons/icon_image to file.png"));
    QMenu *imageMenu = new QMenu(this);
    imageMenu->addAction(actionImageToFile);
    myImageButton->setMenu(imageMenu);
    solutionToolBar->addWidget(myImageButton);

    //! colorbox button
    solutionToolBar->addAction(actionColorBox);

    //! -------------------------------------------------
    //! transformation tool bar - for coordinate systems
    //! -------------------------------------------------
    transformationsToolBar = addToolBar("Transformation tool bar");
    transformationsToolBar->setAllowedAreas(Qt::TopToolBarArea);
    transformationsToolBar->setVisible(true);
    transformationsToolBar->setMovable(true);
    transformationsToolBar->addAction(actionTranslationX);
    transformationsToolBar->addAction(actionTranslationY);
    transformationsToolBar->addAction(actionTranslationZ);
    transformationsToolBar->addAction(actionRotationX);
    transformationsToolBar->addAction(actionRotationY);
    transformationsToolBar->addAction(actionRotationZ);
    transformationsToolBar->addAction(actionDeleteTransformation);

    //! -------------
    //! mesh toolbar
    //! -------------
    meshToolBar = new MeshToolBar("Mesh",this);
    this->addToolBar(meshToolBar);
    meshToolBar->setVisible(true);
    meshToolBar->enableMeshViewButton(false);
    meshToolBar->enableClearMesh(false);
    meshToolBar->enableSurfaceMeshButton(false);
    meshToolBar->enableVolumeMeshButton(false);

    //! ----------------
    //! results toolbar
    //! ----------------
    resultsToolBar = new ResultsToolBar("Results",this);
    this->addToolBar(resultsToolBar);
    resultsToolBar->setVisible(false);

    //! for changing the status variable of the main viewer
    connect(resultsToolBar,SIGNAL(requestUpdateViewerStatus()),myMainOCCViewer,SLOT(updateViewerStatus()));
    connect(meshToolBar,SIGNAL(requestUpdateViewerStatus()),myMainOCCViewer,SLOT(updateViewerStatus()));

    connect(myMainOCCViewer,SIGNAL(resultsPresentationChanged()),mySimulationManager,SLOT(updateResultsPresentation()));
}

//! -------------------------------
//! function: createItemMeshMethod
//! details:
//! -------------------------------
void MainWindow::createItemMeshMethod()
{
    mySimulationManager->createSimulationNode(SimulationNodeClass::nodeType_meshMethod);
}

//! -----------------------------------
//! function: createItemMeshBodySizing
//! function:
//! -----------------------------------
void MainWindow::createItemMeshBodySizing()
{
    mySimulationManager->createSimulationNode(SimulationNodeClass::nodeType_meshBodyMeshControl);
}

//! -----------------------------------
//! function: createItemMeshFaceSizing
//! details:
//! -----------------------------------
void MainWindow::createItemMeshFaceSizing()
{
    mySimulationManager->createSimulationNode(SimulationNodeClass::nodeType_meshFaceSize);
}

//! -----------------------------------
//! function: createItemMeshEdgeSizing
//! details:
//! -----------------------------------
void MainWindow::createItemMeshEdgeSizing()
{
    mySimulationManager->createSimulationNode(SimulationNodeClass::nodeType_meshEdgeSize);
}

//! -----------------------------------
//! function: createItemMeshVertexSizing
//! details:
//! -----------------------------------
void MainWindow::createItemMeshVertexSizing()
{
    mySimulationManager->createSimulationNode(SimulationNodeClass::nodeType_meshVertexSize);
}

//! ---------------------------------------
//! function: createItemMeshPrismaticLayer
//! details:
//! ---------------------------------------
void MainWindow::createItemMeshPrismaticLayer()
{
    mySimulationManager->createSimulationNode(SimulationNodeClass::nodeType_meshPrismaticLayer);
}

//! ---------------------------------
//! function: setSelectionModeSingle
//! details:
//! ---------------------------------
void MainWindow::setSelectionModeSingle()
{
    myMainOCCViewer->setGlobalCurSelectionMode(0);
    myPushButtonTypeOfSelection->setIcon(QIcon(":/icons/icon_single select.png"));
    statusBar()->showMessage("Single selection",TRANSIENT_MESSAGE_TIMEOUT);
}

//! -----------------------------------
//! function: setSelectionModeGeometry
//! details:
//! -----------------------------------
void MainWindow::setSelectionModeGeometry()
{
    cout<<"setSelectionModeGeometry()->____function called____"<<endl;
    myMainOCCViewer->setGeometrySelectionMode();
}

//! -------------------------------
//! function: setSelectionModeMesh
//! details:
//! -------------------------------
void MainWindow::setSelectionModeMesh()
{
    cout<<"setSelectionModeMesh()->____function called____"<<endl;
    myMainOCCViewer->setMeshSelectionMode();
}

//! ------------------------------
//! function: setSelectionModeBox
//! details:
//! ------------------------------
void MainWindow::setSelectionModeBox()
{
    myMainOCCViewer->setGlobalCurSelectionMode(1);
    myPushButtonTypeOfSelection->setIcon(QIcon(":/icons/icon_multiple select.png"));
    statusBar()->showMessage("Multiple selection",TRANSIENT_MESSAGE_TIMEOUT);
}

//! ---------------------------------------
//! function: close the model
//! details:  equivalent to a full restart
//! ---------------------------------------
void MainWindow::closeTheModel()
{
    //! message box
    int ret = QMessageBox::question(this, tr(APPNAME),tr("Close and restart?"),QMessageBox::Close | QMessageBox::Discard);

    switch (ret)
    {
    case (QMessageBox::Discard):
    {
        // do nothing
    }
        break;

    case (QMessageBox::Close):
    {
        //! ---------------
        //! write settings
        //! ---------------
        this->updateSettings();

        //! --------------------------------------------------
        //! status variable: unlock a future import operation
        //! --------------------------------------------------
        myGeometryImportStatus = geometryImportStatus_NotLoaded;

        //! ---------------------------------
        //! disable the mesh toolbar buttons
        //! ---------------------------------
        meshToolBar->enableMeshViewButton(false);
        meshToolBar->enableClearMesh(false);
        meshToolBar->enableSurfaceMeshButton(false);
        meshToolBar->enableVolumeMeshButton(false);

        //! Resets the status of the "view and selection toolbar"
        {
            QList<QAction *> listOfActions = viewAndSelectionToolbar->actions();
            for(int i=0;i<listOfActions.size();i++)
            {
                // uncheck and enable the buttons
                listOfActions.at(i)->setChecked(false);
                listOfActions.at(i)->setEnabled(true);
            }
        }

        //! Reset the main viewers and the contact viewports
        myMainOCCViewer->reset();
        myDockableMasterViewPort->getViewPort()->reset();
        myDockableSlaveViewPort->getViewPort()->reset();

        //! Clean the debug console
        myDebugConsole->clear();

        //! Clean the CCX messages console
        //! to do ...

        //! clear the simulation manager tree
        mySimulationManager->clearTree();

        //! clear the detail viewer tree
        myDetailViewer->clearTree();

        //! re-init the QStandardItemModel
        mySimulationManager->createSimulationDataBaseEmpty();

        //! update the window title
        this->setWindowTitle(APPNAME);
    }
        break;

    default:

          break;
    }
    //! delete <project name>_files
}

//! ---------------------
//! function: importFile
//! details:
//! ---------------------
void MainWindow::importFile(QString &fileName)
{
    //!cout<<"MainWindow::importFile()->____function called: "<<fileName.toStdString()<<"____"<<endl;
    if(myGeometryImportStatus == geometryImportStatus_NotLoaded)
    {
        //! -----------------------------------------------------
        //! Reset the status of the "view and selection toolbar"
        //! -----------------------------------------------------
        QList<QAction*> listOfActions = viewAndSelectionToolbar->actions();
        for(int i=0;i<listOfActions.size();i++)
        {
            //! uncheck and enable the buttons
            listOfActions.at(i)->setChecked(false);
            listOfActions.at(i)->setEnabled(true);
        }

        //! ----------------------------------------------------------------------
        //! Open a modal dialog and return the full name of the file
        //! The selection starts from the GEOMETRY_DIR
        //! The filters are CAD formats (STEP, IGES, ...) which OCC can translate
        //! ----------------------------------------------------------------------
        QString selectedFilter;
        if(fileName=="")
        {
            fileName = QFileDialog::getOpenFileName(this,"Import external geometry file",GEOMETRY_DIR,
                                                    "*.step *.stp *.STEP *STP;;"
                                                    "*.iges *.igs *.IGES *.IGS;;"
                                                    "*.brep *.BREP",&selectedFilter,0);
        }

        cout<<"MainWindow::importFile()->____importing file: "<<fileName.toStdString()<<"____"<<endl;

        if(fileName!="")  //! fileName is absolute
        {
            QList<QString> listOfNames;
            TopoDS_Compound shapeFromReader;

            //! ---------------------------------
            //! create an OCC progress indicator
            //! ---------------------------------
            occHandle(QOccProgressIndicator) aProgressIndicator = new QOccProgressIndicator(0,100,this);

            //! -------------------------------------
            //! load the CAD model
            //! input:   fileName
            //! outputs: shapeFromReader,listOfNames
            //! -------------------------------------
            bool isDone = mySimulationManager->loadCADModel(fileName, shapeFromReader, listOfNames, aProgressIndicator);
            if(isDone)
            {
                //! --------------
                //! mesh tool bar
                //! --------------
                meshToolBar->enableMeshViewButton(true);
                meshToolBar->enableClearMesh(true);
                meshToolBar->enableSurfaceMeshButton(true);
                meshToolBar->enableVolumeMeshButton(true);

                //! ----------------
                //! status variable
                //! ----------------
                myGeometryImportStatus = geometryImportStatus_Loaded;

                //! --------------------------------
                //! create the simulation data base
                //! --------------------------------
                aProgressIndicator->NewScope(10,"Creating simulation database");
                aProgressIndicator->Show();
                mySimulationManager->createSimulationDataBase(shapeFromReader,fileName,listOfNames);
                aProgressIndicator->EndScope();

                //! -------------------------------------------------
                //! feed the clip plane tool with coordinate systems
                //! -------------------------------------------------
                QStandardItem *itemCSRoot = mySimulationManager->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
                myClipTool->setCoordinateSystemRoot(itemCSRoot);

                //! -------------------------------------------------
                //! feed the clip plane tool with the viewer pointer
                //! -------------------------------------------------
                myClipTool->setViewer(myMainOCCViewer);
                myClipTool->setMeshDataBase(mySimulationManager->getDataBase());

                //! -------------------------------------------------------
                //! share the simulation data dase pointer with the viewer
                //! -------------------------------------------------------
                myMainOCCViewer->setSimulationdataBase(mySimulationManager->getDataBase());
                myDockableMasterViewPort->setSimulationdataBase(mySimulationManager->getDataBase());
                myDockableSlaveViewPort->setSimulationdataBase(mySimulationManager->getDataBase());

                //! ---------------------------
                //! the CAD model is displayed
                //! ---------------------------
                aProgressIndicator->NewScope(5,"Building objects");
                aProgressIndicator->Show();

                if(!shapeFromReader.IsNull())
                {
                    myMainOCCViewer->createInteractiveShapes();
                    myDockableMasterViewPort->createInteractiveShapes();
                    myDockableSlaveViewPort->createInteractiveShapes();
                }

                aProgressIndicator->EndScope();
                aProgressIndicator->NewScope(5,"Displaying ...");
                aProgressIndicator->Show();

                if(!shapeFromReader.IsNull())
                {
                    myMainOCCViewer->displayCAD();
                    myDockableMasterViewPort->getViewPort()->displayCAD();
                    myDockableSlaveViewPort->getViewPort()->displayCAD();
                }
                aProgressIndicator->EndScope();

                if(shapeFromReader.IsNull()==false)
                {
                    //! set the selection model - it can be done only when something has been loaded
                    mySimulationManager->setSelectionModel();

                    const occHandle(AIS_InteractiveContext) &aCTX = myMainOCCViewer->getContext();
                    if(aCTX.IsNull())
                    {
                        cout<<"MainWindow::importFile()->____error: context is NULL____"<<endl;
                        exit(1);
                    }

                    //! set the interactive context for simulation manager
                    cout<<"MainWindow::importFile()->____set the context for the simulation manager____"<<endl;
                    mySimulationManager->setContext(myMainOCCViewer->getContext());

                    //! set the context for the Detail viewer
                    cout<<"MainWindow::importFile()->____set the context for the detail viewer____"<<endl;
                    myDetailViewer->setContext(myMainOCCViewer->getContext());

                    cout<<"MainWindow::importFile()->____set the mesh context for the detail viewer____"<<endl;
                    myDetailViewer->setMeshContext(myMainOCCViewer->getMeshContext());
                }

                //! ------------------------------------
                //! change the title of the main window
                //! ------------------------------------
                QString curDirectoryPath = fileName;
                QString relativeName = fileName.split("/").last();
                curDirectoryPath.chop(relativeName.length()+1);
                this->setWindowTitle(curDirectoryPath+"/Untitled.gil");
            }
        }
        else
        {
            //! --------------
            //! mesh tool bar
            //! --------------
            meshToolBar->enableMeshViewButton(false);
            meshToolBar->enableClearMesh(false);
            meshToolBar->enableSurfaceMeshButton(false);
            meshToolBar->enableVolumeMeshButton(false);

            //! ----------------
            //! status variable
            //! ----------------
            myGeometryImportStatus = geometryImportStatus_NotLoaded;
            this->setWindowTitle(APPNAME);
        }
    }
}

//! -------------------------
//! function: updateViewport
//! details:
//! -------------------------
void MainWindow::updateViewport()
{
    cout<<"MainWindow::updateViewport()->____function called____"<<endl;

    //! --------------------------------------------------------------
    //! reset the maps of the interactive objects (shapes and meshes)
    //! contained in the viewport classes. Erase also the interactive
    //! object from the display
    //! --------------------------------------------------------------
    myMainOCCViewer->reset();
    myDockableSlaveViewPort->getViewPort()->reset();
    myDockableMasterViewPort->getViewPort()->reset();

    //! -------------------------------------------------------
    //! share the simulation data dase pointer with the viewer
    //! -------------------------------------------------------
    myMainOCCViewer->setSimulationdataBase(mySimulationManager->getDataBase());
    myDockableMasterViewPort->setSimulationdataBase(mySimulationManager->getDataBase());
    myDockableSlaveViewPort->setSimulationdataBase(mySimulationManager->getDataBase());

    //! -------------------------------
    //! create the interactive objects
    //! -------------------------------
    myMainOCCViewer->createInteractiveShapes();
    myDockableMasterViewPort->createInteractiveShapes();
    myDockableSlaveViewPort->createInteractiveShapes();

    //! ------------------
    //! display functions
    //! ------------------
    myMainOCCViewer->displayCAD();
    myDockableMasterViewPort->getViewPort()->displayCAD();
    myDockableSlaveViewPort->getViewPort()->displayCAD();

    //! -----------------------------------------------------------------------------
    //! set the selection model - it can be done only when something has been loaded
    //! -----------------------------------------------------------------------------
    mySimulationManager->setSelectionModel();

    //! ----------------------------------------------------------------------------
    //! set the interactive context (communication with the display)
    //! ----------------------------------------------------------------------------
    mySimulationManager->setContext(myMainOCCViewer->getContext());

    //! -----------
    //! diagnostic
    //! -----------
    if(myMainOCCViewer->getContext().IsNull()) cout<<"MainWindow::importFile()->____error: context is NULL____"<<endl;
    if(myMainOCCViewer->getMeshContext().IsNull()) cout<<"MainWindow::importFile()->____error: context is NULL____"<<endl;

    //! ---------------------------------------------------------------------------------
    //! set the context for the Detail viewer - could be removed (do not understand why)
    //! ---------------------------------------------------------------------------------
    myDetailViewer->setContext(myMainOCCViewer->getContext());
    myDetailViewer->setMeshContext(myMainOCCViewer->getMeshContext());

    //! ------------------------------------
    //! change the title of the main window
    //! ------------------------------------
    //QString curDirectoryPath = fileName;
    //QString relativeName = fileName.split("/").last();
    //curDirectoryPath.chop(relativeName.length()+1);
    //this->setWindowTitle(curDirectoryPath+"/Untitled.gil");
}

//! ----------------------
//! function: openProject
//! details:
//! ----------------------
void MainWindow::openProject()
{
    //! -----------------------------------------------------
    //! Reset the status of the "view and selection toolbar"
    //! uncheck and enable the buttons
    //! -----------------------------------------------------
    QList<QAction *> listOfActions = viewAndSelectionToolbar->actions();
    for(int i=0;i<listOfActions.size();i++)
    {
        listOfActions.at(i)->setChecked(false);
        listOfActions.at(i)->setEnabled(true);
    }

    //! --------------------------------------------
    //! open a modal dialog: "fileName" is absolute
    //! --------------------------------------------
    QString selectedFilter;
    QString fileName = QFileDialog::getOpenFileName(this,"Select project ", QDir::current().absolutePath(),"*.gil",&selectedFilter);

    if(!fileName.isEmpty())
    {
        //! -----------------------------------
        //! fileName is the full absolute path
        //! -----------------------------------
        myCurrentProjectName = fileName;

        //! ---------------------------------------------------------
        //! create the name of the directory for uncompressing files
        //! ---------------------------------------------------------
        QString deflatedDirName;
        deflatedDirName = fileName;
        deflatedDirName.chop(4) ;
        deflatedDirName.append("_files");

        //! -------------------------------------------
        //! set the current working directory - Simone
        //! -------------------------------------------
        QString relativeName = fileName.split("/").last();
        QString startingDirectory = fileName;
        startingDirectory.chop(relativeName.length()+1);

        myWorkingDir = startingDirectory;
        cout<<"MainWindow::openProject()->____working Directory____"<<startingDirectory.toStdString()<<endl;

        //! ----------------------------------
        //! set the current project directory
        //! ----------------------------------
        //mySimulationManager->setCurrentProjectDir(deflatedDirName);

        //! ------------------------------------------------------
        //! check if the directory <project>_files already exists
        //! ------------------------------------------------------
        if(QDir(deflatedDirName).exists())
        {
            cout<<"MainWindow::openProject()->____the <project>_files directory already exixts____"<<endl;
            QMessageBox::warning(this,APPNAME,"A simulation data base with the same name\n"
                                              "already exists, probably due to a previous crash.", QMessageBox::Ok);

            cout<<"MainWindow::openProject()->____deleting the old <project>_files____";
            tools::clearDir(deflatedDirName);
        }

        cout<<"MainWindow::openProject()->____Deflating file"<<fileName.toStdString()<<" into directory: "<<deflatedDirName.toStdString()<<"____"<<endl;
        cout<<"MainWindow::openProject()->____application current directory: "<<QDir::current().absolutePath().toStdString()<<"____"<<endl;
        cout<<"____"<<myWorkingDir.toStdString()<<"____"<<endl;

        int NbFiles = JlCompress::extractDir(fileName,deflatedDirName).size();
        if(NbFiles!=0)
        {
            //! ----------------
            //! close the model
            //! ----------------
            //this->closeTheModel();

            //! -----------------------------------------
            //! re-build the analysis database from disk
            //! -----------------------------------------
            mySimulationManager->buildDataBaseFromDisk(fileName);

            //! ----------------------------
            //! set the interactive context
            //! ----------------------------            
            const occHandle(AIS_InteractiveContext) &aCTX = myMainOCCViewer->getContext();
            if(aCTX.IsNull())
            {
                cerr<<"MainWindow::openProject()->____cannot get context from viewer____"<<endl;
                exit(9999);
            }

            mySimulationManager->setContext(myMainOCCViewer->getContext());

            myMainOCCViewer->setSimulationdataBase(mySimulationManager->getDataBase());
            myMainOCCViewer->createInteractiveShapes();
            myMainOCCViewer->displayCAD();
            myMainOCCViewer->buildMeshIOs();
            meshToolBar->enableMeshViewButton(true);
            meshToolBar->enableClearMesh(true);

            myDockableMasterViewPort->setSimulationdataBase(mySimulationManager->getDataBase());
            myDockableMasterViewPort->createInteractiveShapes();
            myDockableMasterViewPort->displayCAD();

            myDockableSlaveViewPort->setSimulationdataBase(mySimulationManager->getDataBase());
            myDockableSlaveViewPort->createInteractiveShapes();
            myDockableSlaveViewPort->displayCAD();

            //! -----------------------------------------------------------------
            //! set the interactive shape and mesh context for the Detail viewer
            //! -----------------------------------------------------------------
            myDetailViewer->setContext(myMainOCCViewer->getContext());
            myDetailViewer->setMeshContext(myMainOCCViewer->getMeshContext());

            //! -------------------------------------------------
            //! feed the clip plane tool with coordinate systems
            //! -------------------------------------------------
            QStandardItem *itemCSRoot = mySimulationManager->getTreeItem(SimulationNodeClass::nodeType_coordinateSystems);
            myClipTool->setCoordinateSystemRoot(itemCSRoot);

            //! -------------------------------------------------
            //! feed the clip plane tool with the viewer pointer
            //! and with the geometry/mesh database
            //! -------------------------------------------------
            myClipTool->setViewer(myMainOCCViewer);
            myClipTool->setMeshDataBase(mySimulationManager->getDataBase());

            //! -------------------------
            //! update main window title
            //! -------------------------
            this->setWindowTitle(fileName);

            //! -----------------------------------------------------------------------------------
            //! inform the system the model has been loaded (it avoids multiple import operations)
            //! -----------------------------------------------------------------------------------
            myGeometryImportStatus = geometryImportStatus_Loaded;
        }
        else
        {
            QMessageBox::critical(this, tr(APPNAME),tr("The selected file is empty"));
            myGeometryImportStatus = geometryImportStatus_NotLoaded;

            //! -------------------------
            //! update main window title
            //! -------------------------
            this->setWindowTitle(APPNAME);
        }
    }
}

//! ---------------------------
//! function: toggle3Drotation
//! details:
//! ---------------------------
void MainWindow::toggle3Drotation(bool isActivated)
{
    if(isActivated)
    {
        //!cout<<"MainWindow::toggle3Drotation()->____ROTATE ACTIVATED____"<<endl;
        handleViewAndSelectionButtons(actionToggle3Drotation);
        myMainOCCViewer->setAction3D_Rotation();
        myDockableMasterViewPort->setAction3D_Rotation();
        myDockableSlaveViewPort->setAction3D_Rotation();
        statusBar()->showMessage("Rotating",TRANSIENT_MESSAGE_TIMEOUT);
    }
}

//! ----------------------
//! function: toggle zoom
//! details:
//! ----------------------
void MainWindow::toggleWindowZoom(bool isActivated)
{
    if(isActivated)
    {
        handleViewAndSelectionButtons(actionToggleWindowZoom);
        myMainOCCViewer->setAction3D_WindowZooming();
        myDockableMasterViewPort->setAction3D_WindowZooming();
        myDockableSlaveViewPort->setAction3D_WindowZooming();
        statusBar()->showMessage("Window zooming",TRANSIENT_MESSAGE_TIMEOUT);
    }
}

//! --------------------
//! function: togglePan
//! details:
//! --------------------
void MainWindow::togglePan(bool isActivated)
{
    if(isActivated)
    {
        handleViewAndSelectionButtons(actionTogglePan);
        myMainOCCViewer->setAction3D_Pan();
        myDockableMasterViewPort->setAction3D_Pan();
        myDockableSlaveViewPort->setAction3D_Pan();
        statusBar()->showMessage("Panning",TRANSIENT_MESSAGE_TIMEOUT);
    }
}

//! -----------------------
//! function: toggleFitAll
//! details:
//! -----------------------
void MainWindow::fitAll(bool isActivated)
{
    if(isActivated)
    {
        handleViewAndSelectionButtons(actionFitAll);
        myMainOCCViewer->FitAll();
        myDockableMasterViewPort->FitAll();
        myDockableSlaveViewPort->FitAll();
    }
}

//! -----------------------------------
//! function: toggleSolidSelectionMode
//! details:
//! -----------------------------------
void MainWindow::toggleSolidSelectionMode(bool isActivated)
{
    if(isActivated)
    {
        handleViewAndSelectionButtons(actionToggleSolidSelect);
        myMainOCCViewer->setSelectionMode(CurSelection_Solid);

        //! same selection mode in additional viewports
        myDockableMasterViewPort->setSelectionMode(CurSelection_Solid);
        myDockableSlaveViewPort->setSelectionMode(CurSelection_Solid);
        statusBar()->showMessage("Select Solids",TRANSIENT_MESSAGE_TIMEOUT);
    }
}

//! ----------------------------------
//! function: toggleFaceSelectionMode
//! details:
//! ----------------------------------
void MainWindow::toggleFaceSelectionMode(bool isActivated)
{
    if(isActivated)
    {
        handleViewAndSelectionButtons(actionToggleFaceSelect);
        myMainOCCViewer->setSelectionMode(CurSelection_Face);

        //! same selection mode in additional viewports
        myDockableMasterViewPort->setSelectionMode(CurSelection_Face);
        myDockableSlaveViewPort->setSelectionMode(CurSelection_Face);
        statusBar()->showMessage("Select Faces",TRANSIENT_MESSAGE_TIMEOUT);
    }
}

//! ----------------------------------
//! function: toggleEdgeSelectionMode
//! details:
//! ----------------------------------
void MainWindow::toggleEdgeSelectionMode(bool isActivated)
{
    if(isActivated)
    {
        handleViewAndSelectionButtons(actionToggleEdgeSelect);
        myMainOCCViewer->setSelectionMode(CurSelection_Edge);

        //! same selection mode in additional viewports
        myDockableMasterViewPort->setSelectionMode(CurSelection_Edge);
        myDockableSlaveViewPort->setSelectionMode(CurSelection_Edge);
        statusBar()->showMessage("Select Edges",TRANSIENT_MESSAGE_TIMEOUT);
    }
}

//! ------------------------------------
//! function: toggleVertexSelectionMode
//! details:
//! ------------------------------------
void MainWindow::toggleVertexSelectionMode(bool isActivated)
{
    if(isActivated)
    {
        handleViewAndSelectionButtons(actionToggleVertexSelect);
        myMainOCCViewer->setSelectionMode(CurSelection_Vertex);

        //! same selection mode in additional viewports
        myDockableMasterViewPort->setSelectionMode(CurSelection_Vertex);
        myDockableSlaveViewPort->setSelectionMode(CurSelection_Vertex);
        statusBar()->showMessage("Select Points",TRANSIENT_MESSAGE_TIMEOUT);
    }
}

//! --------------------------------------------
//! function: togglePointCoordinatesPickingMode
//! details:
//! --------------------------------------------
void MainWindow::togglePointCoordinatesPickingMode(bool isActivated)
{
    if(isActivated)
    {
        handleViewAndSelectionButtons(actionTogglePickPointCoordinates);
        myMainOCCViewer->setSelectionMode(CurSelection_PointCoordinatesPicking);

        //! same selection mode in additional viewports
        myDockableMasterViewPort->setSelectionMode(CurSelection_PointCoordinatesPicking);
        myDockableSlaveViewPort->setSelectionMode(CurSelection_PointCoordinatesPicking);
        statusBar()->showMessage("Pick point coordinates",TRANSIENT_MESSAGE_TIMEOUT);
    }
}

//! ----------------------------------------------------
//! function: handleViewAndSelectionButtons toolbar
//! details:  when a selection button is pressed it is
//!           deactivated while the other are activated
//! ----------------------------------------------------
void MainWindow::handleViewAndSelectionButtons(QAction *action)
{
    const QList<QAction*> &listOfToolBarActions= viewAndSelectionToolbar->actions();
    for(int i=0;i<listOfToolBarActions.size();i++)
    {
        if(listOfToolBarActions.value(i)!=action)
        {
            listOfToolBarActions.value(i)->setEnabled(true);
            listOfToolBarActions.value(i)->setChecked(false);
        }
        else
        {
            listOfToolBarActions.value(i)->setEnabled(false);
        }
    }
}

//! ---------------------------------------
//! function: set the high quality display
//! details:
//! ---------------------------------------
void MainWindow::setHighQualityDisplay()
{
    myMainOCCViewer->setDisplayQuality(displayQuality_High);
}

//! ----------------------------------------------
//! function - slot: set the high quality display
//! ----------------------------------------------
void MainWindow::setMediumQualityDisplay()
{
    myMainOCCViewer->setDisplayQuality(displayQuality_Medium);
}

//! ---------------------------------------
//! function: setLowQualityDisplay
//! details:  set the high quality display
//! ---------------------------------------
void MainWindow::setLowQualityDisplay()
{
    myMainOCCViewer->setDisplayQuality(displayQuality_Low);
}

//! -------------------------------
//! function: writeSolverInputFile
//! details:  slot
//! -------------------------------
void MainWindow::writeSolverInputFile()
{
    ccout("MainWindow::writeSolverInputFile()->____function called____");
    mySimulationManager->writeSolverInputFile();
}

//! ------------------------
//! function: startAnalysis
//! details:  slot
//! ------------------------
void MainWindow::startAnalysis()
{
    cout<<"MainWindow::startAnalysis()->____function called____"<<endl;

    //! ---------------------------------------------------------------------
    //! preliminary check: before starting, the project must be saved, since
    //! what follows works if myCurrentProjectName is not empty
    //! ---------------------------------------------------------------------
    if(myCurrentProjectName.isEmpty())
    {
        QMessageBox::information(this, APPNAME,"Please, save project before run",QMessageBox::Ok);
        return;
    }

    //! ---------------------------------------------------------
    //! retrieve the absolute path of the SolutionData directory
    //! ---------------------------------------------------------
    QString project_files = myCurrentProjectName;
    project_files.chop(4);
    project_files.append("_files");
    cout<<"MainWindow::startAnalysis()->____the current saved file has name: "<<myCurrentProjectName.toStdString()<<"____"<<endl;
    cout<<"MainWindow::startAnalysis()->____path of the solution dir: "<<project_files.toStdString()<<"____"<<endl;

    //! -----------------------------
    //! clear the CCX output console
    //! -----------------------------
    mySimulationMonitor->clear();

    //! --------------------------------------------------------------------
    //! tell the SimulationMonitor where to write the "RawSolverOutput.txt"
    //! file and activate the continuos logging into that file
    //! --------------------------------------------------------------------
    SimulationNodeClass *node = mySimulationManager->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    QString timeTag = node->getPropertyValue<QString>("Time tag");

    QString rawSolverOutput = project_files+QString("/SolutionData_"+timeTag+"/RawSolverOutput.txt");
    mySimulationMonitor->setContinuoslyLogging(true, rawSolverOutput);
    cout<<"MainWindow::startAnalysis()->____path of \"RawSolverOutput.txt\" file: "<<rawSolverOutput.toStdString()<<"____"<<endl;

    //! -------------------------------------------------------------------------
    //! tell the SimulationMonitor where to write the "RawSolverOutput.txt" file
    //! and activate the continuos logging into that file
    //! -------------------------------------------------------------------------
    //QString rawSolverOutput = project_files+QString("/SolutionData/RawSolverOutput.txt");
    //mySimulationMonitor->setContinuoslyLogging(true, rawSolverOutput);
    //cout<<"MainWindow::startAnalysis()->____path of \"RawSolverOutput.txt\" file: "<<rawSolverOutput.toStdString()<<"____"<<endl;

    //! -------------------
    //! start the analysis
    //! -------------------
    bool hasStarted = mySimulationManager->startAnalysis(project_files);

    //! -------------------------------------
    //! reset the ConvergenceDataChart panel
    //! -------------------------------------
    if(hasStarted) myConvergenceDataChart->clearViewer();
}

//! -----------------------------------
//! function: setDoubleViewPortVisible
//! details:  slot
//! -----------------------------------
void MainWindow::setDoubleViewPortVisible(bool isVisible)
{
    myDockableMasterViewPort->setVisible(isVisible);
    myDockableSlaveViewPort->setVisible(isVisible);
}

//! ----------------------------------------------------------
//! function: setWorkingMode
//! details:  0 - mesh; 1 - contacts; 2 - model; 3 - solution
//! ----------------------------------------------------------
void MainWindow::setWorkingMode(int workingModeNumber)
{
    cout<<"MainWindow::setWorkingMode()->____working mode: "<<workingModeNumber<<"____"<<endl;
    switch(workingModeNumber)
    {
    case 0:
    {
        //! ------------------
        //! working mode mesh
        //! ------------------
        myMainOCCViewer->setWorkingMode_Mesh();

        //! -----------------------------------------
        //! show the tabular data - messages panel -
        //! and hide the tabular data viewer
        //! -----------------------------------------
        if(!myTabularDataDock->isVisible()) myTabularDataDock->setVisible(true);
        myTabularDataView->setCurrentIndex(1);
        if(myTabularDataViewerDock->isVisible()) myTabularDataViewerDock->setVisible(false);

        myDockableMasterViewPort->setVisible(false);
        myDockableSlaveViewPort->setVisible(false);
    }
        break;

    case 1:
    {
        //! ---------------------
        //! working mode contact
        //! ---------------------
        myMainOCCViewer->setWorkingMode_Contacts();
        /*
        SimulationNodeClass* theCurNode = mySimulationManager->getCurrentNode();
        SimulationNodeClass::nodeType nodeType = theCurNode->getType();
        switch(nodeType)
        {
        case SimulationNodeClass::nodeType_connectionPair:
            myDockableMasterViewPort->setVisible(true);
            myDockableSlaveViewPort->setVisible(true);
        break;

        default:
            //myDockableMasterViewPort->setVisible(false);
            myDockableSlaveViewPort->setVisible(false);
            break;
        }
        */
    }
        break;

    case 2:
    {
        //! -------------------
        //! working mode model
        //! -------------------
        //myDockableMasterViewPort->setVisible(false);
        //myDockableSlaveViewPort->setVisible(false);
        myMainOCCViewer->setWorkingMode_Model();

        //! --------------------------------------------------
        //! show the tabular data and the tabular data viewer
        //! --------------------------------------------------
        if(!myTabularDataDock->isVisible()) myTabularDataDock->setVisible(true);
        if(!myTabularDataViewerDock->isVisible()) myTabularDataViewerDock->setVisible(true);
    }
        break;

    case 3:
    {
        //! ----------------------
        //! working mode solution
        //! ----------------------
        myDockableMasterViewPort->setVisible(false);
        myDockableSlaveViewPort->setVisible(false);
        myMainOCCViewer->setWorkingMode_Solution();
        myTabularDataDock->setVisible(false);
        myTabularDataViewerDock->setVisible(false);
    }
        break;
    }

    //! ---------------------------------------
    //! set the working mode for the clip tool
    //! ---------------------------------------
    if(myClipTool!=Q_NULLPTR) myClipTool->setWorkingMode(workingModeNumber);
}

//! ------------------------------------------------------------
//! function - slot: check the current working mode in the menu
//! details:
//! ------------------------------------------------------------
void MainWindow::checkViewModeItem(CurDisplayMode theCurDisplayMode)
{
    switch(theCurDisplayMode)
    {
    case CurDisplayMode_Wireframe:
        actionWireframe->setChecked(true);
        break;
    case CurDisplayMode_ShadedExterior:
        actionShadedExterior->setChecked(true);
        break;
    case CurDisplayMode_ShadedExteriorAndEdges:
        actionShadedExteriorAndEdges->setChecked(true);
        break;
    }
}

//! --------------------------------------------------------
//! function: getContext
//! details:  get the context of the viewer (shape context)
//! --------------------------------------------------------
const occHandle(AIS_InteractiveContext)& MainWindow::getContext() const
{
    return (myMainOCCViewer->getContext());
}

//! ----------------------------------
//! function: disableSelectionButtons
//! details:
//! ----------------------------------
void MainWindow::disableSelectionButtons()
{
    QList<QString> names;
    names<<"Toggle solid selection"<<"Toggle face selection"<<"Toggle egde selection"<<"Toggle vertex selection";
    QList<QAction*> listOfActions = viewAndSelectionToolbar->actions();
    for(int i=0; i<listOfActions.length(); i++)
    {
        QAction* anAction = listOfActions.at(i);
        if(names.contains(anAction->text()) && anAction->isEnabled()) anAction->setEnabled(false);
    }
}

//! ---------------------------
//! function: disable a button
//! details:
//! ---------------------------
void MainWindow::disableButton(const QString &buttonName, bool disable)
{
    QList<QAction*> listOfActions = viewAndSelectionToolbar->actions();
    for(int i=0; i<listOfActions.length(); i++)
    {
        QAction* anAction = listOfActions.at(i);
        if(anAction->text() == buttonName)
        {
            if(anAction->isEnabled() && disable==true) anAction->setEnabled(false);
            if(!anAction->isEnabled() && disable==false) anAction->setEnabled(true);
        }
    }
}

//! ---------------------------------
//! function: enableSelectionButtons
//! details:
//! ---------------------------------
void MainWindow::enableSelectionButtons()
{
    QList<QString> names;
    names<<"Toggle solid selection"<<"Toggle face selection"<<"Toggle egde selection"<<"Toggle vertex selection";
    QList<QAction*> listOfActions = viewAndSelectionToolbar->actions();

    //! --------------------------------------------------------------------------------
    //! scan all the menu bar actions and enable the buttons for the geometry selection
    //! at the same time set the status of the buttons: the button for the active
    //! selection mode is down, the other are up
    //! --------------------------------------------------------------------------------
    for(int i=0; i<listOfActions.length(); i++)
    {
        QAction* anAction = listOfActions.at(i);
        if(names.contains(anAction->text()) && !anAction->isEnabled())
        {
            anAction->setEnabled(true);
        }
    }
    //! ------------------------------------
    //! reactive the current selection mode
    //! ------------------------------------
    myMainOCCViewer->reactivateSelectionMode();
}

//! --------------------------------------------------------------
//! function: startEditingScope
//! details:  simulate a double click event on the shape selector
//! --------------------------------------------------------------
void MainWindow::startEditingScope()
{
    QModelIndex index = mySimulationManager->myTreeView->currentIndex();
    SimulationNodeClass* curNode = index.data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem *itemGeometry = curNode->getPropertyItem("Geometry");
    if(itemGeometry!=NULL)
    {
        index = itemGeometry->index();
        myDetailViewer->setCurrentIndex(index);
        myDetailViewer->edit(index);
    }
}

//! --------------------------------
//! function: startEditingMagnitude
//! details:
//! --------------------------------
void MainWindow::startEditingMagnitude()
{
    //cout<<"MainWindow::startEditingMagnitude()->____function called____"<<endl;
    QModelIndex index = mySimulationManager->myTreeView->currentIndex();
    SimulationNodeClass* curNode = index.data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem *itemMagnitude = curNode->getPropertyItem("Magnitude");
    if(itemMagnitude!=NULL)
    {
        //cout<<"MainWindow::startEditingMagnitude()->____start editing Magnitude____"<<endl;
        index = itemMagnitude->index();
        myDetailViewer->edit(index);
    }
}

//! ---------------------------------
//! function: startEditingXcomponent
//! details:
//! ---------------------------------
void MainWindow::startEditingXcomponent()
{
    //cout<<"MainWindow::startEditingXcomponent()->____function called____"<<endl;
    QModelIndex index = mySimulationManager->myTreeView->currentIndex();
    SimulationNodeClass* curNode = index.data(Qt::UserRole).value<SimulationNodeClass*>();
    QExtendedStandardItem *itemXcomponent = curNode->getPropertyItem("X component");
    if(itemXcomponent!=Q_NULLPTR)
    {
        cout<<"MainWindow::startEditingMagnitude()->____start editing X component____"<<endl;
        index = itemXcomponent->index();
        myDetailViewer->edit(index);
    }
}

//! -------------------------------------------------------
//! function: addCoordinateSystemTranslationX
//! details:  ask the detail viewer to add a X translation
//! -------------------------------------------------------
void MainWindow::addCoordinateSystemTranslationX()
{
    ccout("MainWindow::addCoordinateSystemTranslationX()->____function called____");
    myDetailViewer->addPropertyToNode("Offset X",Property::PropertyGroup_Transformations);
}

//! -------------------------------------------------------
//! function: addCoordinateSystemTranslationY
//! details:  ask the detail viewer to add a Y translation
//! -------------------------------------------------------
void MainWindow::addCoordinateSystemTranslationY()
{
    ccout("MainWindow::addCoordinateSystemTranslationY()->____function called____");
    myDetailViewer->addPropertyToNode("Offset Y",Property::PropertyGroup_Transformations);
}

//! ---------------------------------------------------------------- //
//! function: addCoordinateSystemTranslationZ()                      //
//! details:  ask the detail viewer to add a Z translation           //
//! ---------------------------------------------------------------- //
void MainWindow::addCoordinateSystemTranslationZ()
{
    ccout("MainWindow::addCoordinateSystemTranslationZ()->____function called____");
    myDetailViewer->addPropertyToNode("Offset Z",Property::PropertyGroup_Transformations);
}

//! ---------------------------------------------------------------- //
//! function: addCoordinateSystemRotationX()                         //
//! details:  ask the detail viewer to add a X translation           //
//! ---------------------------------------------------------------- //
void MainWindow::addCoordinateSystemRotationX()
{
    ccout("MainWindow::addCoordinateSystemRotationX()->____function called____");
    myDetailViewer->addPropertyToNode("Rotation X",Property::PropertyGroup_Transformations);
}

//! ---------------------------------------------------------------- //
//! function: addCoordinateSystemRotationY()                         //
//! details:  ask the detail viewer to add a Y translation           //
//! ---------------------------------------------------------------- //
void MainWindow::addCoordinateSystemRotationY()
{
    ccout("MainWindow::addCoordinateSystemRotationY()->____function called____");
    myDetailViewer->addPropertyToNode("Rotation Y",Property::PropertyGroup_Transformations);
}

//! -------------------------------------------------------
//! function: addCoordinateSystemRotationZ()
//! details:  ask the detail viewer to add a Z translation
//! -------------------------------------------------------
void MainWindow::addCoordinateSystemRotationZ()
{
    ccout("MainWindow::addCoordinateSystemRotationZ()->____function called____");
    myDetailViewer->addPropertyToNode("Rotation Z",Property::PropertyGroup_Transformations);
}

//! ------------------------------------------
//! function: removeTransformation
//! details:  remove a transformation of a CS
//! ------------------------------------------
void MainWindow::removeTransformation()
{
    myDetailViewer->deleteTransformationFromNode();
}

//! ------------------------
//! function: saveProjectAs
//! details:
//! ------------------------
void MainWindow::saveProjectAs()
{
    cout<<"MainWindow::saveProjectAs()->____function called____"<<endl;

    //!--------------------------------------------------
    //! check if it is possible to save
    //! use the number of shapes loaded into the context
    //! -------------------------------------------------
    AIS_ListOfInteractive objectsInside;
    this->myMainOCCViewer->getContext()->ObjectsInside(objectsInside, AIS_KOI_Shape,-1);
    if(objectsInside.Extent()==0)
    {
        cout<<"MainWindow::saveProjectAs()->____nothing to save____"<<endl;
        return;
    }

    //! -------------------------------------
    //! "fileName" is the absolute file name
    //! open a modal dialog
    //! -------------------------------------
    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(this,"Save project as ",myWorkingDir,GIL_FILES,&selectedFilter,0);
    if(fileName.isEmpty())
    {
        cout<<"MainWindow::saveProjectAs()->____project not saved____"<<endl;
        return;
    }

    cout<<"MainWindow::saveProjectAs->____saving as: "<<fileName.toStdString()<<"____"<<endl;;

    //! -------------------
    //! relative file name
    //! -------------------
    QString relativeName = fileName.split("/").last();

    //! ------------------------------------------------------------------
    //! an already saved simulation is open, and "Save as" is called with
    //! a different name
    //! ------------------------------------------------------------------
    if(fileName!=myCurrentProjectName && !myCurrentProjectName.isEmpty())
    {
        QString myCurrentProjectName_copy = myCurrentProjectName;
        QString _filesDir = myCurrentProjectName_copy;
        _filesDir.chop(4);
        _filesDir.append("_files");
        //cout<<"____"<<_filesDir.toStdString()<<"____"<<endl;
        QDir curDir(_filesDir);
        if(curDir.exists())
        {
            QDir SDB_dir(_filesDir+("/SDB"));
            QDir MDB_dir(_filesDir+("/MDB"));
            SDB_dir.removeRecursively();
            MDB_dir.removeRecursively();
        }
    }

    //! ------------------------
    //! set the status variable
    //! ------------------------
    myCurrentProjectName = fileName;

    //! starting directory
    QString startingDirectory = fileName;
    startingDirectory.chop(relativeName.length()+1);

    QDir curDir;
    curDir.cd(startingDirectory);

    myWorkingDir = startingDirectory;
    //cout<<"MainWindow::saveProjectAs()->____start saving from dir: "<<startingDirectory.toStdString()<<"____"<<endl;

    //! -------------------------------------------
    //! build the "<relativeName>_files" directory
    //! chop the suffix
    //! -------------------------------------------
    QString _filesDir = relativeName;
    _filesDir.chop(4);
    _filesDir.append("_files");

    //! -------------------------------------------------------------
    //! check if directory "_files" already exists: if yes delete it
    //! -------------------------------------------------------------
    QFile checkDir(_filesDir);

    //! --------------------------------------------------
    //! clear the directory "_files" and create a new one
    //! --------------------------------------------------
    if(checkDir.exists()) tools::clearDir(_filesDir);

    curDir.mkdir(_filesDir);
    curDir.cd(_filesDir);

    curDir.mkdir("SDB");
    curDir.cd("SDB");
    mySimulationManager->saveSimulationDataBase(curDir.absolutePath(), relativeName);
    //mySimulationManager->setCurrentProjectDir(myWorkingDir+QString("/")+_filesDir);

    curDir.cdUp();
    curDir.cdUp();

    //cout<<"____"<<myWorkingDir.toStdString()<<"____"<<endl;
    //! ------------------------------------------------------------------------------------
    //! set the current position. At next program start "Open project" will start from here
    //! ------------------------------------------------------------------------------------
    this->updateSettings();
    //cout<<"____"<<tools::getWorkingDir().toStdString()<<"____"<<endl;
    //exit(1);

    //! -----------------------------
    //! update the main window title
    //! -----------------------------
    this->setWindowTitle(fileName);

    //cout<<"MainWindow::saveProjectAs()->____project successfully saved____"<<endl;
}

//! ----------------------
//! function: saveProject
//! details:
//! ----------------------
void MainWindow::saveProject()
{
    cout<<"MainWindow::saveProject()->____function called____"<<endl;

    //! ------------------------------------------------
    //! check if there is something to save through the
    //! number of shapes loaded into the context
    //! ------------------------------------------------
    AIS_ListOfInteractive objectsInside;
    this->myMainOCCViewer->getContext()->ObjectsInside(objectsInside, AIS_KOI_Shape,-1);
    if(objectsInside.Extent()==0)
    {
        cout<<"MainWindow::saveProject()->____nothing to save____"<<endl;
        return;
    }

    QString savingDir = myWorkingDir;
    QDir curDir;
    curDir.cd(savingDir);                                           //! moving in <myWorkingDir>

    QString fileName, fullFileName, _filesDir;
    if(myCurrentProjectName.isEmpty())
    {
        //cout<<"____"<<myCurrentProjectName.toStdString()<<"____"<<endl;
        //exit(1);
        this->saveProjectAs();
        return;
    }
    else
    {
        //! --------------------------------------------------
        //! save the current without asking: this is a policy
        //! --------------------------------------------------
        fullFileName = myCurrentProjectName;
        fileName = fullFileName.split("/").last();

        //! -----------------------------------------
        //! build the name of the "_files" directory
        //! -----------------------------------------
        _filesDir = fileName;
    }
    _filesDir.chop(4);
    _filesDir.append("_files");

    //! --------------------------------------------
    //! check if directory "_files" already exists:
    //! if yes enter it and delete SDB and MDB
    //! --------------------------------------------
    QFile checkDir(_filesDir);
    if(checkDir.exists())
    {
        QDir SDB_dir(_filesDir+("/SDB"));
        QDir MDB_dir(_filesDir+("/MDB"));
        SDB_dir.removeRecursively();
        MDB_dir.removeRecursively();
    }

    //! -------------------------------------
    //! move in <myWorkdir>/<fileName_files>
    //! -------------------------------------
    curDir.cd(_filesDir);

    //! ---------------------------------------------------------------------------------------------
    //! create a subdirectory "<myWorkingDir>/<fileName_files>/SDB" for storing the simulation nodes
    //! ---------------------------------------------------------------------------------------------
    cout<<"MainWindow::saveProject()->____creating directory SDB____"<<endl;

    curDir.mkdir("SDB");
    curDir.cd("SDB");                                                //! in "../../untitled_files/SDB"
    QString SDBdir = curDir.absolutePath();                          //! SDBdir = "../../untitled_files/SDB"
    curDir.cdUp();                                                   //! in "untitled_files"

    //! --------------------------------------------------
    //! store the data in SDBdir
    //! fileName is the relative name of the project file
    //! --------------------------------------------------
    mySimulationManager->saveSimulationDataBase(SDBdir,fileName);    //! writing tree items into "../../untitled_files/SDB"
    curDir.cdUp();                                                   //! in "../../untitled_files"
    curDir.cdUp();                                                   //! in "../.."

    //! -----------------------
    //! update status variable
    //! -----------------------
    myCurrentProjectName = fullFileName;

    //! ------------------------
    //! update the window title
    //! ------------------------
    this->setWindowTitle(fullFileName);

    cout<<"MainWindow::saveProject()->____project successfully saved____"<<endl;
}

//! ----------------------------------------
//! function: callExtendSelectionToAdjacent
//! details:
//! ----------------------------------------
void MainWindow::callExtendSelectionToAdjacent()
{
    //cout<<"MainWindow::callExtendSelectionToAdjacent()->____function called____"<<endl;
    myExtendSelectionButton->setIcon(QIcon(":/icons/icon_extend to adjacent.png"));
    emit requestExtendSelectionToAdjacent();
}

//! -------------------------------------------
//! function: callExtendSelectionToConnections
//! details:
//! -------------------------------------------
void MainWindow::callExtendSelectionToConnections()
{
    //cout<<"MainWindow:.callExtendSelectionToConnections->____function called____"<<endl;
    myExtendSelectionButton->setIcon(QIcon(":/icons/icon_extend to limit.png"));
    emit requestExtendSelectionToConnections();
}

//! -------------------------------------------
//! function: callExtendSelectionToConnections
//! details:
//! -------------------------------------------
void MainWindow::callExtendSelectionToLimits()
{
    //cout<<"MainWindow:.callExtendSelectionToLimits->____function called____"<<endl;
    myExtendSelectionButton->setIcon(QIcon(":/icons/icon_extend to limit.png"));
    emit requestExtendSelectioToLimits();
}

//! -----------------------------------
//! function: handleColorBoxVisibility
//! details:
//! -----------------------------------
void MainWindow::handleColorBoxVisibility(bool isVisible)
{
    cout<<"MainWindow::handleColorBoxVisibility->___function called. Colorbox visible: "<<isVisible<<"____"<<endl;
    myMainOCCViewer->displayColorBox(isVisible);
}

//! ----------------------------
//! function: showOptionsWidget
//! details:
//! ----------------------------
void MainWindow::showOptionsWidget()
{
    optionsWidget *a = new optionsWidget();
    a->setWindowModality(Qt::WindowModal);
    a->show();
}

//! ------------------
//! function: myClose
//! details:
//! ------------------
void MainWindow::myClose()
{
    cout<<"MainWindow::myClose()->____function called____"<<endl;

    //! ------------------------------------------------------------------------
    //! the application working directory ("current directory") must be changed
    //! in order to remove both the <project>_files directory content, and both
    //! the directory itself
    //! ------------------------------------------------------------------------

    QDir::setCurrent(tools::getWorkingDir());
    QString dirToDelete = myCurrentProjectName;
    dirToDelete.chop(4);
    dirToDelete.append("_files");

    QDir curDir(dirToDelete);

    cout<<"MainWindow::myClose()->____application current directory: "<<QDir::current().absolutePath().toStdString()<<"____"<<endl;
    cout<<"MainWindow::myClose()->____directory to be removed: "<<dirToDelete.toStdString()<<"____"<<endl;

    bool isDone = curDir.removeRecursively();

    //! alternative to the previous line
    //tools::clearDir(dirToDelete);

    curDir.cdUp();
    cout<<"MainWindow::myClose()->____going up to dir: "<<curDir.absolutePath().toStdString()<<"____"<<endl;

    //bool isDone = curDir.rmdir(dirToDelete);

    if(isDone) cout<<"MainWindow::myClose()->____directory: "<<dirToDelete.toStdString()<<" successfully removed____"<<endl;
    else cerr<<"MainWindow::myClose()->____directory: "<<dirToDelete.toStdString()<<" not removed____"<<endl;
    this->close();
}

//! ---------------------
//! function: closeEvent
//! details:
//! ---------------------
void MainWindow::closeEvent(QCloseEvent *event)
{
    cout<<"MainWindow::closeEvent()->____function called____"<<endl;
    QMessageBox::StandardButton ret = QMessageBox::question(this, tr(APPNAME),tr("Close application"),QMessageBox::Close | QMessageBox::Cancel);
    switch(ret)
    {
    case QMessageBox::Close:
    {
        cout<<"MainWindow::closeEvent()->____accept close____"<<endl;
        event->accept();

        //! ------------------------------------------------
        //! MainWindow::myClose() delete "<name>_files" dir
        //! ------------------------------------------------
        this->myClose();
    }
        break;

    case QMessageBox::Cancel:
    {
        cout<<"MainWindow::closeEvent()->____ignore close____"<<endl;
        event->ignore();
    }
        break;
    }
}

//! --------------------------
//! function: HandleSelectAll
//! details:
//! --------------------------
void MainWindow::HandleSelectAll()
{
    cout<<"MainWindow::HandleSelectAll()->____function called____"<<endl;
    myMainOCCViewer->selectAll();
}

//! ---------------------------
//! function: setUpConnections
//! details:
//! ---------------------------
void MainWindow::setUpConnections()
{
    connect(myCentralTabWidget,SIGNAL(resized()),myMainOCCViewer,SLOT(repaint()));
    connect(mySimulationManager,SIGNAL(requestSetActiveCentralTab(QString)),myCentralTabWidget,SLOT(setCurrentTab(QString)));
    connect(mySimulationManager, SIGNAL(requestMeshInvalidate(std::vector<int>)),myMainOCCViewer, SLOT(invalidateMeshes(std::vector<int>)));

    connect(mySimulationManager,SIGNAL(requestUpdateConvergenceViewer(const QList<solutionInfo> &)), myConvergenceDataChart,SLOT(plotConvergenceData(const QList<solutionInfo> &)));

    //! ----------------
    //! selection modes
    //! ----------------
    connect(myMainOCCViewer,SIGNAL(selectionModeVertex(bool)),this,SLOT(toggleVertexSelectionMode(bool)));
    connect(myMainOCCViewer,SIGNAL(selectionModeEdge(bool)),this,SLOT(toggleEdgeSelectionMode(bool)));
    connect(myMainOCCViewer,SIGNAL(selectionModeFace(bool)),this,SLOT(toggleFaceSelectionMode(bool)));
    connect(myMainOCCViewer,SIGNAL(selectionModeSolid(bool)),this,SLOT(toggleSolidSelectionMode(bool)));
    connect(myMainOCCViewer,SIGNAL(selectionModePickPointCoordinates(bool)),this,SLOT(togglePointCoordinatesPickingMode(bool)));

    connect(myMainOCCViewer,SIGNAL(statusBarMessage(QString)),statusLabel,SLOT(setText(QString)));
    connect(myMainOCCViewer,SIGNAL(viewModeChanged(CurDisplayMode)),this,SLOT(checkViewModeItem(CurDisplayMode)));
    connect(this,SIGNAL(requestExtendSelectionToAdjacent()),myMainOCCViewer,SLOT(extendSelectionToAjacent()));

    //! ------------------------------------------------------------------------------
    //! synchronize the checkmark of the menu action and the visibility of the widget
    //! ------------------------------------------------------------------------------
    connect(actionShowGraphViewer,SIGNAL(triggered(bool)),myTabularDataViewerDock,SLOT(setVisible(bool)));

    //! ------------------------------------------------------------------------------
    connect(myTabularDataViewerDock,SIGNAL(visibilityChanged(bool)),actionShowGraphViewer,SLOT(setChecked(bool)));

    //! ------------------------------------------------------------------------------
    //! synchronize the checkmark of the menu action and the visibility of the widget
    //! ------------------------------------------------------------------------------
    connect(actionShowTabularData,SIGNAL(triggered(bool)),myTabularDataDock,SLOT(setVisible(bool)));
    connect(myTabularDataDock,SIGNAL(visibilityChanged(bool)),actionShowTabularData,SLOT(setChecked(bool)));

    connect(actionShowSimulationManager,SIGNAL(triggered(bool)),mySimulationManagerDock,SLOT(setVisible(bool)));
    connect(mySimulationManagerDock,SIGNAL(visibilityChanged(bool)),actionShowSimulationManager,SLOT(setChecked(bool)));

    connect(actionShowDebugWindow,SIGNAL(triggered(bool)),myDebugConsoleDock,SLOT(setVisible(bool)));
    connect(myDebugConsoleDock,SIGNAL(visibilityChanged(bool)),actionShowDebugWindow,SLOT(setChecked(bool)));

    connect(mySimulationManager,SIGNAL(requestClearMesh()),myMainOCCViewer,SLOT(clearMeshFromViewer()));
    connect(myMainOCCViewer,SIGNAL(requestGenerateMesh(bool)),mySimulationManager,SLOT(buildMesh(bool)));
    connect(myMainOCCViewer,SIGNAL(requestCreateSimulationNode(SimulationNodeClass::nodeType)),mySimulationManager,SLOT(createSimulationNode(SimulationNodeClass::nodeType)));
    connect(myMainOCCViewer,SIGNAL(requestChangeSuppressionStatus(Property::SuppressionStatus)),mySimulationManager,SLOT(changeNodeSuppressionStatus(Property::SuppressionStatus)));

    //! ------------------------------------------------------------
    //! color for boundary conditions, contacts, mesh controls, ...
    //! ------------------------------------------------------------
    connect(mySimulationManager,SIGNAL(requestCustomColor(TopTools_ListOfShape,Quantity_NameOfColor,bool,bool)),
            myMainOCCViewer,SLOT(displayColoredSubshapes(TopTools_ListOfShape,Quantity_NameOfColor,bool,bool)));

    connect(mySimulationManager,SIGNAL(request3DBodySelectionMode(bool)),this,SLOT(toggleSolidSelectionMode(bool)));
    connect(mySimulationManager,SIGNAL(request2DBodySelectionMode(bool)),this,SLOT(toggleFaceSelectionMode(bool)));
    connect(mySimulationManager,SIGNAL(request1DBodySelectionMode(bool)),this,SLOT(toggleEdgeSelectionMode(bool)));
    connect(mySimulationManager,SIGNAL(request0DBodySelectionMode(bool)),this,SLOT(toggleVertexSelectionMode(bool)));
    connect(mySimulationManager,SIGNAL(requestBuildMeshIOs()),myMainOCCViewer,SLOT(buildMeshIOs()));
    connect(mySimulationManager,SIGNAL(requestHideMeshes()),myMainOCCViewer,SLOT(hideAllMeshes()));
    connect(mySimulationManager,SIGNAL(requestShowMeshes(bool)),myMainOCCViewer,SLOT(displayAllMeshes(bool)));
    connect(mySimulationManager,SIGNAL(requestSetWorkingMode(int)),this,SLOT(setWorkingMode(int)));
    connect(mySimulationManager,SIGNAL(requestHideBody(TColStd_ListOfInteger)),myMainOCCViewer,SLOT(hideBody(TColStd_ListOfInteger)));
    connect(mySimulationManager,SIGNAL(requestRemoveObsoleteMeshes()),myMainOCCViewer,SLOT(removeObsoleteMeshes()));
    connect(mySimulationManager,SIGNAL(requestShowBody(TColStd_ListOfInteger)),myMainOCCViewer,SLOT(showBody(TColStd_ListOfInteger)));
    connect(myMainOCCViewer,SIGNAL(requestSynchVisibility()),mySimulationManager,SLOT(synchVisibility()));

    //! ---------------------------
    //! body highlight/unhighlight
    //! ---------------------------
    connect(mySimulationManager,SIGNAL(requestHighlightBody(std::vector<int>)),myMainOCCViewer,SLOT(highlightBody(std::vector<int>)));
    connect(mySimulationManager,SIGNAL(requestUnhighlightBodies(bool)),myMainOCCViewer,SLOT(unhighlightBody(bool)));


    connect(mySimulationManager,SIGNAL(requestReactivateSelectionMode()),myMainOCCViewer,SLOT(reactivateSelectionMode()));

    connect(mySimulationManager,SIGNAL(requestDisplayShapeCopy(TopTools_ListOfShape,TopTools_ListOfShape,Quantity_NameOfColor,Quantity_NameOfColor,QVariant)),
            myMainOCCViewer,SLOT(displayShapeCopy(TopTools_ListOfShape,TopTools_ListOfShape,Quantity_NameOfColor,Quantity_NameOfColor,QVariant)));

    connect(mySimulationManager,SIGNAL(requestTabularData(QModelIndex)),myTabularDataView,SLOT(setTheModel(QModelIndex)));
    connect(mySimulationManager,SIGNAL(requestShowColumns(QList<int>)),myTabularDataView,SLOT(showColumns(QList<int>)));

    connect(mySimulationManager,SIGNAL(requestShowGraph(CustomTableModel*,QList<int>)),myTabularDataGraphViewer1,SLOT(showData(CustomTableModel*,QList<int>)));
    connect(mySimulationManager,SIGNAL(requestClearGraph()),myTabularDataGraphViewer1,SLOT(clearPanel()));

    connect(mySimulationManager,SIGNAL(requestHideFirstRow()),myTabularDataView,SLOT(hideFirstRow()));
    connect(mySimulationManager,SIGNAL(requestShowFirstRow()),myTabularDataView,SLOT(showFirstRow()));

    connect(mySimulationManager,SIGNAL(requestStartEditingScope()),this,SLOT(startEditingScope()));
    connect(mySimulationManager,SIGNAL(requestStartEditingMagnitude()),this,SLOT(startEditingMagnitude()));
    connect(mySimulationManager,SIGNAL(requestStartEditingXcomponent()),this,SLOT(startEditingXcomponent()));

    connect(mySimulationManager,SIGNAL(requestDisplayTrihedron(QVector<double>,QVector<QVector<double>>,int)),
            myMainOCCViewer,SLOT(displayTrihedron(QVector<double>,QVector<QVector<double>>,int)));

    connect(mySimulationManager,SIGNAL(requestCreateColorBox(double,double,int)),myMainOCCViewer,SLOT(createColorBox(double,double,int)));
    connect(mySimulationManager,SIGNAL(requestShowAllBodies()),myMainOCCViewer,SLOT(showAllBodies()));

    connect(mySimulationManager,SIGNAL(requestDisplayResult(sharedPostObject&)),myMainOCCViewer,SLOT(displayResult(sharedPostObject&)));
    connect(mySimulationManager,SIGNAL(requestHideAllResults()),myMainOCCViewer,SLOT(hideAllResults()));

    //! ---------------------------------------------------
    //! connection for hiding the markers (and the triads)
    //! ---------------------------------------------------
    connect(mySimulationManager,SIGNAL(requestHideAllMarkers(bool)),myMainOCCViewer,SLOT(hideAllMarkers(bool)));
    connect(myDetailViewer,SIGNAL(requestHideAllMarkers(bool)),myMainOCCViewer,SLOT(hideAllMarkers(bool)));

    connect(myMainOCCViewer,SIGNAL(requestShowHealingElements()),mySimulationManager,SLOT(showHealingElements()));

    connect(myMainOCCViewer,SIGNAL(requestDisableSelectionButtons()),this,SLOT(disableSelectionButtons()));
    connect(myMainOCCViewer,SIGNAL(requestEnableSelectionButtons()),this,SLOT(enableSelectionButtons()));

    connect(myMainOCCViewer,SIGNAL(requestChangeSuppressionStatus(Property::SuppressionStatus)),mySimulationManager,SLOT(changeNodeSuppressionStatus(Property::SuppressionStatus)));
    connect(mySimulationManager,SIGNAL(requestUpdateMeshView()),myMainOCCViewer,SLOT(updateMeshContext()));

    connect(actionShowDetailViewer,SIGNAL(triggered(bool)),myDetailViewerDock,SLOT(setVisible(bool)));
    connect(myDetailViewerDock,SIGNAL(visibilityChanged(bool)),actionShowDetailViewer,SLOT(setChecked(bool)));

    connect(actionShowSectionPlanes,SIGNAL(triggered(bool)),myClipToolDock,SLOT(setVisible(bool)));

    connect(myClipToolDock,SIGNAL(visibilityChanged(bool)),actionShowSectionPlanes,SLOT(setChecked(bool)));

    connect(actionShowMemoryProfiler,SIGNAL(triggered(bool)),myMemoryProfiler,SLOT(setVisible(bool)));
    connect(myMemoryProfiler,SIGNAL(visibilityChanged(bool)),actionShowMemoryProfiler,SLOT(setChecked(bool)));

    //! -------------------------------------------------------------------------
    //! click on an item of the simulation manager => activate the detail viewer
    //! -------------------------------------------------------------------------
    connect(mySimulationManager,SIGNAL(updatedetailViewer(QModelIndex)),myDetailViewer,SLOT(setTheModel(QModelIndex)));

    connect(myDetailViewer,SIGNAL(requestChangeNodeSuppressionStatus(Property::SuppressionStatus)),mySimulationManager,SLOT(changeNodeSuppressionStatus(Property::SuppressionStatus)));
    connect(myDetailViewer,SIGNAL(requestHandleVisibilityChange(bool)),mySimulationManager,SLOT(handleVisibilityChange(bool)));
    connect(myDetailViewer,SIGNAL(requestChangeElementControl()),mySimulationManager,SLOT(ChangeElementControl()));

    //! -------------
    //! experimental
    //! -------------
    connect(myDetailViewer,SIGNAL(requestChangeColor()),mySimulationManager,SLOT(changeColor()));

    connect(myDetailViewer,SIGNAL(requestResizeTabularData()),mySimulationManager,SLOT(resizeTabularData()));
    connect(myTabularDataView,SIGNAL(requestUpdateDetailViewer(QModelIndex,QModelIndex,QVector<int>)),myDetailViewer,SLOT(updateDetailViewerFromTabularData(QModelIndex,QModelIndex,QVector<int>)));

    //! --------------------------------------------------------------------
    //! the geometry selection has changed => the DetailViewer is informed:
    //! - when choosing a direction an arrow is shown
    //! --------------------------------------------------------------------
    connect(myMainOCCViewer,SIGNAL(selectionChanged()),myDetailViewer,SLOT(selectionHasChanged()));

    connect(myDetailViewer,SIGNAL(requestDisplayTrihedron(QVector<double>,QVector<QVector<double>>,int)),
            myMainOCCViewer,SLOT(displayTrihedron(QVector<double>,QVector<QVector<double> >,int)));

    connect(myDetailViewer,SIGNAL(requestHandleTabularData()),mySimulationManager,SLOT(HandleTabularData()));
    connect(myDetailViewer,SIGNAL(requestHandleFilmCoefficientLoadDefinitionChanged(QString)),mySimulationManager,SLOT(handleFilmCoefficientLoadDefinitionChanged(QString)));
    connect(myDetailViewer,SIGNAL(requestHandleReferenceTemperatureLoadDefinitionChanged(QString)),mySimulationManager,SLOT(handleReferenceTemperatureLoadDefinitionChanged(QString)));
    connect(myDetailViewer,SIGNAL(requestHandleMagnitudeLoadDefinitionChanged(QString)),mySimulationManager,SLOT(handleLoadMagnitudeDefinitionChanged(QString)));
    connect(myDetailViewer,SIGNAL(requestHandleXLoadDefinitionChanged(QString)),mySimulationManager,SLOT(handleLoadXDefinitionChanged(QString)));
    connect(myDetailViewer,SIGNAL(requestHandleYLoadDefinitionChanged(QString)),mySimulationManager,SLOT(handleLoadYDefinitionChanged(QString)));
    connect(myDetailViewer,SIGNAL(requestHandleZLoadDefinitionChanged(QString)),mySimulationManager,SLOT(handleLoadZDefinitionChanged(QString)));
    connect(myDetailViewer,SIGNAL(requestChangeMeshNodesVisibility(bool)),myMainOCCViewer,SLOT(displayAllMeshes(bool)));
    connect(myDetailViewer,SIGNAL(requestInvalidateAllMeshes()),myMainOCCViewer,SLOT(invalidateAllMeshes()));
    connect(myDetailViewer,SIGNAL(requestGlobalMeshControlChange()),mySimulationManager,SLOT(handleGlobalMeshControlChange()));
    connect(myDetailViewer,SIGNAL(requestHandleSolutionInformationUpdateIntervalChanged()),mySimulationManager,SLOT(updateTimer()));

    // to be removed
    //connect(myDetailViewer,SIGNAL(requestHandleBoltControls()),mySimulationManager,SLOT(handleBoltControls()));

    connect(myDetailViewer,SIGNAL(requestHandleTransparencyChanged(double)),myMainOCCViewer,SLOT(setTransparency_onWorkingModeContact(double)));

    //! ------------------------
    //! 3D interpolator, mapper
    //! ------------------------
    connect(myDetailViewer,SIGNAL(startInterpolator()),mySimulationManager,SLOT(interpolate()));
    connect(myDetailViewer,SIGNAL(startOpenFoamScalarDataTranslator()),mySimulationManager,SLOT(translateOpenFoamScalarData()));

    //! ----------------
    //! exporting tools
    //! ----------------
    connect(myMainOCCViewer,SIGNAL(requestExportStepFile()),mySimulationManager,SLOT(exportSTEPFile()));

    //! -----------------------------------------------------------------------
    //! the following connections activate the double viewport functionalities
    //! -----------------------------------------------------------------------
    //! ----------------------------------------------
    //! connections for updating the delegate context
    //! ----------------------------------------------
    connect(myDockableMasterViewPort,SIGNAL(requestChangeDelegateContext(opencascade::handle<AIS_InteractiveContext>)),myDetailViewer,SLOT(setContext(opencascade::handle<AIS_InteractiveContext>)));
    connect(myDockableSlaveViewPort,SIGNAL(requestChangeDelegateContext(opencascade::handle<AIS_InteractiveContext>)),myDetailViewer,SLOT(setContext(opencascade::handle<AIS_InteractiveContext>)));
    connect(myMainOCCViewer,SIGNAL(requestChangeDelegateContext(opencascade::handle<AIS_InteractiveContext>)),myDetailViewer,SLOT(setContext(opencascade::handle<AIS_InteractiveContext>)));

    //! -----------------------------------------------------------------------------------
    //! connection for displaying a selection of bodies on the master and slave view ports
    //! -----------------------------------------------------------------------------------
    connect(mySimulationManager,SIGNAL(requestShowBodyOnMasterViewPort(TColStd_ListOfInteger, TColStd_ListOfInteger)),this,SLOT(showBodiesOnMasterViewPort(TColStd_ListOfInteger, TColStd_ListOfInteger)));
    connect(mySimulationManager,SIGNAL(requestShowBodyOnSlaveViewPort(TColStd_ListOfInteger, TColStd_ListOfInteger)),this,SLOT(showBodiesOnSlaveViewPort(TColStd_ListOfInteger, TColStd_ListOfInteger)));

    //! ---------------------------------------------------------
    //! connection for show/hide the master and slave view ports
    //! ---------------------------------------------------------
    connect(mySimulationManager,SIGNAL(requestShowDoubleViewPort(bool)),this,SLOT(setDoubleViewPortVisible(bool)));

    //! ----------------------------------------------------------
    //! connection for hide bodies from master and slave vieports
    //! ----------------------------------------------------------
    connect(mySimulationManager,SIGNAL(requestHideBodiesFromViewers(TColStd_ListOfInteger)),this,SLOT(clearMasterSlaveViewPorts(TColStd_ListOfInteger)));

    //! -----------------------------------------------------------------
    //! connections for showing the master surfaces on the master viewer
    //! and the slave surfaces on the slave viewer
    //! -----------------------------------------------------------------
    connect(mySimulationManager,SIGNAL(requestDisplayShapeCopy1(TopTools_ListOfShape,Quantity_NameOfColor)),myDockableMasterViewPort,SLOT(displayShapeCopy1(TopTools_ListOfShape,Quantity_NameOfColor)));
    connect(mySimulationManager,SIGNAL(requestDisplayShapeCopy2(TopTools_ListOfShape,Quantity_NameOfColor)),myDockableSlaveViewPort,SLOT(displayShapeCopy2(TopTools_ListOfShape,Quantity_NameOfColor)));

    //! --------------------------------
    //! connection for showing a marker
    //! --------------------------------
    //connect(mySimulationManager,SIGNAL(requestDisplaySphericalMarker(gp_Pnt)),myMainOCCViewer,SLOT(displaySphericalMarker(gp_Pnt)));

    //! --------------------------
    //! mesh tool bar connections
    //! --------------------------
    connect(meshToolBar,SIGNAL(showExteriorMeshRequest(bool)),myMainOCCViewer,SLOT(refreshMeshView(bool)));
    connect(meshToolBar,SIGNAL(requestClearMesh()),myMainOCCViewer,SLOT(clearMeshFromViewer()));
    connect(meshToolBar,SIGNAL(requestGenerateVolumeMesh()),mySimulationManager,SLOT(buildVolumeMesh()));
    connect(meshToolBar,SIGNAL(requestGenerateSurfaceMesh()),mySimulationManager,SLOT(buildSurfaceMesh()));

    //! -------------------------------------------------------
    //! view port update when changing the simulation database
    //! -------------------------------------------------------
    connect(mySimulationManager,SIGNAL(requestUpdateViewport()),this,SLOT(updateViewport()));


#ifdef COSTAMP_VERSION
    connect(myDetailViewer,SIGNAL(requestStartTimeStepBuilder()),mySimulationManager,SLOT(COSTAMP_startTimeStepBuilder()));
#endif

    //! -----------------------------
    //! experimental - custom mesher
    //! -----------------------------
    connect(myMainOCCViewer,SIGNAL(requestBuildFaceMesh(TopoDS_Face)),mySimulationManager,SLOT(customMesherBuildFaceMesh(TopoDS_Face)));

    //! -------------------
    //! refresh the viewer
    //! -------------------
    connect(mySimulationManager,SIGNAL(requestRefreshViewer()),myMainOCCViewer,SLOT(updateViewerAfterDataBaseChange()));

    //! -----------------------------------------
    //! "Show underlying mesh" from context menu
    //! -----------------------------------------
    connect(mySimulationManager,SIGNAL(requestShowMeshDataSources(IndexedMapOfMeshDataSources)),myMainOCCViewer,SLOT(showMeshDataSources(IndexedMapOfMeshDataSources)));

    //! --------------------------------------------------
    //! set the geometry/mesh data base for the clip tool
    //! --------------------------------------------------
    connect(mySimulationManager,SIGNAL(requestSetClipperDataBase(meshDataBase*)),myClipTool,SLOT(setMeshDataBase(meshDataBase*)));

    //! -------------------
    //! mesh quality graph
    //! -------------------
    connect(myDetailViewer,SIGNAL(requestHandleMeshMetricChanged()),mySimulationManager,SLOT(computeAndDisplayMeshMetric()));
    connect(mySimulationManager,SIGNAL(requestDisplayMeshMetricHystogram(const histogramData&)),this,SLOT(displayMeshMetric(const histogramData&)));

    //! -----------------------------------
    //! on click on "Solution information"
    //! -----------------------------------
    connect(mySimulationManager,SIGNAL(requestDisplayTextOnConsole(QString)),mySimulationMonitor,SLOT(setText(QString)));

    //! --------------------------------
    //! experimental - new change color
    //! --------------------------------
    connect(mySimulationManager,SIGNAL(requestResetCustomColors(bool)),myMainOCCViewer,SLOT(resetCustomColors(bool)));
    connect(mySimulationManager,SIGNAL(requestApplyCustomColor(QMap<GeometryTag,TopoDS_Shape>,Quantity_NameOfColor,bool)),myMainOCCViewer,SLOT(applyCustomColors(QMap<GeometryTag,TopoDS_Shape>,Quantity_NameOfColor,bool)));
}

//! -------------------------------------------------------------
//! function: readResultsFile
//! details:  the results of the read operation will be put into
//!           "<theCurrentProjectFileName>_files\SolutionData\"
//! -------------------------------------------------------------
void MainWindow::readResultsFile()
{
    cout<<"MainWindow::readResultsFile()->____function called____"<<endl;

    SimulationNodeClass *curNode = mySimulationManager->myTreeView->currentIndex().data(Qt::UserRole).value<SimulationNodeClass*>();
    if(curNode->isAnalysisRoot()==false)
    {
        QMessageBox::information(this,APPNAME,"Select an analysis root item",QMessageBox::Ok);
        return;
    }

    //if(myGeometryImportStatus == geometryImportStatus_Loaded)
    if(0==0)
    {
        QString selectedFilter;
        QString fileName = QFileDialog::getOpenFileName(this,"Select file to read",tools::getWorkingDir(),CCX_FILES,&selectedFilter,0);

        if(fileName.isEmpty()) return;

        QString projectName = myCurrentProjectName;
        projectName.chop(4);
        QString files_dir = projectName+"_files";

        cout<<"____<_files> directory: "<<files_dir.toStdString()<<"____"<<endl;

        //! ---------------------------------------------
        //! check if the "ResultsData" directory exists
        //! ---------------------------------------------
        QDir dir(files_dir);

        QString timeTag = curNode->getPropertyValue<QString>("Time tag");
        QString solutionDataDir = files_dir+"/"+QString("SolutionData_")+timeTag;
        QString resultsDataDir = files_dir+"/"+QString("SolutionData_")+timeTag+"/"+QString("ResultsData");

        cout<<"____resultsDataDir: "<<resultsDataDir.toStdString()<<"____"<<endl;

        QDir d(resultsDataDir);
        if(d.exists())
        {
            int res = QMessageBox::warning(this,APPNAME,"Previous results data for this simulation found.\n"
                                                        "Do you want to overwrite?",QMessageBox::Ok,QMessageBox::Cancel);
            if(res == QMessageBox::Cancel) return;
        }

        //! ---------------------
        //! delete the directory
        //! ---------------------
        d.removeRecursively();

        //! --------------------------------------------------
        //! create a new directory within the <>_files folder
        //! --------------------------------------------------
        if(!dir.cd(resultsDataDir))
        {
            cout<<"MainWindow::readResultsFile()->____the \"SolutionData\" does not exist: creating it____"<<endl;
            dir.mkdir(resultsDataDir);
            dir.cd(resultsDataDir);
        }

        mySimulationManager->readResultsFile(fileName, solutionDataDir);
    }
    else
    {
        QMessageBox::information(this,"Information","A valid model has not been created or imported yet",QMessageBox::Ok);
    }
}

//! ----------------------------
//! function: displayMeshMetric
//! details:
//! ----------------------------
void MainWindow::displayMeshMetric(const histogramData &hysData)
{
    cout<<"MainWindow::displayMeshMetric()->____function called____"<<endl;
    myHistogramDockContainer->setVisible(true);
    myHistogram->setData(hysData,true);

    //! --------------------------------
    //! position and size on the screen
    //! --------------------------------
    const int POSX = 100;
    const int POSY = 100;
    const int W = 960;
    const int H = 540;

    QPoint tl(POSX,POSY);
    QPoint br(POSX+W,POSY+H);
    QRect geo(tl,br);
    myHistogramDockContainer->setGeometry(geo);
}
