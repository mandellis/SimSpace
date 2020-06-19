#include "timestepbuilder.h"
#include "qcustomplotextended.h"

//! Qt
#include <QMenu>
#include <QStatusBar>
#include <QMenuBar>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QAction>
#include <QFileDialog>
#include <QDir>

//! C++
#include <iostream>
using namespace std;

//! --------------------------------
//! position and size on the screen
//! --------------------------------
#define POSX 100
#define POSY 100
#define W 960
#define H 540

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TimeStepBuilder::TimeStepBuilder(QWidget *parent, Qt::WindowFlags f):QDialog(parent,f),myWorkingDir(QDir::currentPath())
{
    this->createActions();
    this->createContent();    
    connect(myPlot,SIGNAL(tracerInterceptsTime(double)),this,SLOT(updateStatusBar(double)));
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TimeStepBuilder::TimeStepBuilder(const QString &currentDirectory, QWidget *parent,Qt::WindowFlags f):QDialog(parent, f),myWorkingDir(currentDirectory)
{
    this->createActions();
    this->createContent();
    connect(myPlot,SIGNAL(tracerInterceptsTime(double)),this,SLOT(updateStatusBar(double)));
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
TimeStepBuilder::~TimeStepBuilder()
{
    //cout<<"____destructor called____"<<endl;
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void TimeStepBuilder::createContent()
{
    cout<<"TimeStepBuilder::createContent()->____function called____"<<endl;
    //! ---------------
    //! set the layout
    //! ---------------
    //QVBoxLayout *mainLayout = new QVBoxLayout(this);
    QVBoxLayout *mainLayout = new QVBoxLayout();

    mainLayout->setContentsMargins(0,0,0,0);
    mainLayout->setMargin(0);
    mainLayout->setSpacing(0);
    this->setLayout(mainLayout);
    this->setSizeGripEnabled(true);

    //! ---------
    //! menu bar
    //! ---------
    myMenuBar = new QMenuBar(this);
    mainLayout->setMenuBar(myMenuBar);

    menuFile = myMenuBar->addMenu("File");

    QList<QAction*> actionList;
    actionList<<actionOpen<<actionClear<<actionAccept<<actionExit;

    //! -------------------
    //! create the toolbar
    //! -------------------
    myToolBar = new QToolBar("Tool bar",this);
    myToolBar->setVisible(true);
    myToolBar->setMovable(true);

    myToolBar->addActions(actionList);

    actionList.clear();
    actionList<<actionZoom<<actionPan<<actionFit;
    myToolBar->addActions(actionList);

    myToolBar->addSeparator();

    actionList.clear();
    actionList<<actionSaveImage<<actionFitRange;

    myToolBar->addActions(actionList);

    mainLayout->addWidget(myToolBar,0,Qt::AlignTop);

    //! ----------------
    //! create the plot
    //! ----------------
    cout<<"____tag00____"<<endl;
    myPlot = new QCustomPlotExtended(this);
    cout<<"____tag01____"<<endl;

    myPlot->setMouseTracking(true);
    myPlot->axisRect()->setupFullAxesBox();
    myPlot->setVisible(false);
    myPlot->setContentsMargins(0,0,0,0);

    mainLayout->addWidget(myPlot,Qt::AlignTop);

    //! ----------------------------
    //! add the actions to the menu
    //! ----------------------------
    menuFile->addAction(actionOpen);
    menuFile->addAction(actionClear);
    menuFile->addAction(actionAccept);
    menuFile->addAction(actionExit);

    actionList.clear();
    actionList<<actionOpen<<actionClear<<actionAccept<<actionExit;

    //! -----------
    //! status bar
    //! -----------
    myStatusBar = new QStatusBar(this);
    QSizePolicy sizePolicy = myStatusBar->sizePolicy();
    sizePolicy.setVerticalPolicy(QSizePolicy::Fixed);
    myStatusBar->setSizePolicy(sizePolicy);

    mainLayout->addWidget(myStatusBar);

    //! ------------------------
    //! create the context menu
    //! ------------------------
    myContextMenu = new QMenu(this);
    myContextMenu->addAction(actionSaveImage);
    connect(this,SIGNAL(customContextMenuRequested(QPoint)),this,SLOT(showContextMenu(QPoint)));

    QPoint tl(POSX,POSY);
    QPoint br(POSX+W,POSY+H);
    QRect geo(tl,br);
    this->setGeometry(geo);
}

//! -------------------
//! function: openFile
//! details:
//! -------------------
bool TimeStepBuilder::openFile()
{
    cout<<"____openFile____"<<endl;
    QString selectedFilter;
    timeHistoryFile = QFileDialog::getOpenFileName(this,"Select time history file", QDir::current().absolutePath(),"*.txt",&selectedFilter);
    if(timeHistoryFile.isNull() || timeHistoryFile.isEmpty())
    {
        //cout<<"____no file opened____"<<endl;
        return false;
    }
    FILE *timeHistoryFilePointer = fopen(timeHistoryFile.toStdString().c_str(),"r");
    if(timeHistoryFilePointer==NULL)
    {
        //cout<<"____the time history file is empty____";
        return false;
    }

    //! ----------------
    //! clear plot data
    //! ----------------
    timeList.clear();
    valueList.clear();

    //! ----------------------------
    //! read the header of the file
    //! ----------------------------
    firstColumnName = QString("Time");
    secondColumnName = QString("Parameter");
    char firstCol[128], secondCol[128];
    if(2==fscanf(timeHistoryFilePointer,"%s%s",firstCol,secondCol))
    {
        firstColumnName = QString::fromLatin1(firstCol);
        secondColumnName = QString::fromLatin1(secondCol);
    }

    double time,value;
    for(;feof(timeHistoryFilePointer)==0;)
    {
        if(2==fscanf(timeHistoryFilePointer,"%lf%lf",&time,&value))
        {
            //cout<<"____("<<time<<", "<<value<<")____"<<endl;
            timeList<<time;
            valueList<<value;
        }
        else
        {
            //cerr<<"____jumping over an invalid line____"<<endl;
        }
    }

    //! -----------------------------------------------
    //! immediately create the plot after reading data
    //! -----------------------------------------------
    if(!timeList.isEmpty())
    {
        this->clearPanel();
        this->createPlot();
        return true;
    }
    return false;
}

//! ---------------------
//! function: createPlot
//! details:
//! ---------------------
void TimeStepBuilder::createPlot()
{
    myPlot->setVisible(true);

    double ymin = 1e80;
    double ymax = -1e80;
    double xmin = 1e80;
    double xmax = -1e80;

    for(int i=0; i<valueList.length(); i++)
    {
        double curTime = timeList.at(i);
        double curValue = valueList.at(i);
        if(curValue<=ymin) ymin = curValue;
        if(curValue>=ymax) ymax = curValue;
        if(curTime<=xmin) xmin = curTime;
        if(curTime>=xmax) xmax = curTime;
    }
    double deltay = fabs(ymax-ymin)*0.025;
    double deltax = fabs(xmax-xmin)*0.025;

    //! scale y +/- 2.5%
    myPlot->yAxis->setRangeLower(ymin-deltay);
    myPlot->yAxis->setRangeUpper(ymax+deltay);

    //! scale x +/- 2.5%
    myPlot->xAxis->setRangeLower(xmin-deltax);
    myPlot->xAxis->setRangeUpper(xmax+deltax);

    //! zoom out a bit
    myPlot->yAxis->scaleRange(1.1, myPlot->yAxis->range().center());
    myPlot->xAxis->scaleRange(1.1, myPlot->xAxis->range().center());

    //! axis font
    myPlot->xAxis->setLabelFont(QFont("Helvetica",14,QFont::Bold));
    myPlot->yAxis->setLabelFont(QFont("Helvetica",14,QFont::Bold));

    myPlot->xAxis->setTickLabelFont(QFont("Helvetica",12,QFont::Bold));
    myPlot->yAxis->setTickLabelFont(QFont("Helvetica",12,QFont::Bold));

    myPlot->addGraph();
    myPlot->xAxis->setLabel(firstColumnName);
    myPlot->yAxis->setLabel(secondColumnName);

    myPlot->graph(0)->setData(timeList,valueList);
    myPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 10));

    //myPlot->rescaleAxes();
    myPlot->replot();
}

//! --------------
//! function: pan
//! details:
//! --------------
void TimeStepBuilder::pan()
{
    actionZoom->setChecked(false);
    myPlot->pan();
}

//! ---------------
//! function: zoom
//! details:
//! ---------------
void TimeStepBuilder::zoom()
{
    actionPan->setChecked(false);
    myPlot->zoom();
}

//! --------------
//! function: fit
//! details:
//! --------------
void TimeStepBuilder::fit()
{
    actionPan->setChecked(false);
    actionZoom->setChecked(false);
    myPlot->fit();
}

//! -----------------
//! function: accept
//! details:
//! -----------------
void TimeStepBuilder::acceptTimeDivision()
{
    if(!myPlot->getTimeStepMarkerList().isEmpty())
    {
        //! ---------------------------
        //! write settings into a file
        //! ---------------------------
        QString timeStepFile = myWorkingDir+"/timeStepFile.txt";
        //cout<<"____time step file: "<<timeStepFile.toStdString()<<"____"<<endl;
        QFile file(timeStepFile);
        if(file.exists())
        {
            int res = QMessageBox::warning(this,"Time step builder","Overwrite previous settings?", QMessageBox::Ok, QMessageBox::Cancel);
            //cout<<"____res: "<<res<<"____"<<endl;
            switch(res)
            {
            case QMessageBox::Ok:
            {
                //! write time step file
                this->writeTimeStepFile(timeStepFile);
            }
                break;

            case QMessageBox::Cancel:
                return;
                break;
            }
        }
        else
        {
            //! write time step file
            this->writeTimeStepFile(timeStepFile);
        }
        this->close();
    }
    else
    {
        QMessageBox::warning(this,"Time step builder","No time step has been configured yet", QMessageBox::Ok);
    }
}

//! ----------------------------
//! function: writeTimeStepFile
//! details:
//! ----------------------------
void TimeStepBuilder::writeTimeStepFile(const QString &fileName)
{
    FILE *f = fopen(fileName.toStdString().c_str(),"w");
    if(f==NULL)
    {
        QMessageBox::warning(this,"Time step builder","Cannot write time step file", QMessageBox::Ok);
        return;
    }

    const QList<TimeStepMarker> &timeStepMarkerList = myPlot->getTimeStepMarkerList();
    int NbSteps = timeStepMarkerList.length();
    double prevTime = 0;
    for(int i=0; i<NbSteps; i++)
    {
        const TimeStepMarker &curTimeStep = timeStepMarkerList.at(i);
        int timeStepNr = i+1;
        double curTime = curTimeStep.time;
        TimeStepType type =curTimeStep.tsType;
        fprintf(f,"%d\t%d\t%.4f\t%.4f\n",timeStepNr,type, prevTime,curTime);
        prevTime = curTime;
    }
    fclose(f);
}

//! --------------------------
//! function: showContextMenu
//! details:
//! --------------------------
void TimeStepBuilder::showContextMenu(QPoint pos)
{
    QPoint P = this->mapToGlobal(pos);
    myContextMenu->exec(P);
}

//! --------------------------
//! function: updateStatusBar
//! details:
//! --------------------------
void TimeStepBuilder::updateStatusBar(double time)
{
    myStatusBar->showMessage(QString("Time: %1").arg(time),1000);
}

//! --------------------
//! function: saveImage
//! details:
//! --------------------
void TimeStepBuilder::saveImage()
{
    myPlot->saveImage();
}

//! -------------------
//! function: fitRange
//! details:
//! -------------------
void TimeStepBuilder::fitRange()
{
    myPlot->fitRange();
}

//! ---------------------
//! function: clearPanel
//! details:
//! ---------------------
void TimeStepBuilder::clearPanel()
{
    myPlot->clearPanel();
}

//! ------------------------
//! function: createActions
//! details:
//! ------------------------
void TimeStepBuilder::createActions()
{
    cout<<"TimeStepBuilder::createActions()->____function called____"<<endl;
    actionOpen = new QAction("Open",this);
    actionOpen->setIcon(QIcon(":/icons/icon_open.png"));
    actionOpen->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_O));
    actionOpen->setToolTip("Open time history file");
    connect(actionOpen,SIGNAL(triggered(bool)),this,SLOT(openFile()));

    actionClear = new QAction("Clear",this);
    actionClear->setIcon(QIcon(":/icons/icon_clear.png"));
    actionClear->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_X));
    actionClear->setToolTip("Clear all settings");
    connect(actionClear,SIGNAL(triggered(bool)),this,SLOT(clearPanel()));

    actionAccept = new QAction("Accept",this);
    actionAccept->setIcon(QIcon(":/icons/icon_accept.png"));
    actionAccept->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_A));
    actionAccept->setToolTip("Accept current settings and exit");
    connect(actionAccept,SIGNAL(triggered(bool)),this,SLOT(acceptTimeDivision()));

    actionExit = new QAction("Exit",this);
    actionExit->setIcon(QIcon(":/icons/icon_exit.png"));
    actionExit->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_Q));
    connect(actionExit,SIGNAL(triggered(bool)),this,SLOT(close()));

    actionZoom = new QAction("Zoom",this);
    actionZoom->setIcon(QIcon(":/icons/icon_zoom.png"));
    actionZoom->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_Z));
    actionZoom->setCheckable(true);
    actionZoom->setToolTip("Window zoom");
    connect(actionZoom,SIGNAL(triggered(bool)),this,SLOT(zoom()));

    actionPan = new QAction("Pan",this);
    actionPan->setIcon(QIcon(":/icons/icon_pan.png"));
    actionPan->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_P));
    actionPan->setCheckable(true);
    actionPan->setToolTip("Pan");
    connect(actionPan,SIGNAL(triggered(bool)),this,SLOT(pan()));

    actionFit = new QAction("Fit",this);
    actionFit->setIcon(QIcon(":/icons/icon_fit all.png"));
    actionFit->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_O));
    actionFit->setToolTip("Fit all");
    connect(actionFit,SIGNAL(triggered(bool)),this,SLOT(fit()));

    actionSaveImage = new QAction("Save image",this);
    actionSaveImage->setIcon(QIcon(":/icons/icon_save image.png"));
    actionSaveImage->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_I));
    actionSaveImage->setToolTip("Export panel to file");
    connect(actionSaveImage,SIGNAL(triggered(bool)),this,SLOT(saveImage()));

    actionFitRange = new QAction("Fit range",this);
    actionFitRange->setIcon(QIcon(":/icons/icon_expand.png"));
    actionFitRange->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_E));
    actionFitRange->setToolTip("Fit range");
    connect(actionFitRange,SIGNAL(triggered(bool)),this,SLOT(fitRange()));
}
