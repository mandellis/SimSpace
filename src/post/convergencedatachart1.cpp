//! ----------------
//! custom includes
//! ----------------
#include "convergencedatachart1.h"
#include "tools.h"
#include "solutioninfo.h"

//! ---
//! Qt
//! ---
#include <QVBoxLayout>
#include <QMenu>
#include <QFileDialog>
#include <QAction>

//! ------------
//! QCustomPlot
//! ------------
#include <QCustomPlot/qcp/qcustomplot.h>

//! --------------------------------
//! function: ConvergenceDataChart1
//! details:
//! --------------------------------
ConvergenceDataChart1::ConvergenceDataChart1(QWidget *parent):QWidget(parent)
{
    cout<<"ConvergenceDataChart1::ConvergenceDataChart1()->____contructor called____"<<endl;

    this->setMouseTracking(true);

    //! --------
    //! margins
    //! --------
    this->setContentsMargins(0, 0, 0, 0);

    //! -----------------
    //! scale type chart
    //! -----------------
    myScaleType = scaleType::linear;
    myScaleType1 = scaleType::linear;

    myLabelsFont = QFont("Helvetica",10);

    //! --------
    //! actions
    //! --------
    myActionHideLegend = new QAction("Hide legend",this);
    myActionHideLegend->setIcon(QIcon(":/icons/icon_lamp OFF.png"));
    connect(myActionHideLegend,SIGNAL(triggered(bool)),this,SLOT(hideLegend()));

    myActionShowLegend = new QAction("Show legend",this);
    myActionShowLegend->setIcon(QIcon(":/icons/icon_lamp ON.png"));
    connect(myActionShowLegend,SIGNAL(triggered(bool)),this,SLOT(showLegend()));

    myActionSaveImage = new QAction("Save image",this);
    myActionSaveImage->setIcon(QIcon(":/icons/icon_image to file.png"));
    connect(myActionSaveImage,SIGNAL(triggered(bool)),this,SLOT(saveImage()));

    myActionToggleYScale = new QAction("Toggle linear/log scale",this);
    myActionToggleYScale->setIcon(QIcon(":/icons/icon_graph scale.png"));
    connect(myActionToggleYScale,SIGNAL(triggered(bool)),this,SLOT(changeYScale()));

    //! -----------------------
    //! pen for vertical lines
    //! -----------------------
    QPen penVertical;
    penVertical.setColor(Qt::green);
    penVertical.setStyle(Qt::DashDotLine);
    penVertical.setWidth(1);

    //! ----------------------------------
    //! pen for the first and second plot
    //! ----------------------------------
    QPen pen1,pen2;
    pen1.setWidth(2);
    pen1.setColor(Qt::red);
    pen2.setWidth(2);
    pen2.setColor(Qt::blue);

    //! -----------
    //! containers
    //! -----------
    myChartView = new QCustomPlot(this);
    myChartView->addGraph();
    myChartView->addGraph();
    myChartView->setObjectName("top chart");
    myChartView1 = new QCustomPlot(this);
    myChartView1->addGraph();
    myChartView1->setObjectName("bottom chart");

    //! ------------------
    //! legend visibility
    //! ------------------
    myLegendIsVisible = false;
    myChartView->legend->setVisible(false);
    myChartView->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignBottom);
    myChartView1->legend->setVisible(false);
    myChartView1->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignBottom);
    //myChartView1->graph()->setLineStyle(QCPGraph::lineStyle());
    //myChartView1->graph(1)->setLineStyle(QCPGraph::lineStyle());

    //! --------------
    //! customization
    //! --------------
    myChartView->axisRect()->setupFullAxesBox();
    myChartView->rescaleAxes(true);
    myChartView->graph(0)->setPen(pen1);
    myChartView->graph(1)->setPen(pen2);
    myChartView->graph(0)->setLineStyle(QCPGraph::lsLine);
    myChartView->graph(1)->setLineStyle(QCPGraph::lsLine);
    myChartView->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    myChartView->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    myChartView->xAxis->setLabelFont(myLabelsFont);
    myChartView->yAxis->setLabelFont(myLabelsFont);

    if(myScaleType == scaleType::linear) myChartView->yAxis->setScaleType(QCPAxis::stLinear);
        else myChartView->yAxis->setScaleType(QCPAxis::stLogarithmic);

    myChartView1->axisRect()->setupFullAxesBox();
    myChartView1->rescaleAxes(true);
    myChartView1->graph()->setPen(pen2);
    myChartView1->xAxis->setLabelFont(myLabelsFont);
    myChartView1->yAxis->setLabelFont(myLabelsFont);

    //! ------------------------
    //! font of the thick label
    //! ------------------------
    myChartView->xAxis->setTickLabelFont(myLabelsFont);
    myChartView->yAxis->setTickLabelFont(myLabelsFont);
    myChartView1->xAxis->setTickLabelFont(myLabelsFont);
    myChartView1->yAxis->setTickLabelFont(myLabelsFont);
    myChartView1->graph(0)->setLineStyle(QCPGraph::lsLine);
    myChartView1->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));

    //! --------------
    //! hide the axes
    //! --------------
    this->hideAxes();

    //! ----------
    //! axes font
    //! ----------
    myLabelsFont.setPixelSize(14);
    myLabelsFont.setBold(true);

    //! -------
    //! layout
    //! -------
    QVBoxLayout *vb = new QVBoxLayout(this);
    vb->setContentsMargins(0,0,0,0);
    vb->addWidget(myChartView);
    vb->addWidget(myChartView1);

    //! ------------------
    //! mouse tracking on
    //! ------------------
    myChartView->setMouseTracking(true);
    myChartView1->setMouseTracking(true);

    //! -------------
    //! context menu
    //! -------------
    myChartView->setContextMenuPolicy(Qt::CustomContextMenu);
    myChartView->setInteraction(QCP::iSelectItems);
    connect(myChartView,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(showContextMenu(const QPoint&)));

    myChartView1->setContextMenuPolicy(Qt::CustomContextMenu);
    myChartView1->setInteraction(QCP::iSelectItems);
    connect(myChartView1,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(showContextMenu1(const QPoint&)));
}

//! ---------------------
//! function: destructor
//! details:
//! ---------------------
ConvergenceDataChart1::~ConvergenceDataChart1()
{
    cout<<"ConvergenceDataChart1::~ConvergenceDataChart1()->____DESTRUCTOR CALLED____"<<endl;
}

//! -------------------
//! function: showData
//! details:
//! -------------------
void ConvergenceDataChart1::showData()
{
    cout<<"ConvergenceDataChart1::showData()->____function called____"<<endl;
}

//! ------------------------------
//! function: handleMarkerClicked
//! details:
//! ------------------------------
void ConvergenceDataChart1::handleMarkerClicked()
{
    ;
}

//! -------------------------
//! function: connectMarkers
//! details:
//! -------------------------
void ConvergenceDataChart1::connectMarkers()
{
    ;
}

//! ----------------------
//! function: clearViewer
//! details:
//! ----------------------
void ConvergenceDataChart1::clearViewer()
{
    for(int i=0; i<myChartView->graphCount(); i++)
    {
        myChartView->graph(i)->data()->clear();
        myChartView->replot();
    }
    for(int i=0; i<myChartView1->graphCount(); i++)
    {
        myChartView1->graph(i)->data()->clear();
        myChartView->replot();
    }

    //! ------------------------------------------------------
    //! hide the legend without changing the status variables
    //! ------------------------------------------------------
    //this->legend->setVisible(false);
    //this->replot();
}

//! --------------------------
//! function: showContextMenu
//! details:
//! --------------------------
void ConvergenceDataChart1::showContextMenu(const QPoint& p)
{
    cout<<"ConvergenceDataChart1::showContextMenu()->____function called____"<<endl;

    if(myChartView->graph()->data()->isEmpty()) return;

    QPoint globalPos = myChartView->mapToGlobal(p);
    QMenu *contextMenu = new QMenu(this);
    contextMenu->setTitle("Convergence info");

    if(myLegendIsVisible==true) contextMenu->addAction(myActionHideLegend);
    else contextMenu->addAction(myActionShowLegend);
    contextMenu->addSeparator();
    contextMenu->addAction(myActionSaveImage);
    contextMenu->addSeparator();
    contextMenu->addAction(myActionToggleYScale);
    contextMenu->exec(globalPos);
}

//! --------------------------
//! function: showContextMenu
//! details:
//! --------------------------
void ConvergenceDataChart1::showContextMenu1(const QPoint& p)
{
    cout<<"ConvergenceDataChart1::showContextMenu1()->____function called____"<<endl;

    if(myChartView1->graph()->data()->isEmpty()) return;

    QPoint globalPos = myChartView1->mapToGlobal(p);
    QMenu *contextMenu = new QMenu(this);
    contextMenu->setTitle("Timing info");

    if(myLegendIsVisible==true) contextMenu->addAction(myActionHideLegend);
    else contextMenu->addAction(myActionShowLegend);
    contextMenu->addSeparator();

    contextMenu->addAction(myActionSaveImage);
    contextMenu->exec(globalPos);
}

//! ------------------------------------
//! function: plotConvergenceData
//! details:  read and display the data
//! ------------------------------------
void ConvergenceDataChart1::plotConvergenceData(const QList<solutionInfo> &solutionInfoList)
{
    cout<<"ConvergenceDataChart1::plotConvergenceData()->____function called: updating convergence viewer____"<<endl;

    //! --------------
    //! show the axes
    //! --------------
    if(!myChartView->xAxis->visible()) this->showAxes();

    //! ---------------------------------------
    //! 1) global iteration number
    //! 2) step number
    //! 3) increment
    //! 4) attempt
    //! 5) actual step time
    //! 6) actual total time
    //! 7) average force
    //! 8) largest residual force
    //! 9) largest increment of displacement
    //! 10) largest correction to displacement
    //! 11) substep converged
    //! ---------------------------------------
    if(solutionInfoList.isEmpty()) return;
    int NbData = solutionInfoList.length();

    //! ---------------
    //! top panel data
    //! ---------------
    QVector<double> X_top;
    QVector<double> Y1_top, Y2_top;

    //! ------------------
    //! bottom panel data
    //! ------------------
    QVector<double> X_bottom;
    QVector<double> Y1_bottom;

    //! ---------------------------------
    //! display flux (force) convergence
    //! ---------------------------------
    for(int i=0; i<NbData; i++)
    {
        X_top.push_back(solutionInfoList.at(i).globalIterationNb);
        Y1_top.push_back(solutionInfoList.at(i).largestResidual);
        Y2_top.push_back(solutionInfoList.at(i).average);
    }

    //! -------
    //! labels
    //! -------
    QString xLabel,yLabel;
    xLabel = "Global iteration number";
    yLabel = "Largest residual force/Average force";

    myChartView->xAxis->setLabel(xLabel);
    myChartView->yAxis->setLabel(yLabel);

    myChartView->graph(0)->setData(X_top,Y1_top);
    myChartView->graph(1)->setData(X_top,Y2_top);
    myChartView->rescaleAxes(true);

    //! --------
    //! scale Y
    //! --------
    if(myScaleType==scaleType::linear)
        myChartView->yAxis->setScaleType(QCPAxis::stLinear);
    else myChartView->yAxis->setScaleType(QCPAxis::stLogarithmic);

    myChartView->replot();

    //! ---------------------------------------
    //! display total time vs iteration number
    //! ---------------------------------------
    for(int i=0; i<NbData; i++)
    {
        X_bottom.push_back(solutionInfoList.at(i).globalIterationNb);
        Y1_bottom.push_back(solutionInfoList.at(i).time);
    }

    //! -------
    //! labels
    //! -------
    QString xLabel1,yLabel1;
    xLabel1 = "Global iteration number";
    yLabel1 = "Total time";

    myChartView1->graph(0)->setData(X_bottom,Y1_bottom);
    myChartView1->rescaleAxes(true);
    myChartView1->xAxis->setLabel(xLabel1);
    myChartView1->yAxis->setLabel(yLabel1);
    myChartView1->replot();
}

//! ---------------------
//! function: initMinMax
//! details:
//! ---------------------
void ConvergenceDataChart1::initMinMax()
{
    //! x axis max value - initial value (@ start)
    myXmax = 5;

    //! y range - initial (@ start)
    myYmin = -1.0;
    myYmax = 1.0;
}

//! ---------------------
//! function: hideLegend
//! details:
//! ---------------------
void ConvergenceDataChart1::hideLegend()
{
    ;
}

//! ---------------------
//! function: showLegend
//! details:
//! ---------------------
void ConvergenceDataChart1::showLegend()
{

}

//! --------------------
//! function: saveImage
//! details:
//! --------------------
void ConvergenceDataChart1::saveImage()
{
    QString filters("PNG files (*.png);;JPG files (*.jpg);; BMP files (*.bmp);; PDF files (*.pdf)");
    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(0,"Save image as ",QDir::currentPath(),filters,&selectedFilter,0);
    if(!fileName.isEmpty())
    {
        if(selectedFilter=="JPG files (*.jpg)") myChartView->saveJpg(fileName);
        if(selectedFilter=="BMP files (*.bmp)") myChartView->saveBmp(fileName);
        if(selectedFilter=="PNG files (*.png)") myChartView->savePng(fileName);
        if(selectedFilter=="PDF files (*.pdf)") myChartView->savePdf(fileName);
    }
}

//! -------------------
//! function: hideAxes
//! details:
//! -------------------
void ConvergenceDataChart1::hideAxes()
{
    myChartView->xAxis->setVisible(false);
    myChartView->xAxis2->setVisible(false);
    myChartView->yAxis->setVisible(false);
    myChartView->yAxis2->setVisible(false);
    myChartView1->xAxis->setVisible(false);
    myChartView1->xAxis2->setVisible(false);
    myChartView1->yAxis->setVisible(false);
    myChartView1->yAxis2->setVisible(false);
}

//! -------------------
//! function: showAxes
//! details:
//! -------------------
void ConvergenceDataChart1::showAxes()
{
    myChartView->xAxis->setVisible(true);
    myChartView->xAxis2->setVisible(true);
    myChartView->yAxis->setVisible(true);
    myChartView->yAxis2->setVisible(true);
    myChartView1->xAxis->setVisible(true);
    myChartView1->xAxis2->setVisible(true);
    myChartView1->yAxis->setVisible(true);
    myChartView1->yAxis2->setVisible(true);
}


//! -----------------------
//! function: changeYScale
//! details:
//! -----------------------
void ConvergenceDataChart1::changeYScale()
{
    if(myScaleType == scaleType::linear)
    {
        myScaleType = scaleType::log;
        myChartView->yAxis->setScaleType(QCPAxis::stLogarithmic);
    }
    else
    {
        myScaleType = scaleType::log;
        myChartView->yAxis->setScaleType(QCPAxis::stLinear);
    }
    myChartView->replot();
}
