//! ----------------
//! custom includes
//! ----------------
#include <tabularData/tabulardataviewerclass.h>

//! ---
//! Qt
//! ---
#include <QFileDialog>
#include <QDir>
#include <QRect>
#include <QAction>
#include <QVariant>

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TabularDataViewerClass::TabularDataViewerClass(QWidget *parent):QCustomPlot(parent)
{
    //! ------------------
    //! mouse tracking on
    //! ------------------
    this->setMouseTracking(true);

    //! --------------
    //! hide the axes
    //! --------------
    this->hideAxes();

    //! -------
    //! legend
    //! -------
    myLegendIsVisible = true;
    this->legend->setVisible(false);
    this->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignBottom);

    //! -------------
    //! context menu
    //! -------------
    this->setContextMenuPolicy(Qt::CustomContextMenu);
    this->setInteraction(QCP::iSelectItems);
    connect(this,SIGNAL(customContextMenuRequested(const QPoint&)),this,SLOT(ShowContextMenu(const QPoint&)));

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

    //! -----------
    //! set colors
    //! -----------
    myColors.push_back(QColor(255,0,0,255));
    myColors.push_back(QColor(255,0,255,255));
    myColors.push_back(QColor(0,0,255,255));
    myColors.push_back(QColor(0,255,255,255));
    myColors.push_back(QColor(0,255,0,255));
    myColors.push_back(QColor(0,0,255,255));
    myColors.push_back(QColor(255,255,0,255));
}

//! ---------------------
//! function: showLegend
//! details:
//! ---------------------
void TabularDataViewerClass::showLegend()
{
    myLegendIsVisible = true;
    this->legend->setVisible(true);
    this->replot();
}

//! ---------------------
//! function: hideLegend
//! details:
//! ---------------------
void TabularDataViewerClass::hideLegend()
{
    myLegendIsVisible = false;
    this->legend->setVisible(false);
    this->replot();
}

//! --------------------------
//! function: ShowContextMenu
//! details:
//! --------------------------
void TabularDataViewerClass::ShowContextMenu(QPoint pos)
{
    QPoint globalPos = this->mapToGlobal(pos);

    QMenu *contextMenu = new QMenu(this);
    if(myLegendIsVisible==true) contextMenu->addAction(myActionHideLegend);
    else contextMenu->addAction(myActionShowLegend);

    contextMenu->addSeparator();
    contextMenu->addAction(myActionSaveImage);

    contextMenu->exec(globalPos);
}

//! -------------------
//! function: showAxes
//! details:
//! -------------------
void TabularDataViewerClass::showAxes()
{
    this->xAxis->setVisible(true);
    this->xAxis2->setVisible(true);
    this->yAxis->setVisible(true);
    this->yAxis2->setVisible(true);
}

//! -------------------
//! function: hideAxes
//! details:
//! -------------------
void TabularDataViewerClass::hideAxes()
{
    this->xAxis->setVisible(false);
    this->xAxis2->setVisible(false);
    this->yAxis->setVisible(false);
    this->yAxis2->setVisible(false);
}

//! -------------------------------------------------------------
//! function: plotXYs
//! details:  the first colum data is always put onto the X axis
//! -------------------------------------------------------------
void TabularDataViewerClass::plotXYs()
{
    //! -----------------
    //! clear all graphs
    //! -----------------
    this->clearGraphs();

    //! ------------------
    //! activate the axes
    //! ------------------
    this->axisRect()->setupFullAxesBox();

    //! -------------------
    //! display the legend
    //! -------------------
    if(myLegendIsVisible) this->legend->setVisible(true);

    //! -----------------
    //! get the X series
    //! -----------------
    int firstColumn = myColumns.at(0);

    const load &aLoad = myData->getColumn(firstColumn);
    QVector<QVariant> XData = aLoad.values();

    this->xAxis->setLabel("Time [s]");

    double xmin, xmax;
    xmin = 1e80; xmax = -1e80;
    QVector<double> XSeries;
    for(int n = 0; n<XData.size(); n++)
    {
        double xval = XData.at(n).toDouble();
        if(xval<xmin) xmin=xval;
        if(xval>xmax) xmax=xval;
        XSeries.push_back(xval);
    }
    this->xAxis->setRange(xmin,xmax);

    //! ----------------------
    //! get the Y data series
    //! ----------------------
    double ymin, ymax;
    ymin = 1e80; ymax = -1e80;

    int NbYSeries = myColumns.length()-1;
    for(int n=1; n<=NbYSeries; n++)
    {
        QCPGraph *aGraph = this->addGraph();
        const load &aLoadY = myData->getColumn(myColumns.at(n));
        QVector<QVariant> YData = aLoadY.values();

        QVector<double> YSeries;
        for(int j=0; j<YData.size(); j++)
        {
            double yval = YData.at(j).toDouble();
            if(yval<ymin) ymin = yval;
            if(yval>ymax) ymax = yval;
            YSeries.push_back(yval);
        }
        aGraph->setData(XSeries,YSeries);
        aGraph->setName(myData->getHeaderString(myColumns.at(n)));
        aGraph->setPen(QPen(myColors.at((n-1)%myColors.length())));
    }

    this->yAxis->setLabel("Y Values");
    double dy = ymax-ymin;
    this->yAxis->setRange(ymin-dy*0.025,ymax+dy*0.025);
    this->replot();
}

//! --------------------
//! function: show data
//! details:
//! --------------------
void TabularDataViewerClass::showData(CustomTableModel *tabData, const QList<int> &columnsToShow)
{
    //cout<<"TabularDataViewerClass::showData1()->____function called____"<<endl;

    myData = tabData;
    myColumns = columnsToShow;

    //! -----------------------------------------------------------------------------------
    //! actually the effect of this connection consists in updating of the axes and labels
    //! ------------------------------------------------------------------------------------
    disconnect(myData,SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),this,SLOT(updateViewer()));
    connect(myData,SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),this,SLOT(updateViewer()));

    this->plotXYs();

    //cout<<"TabularDataViewerClass::showData1()->____exiting function____"<<endl;
}

//! -----------------------
//! function: updateViewer
//! detail:
//! -----------------------
void TabularDataViewerClass::updateViewer()
{
    //cout<<"TabularDataViewerClass::updateViewer()->____function called____"<<endl;
    this->plotXYs();
}

//! ---------------------
//! function: clearPanel
//! details:
//! ---------------------
void TabularDataViewerClass::clearPanel()
{
    //cout<<"TabularDataViewerClass::clearPanel()->____function called____"<<endl;
    for(int i=0; i<this->graphCount(); i++) this->graph(i)->data()->clear();
    this->hideAxes();
    //! ------------------------------------------------------
    //! hide the legend without changing the status variables
    //! ------------------------------------------------------
    this->legend->setVisible(false);
    this->replot();
}

//! -------------------------
//! function: mouseMoveEvent
//! details:
//! -------------------------
void TabularDataViewerClass::mouseMoveEvent(QMouseEvent *event)
{
    //cout<<"TabularDataViewerClass::mousePressEvent()->____MOUSE MOVE____"<<endl;
    myX = event->pos().x();
    myY = event->pos().y();
}

//! --------------------------
//! function: mousePressEvent
//! details:
//! --------------------------
void TabularDataViewerClass::mousePressEvent(QMouseEvent *event)
{
    cout<<"TabularDataViewerClass::mousePressEvent()->____MOUSE PRESS____"<<endl;
    myX = event->pos().x();
    myY = event->pos().y();
}

//! --------------------
//! function: saveImage
//! details:
//! --------------------
void TabularDataViewerClass::saveImage()
{
    QString filters("PNG files (*.png);;JPG files (*.jpg);; BMP files (*.bmp);; PDF files (*.pdf)");
    QString selectedFilter;
    QString fileName = QFileDialog::getSaveFileName(0,"Save image as ",QDir::currentPath(),filters,&selectedFilter,0);

    //cout<<"____"<<fileName.toStdString()<<"____"<<endl;

    if(!fileName.isEmpty())
    {
        cout<<"____"<<selectedFilter.toStdString()<<"____"<<endl;
        if(selectedFilter=="JPG files (*.jpg)") this->saveJpg(fileName);
        if(selectedFilter=="BMP files (*.bmp)") this->saveBmp(fileName);
        if(selectedFilter=="PNG files (*.png)") this->savePng(fileName);
        if(selectedFilter=="PDF files (*.pdf)") this->savePdf(fileName);
    }
}
