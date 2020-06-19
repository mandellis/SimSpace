//! ----------------
//! custom includes
//! ----------------
#include "qhistogram.h"

//! ------------
//! QCustomPlot
//! ------------
#include <qcustomplot.h>

//! ---
//! Qt
//! ---
#include <QVBoxLayout>
#include <QPoint>
#include <QRect>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
QHistogram::QHistogram(QWidget *parent):QWidget(parent)
{
    cout<<"QHistogram::QHistogram()->____constructor I called____"<<endl;

    //! ----------------------------------
    //! a single bin [0,1] with 0 samples
    //! ----------------------------------
    //std::pair<hbin,int> element;
    //element.first = hbin(0,1);
    //element.second = 0;
    //myhistogramData.insert(element);

    //cout<<"QHistogram::QHistogram()->____make fake data____"<<endl;
    //myhistogramData = fakeData();

    this->createContent();
    this->updatePlot();
}

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
QHistogram::QHistogram(const histogramData &hdata, QWidget *parent):QWidget(parent)
{
    cout<<"QHistogram::QHistogram()->____constructor II called____"<<endl;

    //! ------------------------
    //! init data with external
    //! ------------------------
    std::map<hbin,double> hisData;
    hdata.getData(hisData);
    for(std::map<hbin,double>::iterator it = hisData.begin(); it!=hisData.end(); ++it)
    {
        const std::pair<hbin,int> &element = *it;
        myHistogramData.insert(element);
    }

    this->createContent();
    this->updatePlot();
}

//! ------------------
//! function: setData
//! details:
//! ------------------
void QHistogram::setData(const histogramData &hdata, bool updateViewer)
{
    cout<<"QHistogram::setData()->____function called____"<<endl;

    //! ------------------------
    //! init data with external
    //! ------------------------
    std::map<hbin,double> hisData;
    hdata.getData(hisData);
    for(std::map<hbin,double>::iterator it = hisData.begin(); it!=hisData.end(); ++it)
    {
        const std::pair<hbin,int> &element = *it;
        cout<<"QHistogram::setData()->____("<<element.first.xmin<<", "<<element.first.xmax<<") => "<<element.second<<")____"<<endl;
        myHistogramData.insert(element);
    }

    //! ----------------
    //! update the plot
    //! ----------------
    if(updateViewer) this->updatePlot();
}

//! ------------------------
//! function: createContent
//! details:
//! ------------------------
void QHistogram::createContent()
{
    myPlot = new QCustomPlot(this);

    //! ----------------
    //! vertical layout
    //! ----------------
    QVBoxLayout *mainLayout = new QVBoxLayout();
    mainLayout->setContentsMargins(0,0,0,0);
    mainLayout->setMargin(0);
    mainLayout->setSpacing(0);
    this->setLayout(mainLayout);

    //! ---------------
    //! add the widget
    //! ---------------
    mainLayout->addWidget(myPlot,Qt::AlignTop);
}

//! ---------------------
//! function: updatePlot
//! details:
//! ---------------------
void QHistogram::updatePlot()
{
    cout<<"QHistogram::updatePlot()->____function called____"<<endl;

    //! -------------------------------
    //! clear previous data and replot
    //! -------------------------------
    this->clear();

    //! -------------------------------
    //! show plot if previously hidden
    //! -------------------------------
    QVector<double> ticks;
    QVector<QString> labels;

    //! -----------------
    //! xtics and labels
    //! -----------------
    QVector<double> barData;
    QVector<double> barTics;
    double xMax = -1e80;
    double xMin = 1e80;
    double yMax = -1e80;
    double w;

    //! -------------------------------------
    //! retrieve the data from the container
    //! -------------------------------------
    std::map<hbin,double> hisData;
    myHistogramData.getData(hisData);

    for(std::map<hbin,double>::iterator it = hisData.begin(); it != hisData.end(); it++)
    {
        const std::pair<hbin,double> &element = *it;

        //cout<<"QHistogram::updatePlot()->____("<<element.first.xmin<<", "<<element.first.xmax<<") => "<<element.second<<"____"<<endl;

        ticks<<element.first.xmin;
        ticks<<element.first.xmax;

        labels<<QString("%1").arg(element.first.xmin);
        labels<<QString("%1").arg(element.first.xmax);

        if(element.first.xmax>xMax) xMax = element.first.xmax;
        if(element.first.xmin<xMin) xMin = element.first.xmin;
        if(element.second>yMax) yMax = element.second;

        barTics<<double(element.first.intervalCenter());
        barData<<double(element.second);
        w = element.first.xmax-element.first.xmin;
    }
    //QCPBars *aBar = new QCPBars(myPlot->xAxis, myPlot->yAxis);
    myBars = new QCPBars(myPlot->xAxis, myPlot->yAxis);

    //aBar->setWidthType(QCPBars::wtPlotCoords);
    //aBar->setWidth(w);
    //bool alreadySorted = true;
    //aBar->setData(barTics,barData,alreadySorted);
    //aBar->setPen(QPen(QColor(0, 168, 140).lighter(130)));
    //aBar->setBrush(QColor(0, 168, 140));

    myBars->setWidthType(QCPBars::wtPlotCoords);
    myBars->setWidth(w);
    bool alreadySorted = false;
    myBars->setData(barTics,barData,alreadySorted);
    myBars->setPen(QPen(QColor(0, 168, 140).lighter(130)));
    myBars->setBrush(QColor(0, 168, 140));

    //! -----------------------------------------------
    //! make top right axes clones of bottom left axes
    //! -----------------------------------------------
    myPlot->axisRect()->setupFullAxesBox();

    //! -----------------
    //! configure x axis
    //! -----------------
    QSharedPointer<QCPAxisTickerText> textTicker(new QCPAxisTickerText);
    textTicker->addTicks(ticks,labels);
    myPlot->xAxis->setTicker(textTicker);
    myPlot->xAxis->setSubTicks(false);
    myPlot->xAxis->setTickLength(0, 2);
    myPlot->xAxis->setRange(xMin, xMax);
    myPlot->xAxis->setLabel("Quality parameter");
    myPlot->xAxis->setBasePen(QPen(Qt::black));
    myPlot->xAxis->setTickPen(QPen(Qt::black));
    myPlot->xAxis->grid()->setVisible(true);
    myPlot->xAxis->grid()->setPen(QPen(QColor(130, 130, 130), 0, Qt::DotLine));
    myPlot->xAxis->setTickLabelColor(Qt::black);
    myPlot->xAxis->setLabelColor(Qt::black);
    myPlot->xAxis->setTickLabelRotation(60);
    myPlot->xAxis->setLabelFont(QFont("Helvetica",12,QFont::Bold));

    //! -----------------
    //! configure y axis
    //! -----------------
    myPlot->yAxis->setRange(0, 1.05*yMax);
    myPlot->yAxis->setLabel("Samples");
    myPlot->yAxis->setBasePen(QPen(Qt::black));
    myPlot->yAxis->setTickPen(QPen(Qt::black));
    myPlot->yAxis->setSubTickPen(QPen(Qt::black));
    myPlot->yAxis->grid()->setSubGridVisible(true);
    myPlot->yAxis->setTickLabelColor(Qt::black);
    myPlot->yAxis->setLabelColor(Qt::black);
    myPlot->yAxis->grid()->setPen(QPen(QColor(130, 130, 130), 0, Qt::SolidLine));
    myPlot->yAxis->grid()->setSubGridPen(QPen(QColor(130, 130, 130), 0, Qt::DotLine));
    myPlot->yAxis->setLabelFont(QFont("Helvetica",12,QFont::Bold));

    //! -------
    //! replot
    //! -------
    myPlot->replot();
}

//! -------------------
//! function: fakeData
//! details:
//! -------------------
histogramData QHistogram::fakeData()
{
    histogramData data;
    const double dx = 0.2;
    for(int n=0; n<10; n++)
    {
        double x = n*dx;
        std::pair<hbin,double> element;
        element.first.xmin=x;
        element.first.xmax=x+dx;
        element.second = 10.0*double(n);
        data.insert(element);
    }
    return data;
}

//! ----------------
//! function: clear
//! details:
//! ----------------
void QHistogram::clear()
{
    cout<<"QHistogram::clear()->____function called____"<<endl;
    if(myBars!=Q_NULLPTR) myPlot->removePlottable(myBars);
    //myPlot->graph()->data().clear();
}
