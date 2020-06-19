//! ----------------
//! custom includes
//! ----------------
#include "tabulardataviewerclass.h"
#include "mydefines.h"

//! ---
//! Qt
//! ---
#include <QGraphicsLayout>
#include <QLegendMarker>
#include <QMessageBox>
#include <QValueAxis>

//! ----
//! C++
//! ----
#include <iostream>
using namespace std;

//! ----------------------
//! function: constructor
//! details:
//! ----------------------
TabularDataViewerClass::TabularDataViewerClass(QWidget *parent):QChartView(parent)
{
    cout<<"TabularDataViewerClass::TabularDataViewerClass()->____contructor called____"<<endl;

    myChart = new QChart;
    myChart->setBackgroundRoundness(0);
    myChart->layout()->setContentsMargins(0,0,0,0);
    this->setChart(myChart);
    this->setRenderHint(QPainter::Antialiasing);

    //! set object name
    this->setObjectName("tabularDataViewer");
}

//! -----------------------
//! function: updateViewer
//! details:
//! -----------------------
void TabularDataViewerClass::updateViewer()
{
    cout<<"TabularDataViewerClass::updateViewer()->____function called____"<<endl;
    if(myData->columnCount()>=NUMBER_OF_COLUMNS_BEFORE_BC_DATA)
    {
        this->showData(myData,myYseries);
    }
}

//! -------------------
//! function: showData
//! details:
//! -------------------
void TabularDataViewerClass::showData(CustomTableModel *tabData, const QList<int> &columnsToShow)
{
    cout<<"TabularDataViewerClass::showData()->____function called____"<<endl;

    myYseries = columnsToShow;
    myData = tabData;
    myChart->removeAllSeries();

    //! -----------------------------------------------------------------------------------
    //! actually the effect of this connection consists in updating of the axes and labels
    //! ------------------------------------------------------------------------------------
    disconnect(myData,SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),this,SLOT(updateViewer()));
    connect(myData,SIGNAL(dataChanged(QModelIndex,QModelIndex,QVector<int>)),this,SLOT(updateViewer()));

    myChart->legend()->setAlignment(Qt::AlignBottom);
    for(int k=1; k<columnsToShow.length(); k++)
    {
        cout<<"TabularDataViewerClass::showData()->____"<<columnsToShow.at(k)<<"____"<<endl;

        myMapper = new QVXYModelMapper(this);
        myMapper->setModel(myData);
        mySeries = new QLineSeries;

        QVariant headerString = myData->headerData(columnsToShow.at(k),Qt::Horizontal,Qt::DisplayRole);
        mySeries->setName(headerString.toString());

        //! ------------------------------------------------------
        //! first column: times => they go on the horizontal axis
        //! ------------------------------------------------------
        myMapper->setXColumn(columnsToShow.at(0));

        //! ------------------------------------------------------
        //! second column: values => they go on the vertical axis
        //! ------------------------------------------------------
        myMapper->setYColumn(columnsToShow.at(k));
        myMapper->setSeries(mySeries);

        myChart->addSeries(mySeries);
        this->configureAxis();
        this->connectMarkers();
    }
    cout<<"TabularDataViewerClass::showData()->____exiting function____"<<endl;
}

//! --------------------------
//! function: configureAxis()
//! details:
//! --------------------------
void TabularDataViewerClass::configureAxis()
{
    myChart->createDefaultAxes();
    static_cast<QValueAxis>(myChart->axisY()).setLabelFormat("%g");
}

//! ------------------------------
//! function: handleMarkerClicked
//! details:
//! ------------------------------
void TabularDataViewerClass::handleMarkerClicked()
{
    //cout<<"TableWidget::handleMarkerClicked()->____handle marker clicked____"<<endl;
    QLegendMarker* marker = qobject_cast<QLegendMarker*> (sender());

    //! Toggle visibility of series
    marker->series()->setVisible(!marker->series()->isVisible());

    //! Turn legend marker back to visible, since hiding series also hides the marker
    marker->setVisible(true);

    //! Dim the marker, if series is not visible
    qreal alpha = 1.0;

    if (!marker->series()->isVisible()) alpha = 0.5;

    QColor color;
    QBrush brush = marker->labelBrush();
    color = brush.color();
    color.setAlphaF(alpha);
    brush.setColor(color);
    marker->setLabelBrush(brush);

    brush = marker->brush();
    color = brush.color();
    color.setAlphaF(alpha);
    brush.setColor(color);
    marker->setBrush(brush);

    QPen pen = marker->pen();
    color = pen.color();
    color.setAlphaF(alpha);
    pen.setColor(color);
    marker->setPen(pen);
}

//! -------------------------
//! function: connectMarkers
//! details:
//! -------------------------
void TabularDataViewerClass::connectMarkers()
{
    foreach (QLegendMarker* marker, myChart->legend()->markers())
    {
        disconnect(marker, SIGNAL(clicked()), this, SLOT(handleMarkerClicked()));
        connect(marker, SIGNAL(clicked()), this, SLOT(handleMarkerClicked()));
    }
}

//! ----------------------
//! function: clearViewer
//! details:
//! ----------------------
void TabularDataViewerClass::clearGraphViewer()
{
    if(myChart==Q_NULLPTR) return;
    if(myChart->series().length()!=0) myChart->removeAllSeries();
    if(myChart->axisY()!=Q_NULLPTR) myChart->removeAxis(myChart->axisY());
    if(myChart->axisX()!=Q_NULLPTR) myChart->removeAxis(myChart->axisX());
}
