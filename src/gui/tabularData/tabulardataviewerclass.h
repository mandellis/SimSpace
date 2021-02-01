#ifndef TABULARDATAVIEWERCLASS_H
#define TABULARDATAVIEWERCLASS_H

//! ----------------
//! custom includes
//! ----------------
#include "customtablemodel.h"
#include "tabulardatacolumns.h"

//! ---
//! Qt
//! ---
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QVXYModelMapper>
#include <QtCharts/QValueAxis>
#include <QHeaderView>
#include <QWidget>

using namespace QtCharts;

class TabularDataViewerClass: public QChartView
{
    Q_OBJECT

private:

    QVXYModelMapper *myMapper;

    QChart *myChart;
    QValueAxis *myAxisX;
    QValueAxis *myAxisY;
    QLineSeries *mySeries;
    CustomTableModel *myData;

    QList<int> myYseries;

private:

    void connectMarkers();

private slots:

    void configureAxis();
    void handleMarkerClicked();

public:

    TabularDataViewerClass(QWidget *parent=0);
    virtual ~TabularDataViewerClass()
    {
        cout<<"TabularDataViewerClass::~TabularDataViewerClass()->____DESTRUCTOR CALLED____"<<endl;
    }

    CustomTableModel* getTabularDataModel() { return myData; }
    void setData(CustomTableModel *tabularData){ myData = tabularData; }

public slots:

    void showData(CustomTableModel *tabData, const QList<int> &columnsToShow);
    void clearGraphViewer();

    void updateViewer();
};

#endif // TABULARDATAVIEWERCLASS_H
