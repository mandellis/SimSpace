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
#include <QWidget>
#include <QList>
#include <QVector>
#include <QPoint>
#include <QMenu>
#include <QColor>
#include <ext/QCustomPlot/qcp/qcustomplot.h>
#include <QMouseEvent>

//! ----
//! C++
//! ----
#include <vector>
#include <iostream>
using namespace std;

class QAction;

class TabularDataViewerClass: public QCustomPlot
{
    Q_OBJECT

private:

    CustomTableModel *myData;
    QList<int> myColumns;

    int myX,myY;

private:

    QCPLegend *myLegend;
    bool myLegendIsVisible;

private:

    void connectMarkers(){;}

    //! --------------------
    //! context menu action
    //! --------------------
    QAction *myActionHideLegend;
    QAction *myActionShowLegend;
    QAction *myActionSaveImage;
    QVector<QColor> myColors;

private slots:

    void configureAxis(){;}
    void handleMarkerClicked(){;}
    void ShowContextMenu(QPoint pos);
    void showLegend();
    void hideLegend();
    void saveImage();

public:

    //! ------------
    //! constructor
    //! ------------
    TabularDataViewerClass(QWidget *parent=0);

    //! -----------
    //! destructor
    //! -----------
    virtual ~TabularDataViewerClass()
    {
        cout<<"TabularDataViewerClass::~TabularDataViewerClass()->____DESTRUCTOR CALLED____"<<endl;
    }

    //! ------------------------------
    //! function: getTabularDataModel
    //! details:
    //! ------------------------------
    CustomTableModel* getTabularDataModel() { return myData; }

    //! -------------------
    //! function: setData
    //! details:
    //! -------------------
    void setData(CustomTableModel *tabularData){ myData = tabularData; }

public slots:

    //! ----------
    //! show data
    //! ----------
    void showData(CustomTableModel *tabData, const QList<int> &columnsToShow);

    //! -----------
    //! clearPanel
    //! -----------
    void clearPanel();

    //! ------------
    //! updateViewer
    //! ------------
    void updateViewer();

private:

    //! --------
    //! plotXYs
    //! --------
    void plotXYs();

    void showAxes();
    void hideAxes();

protected:

    //! -------------
    //! mouse events
    //! -------------
    void mouseMoveEvent(QMouseEvent* event);
    void mousePressEvent(QMouseEvent *event);

};

#endif // TabularDataViewerClass_H
