#ifndef CONVERGENCEDATACHART1_H
#define CONVERGENCEDATACHART1_H

//! ---
//! Qt
//! ---
#include <QHeaderView>
#include <QWidget>
#include <QMouseEvent>

//! ----------------
//! custom includes
//! ----------------
#include "solutioninfo.h"

//! ----
//! C++
//! ----
#include <iostream>
using namespace  std;

class QAction;
class QCustomPlot;

class ConvergenceDataChart1: public QWidget
{
    Q_OBJECT

private:

    enum scaleType
    {
        linear,
        log
    };

    //! -------------
    //! custom plots
    //! -------------
    QCustomPlot *myChartView;
    QCustomPlot *myChartView1;

    QList<solutionInfo> myData;
    QFont myLabelsFont;

    double myXmax;
    double myYmin;
    double myYmax;

    //! --------
    //! Actions
    //! --------
    QAction* myActionHideLegend;
    QAction* myActionShowLegend;
    QAction* myActionSaveImage;
    QAction *myActionToggleYScale;

    bool myLegendIsVisible;

    scaleType myScaleType;
    scaleType myScaleType1;

private:

    void connectMarkers();
    void hideAxes();
    void showAxes();

private slots:

    void initMinMax();
    void showLegend();
    void hideLegend();
    void handleMarkerClicked();
    void saveImage();

    void changeYScale();

public:

    ConvergenceDataChart1(QWidget *parent=0);
    virtual ~ConvergenceDataChart1();

public slots:

    //! update the data from QList<solutionInfo>
    void plotConvergenceData(const QList<solutionInfo> &solutionInfoList);

    //! display the data
    void showData();

    //! clear the viewer (two graphs)
    void clearViewer();

private slots:

    //! show context menu
    void showContextMenu(const QPoint &p);
    void showContextMenu1(const QPoint &p);

protected:

    //! -------------
    //! mouse events
    //! -------------
    //void mousePressEvent(QMouseEvent* event);
};

#endif // CONVERGENCEDATACHART1_H
