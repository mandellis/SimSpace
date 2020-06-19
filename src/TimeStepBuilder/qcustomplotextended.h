#ifndef QCUSTOMPLOTEXTENDED_H
#define QCUSTOMPLOTEXTENDED_H

//! custom includes
#include "qcustomplot.h"
#include "Costamp/timesteptype.h"

//! Qt
#include <QWidget>
#include <QMouseEvent>
#include <QPoint>
#include <QMenu>
#include <QMap>

class QRubberBand;

//! C++
#include <iostream>
using namespace std;

class QCustomPlotExtended: public QCustomPlot
{
   Q_OBJECT

public:

    QCustomPlotExtended(QWidget *parent=0);

private:

    enum action
    {
        action_none,
        action_drag,
        action_zoom,
        action_fit,
        action_pan
    };

    //! current action
    action myCurAction;

    //! rubber band
    QRubberBand *myRectBand;

    //! current mouse pointer coordinates
    int myX, myY;

    //! start point after mouse lb click
    int myX_start;
    int myY_start;

    //! ranges at the creation of a graph
    double xinitmin,xinitmax,yinitmin,yinitmax;

    //! context menu
    QMenu *myContextMenu;

    //! context menu actions
    QAction *actionSaveImage;

    //! a very custom context menu - Costamp stuff
    QMenu *myVeryCustomMenu;

    //! tracer
    QCPItemTracer *myTracer;

    //! tracing enabled/disabled
    bool isTracingEnabled;

    //! item for selection
    QCPAbstractItem *myCurSelectedItem;

    QList<TimeStepMarker> TimeStepMarkerList;
    QMap<TimeStepMarker,QCPItemLine*> vLineMap;
    QMap<TimeStepMarker,QCPItemRect*> rectMap;

private:

    //! --------------
    //! Costamp stuff
    //! --------------
    void fillVeryCustomMenu(QCPItemTracer *aTimeTracer);
    QBrush getBrush(const TimeStepMarker &aTimeStepMarker);
    void updateBoxes();

private slots:

    void drawRubberBand(int tlx, int tly, int brx, int bry);

protected:

    //! mouse events
    void mouseMoveEvent(QMouseEvent* event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

    void showPointValue(QMouseEvent* event);

protected slots:

    void ShowContextMenu(QPoint pos);

public slots:

    void setTracerVisible(bool isVisible) { myTracer->setVisible(isVisible); }
    void fit();
    void saveImage();
    void fitRange();
    void zoom();
    void pan();
    void clearPanel();

    QList<TimeStepMarker> getTimeStepMarkerList() const { return TimeStepMarkerList; }

signals:

    void tracerInterceptsTime(double time);
};

#endif // QCUSTOMPLOTEXTENDED_H
